#!/usr/bin/env perl

=head 1 SUMMARY

 3.4-cluster_tree.pl
 
 This script cluster sequences in a phylgentic tree so that all clusters are 
       monophyletic clades (in addition to meeting the sequence ID threshold).
 
 Usage: 3.4-cluster_tree.pl --nats nativeNames.txt --tree project.tree
                            --seq project-collected.fa --outtree clustered.tree
                            --outtable major_branches.txt --minID 90
                            --minSize 3 --format <list|table|fasta>
                            --label time1 --label time2 --label time3...
                            -h 

 Invoke with -h or --help to print this documentation.

 Required parameters:
     nats     - File with the native/known sequence IDs. These will be
                    preserved in the condensed tree without being clustered
                    into groups with NGS reads. It is recommended to include
                    the germline/outgroup on this list, as well. Sequence IDs
                    in the file that do not appear in the tree will be ignored.
                    Can be a list (default expectation), the first column in
                    a tab-delimited table, or a sequence (fasta) file.
                    Specify using the -format option (see below).

 Optional parameters:
     tree     - Newick-formatted representation of the tree to be clustered.
                    By default, pulls pipeline output (output/project.tree).
     seq      - Fasta file with CDR3 sequences for the sequences in the tree.
                    Recommended usage is manual specification with nucleotide
                    sequences. By default, pulls amino acide sequences from 
                    pipeline output (output/sequences/amino_acid/project-CDR3.fa).
     outtree  - Newick-formatted output file with clustered tree. Defaults
                    to output/project_clustered.tree.
     outtable - Tab-delimited text file with summary of clusters. Defaults
                    to output/tables/major_brances.txt.
     minID    - Clustering threshold, as percent sequence ID. Defaults to 90.
     minSize  - The minimum number of sequences to consider a cluster 
                    'significant' and to display in the condesed tree.
                    Defaults to 3.
     format   - Indicates format of the input file containing a list of
                    native/known sequences. Possible options are:
                         'list' (default) - a single column text file
                         'table' - multicolumn, tab-delimited text file
                                   with names in first column
                         'fasta' - a fasta-formatted sequence file with
                                   the names as sequence IDs
     label    - Sequence ID prefixes indicating the time point from which a
                    particular NGS read originated. Used to calculate the
                    temporal persistence of a major cluster in the condensed
                    tree. Must be in proper chronological order. Missing
                    labels will cause program to crash. If *no* labels are
                    provided, calculation will be skipped and "nd" will be
                    output in persistence column.

 Created by Chaim A. Schramm 2013-04-18
 Edited and commented for publication by Chaim A. Schramm on 2015-08-13

 Copyright (c) 2013 Columbia University and Vaccine Research Center, National
                     Institutes of Health, USA. All rights reserved.

=cut

use warnings;
use strict;
use diagnostics;
use Pod::Usage;
use Getopt::Long;
use Cwd;
use File::Basename;
use Switch;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Align::PairwiseStatistics;
use Bio::TreeIO;
use List::Util qw(min sum);
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);
use version qw/logCmdLine/;

if ($#ARGV<0) { pod2usage(1); }

my $prj_name = basename(getcwd);

#set defaults
my @saveArgs = @ARGV;
my $natFile   = "";
my $treeFile  = "output/$prj_name.tree";
my $seqFile   = "output/sequences/amino_acid/$prj_name-CDR3.fa";
my $outTree   = "output/$prj_name\_clustered.tree";
my $outFile   = "output/tables/major_branches.txt";
my $minID     = 90;
my $minSize   = 3;
my $natFormat = "list";
my @labels;
my $help = 0;

GetOptions("tree=s"     => \$treeFile,
	   "seq=s"      => \$seqFile,
	   "nats=s"     => \$natFile,
	   "outtable=s" => \$outFile,
	   "outtree=s"  => \$outTree,
	   "minID=i"    => \$minID,
	   "minSize=i"  => \$minSize,
	   "format=s"   => \$natFormat,
	   "label=s"    => \@labels,
	   "help!"      => \$help
    );


if ( $help || ! -e($treeFile) || ! -e($seqFile) ) { pod2usage(1); }

&logCmdLine($0,@saveArgs);


#invert labels
my %labelOrder;
for my $i (0 ..$#labels) { $labelOrder{$labels[$i]} = $i; }

#start by loading natives so we don't count them
my %natives;
if ( $natFile ne "" ) {
    switch ($natFormat) {
	case "list"  { open IN, $natFile or die "Can't read $natFile: $!";
		       while (<IN>) {
			   chomp;
			   $natives{$_} = 1;
		       }
		       close IN;
	}
	case "table" { open IN, $natFile or die "Can't read $natFile: $!";
		       while (<IN>) {
			   my @a = split;
			   $natives{$a[0]} = 1;
		       }
		       close IN;
	}
	case "fasta" { my $in = Bio::SeqIO->new( -file=>$natFile );
		       while (my $seq=$in->next_seq) {
			   $natives{$seq->id} = 1;
		       }
	}
	else         { die "Cannot recognize specified format '$natFormat'.\nAllowed values are 'list' (default), 'table', or 'fasta'." }
    }
}


#loading sequences and score them
my $in = Bio::SeqIO->new(-file=>$seqFile);
my @cdr3seqs;
while (my $seq = $in->next_seq) {
    push @cdr3seqs, $seq;
}
my %pairwiseScores;
for my $i (0 .. $#cdr3seqs) {

    #fill in diagonal for later ease
    my $seq1 = $cdr3seqs[$i];
    $pairwiseScores{$seq1->id}{$seq1->id} = 100;

    for my $j ($i+1 .. $#cdr3seqs) {
	my $seq2 = $cdr3seqs[$j];
	my $seqObjs = [ $seq1, $seq2 ];
	my $factory = Bio::Tools::Run::Alignment::Muscle->new();
	my $align = $factory->align($seqObjs);
	$pairwiseScores{$seq1->id}{$seq2->id} = $align->percentage_identity;
	$pairwiseScores{$seq2->id}{$seq1->id} = $align->percentage_identity;
    }

}



#now do the actual clustering
my @clusters;
my $input = new Bio::TreeIO(-file=>$treeFile, -format=>'newick');
my $tree = $input->next_tree;
my $rnode = $tree->get_root_node;
checkNode($rnode, \@clusters, \%pairwiseScores, \%natives);


#print the output!
open OUT, ">$outFile" or die "Can't write to $outFile: $!\n";
print OUT "Centroid\tSize\tPersistence\n";
for my $cluster ( sort { $b->{'size'} <=> $a->{'size'} } @clusters ) {
    my $persist = "nd";
    if ($#labels > 0) {
        my $date = findLastTimepoint($cluster->{'members'}, \%labelOrder);
        $persist = "$labels[$date->[0]]-$labels[$date->[1]]";
    }
    print OUT "$cluster->{'centroid'}\t$cluster->{'size'}\t$persist\n";
}
close OUT;


#fix up tree and save condensed version
$tree = contract_paths($tree);
my $out = new Bio::TreeIO(-file=>">$outTree", -format=>'newick');
$out->write_tree($tree);



###########################################################


sub checkNode {
    my ($node, $clustersRef, $scoresRef, $natsRef) = @_;

    if ($node->is_Leaf) {
	# singleton; delete unless it is a native
	&removeNode($node->id, 1) unless -defined($natsRef->{$node->id});
	return;
    }

    
    my @leaves;
    for my $desc ($node->get_all_Descendents) {
	if ($desc->is_Leaf) {
	    next if ( -defined($natsRef->{$desc->id}) ); #natives go in the reduced tree no matter what so don't care
	    push @leaves, $desc->id;
	}
    }

    my $thresholdMet = 1;
    for my $first (@leaves) {
	for my $second (@leaves) {
	    if ($scoresRef->{$first}{$second} < $minID) {
		$thresholdMet = 0;
		last;
	    }
	}
	last unless $thresholdMet;
    }

    if ($thresholdMet && scalar(@leaves)>$minSize) {

	#This sorts the sequences in the cluster by sum of the identity to *all* leaves in the clade.
	# The highest result represents the read that is most similar, on average, to all others, so will be the centroid.
	my @pickRep = sort { sum(@{$scoresRef->{$b}}{@leaves}) <=> sum(@{$scoresRef->{$a}}{@leaves}) } @leaves;

	#save centroid
	push @{$clustersRef}, { 'centroid'=>$pickRep[0], 'size'=>scalar(@leaves), 'members'=>join(",",@leaves) };

	#mark other members of cluster for removal
	for my $n (1 .. $#pickRep) { &removeNode($pickRep[$n], 1); }

    } else {
	for my $child ($node->each_Descendent) { checkNode($child, $clustersRef, $scoresRef, $natsRef); }
    }

}


sub findLastTimepoint {
    my ($list, $dictRef) = @_;

    my @members = split(/,/,$list);
    my $earliest = 999;
    my $latest   = 0;
    for my $read (@members) {
	my ($time) = $read =~ /^(.*)-/;
	if ($dictRef->{$time} > $latest)   { $latest   = $dictRef->{$time}; }
	if ($dictRef->{$time} < $earliest) { $earliest = $dictRef->{$time}; }
    }

    return [$earliest, $latest];
}
	    

=head2

this subroutine calls Bio::Tree's own delete-node routine, but it also checks ancestor nodes to remove those that no longer have any descendents

=cut
sub removeNode {
    my ($id, $leaf) = @_;
    my $node = $leaf ? $tree->find_node(id=>$id) : $tree->find_node(internal_id=>$id);

    #can't find node, it's already gone...
    unless (-defined($node)) { print "Node $id is already gone\n"; return; }

    #never remove the root node
    return unless -defined($node->ancestor);

    my $ancestor = $node->ancestor;
    
    $tree->remove_Node($node);
    if (scalar $ancestor->each_Descendent == 0) { removeNode($ancestor->internal_id, 0); }
}


sub contract_paths {
    my $tree = shift;
    my @remove;
    foreach my $node ($tree->get_nodes) {
        if ($node->ancestor && $node->each_Descendent == 1) {
	    $node->id("myUnnamedNode".$node->internal_id) unless -defined($node->id);
            push(@remove, $node->id);
        }
    }
    $tree->splice(-remove_id => \@remove, -preserve_lengths => 1) if @remove;
    return $tree

}
