#!/usr/bin/perl -w

=head 1 SUMMARY

 clusterByClade.pl
 Usage: clusterByClade.pl --tree <file.tree> --aln <file.sc> --orig <reduced.fa.clstr> [--minID <95>]
 
 This script cluster sequences in a phylgentic tree so that all clusters are monophyletic clades (in addition to meeting the sequence ID threshold).
 
 Parameters:
     tree  - File with the tree in it (Newick format)
     aln   - File with the alignment scores produced by clustalw
     orig  - File with the original clusters (CD-Hit) that were used to produce the tree
     minID - Sequence ID threshold for clustering (default = 95%)

 Created by Chaim A. Schramm 2013-04-18
 Copyright (c) 2013 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
=cut

use strict;
use diagnostics;
use Pod::Usage;
use Getopt::Long;
use Bio::TreeIO;
use List::Util qw(min sum);

if ($#ARGV<0) { pod2usage(1); }

my ($treeFile, $scoreFile, $cdhitFile, $minID, $help) = ("", "", "", 95, 0);
GetOptions("tree=s"  => \$treeFile,
	   "aln=s"   => \$scoreFile,
	   "orig=s"  => \$cdhitFile,
	   "minID=i" => \$minID,
	   "help!"   => \$help
    );

if ($help || ! -e($treeFile) || ! -e($scoreFile) || ! -e($cdhitFile)) { pod2usage(1); }


#start by loading pairwise scores
my @idTable;
my %pairwiseScores;
open SC, $scoreFile or die "Can't read $scoreFile: $!\n";
while (<SC>) {
    if ($_ =~ /Sequence (\d+): (\S*)/) { 
	$idTable[$1] = $2;
	$pairwiseScores{$2}{$2} = 100; #by definition; makes cluster checking later a bit easier
    } elsif ($_ =~ /^Sequences \((\d+):(\d+)\) Aligned\. Score:\s+(\d+)$/) {
	$pairwiseScores{$idTable[$1]}{$idTable[$2]} = $3;
	$pairwiseScores{$idTable[$2]}{$idTable[$1]} = $3;
    }
}
close SC;


#now do the actual clustering
my @clusters;
my $input = new Bio::TreeIO(-file=>$treeFile, -format=>'newick');
my $tree = $input->next_tree;
my $rnode = $tree->get_root_node;
checkNode($rnode, \@clusters, \%pairwiseScores);


#load the original CD-Hit clusters
my %cdHitClusters;
my %currentCluster;
my ($rep, $repTime, $total) = ("","",0);
open CDH, "$cdhitFile" or die "Can't read $cdhitFile: $!\n";
while (<CDH>) {
    if ($_ =~ />Cluster/) {
	my @times = sort keys %currentCluster;
	next if $#times < 0; #takes care of first ">Cluster" line
#	if ($repTime ne $times[0]) {
#	    my @best = sort { $currentCluster{$times[0]}{$b} <=> $currentCluster{$times[0]}{$a} } keys %{$currentCluster{$times[0]}};
#	    my $rep = $best[0];
#	}
	unless ($rep eq "") { $cdHitClusters{$rep} = { 'last'=>$times[$#times], 'numReads'=>$total }; } #no rep just means not in subtree
        %currentCluster=();
        ($rep, $repTime, $total) = ("","",0);
    } elsif ($_ =~ />((\d\d\d).*)... (at (\d+\.\d+)%|\*)/) {
        #if ($3 eq "*") { $rep = $1; $repTime = $2; }
	if (-defined($pairwiseScores{$1})) { 
	    if ($rep ne "") { warn ("Two reps ($rep and $1)??\n"); }
	    else { $rep = $1; }
	}
        $currentCluster{$2}{$1} = $4;
	$total++;
    }

}
#do last cluster
my @times = sort keys %currentCluster;
if ($repTime ne $times[0]) {
    my @best = sort { $currentCluster{$times[0]}{$b} <=> $currentCluster{$times[0]}{$a} } keys %{$currentCluster{$times[0]}};
    my $rep = $best[0];
}
$cdHitClusters{$rep} = { 'last'=>$times[$#times], 'numReads'=>$total };
close CDH;


#print the output!
print "#Centroid\tSize\tPersistance\tTotalReads\n";
for my $cluster (sort { $b->{'size'} <=> $a->{'size'} } @clusters) {
    my ($latest, $numReads) = @{ findLastTimepoint($cluster->{'members'}, \%cdHitClusters) };
    print "$cluster->{'centroid'}\t$cluster->{'size'}\tweek$latest\t$numReads\n";
}


###########################################################


sub checkNode {
    my ($node, $clustersRef, $scoresRef) = @_;
    
    return if $node->is_Leaf; #don't care about singletons
    
    my @leaves;
    for my $desc ($node->get_all_Descendents) {
	if ($desc->is_Leaf) {
	    next if ($desc->id =~/(VRC|IG)/); #natives go in the reduced tree no matter what so don't care
	    push @leaves, $desc->id;
	}
    }

    my $thresholdMet = 1;
    for my $first (@leaves) {
	for my $second (@leaves) {
	    if ($scoresRef->{$first}{$second} < $minID) {
		#print "$first, $second: $scoresRef->{$first}{$second}\n";
		$thresholdMet = 0;
		last;
	    }
	}
	last unless $thresholdMet;
    }

    if ($thresholdMet && scalar(@leaves)>0) {

	#find the reads from the earliest possible time point
	my $earliest = min(map {substr($_,0,3)} @leaves);
	my @candidates = grep(/^$earliest/, @leaves);

	#This sorts the 'candidates' above by sum of the identity to *all* leaves in the clade.
	# The highest result represents the read that is most similar, on average, to all others, so will be the centroid.
	my @pickRep = sort { sum(@{$scoresRef->{$b}}{@leaves}) <=> sum(@{$scoresRef->{$a}}{@leaves}) } @candidates;

	push @{$clustersRef}, { 'centroid'=>$pickRep[0], 'size'=>scalar(@leaves), 'members'=>join(",",@leaves) };

    } else {
	for my $child ($node->each_Descendent) { checkNode($child, $clustersRef, $scoresRef); }
    }

}


sub findLastTimepoint {
    my ($list, $dictRef) = @_;

    my @members = split(/,/,$list);
    my $latest = '000';
    my $total = 0;
    for my $read (@members) {
	if ($dictRef->{$read}{'last'} > $latest) { $latest = $dictRef->{$read}{'last'}; }
	$total += $dictRef->{$read}{'numReads'};
    }

    return [$latest, $total];
}
	    
