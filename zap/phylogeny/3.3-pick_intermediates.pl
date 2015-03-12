#!/usr/bin/perl -w

=head 1 SUMMARY

 findIntermediates.pl
 Usage: findIntermediates.pl dnaml.out native 
 
 This script processes the output from DNAML to select developmental intermediates for the specified native mAb based on numbers of AA changes between nodes.
 
 Parameters:
     dnaml.out - output from dnaml (including inferred sequences)
     native    - name of native mAb to trace

 Created by Chaim A. Schramm 2013-05-27
 Copyright (c) 2013 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
=cut

use strict;
use diagnostics;
use Pod::Usage;
use Bio::Seq;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Align::PairwiseStatistics;
use POSIX;


if ($#ARGV<0 || $ARGV[0] =~ /-h/) { pod2usage(1); }

my %parent;
my %sequence;
open IN, $ARGV[0] or die "Can't read from $ARGV[0]: $!\n";

while (<IN>) {
#list of connections
#    45          IGHV3-30*1        0.00003     (     zero,    infinity)
#    45            39              0.00001     (     zero,     0.00585)
    if ($_ =~ /^\s+(\d+)\s+(\S+)\s+0/) {
	$parent{$2} = $1;
	#print "assigning $1 as parent of $2\n";

#sequence format
#   45        CAGGTGCAGC TGGTGGAGTC TGGGGGAGGC GTGGTCCAGC CTGGGAGGTC CCTGAGACTC
#   39        CAGGTGCAGC TGGTGGAGTC TGGGGGAGGC GTGGTCCAGC CTGGGAGGTC CCTGAGACTC
    } elsif ($_ =~ /^\s+(\S+)\s+([ACGTacgtyrwskm\- ]+)/) {
	my $id = $1;
	my $seq = $2; $seq =~ s/[ \-]//g; $seq =~ s/kkk//g;
	if (-defined($sequence{$id})) {
	    $sequence{$id} .= $seq;
	} else {
	    $sequence{$id} = $seq;
	}

    }
}
close IN;

my $current = $ARGV[1];
my @seqObjs;
#print "$current has parent $parent{$current}!\n\n";
#print join("*",%parent), "\n";
while ( -defined($parent{$current}) ) {
    push @seqObjs, Bio::Seq->new(-id=>$current, -seq=>$sequence{$current})->translate();
    $current = $parent{$current};
}
#for uca (so skip if we are doing 26a/b where it's not really UCA
#unless ($ARGV[1] =~ /VRC26b?$/) { 
push @seqObjs, Bio::Seq->new(-id=>$current, -seq=>$sequence{$current})->translate(); 
#}
#my $foo=Bio::SeqIO->new(-file=>">test.fa", -format=>'fasta'); for my $s (@seqObjs) { $foo->write_seq($s); }

#get alignments to count mutations
my @clustalParams = ('matrix'=>'BLOSUM', 'outorder'=>'input');
my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@clustalParams);
my $aln = $factory->align(\@seqObjs);

my $stats = Bio::Align::PairwiseStatistics->new();
my @intermediates = ($ARGV[1]);

my $native_to_uca = $aln->select_noncont(1,scalar(@seqObjs));
my $totalDist = $stats->number_of_differences($native_to_uca);
#for my $s ($native_to_uca->each_seq) { print sprintf(">%s\n%s\n", $s->id, $s->seq); }
print "$totalDist mutations from UCA to $ARGV[1]\n";

my $i = 1;
if ($totalDist >= 10) {
    my $numSteps = ceil($totalDist/10);
    my $stepDist = floor($totalDist/$numSteps);

    while ($i < $#seqObjs) {
	my $found = 0;
	for my $j ($i+1 .. $#seqObjs) {
	    my $pairwise = $aln->select_noncont($i, $j);
	    my $distance = $stats->number_of_differences($pairwise);
	    if ($distance >= $stepDist) {
		push @intermediates, $seqObjs[$j-1]->id;
		my $skip = $j-$i; print "Previous intermediate is $skip steps up ($distance AA mutations)\n";
		$i = $j;
		$found = 1;
		last;
	    }
	}
	last if !$found;
    }
}

if ($i < $#seqObjs) {
    my $uca_steps = $#seqObjs - $i; my $pairwise = $aln->select_noncont($i, $#seqObjs); my $distance = $stats->number_of_differences($pairwise); print "UCA is $uca_steps steps up ($distance AA mutations)\n";
}

push @intermediates, $seqObjs[$#seqObjs]->id;

for my $node (reverse(@intermediates)) {
    print ">$node\n$sequence{$node}\n";
}

