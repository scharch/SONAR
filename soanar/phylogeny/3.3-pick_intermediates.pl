#!/usr/bin/env perl -w

=head 1 SUMMARY

 3.3-pick_intermediates.pl
 
 This script processes the output from DNAML to select developmental intermediates 
       for the specified native mAb based on numbers of AA changes between nodes.
 
 Usage: 3.3-pick_intermediates.pl --target mAb1 [ --target mAb2 --target mAb3....
                                                  --numSteps 3 --stepDist 10
                                                  --input dnaml.out 
                                                  --output intermediates.fa
                                                  -h ]

 Invoke with -h or --help to print this documentation.

 Parameters:
     dnaml.out - output from dnaml (including inferred sequences)
     native    - name of native mAb to trace

 Created by Chaim A. Schramm 2013-05-27
 Edited and commented for publication with extensive modifications by Chaim A
         Schramm on 2015-07-16

 Copyright (c) 2011-2015 Columbia University and Vaccine Research Center, National
                          Institutes of Health, USA. All rights reserved.

=cut

use strict;
use diagnostics;
use Pod::Usage;
use Getopt::Long;
use Bio::Seq;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Align::PairwiseStatistics;
use Algorithm::Combinatorics qw/combinations/;
use Statistics::Basic qw/stddev/;
use List::Util qw/min max sum/;
use POSIX qw/floor getcwd/;
use File::Basename;



if ($#ARGV < 0) { pod2usage(1); }

my $prj_name = basename(getcwd);
my @targets;
my ($dnaml, $numSteps, $stepDist, $output, $help) = ( "output/logs/$prj_name.dnaml.out", 0, 0, "", 0 );
GetOptions("input=s"    => \$dnaml,
	   "target=s"   => \@targets,
	   "numSteps=i" => \$numSteps,
           "stepDist=i" => \$stepDist,
	   "output=s"   => \$output,
           "help!"      => \$help
    );


#check inputs
if ($help) { pod2usage(1); }

if ($#targets < 0) {
    die( "Error: please specify at least one target (run \"3.3-pick_intermediates.pl -h\" for help)\n\n" );
}
if (! -e $dnaml) {
    die( "Error: cannot find input file $dnaml!\n\n" );
}
$stepDist = 10 if $numSteps == 0 && $stepDist == 0;
if ($numSteps > 0 && $stepDist > 0) {
    die( "Error: \"numSteps\" and \"stepDist\" are mutually exclusive options (run \"3.3-pick_intermediates.pl -h\" for help)\n\n" );
}
if ($output eq "") {
    if ($#targets == 0) { $output = "output/sequences/nucleotide/$targets[0]-intermediates.fa"; }
    else { $output = "output/sequences/nucleotide/" . join("-", @targets) . "_MRCA.fa"; }
}    
open OUT, ">$output" or die "Can't write to $output: $!\n";


#set up important variables
my %nodes;
my %path;
my $root = "";
my $germline = "";
my $uca = "";


open IN, $dnaml or die "Can't read from $dnaml: $!\n";
while (<IN>) {

# format for list of connections in DNAML output
#    45          IGHV3-30*1        0.00003     (     zero,    infinity)
#    45            39              0.00001     (     zero,     0.00585)

    if ($_ =~ /^\s+(\d+)\s+(\S+)\s+0/) {

	my $pid = $1;
	my $id  = $2;

	#create root
	if ($root eq "") {
	    $root = $pid;
	    $nodes{$root} = { 'id'=>$root, 'sequence'=>"" }; #id is for accessing a node through a parent or child
	}

	$nodes{$id} = { 'id'=>$id, 'sequence'=>"", 'parent'=>$nodes{$pid} };
	push @{$nodes{$pid}{'children'}}, $nodes{$id};


	if ($id =~ /(IG|VH|VK|VL|HV|KV|LV)/) { $germline= $id; }
	elsif ($pid eq $root) { $uca = $id; }

# format for input or inferred sequences in DNAML output
#   45        CAGGTGCAGC TGGTGGAGTC TGGGGGAGGC GTGGTCCAGC CTGGGAGGTC CCTGAGACTC
#   39        CAGGTGCAGC TGGTGGAGTC TGGGGGAGGC GTGGTCCAGC CTGGGAGGTC CCTGAGACTC
# at least some IUPAC ambiguity codes are allowed, so we have to account for that
# also, lower confidence nucleotides are represented in lower case, so use /i

    } elsif ($_ =~ /^\s+(\S+)\s+([acgtyrwskmbdhvn\- ]+)/i) {
	$nodes{$1}{'sequence'} .= $2;
    }


}
close IN;


#check that the mAbs we want are present
for my $t (@targets) {
    if (! -defined($nodes{$t})) { 
	die("Can't find desired mAb $t in tree from file $dnaml\n")
    }
}


#check to make sure tree is rooted correctly and fix if necessary
if ( -defined( $nodes{$germline}{'parent'}{'parent'} ) ) {
    my $oldParent = $nodes{$germline}{'parent'};
    $nodes{'new root'} = { 'id'       => 'new root',
			   'sequence' => $oldParent->{'sequence'}, # use sequence from what will become the new UCA.
			                                          # shouldn't be relevant, but just to prevent possible errors
			   'children' => [ $nodes{$germline}, $oldParent ] };
    $nodes{$germline}{'parent'} = $nodes{"new root"};
    $oldParent->{'newParent'}   = $nodes{"new root"};
    &reroot($oldParent);
    delete( $nodes{$root} );
    $root = "new root";
    $uca  = $oldParent->{'id'};
}


&assignGaps($nodes{$root});


#get the all intermediates on the pathway to first target
#  file output is to reorder correctly after muscle alignment below
my $order = 0;
my @seqObjs;
my $current = $nodes{$targets[0]};
open TMP, ">work/phylo/order.txt" or die "Can't write to work/phylo/order.txt: $!\n\n";
while ( -defined($current->{'parent'}) ) {
    print TMP $current->{'id'} . "\n";
    $path{$current->{'id'}} = $order++;
    my $sequence = &degapSequence($current);
    push @seqObjs, Bio::Seq->new(-id=>$current->{'id'}, -seq=>$sequence)->translate();
    $current = $current->{'parent'};
}
close TMP;


if ($#targets > 0) {

    #multiple targets, find MRCA
    my @mrca;
    for my $t ( 1 .. $#targets ) {
	my $current = $nodes{$targets[$t]};
	while ( -defined($current->{'parent'}) ) {
	    if (-defined($path{$current->{'id'}})) {
		push @mrca, $current->{'id'}; #mrca of this pair
		last;
	    }
	    $current = $current->{'parent'};
	}
    }

    my @sorted = sort { $path{$b} <=> $path{$a} } @mrca; #pairwise mrca that is highest on the tree is the overall mrca
    my $sequence = &degapSequence($nodes{$sorted[0]});

    print OUT ">Node$sorted[0] MRCA of " . join(",",@targets) . "\n$sequence\n";

    close OUT; #end of program (rest in else statement for single target)

} else {

    #get alignments to count mutations 
    # (muscle no longer has an option to preserve input order, so use tmp file to re-sort)
    my $factory = Bio::Tools::Run::Alignment::Muscle->new();
    my $preAlign = $factory->align(\@seqObjs)->sort_by_list("work/phylo/order.txt");

    #eliminate consecutive nodes with no differences
    my ($aln, $uniqRef) = $preAlign->uniq_seq;

    my @diffCounts;
    my $stats = Bio::Align::PairwiseStatistics->new();
    for my $i (1 .. $aln->num_sequences-1) {
	for my $j ($i+1 .. $aln->num_sequences) {
	    my $pairwise = $aln->select_noncont($i, $j);
	    my $distance = $stats->number_of_differences($pairwise);
	    $diffCounts[$i][$j] = $distance;
	}
    }

    my @intermediates;

    if ($stepDist) {
	# easiest way to do this is to approximate the number of steps involved
	# start by counting the maximum possible distance if we step to every intermediate
	my $maxDist = sum( map { $diffCounts[$_][$_+1] } (1..$#diffCounts) );
	$numSteps = floor( $maxDist / $stepDist );
    }


    #error checking
    #make sure there is at least 1 intermediate and not more than there are on the tree
    $numSteps = 1 if $numSteps < 1;
    $numSteps = scalar(@diffCounts) if $numSteps > scalar(@diffCounts);


    # Now go through all possible paths of N steps
    # Choose the one with the most equal step sizes (smallest standard deviation)
    my $minStDev = 9999;

    #it's n-1 because starting (target) and ending (uca) points are fixed
    my $iter = combinations( [2 .. $#diffCounts], $numSteps-1 );
    while (my $combo = $iter->next) {
	#add the endpoints back in
	my @steps = ( 1, sort {$a<=>$b} @$combo, scalar(@diffCounts) );
	
	#now get the distance of each step
	my @stepSizes;
	for my $first (0 .. $#steps-1) {
	    push @stepSizes, $diffCounts[$steps[$first]][$steps[$first+1]];
	}
	my $stdev = stddev(@stepSizes);
	if ($stdev < $minStDev) {
	    $minStDev = $stdev;
	    @intermediates = @steps;
	}
    }

    my $prevID = "";
    for my $ind (reverse(0 .. $#intermediates)) {

	# NOTE
	# @intermediates is just a list of numbers pointing to postions in @diffCounts
	#      (ie sequences in the unique-ified alignment)
	# From there, we need to go back to $uniqRef to get the original ids to be
	#      able to look them up in %path
	# Also, the IDs in the alignment are "ST1", "ST2", etc, but in the lookup they
	#      are just "1", "2", ...,  hence the substring

	my $uniqID         = substr( $aln->get_seq_by_pos($intermediates[$ind])->id, 2 );
	my @original       = @{$uniqRef->{$uniqID}};
	my @sortByPosition = sort { $path{$b} <=> $path{$a} } @original;
	my $id             = $sortByPosition[0]; 
	$prevID = $id if $prevID eq ""; #initiate

	#setup vars that will be used for actual output below
	my $intNum = "intermediate" . ($#intermediates - $ind);
	my $desc   = "";
	my $seq    = "";

	my $changes = $ind < $#intermediates ? $diffCounts[$intermediates[$ind]][$intermediates[$ind+1]] : ""; #mutations to UCA is undefined
	my $steps  = "";
	if ($ind < $#intermediates) {  #steps to UCA is undefined
	    $steps = $path{$prevID} - $path{$id};
	    $prevID = $id;
	}


	if ($ind == $#intermediates) {
	    $intNum = "UCA";
	    $desc = "Node$id $diffCounts[1][$#{$diffCounts[1]}] AA changes from $targets[0]"; #would be interesting to include changes from germline
	    $seq = &degapSequence($nodes{$id});
	} elsif ($ind == 0) {
	    $intNum = $targets[0];
	    $desc = "$steps steps on tree and $changes AA changes from previous intermediate";
	    $seq = &degapSequence($nodes{$targets[0]});
	} else {
	    $desc = "Node$id: $steps steps on tree and $changes AA changes from previous intermediate";
	    $seq = &degapSequence($nodes{$id});
	}

	print OUT ">$intNum $desc\n$seq\n";

    }
    close OUT;
    

} #ends if/else for MRCA vs single pathway intermediates.




sub reroot {
    
    my ($node) = @_;

    if ( -defined($node->{'parent'}) ) {
	my $parent = $node->{'parent'};

	#this was passed up (now down) from the previous node
	#  complete the transformation
	$node->{'parent'} = $node->{'newParent'};
	delete( $node->{'newParent'} );

	#add old parent to children list
	push @{$node->{'children'}}, $parent;

	#remove myself from old parent's children list
	my @delIndex = grep { $parent->{'children'}[$_] == $node } ( 0 .. $#{$parent->{'children'}});
	splice( @{$parent->{'children'}}, $delIndex[0], 1 );

	#tell old parent that I am now its parent
	$parent->{'newParent'} = $node;

	#call next round
	reroot( $parent );
    } else {
	#no parent means that this was the original root, which will be removed
	#  transfer other child(ren) to new parent
	for my $child (@{$node->{'children'}}) {
	    $child->{'parent'} = $node->{'newParent'};
	    push @{$node->{'newParent'}{'children'}}, $child;
	}
    }

}



sub assignGaps {

    my ($node) = @_;

    #first remove spaces
    $node->{'sequence'} =~ s/\s+//g;

    #no children means leaf, get gaps from actual sequence
    unless ( -defined($node->{'children'}) ) {

	while ( $node->{'sequence'} =~ /(-+)/g ) {
	    my $gid = sprintf("%03d-%03d", $-[1], $+[1]);

	    #don't count terminal gap on V gene
	    last if $node->{'id'} =~ /(IG|HV|KV|LV|VH|VK|VL)/i && $+[1] == length($node->{'sequence'});

	    $node->{'gaps'}{$gid} = 1; #1 means a real gap ( 0.5 for provisional )
	}
	
    } else {

	#for an internal node, start by processing all children and counting their gaps
	my %possibleGaps;
	for my $child (@{$node->{'children'}}) {
	    assignGaps($child);
	    for my $g (keys %{$child->{'gaps'}}) {

		if ( -defined($possibleGaps{$g}) ) {
		    $possibleGaps{$g} = $possibleGaps{$g} + $child->{'gaps'}{$g};
		} else {
		    $possibleGaps{$g} = $child->{'gaps'}{$g};
		}

	    }
	}

	# sort the gaps, to make detecting overlaps easier
	my @sortedGaps = sort keys %possibleGaps;

	# For each gap, a score >1 means it appears in both children, and at least
	#    one is definitive. In this case, we keep the gap as definitive, which  
	#    has the (desired) side effect of also fixing the gap in any descendants 
	#    in which it was previously marked as provisional. (This could also 
	#    happen for the occasional polytomy, if 3 children each have the gap 
	#    provisionally, but that combination should be rare enough that I'm
	#    just not going to worry about it.)
	# A score of 1 means either provisional in both children or definitive in
	#    one and not present in the other. In either case, we mark it as 
	#    provisional for this node.
	# Overlapping gaps get parsed (VERY KLUDGILY) into shared and unshared 
	#    pieces, each processed as above.

	for (my $ind=0; $ind<scalar(@sortedGaps); $ind++) {
	    my $thisGap = $sortedGaps[$ind]; 
	    my ($start, $stop) = split(/-/,$thisGap);
	    my $pass = 0; #a tracking variable for overlapping gaps

	    if ($possibleGaps{$thisGap} > 1) {

		$node->{'gaps'}{$thisGap} = 1;

	    } else {

		#check for other gaps overlapping with this one
		while ( $ind < $#sortedGaps  &&  (split(/-/,$sortedGaps[$ind+1]))[0] < $stop ) {
		    
		    my ($newStart, $newStop) = split(/-/,$sortedGaps[$ind+1]);

		    #process non-overlapping 5' gap
		    if ($newStart > $start) {
			if ($possibleGaps{$thisGap} == 1) {
			    # it was definitive in child, so keep it here as provisional
			    $node->{'gaps'}{sprintf("%03d-%03d", $start, $newStart)} = 0.5;
			} else {
			    # it was provisional below, so delete it
			    &deleteProvisionalGap($node, $start, $newStart);
			}
		    }

		    #process overlapping gap
		    my $lowerEnd = min( $stop, $newStop );
		    my $overlapScore = $possibleGaps{$thisGap} + $possibleGaps{$sortedGaps[$ind+1]};
		    if ($overlapScore > 1) {
			$node->{'gaps'}{sprintf("%03d-%03d", $newStart, $lowerEnd)} = 1;
		    } else {
			$node->{'gaps'}{sprintf("%03d-%03d", $newStart, $lowerEnd)} = 0.5;
		    }

		    # Reset boundaries and try again
		    if ($stop == $newStop) { 
			$ind++;    #we've taken care of both this one and the next one
			$pass = 1; #record keeping for below so we don't write it twice
			last;      #stop looking for overlaping gaps
		    } else {
			$start = $lowerEnd;
			$stop = max($stop,$newStop);
			$thisGap = $sortedGaps[++$ind];
		    }

		} #no more overlapping gaps

		if ($pass) { 
		    $pass = 0;
		    next;
		}

		if ($possibleGaps{$thisGap} == 1) {
		    # it was definitive in one child, or provisional in both
		    #     so keep it here as provisional
		    $node->{'gaps'}{sprintf("%03d-%03d", $start, $stop)} = 0.5;
		} else {
		    # it was provisional below, so delete it
		    &deleteProvisionalGap($node, $start, $stop);
		}
		
	    }

	} #done iterating through possible gaps

    } #done with unless( children )

}


sub deleteProvisionalGap {

    my ($node,$start,$stop) = @_;

    for my $child (@{$node->{'children'}}) {
	for my $gap (keys %{$child->{'gaps'}}) {

	    my ($currentStart, $currentStop) = split(/-/, $gap);
	    #explicitly check provisional status in addition to boundaries
	    # allows me to call the delete function recursively without worrying
	    if ( $start <= $currentStart && $stop >= $currentStop && $child->{'gaps'}{$gap} == 0.5 ){

		#first, adjust boundaries to save part of current gap, if necessary
		if ($start < $currentStart) {
		    $child->{'gaps'}{"$start-$currentStart"} = 0.5;
		}
		if ($stop > $currentStop) {
		    $child->{'gaps'}{"$currentStart-$stop"} = 0.5;
		}
		 
		#actual delete
		delete $child->{'gaps'}{$gap};

		#and recurse
		deleteProvisionalGap($child, $start, $stop);

	    }

	}
    }

}
	      
	    

sub degapSequence {

    my ($node) = @_;
    my $sequence = $node->{'sequence'};
    my @gaps = keys %{$node->{'gaps'}};

    for my $g (sort { $b cmp $a } @gaps) { #lexical sort but in reverse order so we don't screw up indices splicing out the gaps
	my ($start, $stop) = split /-/, $g;
	my $len = $stop - $start;
	substr($sequence, $start, $len, "");
    }

    my $check  = Bio::Seq->new(-id=>$node->{'id'}, -seq=>$sequence)->translate();
    if ($check->seq =~ /\*/) {
	print "Whoops! I seemed to have mis-spliced the gaps from $node->{'id'} and introduced stop codons!\n\tI will continue, but please be aware of the likely error!n\n";
    }

    return $sequence;

}
