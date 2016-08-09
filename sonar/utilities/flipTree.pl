#!/usr/bin/env perl -w

=head 1 SUMMARY
 
 flipTree.pl

 This program reorders the branches of a phylogenetic tree so that the
       longest branches are toward the bottom.

 Usage: flipTree.pl <in.tree> <out.tree> [<root>]

     Invoke with -h or --help to print this documentation. 
 
 Required Parameters:
     <in.tree>  => A file with the tree to be annotated in Newick format.
     <out.tree> => Where to save the annotated tree.

 Optional Parameters:
     <root>     => Name of sequence at which to re-root the tree
 
 Created by Chaim A. Schramm 2013-02-28.
 Edited and commented for publication by Chaim A Schramm on 2016-04-21 -happy birthday, Mom!

 Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National
                          Institutes of Health, USA. All rights reserved.
=cut

use strict;
use Pod::Usage;
use Bio::TreeIO;
use Bio::Tree::Node;
use Bio::Tree::TreeFunctionsI;
use diagnostics;
use List::Util qw(min sum);
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);
use version qw/logCmdLine/;


if ($#ARGV<2 || $ARGV[0] =~ /-h/) { pod2usage(1); }
&logCmdLine($0,@ARGV);


#load input tree
my $input = new Bio::TreeIO(-file=>$ARGV[0], -format=>'newick');
my $tree = $input->next_tree;


#create custom sort function
my $best3_sort = sub {
    my ($a1, $b1) = @_;
    return (threeDeep(\$a1) <=> threeDeep(\$b1));
};


#re-root if desired
if ($#ARGV >= 2) {
    my @nodes = $tree->find_node("$ARGV[2]");
    if ($#nodes >= 0) { $tree->reroot_at_midpoint($nodes[0]); }
}


#and output
my $out = new Bio::TreeIO(-file=>">$ARGV[1]", -format=>'newick');
$out->get_params->{'order_by'}=$best3_sort;
$out->write_tree($tree);


#this is the actual sorting function
# assign each internal node a value based on the depth of the three
# most distant child leaves
sub threeDeep {
    my ($node) = @_;
    
    #if it's a leaf node, there's nothing to do
    return ${$node}->branch_length if (${$node}->is_Leaf);

    #otherwise
    my @leaf_distances =();
    for my $child (${$node}->get_all_Descendents) {
	next unless $child->is_Leaf;
	push @leaf_distances, $child->depth-${$node}->depth;
    }

    #prevent an error if this node only has 1 or 2 children
    my $num2avg = min($#leaf_distances, 2); #will be 0-indexed

    my @sorted = sort {$b <=> $a} @leaf_distances;
    return sum(@sorted[0..$num2avg])/($num2avg+1) + ${$node}->branch_length;
}


#left over error checking function which prints information about a
# node and its descendants
sub printNode {
    my $thisNode = shift;
    if (! $thisNode->is_Leaf) { $thisNode->id($thisNode->internal_id); }
    print $thisNode->id . ": " . $thisNode->height . " plus " . $thisNode->depth . "\n";
    if (! $thisNode->is_Leaf) {
	my @children = $thisNode->each_Descendent('height');
	for my $otherNode (@children) { &printNode($otherNode); }
    }
}
