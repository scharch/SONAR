#!/usr/bin/env perl -w

=head 1 SUMMARY
 
 pruneTree.pl
 This program prunes unwanted/uninteresting branches from an existing 
       phylogenetic tree.


 Usage: pruneTree.pl --input <in.tree> --output <out.tree> 
                     (--keep <keep.list> | --remove <remove.list>)
                     [--annotate <data.txt>]

    Invoke with -h or --help to print this documentation. 
 
 Parameters:
     --input       => a file with the tree to be pruned
     --output      => where to save the pruned tree
     --keep        => a list of the nodes that should be kept 
                         (ALL others will be deleted)
     --remove      => a list of specifc nodes to be pruned 
                         (all others will be kept)
     --annotate    => a data file with annotations to pass to annotate the pruned
                         tree. Only works with the --keep option and if the 
                         annotations are in the default (2nd) column. Sets the 
                         output of that program to "out.tree_annotated". 
                         See annotateTree.pl for more details
 
 Created by Chaim A. Schramm 2013-02-28.
 
Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National
                          Institutes of Health, USA. All rights reserved.
 
=cut

use strict;
use diagnostics;
use Pod::Usage;
use Getopt::Long;
use Bio::TreeIO;
use Bio::Tree::Node;

if ($#ARGV < 0) { pod2usage(1); }

my ($infile, $outfile, $keep, $remove, $annotate, $help) = ("", "", "", "", "", 0);
GetOptions("input=s"    => \$infile,
           "output=s"   => \$outfile,
           "keep=s"     => \$keep,
           "remove=s"   => \$remove,
	   "annotate=s" => \$annotate,
           "help!"      => \$help
    );

if ($help || ! -e($infile) || (! -e($keep) && ! -e($remove))) { pod2usage(1); }
if ($keep && $remove) { die "Only *ONE* of keep, remove may be specified!\n"; }
if ($remove && $annotate) { die "The 'annotate' option may only be used when inputting a list of nodes to keep.\n"; }

=head2

First, read in the tree

=cut

my $input = new Bio::TreeIO(-file=>$infile, -format=>'newick');
my $tree = $input->next_tree;


=head2

now, get the list of nodes to remove

=cut

my @removeList = ();
if ($remove) {
    open RE, "$remove" or die "Can't read $remove: $!\n";
    while (<RE>) {
	chomp;
	my @a = split;
	push @removeList, $a[0];
    }
    close RE;
} else {
    my %keepList;
    open KEEP, "$keep" or die "Can't read $keep: $!\n";
    while (<KEEP>) {
	chomp;
	my @a = split;
	$keepList{$a[0]} = 1;
    }
    close KEEP;

    #invert keep list to get a remove list
    for my $node ($tree->get_nodes) {
	if ($node->is_Leaf && ! -defined($keepList{$node->id})) { push @removeList, $node->id; }
    }
}


=head2

remove unwanted nodes and save the resulting tree

=cut

for my $id (@removeList) { removeNode($id, 1); }

# remove any nodes with only 1 descendant (without changing total path length)
$tree = contract_paths($tree);

#output results
my $out = new Bio::TreeIO(-file=>">$outfile", -format=>'newick');
$out->write_tree($tree);


=head2

call annotateTree.pl, if desired

=cut

if ($annotate) {
    local @ARGV = ("--list", $keep, "--data", $annotate, "--input", $outfile, "--output", "$outfile\_annotated");
    my $result = do 'annotateTree.pl';
    if (! -defined($result)) { print "Error: " . $@ . "\n"; }
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


=head2

this subroutine is cribbed directly from the BioPerl code for TreeFunctionsI->contract_linear_paths()
 but I've added the flag '-perserve_lengths=>1' in the call to splice.
 also cut out code for rerooting, because I don't need it.

Original code called splice with an array of node refs. To use -preserve_lengths, I've switched to 
 using node->id, but this also requires a check to make sure that the name is not null

=cut

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
