#!/usr/bin/env perl

=head 1 SUMMARY
 
 annotateTree.pl

 This program renames the leaves of a phylogenetic tree to include other 
       information (such as CDR3 sequence or number of raw reads represented).


 Usage: annotateTrees.pl --data <file.txt> --input <in.tree> --output <out.tree>
                         [ --ncol <0> --dcol <1> --[no]natives ]

     Invoke with -h or --help to print this documentation. 
 
 Required Parameters:
     --data        => A tab-delimited file with a list of sequences and the
                         annotation for each one.
     --input       => A file with the tree to be annotated in Newick format.
     --output      => Where to save the annotated tree.

 Optional Parameters:
     --ncol        => The column number that contains the sequence names.
                         Default = 0 (first column).
     --dcol        => The column number that contains the annotations.
                         Default = 1 (second column).
     --[no]native  => Annotate the known antibody sequences, as well. Default = Y.
                         For now, assumes any sequnce with an id that is not
                         strictly numeric is a known anitbody.
 
 Created by Chaim A. Schramm 2013-02-28.
 Edited and commented for publication by Chaim A Schramm on 2015-05-04.

 Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National
                          Institutes of Health, USA. All rights reserved.
=cut


use warnings;
use strict;
use diagnostics;
use Pod::Usage;
use Getopt::Long;
use Bio::TreeIO;
use Bio::Tree::Node;
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);
use version qw/logCmdLine/;

my @saveArgs = @ARGV;
my ($data, $infile, $outfile, $ncol, $dcol, $native, $help) = ("", "", "", 0, 1, 1, 0);
GetOptions("data=s"   => \$data,
	   "input=s"  => \$infile,
	   "output=s" => \$outfile,
	   "ncol=i"   => \$ncol,
	   "dcol=i"   => \$dcol,
	   "native!"  => \$native,
	   "help!"    => \$help
    );


if ($help || ! -e($data) || ! -e($infile) || $ncol<0 || $dcol<0) { pod2usage(1); }


&logCmdLine($0,@saveArgs);


=head2

first, read in the tree

=cut
my $input = new Bio::TreeIO(-file=>$infile, -format=>'newick');
my $tree = $input->next_tree;




=head2

now get the annotations from the data file

=cut
my %annotations;
open DAT, "$data" or die "Can't read $data: $!\n";
while (<DAT>) {
    next if ($_ =~ /^#/);
    chomp;
    my @a = split /\t/,$_;
    if ( $native || $a[$ncol] =~ /^\d+$/ ) {
	$annotations{$a[$ncol]} = $a[$dcol]; #saves the annotation
    }
}
close DAT;



=head2

now do the actual annotation

=cut
for my $id (keys %annotations) {

    my $node = $tree->find_node(id=>$id);
    if (-defined $node) {
	$node->id("$id          $annotations{$id}");
    } else {
	warn("$id: No such node was found!\n");
    }

}



=head2

Finally, write out the annotated tree

=cut
my $out = new Bio::TreeIO(-file=>">$outfile", -format=>'newick');
$out->write_tree($tree);
