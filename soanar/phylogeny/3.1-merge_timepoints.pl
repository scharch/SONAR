#!/usr/bin/perl -w

=head 1 SUMMARY

 3.1-merge_timepoints.pl

 This script uses USearch to recluster selected sequences from multiple timepoints
       to determine "birthdays" and persistence times.

 Usage: 3.1-merge_timepoints.pl --seqs time1.fa --seqs time2.fa ...
                                  [ --labels t1 --labels t2 ... --t 100 -f]

 Invoke with -h or --help to print this documentation.

 Parameters:
     --seqs   => list of files with selected sequences from each timepoint.
                     Must be in proper chronological order.
     --labels => Optional list of timepoint IDs to be used in relabeling
                     input sequences. Must be in proper chronological order.
                     Defaults to "01-", "02-", etc.
     --t      => Optional % threshold for clustering reads across time points.
                     Defaults to 100.

 Created by Chaim A. Schramm 2015-06-12.

 Copyright (c) 2011-2015 Columbia University and Vaccine Research Center, National
                          Institutes of Health, USA. All rights reserved.

=cut

use strict;
use diagnostics;
use Pod::Usage;
use Getopt::Long;
use List::Util qw/max/;
use PPvars qw(ppath);
use Cwd;
use File::Basename;
use Bio::SeqIO;


if ($#ARGV < 0) { pod2usage(1); }

my @seqFiles   = ( );
my @prefix     = ( );
my ($t, $force, $help) = ( 100, 0, 0 );
GetOptions("seqs=s{1,}" => \@seqFiles,
	   "labels=s"   => \@prefix,
	   "t=i"        => \$t,
           "force!"     => \$force,
           "help!"      => \$help
    );

if ($help) { pod2usage(1); }
my $threshold = $t/100;

#check to make sure inputs match
if ($#prefix >= 0) {
    if ($#prefix != $#seqFiles) {
	print "\nCurrent input:\n\tfile\tprefix\n";
	for (my $i=0; $i<=max($#prefix, $#seqFiles); $i++) {
	    my $name = $prefix[$i] || "??";
	    my $file = $seqFiles[$i] || "??";
	    print "\t$file\t$name\n";
	}
	die("Please make sure that the number of custom prefixes matches the number of sequence files!\n");
    }
} else { #no prefixes provided, use defaults
    for (my $i=0; $i<=$#seqFiles; $i++) {
	push @prefix, sprintf("%02d", $i+1);
    }
}

#invert array into a hash so we have the right order
my %order;
for (my $i=0; $i<=$#prefix; $i++) { $order{$prefix[$i]} = $i; }



#setup directory structure
my $prj_name = basename(getcwd);

for my $dir ("output", "output/sequences", "output/sequences/nucleotide", 
	     "output/sequences/amino_acid", "output/tables", "work",
	     "work/phylo") {
    if ( ! -d $dir ) { mkdir $dir; }
}

if ( -e "output/sequences/nucleotide/$prj_name-collected.fa" && ! $force ) {
    die "Output already exists. Please use the -f(orce) option to re-intiate and overwrite.\n";
}


#Read in selected seqs, rename, and write to temporary output
my $all = Bio::SeqIO->new(-file=>">work/phylo/all_seqs.fa", -format=>'fasta');
$all->width(600);
for (my $i=0; $i<=$#seqFiles; $i++) {
    my $in = Bio::SeqIO->new(-file=>$seqFiles[$i]);
    while (my $seq = $in->next_seq) {
	my $newName = $prefix[$i] . "-" . $seq->id;
	$seq->id($newName);
	$all->write_seq($seq);
    }
}
$all->close();



# run USearch
# -maxgaps parameter treats sequences unique except for an indel as distinct
# Since we have order the sequences by timepoint, we are guaranteed that the
#            centroid from each cluster will be from the earliest time point
#            at which that cluster is found, which saves us a re-sort step.

my $cmd = ppath('usearch') . " -cluster_smallmem work/phylo/all_seqs.fa -sortedby other -id $threshold -uc work/phylo/uc -maxgaps 2";
print "$cmd\n";
system($cmd);



#Parse USearch output (sample lines below)
#S       0       348     *       .       *       *       *       00-000154        *
#H       0       348     100.0   .       0       348     =       00-000180        00-000154
#C       0       18258   *       *       *       *       *       00-000154        *

my %cluster;
open UC, "work/phylo/uc" or die "Can't find output from USearch: $!. Please check parameters.\n";
while (<UC>) {
    last if /^C/; #speed processing by skipping summary lines

    my @a = split;
    my ($time) = $a[8] =~ /^(.*)\-/; #greedy operators mean we don't have to worry about dashes in user-supplied prefixes

    if ($a[0] eq "S") {

	#initialize new cluster
	$cluster{$a[8]}{'times'} = 1;
	$cluster{$a[8]}{'seen'}[$order{$time}] = 1;
	$cluster{$a[8]}{'first'} = $order{$time};
	$cluster{$a[8]}{'latest'} = $time;
	$cluster{$a[8]}{'persist'} = 1;
	$cluster{$a[8]}{'count'} = 1;
	next; #allows me to increment the total count for H lines outside the else below

    } elsif (! -defined $cluster{$a[9]}{'seen'}[$order{$time}] ) {

	#We haven't seen this cluster at this time point yet
	$cluster{$a[9]}{'times'}++;
	$cluster{$a[9]}{'seen'}[$order{$time}] = 1;
	if ( $order{$time} > $order{$cluster{$a[9]}{'latest'}} ) {
	    $cluster{$a[9]}{'latest'} = $time;
	    $cluster{$a[9]}{'persist'} = $order{$time} - $cluster{$a[9]}{'first'} + 1;
	}
    }

    #keep a total count
    $cluster{$a[9]}{'count'}++;

}
close UC;




my $in    = Bio::SeqIO->new(-file=>"work/phylo/all_seqs.fa");
my $out   = Bio::SeqIO->new(-file=>">output/sequences/nucleotide/$prj_name-collected.fa", -format=>'fasta');
my $outAA = Bio::SeqIO->new(-file=>">output/sequences/amino_acid/$prj_name-collected.fa", -format=>'fasta');
open TABLE, ">output/tables/$prj_name.txt" or die "Can't write to output/tables/$prj_name.txt: $!\n";
print TABLE "rep\tcount\ttimes\tpersist\tlatest\n";

while (my $seq = $in->next_seq) {

    next unless -defined( $cluster{$seq->id} );

    my $desc = $seq->desc . " long_obs=$cluster{$seq->id}{'count'} long_timepoints=$cluster{$seq->id}{'times'} persist=$cluster{$seq->id}{'persist'} last_timepoint=$cluster{$seq->id}{'latest'}";
    $seq->desc($desc);
    $out->write_seq($seq);
    $outAA->write_seq($seq->translate);
    print TABLE $seq->id."\t$cluster{$seq->id}{'count'}\t$cluster{$seq->id}{'times'}\t$cluster{$seq->id}{'persist'}\t$cluster{$seq->id}{'latest'}\n";

}

close TABLE;

