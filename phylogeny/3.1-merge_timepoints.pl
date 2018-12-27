#!/usr/bin/env perl

=head 1 SUMMARY

 3.1-merge_timepoints.pl

 This script uses vsearch to recluster selected sequences from multiple timepoints
       to determine "birthdays" and persistence times. It is recommended to use
       custom prefixes that correspond to actual time points from longitudinal
       samples, whether relative (weeks or months post infection) or absolute
       (sample collection date).

 Usage: 3.1-merge_timepoints.pl --seqs time1.fa --seqs time2.fa ...
                                  [ --labels t1 --labels t2 ... 
                                    -c --t 100 -m 1 -f]

 Invoke with -h or --help to print this documentation.

 Parameters:
     --seqs   => list of files with selected sequences from each timepoint.
                     Must be in proper chronological order.
     --labels => Optional list of timepoint IDs to be used in relabeling
                     input sequences. Must be in proper chronological order.
                     Defaults to "01", "02", etc and will be prepended with a
                     dash separating it from the orginal sequence name. Dashes
                     are not allowed in the labels themselves, to ease processing.
      -c      => Use vsearch's cluster_fast command instead of derep_fulllength.
                     Allows fragments missing a few AA at either end to be
                     counted as the same as a full-length sequence.
     --t      => Clustering threshold to use. Not setable when using the 
                     derep_fulllength algorithm. Defaults to 100 (%) even when
                     using cluster_fast.
     --m      => Optional minimum number of time points a sequence must be
                     observed in to be saved. Defaults to 1.

 Created by Chaim A. Schramm 2015-06-12.
 Fixed algorithm and added some options 2015-07-22.
 Changed label handling 2017-04-25.
 Updated for AIRR-format compatibility

 Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
                          Institutes of Health, USA. All rights reserved.

=cut

use warnings;
use strict;
use diagnostics;
use Pod::Usage;
use Getopt::Long;
use List::Util qw/max/;
use Cwd;
use File::Basename;
use Bio::SeqIO;
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);
use version qw/logCmdLine/;


if ($#ARGV < 0) { pod2usage(1); }

my @saveArgs = @ARGV;
my @seqFiles   = ( );
my @prefix     = ( );
my ($t, $min, $clustFast, $force, $help) = ( 100, 1, 0, 0, 0 );
GetOptions("seqs=s{1,}" => \@seqFiles,
	   "labels=s"   => \@prefix,
	   "c!"         => \$clustFast,
	   "t=f"        => \$t,
	   "m=i"        => \$min,
           "force!"     => \$force,
           "help!"      => \$help
    );

if ($help) { pod2usage(1); }
if ( $t < 100 && ! $clustFast ) { die("Must specify use of the cluster_fast command to use a custom id threshold.\n"); }
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
    #check for reg-ex compatibility
    for my $label (@prefix) {
	if ($label =~ /-/) { die("Please do not use dash characters in time point labels"); }
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
	     "output/sequences/amino_acid", "output/tables", "output/logs",
	     "output/rates", "output/plots", "work", "work/phylo") {
    if ( ! -d $dir ) { mkdir $dir; }
}

if ( -e "output/sequences/nucleotide/$prj_name-collected.fa" && ! $force ) {
    die "Output already exists. Please use the -f(orce) option to re-intiate and overwrite.\n";
}


#put this way down here so it's after the directory structure has been created
&logCmdLine($0,@saveArgs);


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



# run vsearch
my $cmd =  ppath() . "vsearch -derep_fulllength work/phylo/all_seqs.fa -uc work/phylo/uc";
if ($clustFast){
    # -maxgaps parameter treats sequences unique except for an indel as distinct
    # vsearch automatically sorts by length, accountign for possible fragments
    $cmd = ppath() . "vsearch -cluster_fast work/phylo/all_seqs.fa -id $threshold -uc work/phylo/uc -maxgaps 2";
}
print "$cmd\n";
system($cmd);



#Parse vsearch output (sample lines below)
#S       0       348     *       .       *       *       *       00-000154        *
#H       0       348     100.0   .       0       348     =       00-000180        00-000154
#C       0       18258   *       *       *       *       *       00-000154        *

my %cluster;
open UC, "work/phylo/uc" or die "Can't find output from vsearch: $!. Please check parameters.\n";
while (<UC>) {
    last if /^C/; #speed processing by skipping summary lines

    chomp;
    my @a = split /\t/; #usearch9 now includes the full fasta def line, which may include spaces 
    $a[8] =~ s/\S*\K.*//; #if there are spaces in the fasta def line, this dumps them, for compatibility with Bio::SeqIO below
    my ($time) = $a[8] =~ /^(.+?)\-/; #non-greedy operator means we don't have to worry about dashes in original sequence id

    if ($a[0] eq "S") {

	#initialize new cluster
	$cluster{$a[8]}{'rep'} = $a[8];	
	$cluster{$a[8]}{'seen'}{$time} = 1;
	$cluster{$a[8]}{'first'} = $order{$time};
	$cluster{$a[8]}{'latest'} = $time;
	$cluster{$a[8]}{'persist'} = 1;
	$cluster{$a[8]}{'count'} = 1;

    } else {

	$a[9] =~ s/\S*\K.*//; #if there are spaces in the fasta def line, this dumps them, for compatibility with Bio::SeqIO below
	$cluster{$a[9]}{'seen'}{$time}++;
	if ( $order{$time} < $cluster{$a[9]}{'first'} ) {
	    $cluster{$a[9]}{'first'} = $order{$time};
	    $cluster{$a[9]}{'persist'} = $order{$cluster{$a[9]}{'latest'}} - $order{$time} + 1;
	    $cluster{$a[8]}{'rep'} = $a[8];
	} elsif ( $order{$time} > $order{$cluster{$a[9]}{'latest'}} ) {
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
$out->width(600);
my $outAA = Bio::SeqIO->new(-file=>">output/sequences/amino_acid/$prj_name-collected.fa", -format=>'fasta');
$outAA->width(200);
my @cdr3;
open TABLE, ">output/tables/$prj_name.txt" or die "Can't write to output/tables/$prj_name.txt: $!\n";
print TABLE "rep\tcount\ttimes\tpersist\tlatest\n";

while (my $seq = $in->next_seq) {

    next unless -defined( $cluster{$seq->id} );

    my $numTimes = scalar( keys( %{$cluster{$seq->id}{'seen'}} ) );
    next unless $numTimes >= $min;

    #this sets the "birthday" time correctly
    #for derep, the sequence is by definition the same, anyway
    #for cluster fast, especially at lower thresholds, it's probably
    #  different, but we want to keep the centroid sequence so we'll
    #  just rename it
    $seq->id($cluster{$seq->id}{'rep'}); 

    my $desc = $seq->desc . " total_observations=$cluster{$seq->id}{'count'} num_timepoints=$numTimes persist=$cluster{$seq->id}{'persist'} last_timepoint=$cluster{$seq->id}{'latest'}";
    $seq->desc($desc);
    $out->write_seq($seq);
    $outAA->write_seq($seq->translate);
    print TABLE $seq->id."\t$cluster{$seq->id}{'count'}\t$numTimes\t$cluster{$seq->id}{'persist'}\t$cluster{$seq->id}{'latest'}\n";

    if ( $desc =~ /junction_aa=([A-Z]+)/ ) { 
	push @cdr3, Bio::Seq->new( -display_id => $cluster{$seq->id}{'rep'}, -seq=> $1 );
    }

}

close TABLE;

if ($#cdr3 >= 0) {
    my $outCDR3 = Bio::SeqIO->new(-file=>">output/sequences/amino_acid/$prj_name-CDR3.fa", -format=>'fasta');
    for my $s (@cdr3) { $outCDR3->write_seq($s); }
}
