#!/usr/bin/env perl

=head 1 SUMMARY

 mergeReplicates.pl

 This script uses USearch's dereplication command to find overlapping sequences
       between multiple technical replicates of a sample. Expected usage is with
       processed SOANAR pipeline output, but the algorithm will take any sequence
       files...

 Usage: mergeReplicates.pl --min 1 --seqs rep1.fa --seqs rep2.fa ...

 Invoke with -h or --help to print this documentation.

 Parameters:
     --min    => Minimum number of occurences in each replicate in order to keep
                    a particular read. Default=1.
     --seqs   => list of files with selected sequences from each timepoint.

 Created by Chaim A. Schramm 2015-06-12.
 Modified to current form by CA Schramm, 2015-07-21
 Added to SONAR utilities 2017-02-24
 Added min option 2017-04-23
 Modified to use VSearch by CAS 2018-07-30.

 Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
                          Institutes of Health, USA. All rights reserved.

=cut

use warnings;
use strict;
use diagnostics;
use Pod::Usage;
use Getopt::Long;
use Cwd;
use File::Basename;
use Bio::SeqIO;
use Statistics::Basic qw(correlation);
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);
use version qw/logCmdLine/;

my @saveArgs = @ARGV;


if ($#ARGV < 0) { pod2usage(1); }

my @seqFiles   = ( );
my @prefix     = ( );
my ($min, $force, $help) = ( 1, 0, 0 );
GetOptions("seqs=s{1,}" => \@seqFiles,
	   "min=i"      => \$min,
           "force!"     => \$force,
           "help!"      => \$help
    );

if ($help) { pod2usage(1); }


my $prj_name = basename(getcwd);
if ( -e "$prj_name-merged.fa" && ! $force ) {
    die "Output already exists. Please use the -f(orce) option to re-intiate and overwrite.\n";
}


&logCmdLine($0,@saveArgs);


#Read in selected seqs, rename, and write to temporary output
my $all = Bio::SeqIO->new(-file=>">all_seqs.fa", -format=>'fasta');
$all->width(600);
my %sizes;
my %correlated;
for my $f (@seqFiles) {

    my $base = basename( $f );
    my $prefix = (split(/[_]/, $base))[0];
    $correlated{$prefix} = [];
    
    my $cmd = ppath('usearch') . " -derep_fulllength $f -output derep-$base -sizein -sizeout -minuniquesize $min";
    print "\n$cmd\n\n";
    system( $cmd );

    #this is annoying; it's because vsearch doesn't preserve fasta description lines
    my %unique;
    my $usearch = Bio::SeqIO->new( -file=>"derep-$base" );
    while (my $seq = $usearch->next_seq) { 
	my $seqID = $seq->id;
	my $size  = 0;
	if ($seqID =~ /;size=(\d+)/) {
	    $size = $1;
	    $seqID =~ s/;.*//;
	}
	$unique{$seqID} = 1;
	$sizes{"$prefix-".$seqID} = $size;
    }

    my $in = Bio::SeqIO->new( -file=>$f );
    while (my $seq = $in->next_seq) {
	next unless -defined( $unique{$seq->id} );
	my $newName = "$prefix-" . $seq->id;
	$seq->id($newName);
	$all->write_seq($seq);
    }
}
$all->close();



# run USearch
my $cmd = ppath('usearch') . " -derep_fulllength all_seqs.fa -uc derep_all.uc";
print "\n$cmd\n\n";
system($cmd);



#Parse USearch output (sample lines below)
#S       0       348     *       .       *       *       *       00-000154        *
#H       0       348     100.0   .       0       348     =       00-000180        00-000154
#C       0       18258   *       *       *       *       *       00-000154        *

my %cluster;
my ($total, $both) = (0, 0);
open UC, "derep_all.uc" or die "Can't find output from vsearch: $!. Please check parameters.\n";
while (<UC>) {
    last if /^C/; #speed processing by skipping summary lines

    my @a = split /\t/;
    my @fastaDefLine = split /\s/, $a[8]; #legacy from usearch9, which included fasta definition line
    my ($time) = $fastaDefLine[0] =~ /^(.*)\-/; #greedy operators mean we don't have to worry about dashes in user-supplied prefixes

    if ($a[0] eq "S") {
	#initialize new cluster
	$cluster{$fastaDefLine[0]}{$time} = $fastaDefLine[0];
	$total++;
    } else {
	my @centroid = split /\s/, $a[9]; #legacy from usearch9, which included fasta definition line
	$cluster{$centroid[0]}{$time} = $fastaDefLine[0];
    }
}
close UC;



my $in    = Bio::SeqIO->new(-file=>"all_seqs.fa");
my $out   = Bio::SeqIO->new(-file=>">$prj_name-merged.fa", -format=>'fasta');
$out->width(600);
while (my $seq = $in->next_seq) {

    next unless -defined( $cluster{$seq->id} );

    for my $time (keys %correlated) {
	my $size = -defined($cluster{$seq->id}{$time}) ? $sizes{$cluster{$seq->id}{$time}} : 0;
	push @{$correlated{$time}}, $size;
    }
    
    next unless scalar(keys %{$cluster{$seq->id}}) >= 2; #don't really expect more than 2...

    $both++;
    $out->write_seq($seq);

}

print sprintf("\n\nFound $both/$total unique reads in both replicates (%d%%)\n", 100*$both/$total);
my @times = keys(%correlated);
my $r = correlation( $correlated{$times[0]}, $correlated{$times[1]} );
print sprintf("Pearson's rho of read abundance is %0.3f between the two replicates\n", $r);
