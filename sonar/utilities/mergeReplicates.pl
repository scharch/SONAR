#!/usr/bin/env perl

=head 1 SUMMARY

 mergeReplicates.pl

 This script uses USearch's dereplication command to find overlapping sequences
       between multiple technical replicates of a sample. Expected usage is with
       processed SOANAR pipeline output, but the algorithm will take any sequence
       files...

 Usage: mergeReplicates.pl --seqs rep1.fa --seqs rep2.fa ...

 Invoke with -h or --help to print this documentation.

 Parameters:
     --seqs   => list of files with selected sequences from each timepoint.

 Created by Chaim A. Schramm 2015-06-12.
 Modified to current form by CA Schramm, 2015-07-21
 Added to SONAR utilities 2017-02-24

 Copyright (c) 2011-2017 Columbia University and Vaccine Research Center, National
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
my ($force, $help) = ( 0, 0 );
GetOptions("seqs=s{1,}" => \@seqFiles,
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
    
    my $cmd = "usearch -derep_fulllength $f -fastaout derep-$base -sizein -sizeout";
    print "\n$cmd\n\n";
    system( $cmd );

    #this is annoying; it's because usearch doesn't preserve fasta description lines
    my %unique;
    my $usearch = Bio::SeqIO->new( -file=>"derep-$base" );
    while (my $seq = $usearch->next_seq) { 
	$unique{$seq->id} = 1;
	my ($size) = $seq->description =~ /;size=(\d+);/; #for old usearch it would have to be $seq->id
	$sizes{"$prefix-".$seq->id} = $size;
    }

    my $in = Bio::SeqIO->new( -file=>$f ); #with usearch9, this is an unnecessary step
    while (my $seq = $in->next_seq) {
	next unless -defined( $unique{$seq->id} );
	my $newName = "$prefix-" . $seq->id;
	$seq->id($newName);
	$all->write_seq($seq);
    }
}
$all->close();



# run USearch
my $cmd = "usearch -derep_fulllength all_seqs.fa -uc derep_all.uc";
print "\n$cmd\n\n";
system($cmd);



#Parse USearch output (sample lines below)
#S       0       348     *       .       *       *       *       00-000154        *
#H       0       348     100.0   .       0       348     =       00-000180        00-000154
#C       0       18258   *       *       *       *       *       00-000154        *

my %cluster;
my ($total, $both) = (0, 0);
open UC, "derep_all.uc" or die "Can't find output from USearch: $!. Please check parameters.\n";
while (<UC>) {
    last if /^C/; #speed processing by skipping summary lines

    my @a = split /\t/;
    my @fastaDefLine = split /\s/, $a[8];
    my ($time) = $fastaDefLine[0] =~ /^(.*)\-/; #greedy operators mean we don't have to worry about dashes in user-supplied prefixes

    if ($a[0] eq "S") {
	#initialize new cluster
	$cluster{$fastaDefLine[0]}{$time} = $fastaDefLine[0];
	$total++;
    } else {
	my @centroid = split /\s/, $a[9];
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
