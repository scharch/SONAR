#!/usr/bin/perl -w

=head 1 SUMMARY

 parseCDHit.pl
 Usage: parseCDHit.pl --input <cdhit.fa.clstr> --output <representatives.txt>
                      [--minSize <1> --imgt 1_Summary_imgt.txt --sortOn <\d\d\d> --sortOrder '("DKJ"=>1, "ABC"=>2 ...)']

 This program parses the results of clustering done by CD-Hit by picking a sequence from the earliest possible timepoint as the representative,
    subject to a minimum cluster size.

 Parameters:
     --input       => list of clusters output by CD-Hit
     --output      => where to save the list of representatives
     --minSize     => ignore clusters smaller than some size. Defaults to 1 (keep all clusters)
     --imgt        => when picking representatives, look this IMGT output file (ids must match!) to try and skip bad reads. Will not change the
                         earliest timepoint, but may associate a later sequence with an earlier ID.
     --sortOn      => regex for a leader sequence that will be sorted to check the time point.
     --sortOrder   => string representation of a hash that describes the correct sort order (keys are the leader sequences and values are numeric
                         order. Defaults to Perl's native lexical sort.

 Created by Chaim A. Schramm 2013-02-28.
 Modified 2013-05-06 to check IMGT output for low-quality reads.
 Modified 2013-08-15 to allow a user-supplied sort function
 Modified 2013-09-15 to account for possible left-zero difference between IMGT and CDHit ids. Also now ignores an HX/LX clade designation for donor 45 data.

 Copyright (c) 2013 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

=cut

use strict;
use diagnostics;
use Pod::Usage;
use Getopt::Long;
use Bio::TreeIO;
use Bio::Tree::Node;

if ($#ARGV < 0) { pod2usage(1); }

my ($infile, $outfile, $minSize, $imgt, $sortOn, $sortOrder, $help) = ("", "", 1, "", '\d\d\d', 0);
GetOptions("input=s"     => \$infile,
           "output=s"    => \$outfile,
           "minSize=i"   => \$minSize,
	   "imgt=s"      => \$imgt,
           'sortOn=s'    => \$sortOn,
	   'sortOrder=s' => \$sortOrder,
           "help!"       => \$help
    );

if ($help || ! -e($infile) || $minSize < 1) { pod2usage(1); }


#load imgt, if flagged
my %quality;
if (-e $imgt) {
    open IMGT, $imgt or die "Can't read $imgt: $!\n";
    while (<IMGT>) {
	next if /(Sequence number|No results)/;
	my @imgt_line = split(/\t/, $_);
	my @name = split /\s+/, $imgt_line[1];
	my ($date,$id);
	if ($name[0]=~/([HL]\d-)?($sortOn)-?(\d+)/) { $date = $2; $id = $3; }
	my $gap_count = 0;
	if ($imgt_line[2] =~ /deletion: (\d+)/) { $gap_count += $1; }
	if ($imgt_line[8] =~ /^(\.+)/) { $gap_count += length($1); } # 5'V is missing
	$gap_count += 34-length($imgt_line[17]); #assumption at proper length of J; may be clone specific?
	
	$quality{$date}[$id] = $gap_count;
    }
    close IMGT;
}

open OUT, ">$outfile" or die "Can't write to $outfile: $!\n";

my %seqs = ();
my ($rep, $repTime, $total) = ("","",0);
my %cladeNames; #allows us to sort numerically on IDs and not worry about different numbers of 0s on the left side, while still maintaining the convention used in the CDHit output for final output and downstream matching.
open IN, "$infile" or die "Can't read from $infile: $!\n";
print "#Centroid\tSize\n";
while (<IN>) {
    
    if ($_ =~ />Cluster/) {
	if ($total >= $minSize) {
	    my @times;
	    if ($sortOrder eq "") {
		@times = sort keys %seqs;
	    } else {
		my %order = eval($sortOrder);
		@times = sort { $order{$a} <=> $order{$b} } keys %seqs;
	    }
	    my @best; my $printed = 0;
	    if ($repTime ne $times[0]) {
		print "earliest is $times[0], but $repTime was assigned...\n";
		@best = sort { $seqs{$times[0]}{$b} <=> $seqs{$times[0]}{$a} } keys %{$seqs{$times[0]}}; #grab highest ID from earliest time point as new rep
		$rep = $best[0]; $repTime=$times[0];
	    }

	    if (scalar(keys %quality) > 0) {
		    if ($quality{$repTime}[$rep] > 0) {
			#first check other reads at same time
			@best = sort { $quality{$times[0]}[$a] <=> $quality{$times[0]}[$b] } keys(%{$seqs{$times[0]}});
			#print "Best at this point is $best[0]\n";
			if ($quality{$times[0]}[$best[0]] < $quality{$repTime}[$rep]) { 
			    print "Change rep from $repTime-$rep to $times[0]-$best[0]\n"; 
			    $rep = $best[0]; 
			}

			if ($quality{$repTime}[$rep] > 0) {
			    #now gotta check other time points
			    my $bestQual=9999; my $bestTime=""; my $bestID=0;
			    for my $t (@times) {
				my @all_by_quality = sort { $quality{$t}[$a] <=> $quality{$t}[$b] } keys(%{$seqs{$t}});
				if ($quality{$t}[$all_by_quality[0]] < $bestQual) {
				    $bestQual = $quality{$t}[$all_by_quality[0]];
				    $bestTime = $t; $bestID = $all_by_quality[0];
				}
			    }
			    if ($bestQual < $quality{$repTime}[$rep]) {
				if (! -defined($cladeNames{$bestTime}[$bestID]) ) { print "Where is $bestTime-$bestID?\n"; }
				#print "$rep: $quality{$rep} VS. $all_by_quality[0]: $quality{$all_by_quality[0]}\n";
				print "Use sequence from $cladeNames{$bestTime}[$bestID] with ID $cladeNames{$repTime}[$rep] for sake of quality! ($total)\n";
				print OUT "$cladeNames{$repTime}[$rep] => $cladeNames{$bestTime}[$bestID]\n";
				$printed=1;
			    } 
			} 
		    }
	    }
	    
	    if (!$printed) { 
		print OUT "$cladeNames{$repTime}[$rep]\n"; print "$cladeNames{$repTime}[$rep]\t$total\n";
	    }
	
	%seqs=();
	($rep, $repTime, $total) = ("","",0);
	}
    } elsif ($_ =~ />(([HL]\d-)?($sortOn)-?(\d+))... (at [\d:]*\/?(\d+\.\d+)%|\*)/) {
	$total++;
	if ($5 eq "*") { $rep = $4; $repTime = $3; }
	$seqs{$3}{$4} = $6;
	$cladeNames{$3}[$4] = $1;
	if (! -defined($quality{$3}[$4]) ) { $quality{$3}[$4] = 1000; } #probably due to IMGT "No Results"
    }
	
}
close IN;

#last cluster!
if ($total >= $minSize) {
    my @times;
    if ($sortOrder eq "") {
	@times = sort keys %seqs;
    } else {
	my %order = eval($sortOrder);
	@times = sort { $order{$a} <=> $order{$b} } keys %seqs;
    }
    my @best; my $printed = 0;
    if ($repTime ne $times[0]) {
	print "earliest is $times[0], but $repTime was assigned...\n";
	@best = sort { $seqs{$times[0]}{$b} <=> $seqs{$times[0]}{$a} } keys %{$seqs{$times[0]}}; #grab highest ID from earliest time point as new rep
	$rep = $best[0]; $repTime=$times[0];
    }
    
    if (scalar(keys %quality) > 0) {
	if ($quality{$repTime}[$rep] > 0) {
	    #first check other reads at same time
	    @best = sort { $quality{$times[0]}[$a] <=> $quality{$times[0]}[$b] } keys(%{$seqs{$times[0]}});
	    #print "Best at this point is $best[0]\n";
	    if ($quality{$times[0]}[$best[0]] < $quality{$repTime}[$rep]) { 
		print "Change rep from $repTime-$rep to $times[0]-$best[0]\n"; 
		$rep = $best[0]; 
	    }
	    
	    if ($quality{$repTime}[$rep] > 0) {
		#now gotta check other time points
		my $bestQual=999; my $bestTime=""; my $bestID=0;
		for my $t (@times) {
		    my @all_by_quality = sort { $quality{$t}[$a] <=> $quality{$t}[$b] } keys(%{$seqs{$t}});
		    if ($quality{$t}[$all_by_quality[0]] < $bestQual) {
			$bestQual = $quality{$t}[$all_by_quality[0]];
			$bestTime = $t; $bestID = $all_by_quality[0];
		    }
		}
		if ($bestQual < $quality{$repTime}[$rep]) {
		    #print "$rep: $quality{$rep} VS. $all_by_quality[0]: $quality{$all_by_quality[0]}\n";
		    print "Use sequence from $cladeNames{$bestTime}[$bestID] with ID $cladeNames{$repTime}[$rep] for sake of quality! ($total)\n";
		    print OUT "$cladeNames{$repTime}[$rep] => $cladeNames{$bestTime}[$bestID]\n";
		    $printed=1;
		} 
	    } 
	}
    }
	    
    if (!$printed) { 
	print OUT "$cladeNames{$repTime}[$rep]\n"; print "$cladeNames{$repTime}[$rep]\t$total\n";
    }
    
}    

close OUT;
