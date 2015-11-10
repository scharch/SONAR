#!/usr/bin/env perl

=head 1 SUMMARY

 4.1-setup_plots.pl

 This script parses SONAR output tables and generates small flat files for
       easy plotting with 4.2-generate_plots.R.

 Usage: 4.1-setup_plots.pl

 Invoke with -h or --help to print this documentation.

 Parameters:
     --plot    => Statistic to be plotted. Options [with description]:
                      * ntraw     [raw read length in nucleotides]
                      * aaraw     [raw read length in amino acids]
                      * nttrim    [trimmed read length in nucleotides]
                      * aatrim    [trimmed read length in amino acids]
                      * ntcdr3    [cdr3 length in nucleotides]
                      * aacdr3    [cdr3 length in amino acids]
                      * v         [V gene usage]
                      * j         [J gene usage]
                      * d         [D gene usage, ignores unassigned reads]
                      * c         [subtype usage, ignores unassigned reads]
                      * div       [germline divergence (nt %) from muscle]
                      * blastdiv  [estimated germline divergence from blast]
                      * charge    [cdr3 charge, not including possible TYS]
                      * features  [presence of indels and stop codons]
                      * status    [V/J assignment success, in-frame, ORF, unique]
                      * size      [number of NGS reads per representative]

     --subset  => Which sequence set to plot for. Not relevant for the "status"
                      or "size" plots, which can only be done for 'all' and 
                      'unique', repsectively. Otherwise can be specified multiple
                      time to see how distributions shift.
                  Options are:
                      * all
                      * orf (default)    
                      * unique
                      * manual (use --list to specify)

     --compare => Add another dataset to the same plot. Can specify another plot
                     type, eg nttrim to go with ntraw or div and blastdiv
                     (WARNING: no sanity checks, it'll put nttrim and div on the
                     same plot if you tell it to!); another subtype; or another 
                     sequencing run (by path to project folder). Can be specified
                     multiple times, but user discretion is key.

     --output  => File name for output plot. Filetype (png, pdf, etc) specified
                     automatically by extension. Will be saved in the output/plots
                     directory. Default file name is plot_of_PLOT_for_SUBSET.png.

     --title   => Display title for plot. Default is just the current project name.

     --list    => If --subset manual is specificed, this opition should point to a
                     text file with one sequence ID per line and statistics will be
                     calculated only for the listed sequences.

 Created by Chaim A. Schramm 2015-09-01.

 Copyright (c) 2011-2015 Columbia University and Vaccine Research Center, National
                          Institutes of Health, USA. All rights reserved.

=cut

use strict;
use warnings;
use diagnostics;
use Pod::Usage;
use Getopt::Long;
use Switch;
use File::Basename;
use POSIX;

use List::Util qw/max min/;
use PPvars qw(ppath);
use Bio::SeqIO;


if ($#ARGV < 0) { pod2usage(1); }

my @compare  = ( );
my ($type, $subset, $outFile, $title, $listFile, $help) = ( "ntraw", "orf", "", "", "", 0 );
GetOptions("plot=s"    => \$type,
	   "subset=s"  => \$subset,
	   "compare=s" => \@compare,
	   "output=s"  => \$outFile,
	   "title=s"   => \$title,
	   "list=s"    => \$listFile,
           "help!"     => \$help
    );

if ($help) { pod2usage(1); }


#set up some defaults
$outFile = "output/plots/plot_of_$type\_for_$subset.png" if $outFile eq "";
$title   =  basename( getcwd ) if $title eq "";

#array to store output
my @lines = ();

#basic plot:
&chooseSubset( $subset, $listFile, $type, \@lines, getcwd );

#add overlays
for my $over (@compare) {
    if (-d $over) { &chooseSubset( $subset, $listFile, $type, \@lines, $over );
    } elsif ($over =~ /(all|orf|unique|manual)/i) { &chooseSubset( $over, $listFile, $type, \@lines, getcwd );
    } else { &chooseSubset( $subset, $listFile, $over, \@lines, getcwd ); }
}    


open OUT, ">work/internal/dataForPlotting.txt" or die "Can't write to work/internal/dataForPlotting.txt: $!\n\n";
print OUT join("\n",@lines)."\n";
close OUT;
#call 4.2 with appropriate options
system( "4.2-plot_histograms.R work/internal/dataForPlotting.txt $outFile $title");


sub chooseSubset {

    my ($subset, $listFile, $stat, $linesRef, $project) = @_;

    my $keepVal = 1; #default, corresponds to ORF
    my %lookUp;
    switch ($subset) {
	case (/^a/i) { $keepVal = 0;  }
	case (/^u/i) { $keepVal = 2;  }
	case (/^m/i) {
	    $keepVal = 999;
	    if (-e($listFile)) {
		open LIST, $listFile or die "Can't read from $listFile: $!\n\n";
		while (<LIST>) {
		    chomp;
		    $lookUp{$_} = 1;
		}
		close LIST;
	    } else {
		die "Please specify a valid file with sequence IDs when using --subset manual\n\n";
	    }
	}
    }

    &chooseStatistic($stat, $linesRef, $project, $keepVal, \%lookUp )

}


sub chooseStatistic {

    my ($stat, $linesRef, $project, $keepVal, $listRef) = @_;
 
    my $colNum = 3;
    my $type = "hist"; 
    my $divide = "";

    switch ($stat) {
	case /ntraw/    { ; } #default values work
	case /nttrim/   { $colNum = 4; } 
	case /ntcdr3/   { $colNum = 13; }
	case /aaraw/    { $divide = 3; }
	case /aatrim/   { $colNum = 4; $divide = 3; }
	case /aacdr3/   { $colNum = 14; }
	case /^v$/i     { $colNum = 5; $type = "count"; $divide = "\*"; }
	case /^j$/i     { $colNum = 7; $type = "count"; $divide = "\*"; }
	case /^d$/i     { $colNum = 6; $type = "count"; $divide = "\*"; }
	case /^c$/i     { $colNum = 8; $type = "count"; $divide = "\*"; }
	case /div/      { $colNum = 19; }
	case /blast/    { $colNum = 12; }
	case /charge/   { $colNum = 15; $divide = "charge"; }
	case /features/ { $colNum = 9; $type = "count"; }
	case /status/   { print "Using subset=all to plot read statuses...\n\n";
			  $keepVal = 0; $colNum = 11; $type = "count"; }
	case /size/     { print "Using subset=unique to plot read sizes...\n\n";
			  $keepVal = 2; $colNum = 18; }
	else {
	    print "Sorry, I don't recognize the plot option $stat. Acceptable values are:\n\t";
	    print join("\n\t", ("ntraw",'aaraw','nttrim','aatrim','ntcdr3','aacdr3','v','j','d','c','div','blastdiv','charge','features','status','size'));
	    die"\nPlease see program help for more information.\n\n";
	}
    }

    push @$linesRef, @{&plot($project, $keepVal, $listRef, $colNum, $type, $divide)}; 
    #I think the fact that I am passing by reference from the main program means I don't need to return the modified array

}


sub plot {

    my ($project, $subset, $listRef, $column, $type, $divide) = @_;

    #for results
    my @values;
    my %categories;

    my $prj_name = basename($project);
    open DATA, "$project/output/tables/$prj_name\_all_seq_stats.txt" or die "Can't read $project/output/tables/$prj_name\_all_seq_stats.txt: $!\n\n";
    while (<DATA>) {
	chomp;
	my @arr = split;
	
	#check if this sequence matches the subset
	my $value = 0; $value++ if $arr[12] =~ /T/; $value++ if $arr[17] =~ /T/;
	next unless ($value >= $subset || -defined($listRef->{$arr[0]}));
	
	
	#special cases
	#features
	if ($arr[9] =~ /T/) {
	    if ($arr[10] =~ /T/) {
		$arr[9] = "both";
	    } else {
		$arr[9] = "in-del only";
	    }
	} else {
	    if ($arr[10] =~ /T/) {
		$arr[9] = "stop codon only";
	    } else {
		$arr[9] = "orf";
	    }
	}
	#status
	if ($arr[17] =~ /T/) { $arr[11] = "unique"; }
	
	
	#get data
	if ($divide =~ /\d+/) { 
	    $arr[$column] /= $divide;
	} elsif ($divide == "charge") {
	    my $pos = () = $arr[$column] =~ /[KR]/gi;
	    my $neg = () = $arr[$column] =~ /[DE]/gi;
	    $arr[$column] = $pos - $neg;
	} elsif ($divide ne "") {
	    my @spl = split(/$divide/, $arr[$column]);
	    $arr[$column] = $spl[0];
	}
	
	#put it in the right place
	if ($type == "count") {
	    $categories{$arr[$column]}++;
	} else {
	    push @values, $arr[$column];
	}
	
    }
    close DATA;
    
    if ($type == "hist") {
	
	#get range
	my $min = floor(min(@values)); 
	my $max = ceil( max(@values));
	
	#choose bins
	if ($max > 700) {
	    #size plot, use log-scaled x
	    $min = floor(log10($min));
	    $max = ceil( log10($max));
	    for (my $v=$min; $v<$max+0.2; $v+=0.2) {
		$categories{sprintf("%d",10**$v)} = 0;
	    }
	} elsif ($max > 150) {
	    #nt raw, nt trim, or aa raw read lengths
	    $min = 25 * floor( $min/25 );
	    for (my $v=$min; $v<$max+25; $v+=25) {
		$categories{$v} = 0;
	    }
	} elsif ($max > 50 && $min > 30) {
	    #aa trim lengths
	    $min = 10 * floor( $min/10 );
	    for (my $v=$min; $v<$max+10; $v+=10) {
		$categories{$v} = 0;
	    }
	} else {
	    #nt/aa CDRH3, SHM, or charge
	    for (my $v=$min; $v<=$max; $v++) {
		$categories{$v} = 0;
	    }
	}
	
	#now sort into the bins
	my @bins = sort { $b <=> $a } keys %categories;
	for my $obs (@values) {
	    for my $b (@bins) {
		if ($b <= $obs) {
		    $categories{$b}++;
		    last;
		}
	    }
	}
	
	}
    
    my @lines;
    for my $cat (keys %categories) { #need to sort these
	#line format:
	#group x y
	#need to determine if group is just file or file+subset (eg multiple length comparisons)
	#for status and features, group and x are reversed compared to genes
    }
    
    return \@lines;

}
