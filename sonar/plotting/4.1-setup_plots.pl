#!/usr/bin/env perl

=head 1 SUMMARY

 4.1-setup_plots.pl

 This script parses SONAR output tables and generates small flat files for
       easy plotting with 4.2-generate_plots.R.

 Usage: 4.1-setup_plots.pl --statistic ntraw --subset orf
                           [ --compare /path/to/other_project1 (...) 
                             --output output/plots/plot_of_ntraw_for_orf.png
                             --plotOpt "title=ProjectName" --plotOpt "h=3" (...)
                             --list readsToInclude.txt

 Invoke with -h or --help to print this documentation.

 Parameters:
     --statistic => Statistic to be plotted. Options [with description]:
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

     --subset    => Which sequence set to plot for. Default is 'orf' except for
                        "status" "features" or "size" plots, which default to 'all' 
                        'all' and 'unique', respectively, and can only be overridden
                        with the manual option. Otherwise can be specified multiple
                        multiple time to see how distributions shift.
                    Options are:
                        * all
                        * orf (default)    
                        * unique
                        * manual (use --list to specify)

     --compare   => Add another dataset to the same plot. Can specify another plot
                       type, eg nttrim to go with ntraw or div and blastdiv
                       (WARNING: no sanity checks, it'll put nttrim and div on the
                       same plot if you tell it to!); another subtype; or another 
                       sequencing run (by *ABSOLUTE* path to project folder). Can 
                       be specified multiple times, but user discretion is key.

     --output    => File name for output plot. Filetype (png, pdf, etc) specified
                       automatically by extension. Default file name is 
                       plot_of_PLOT_for_SUBSET.png in the output/plots directory. 

     --plotOpt   => Formatting options to pass to the actual plotting script in the 
                       form of "option=value". Can be specified multiple times.
                    Possible values (with defaults in parentheses) are:
                       * title (project name)
                       * bars [stack or dodge] (dodge except for features and status)
                       * percent (T)
                       * xlab (number OR category)
                       * ylab (percent OR number)
                       * xlim (chosen by system)
                       * ylim (chosen by system)
                       * logx (F except for size plot)
                       * logy (F except for size plot)
                       * mids [ie plot histogram bars at midpoint of bin]
                              ( T for nt raw/trim and aa raw/trim;
                                F for nt/aa CDR3, div/blastdiv, and charge )
                       * magnify (1)
                       * showlegend (T)
                       * legendpos ("right")
                       * h [figure height] (chosen by system)
                       * w [figure width] (chosen by system)
                       * dpi (150)

     --list    => If --subset manual is specificed, this opition should point to a
                     text file with one sequence ID per line and statistics will be
                     calculated only for the listed sequences.

 Created by Chaim A. Schramm 2015-09-01.

 Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National
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
use Bio::SeqIO;
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);
use version qw/logCmdLine/;


if ($#ARGV < 0) { pod2usage(1); }

my @saveArgs = @ARGV;
my @compare     = ( );
my @plotOptions = ( );
my ($type, $subset, $outFile, $listFile, $help) = ( "ntraw", "orf", "", "", 0 );
GetOptions("statistic=s" => \$type,
	   "subset=s"    => \$subset,
	   "compare=s"   => \@compare,
	   "output=s"    => \$outFile,
	   "plotOpt=s"   => \@plotOptions,
	   "list=s"      => \$listFile,
           "help!"       => \$help
    );

if ($help) { pod2usage(1); }

&logCmdLine($0,@saveArgs);


#set up some defaults
# put them at front of options array so they can be overriden if user specified alternates
$outFile = "output/plots/plot_of_$type\_for_$subset.png" if $outFile eq "";
unshift @plotOptions, "title=".basename( getcwd );
unshift @plotOptions, ("logx=T","logy=T") if $type =~ /size/;
unshift @plotOptions, "mids=T" if $type =~ /(raw|trim)/;
unshift @plotOptions, "bars=stack" if $type =~ /(features|status)/;


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
#add single quotes so spaces and parens are handled as expected
system( "4.2-plot_histograms.R work/internal/dataForPlotting.txt $outFile '" . join("' '",@plotOptions) . "'" );


sub chooseSubset {

    my ($subset, $listFile, $stat, $linesRef, $project) = @_;

    my $keepVal = 1; #default, corresponds to ORF
    my %lookUp;
    switch ($subset) {
	case (/^a/i) { $keepVal = 0;  }
	case (/^u/i) { $keepVal = 2;  }
	case (/^m/i) {
	    $keepVal = 3;
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

    &chooseStatistic($stat, $linesRef, $project, $keepVal, \%lookUp, $subset )

}


sub chooseStatistic {

    my ($stat, $linesRef, $project, $keepVal, $listRef, $subset) = @_;
 
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
	case /^v$/i     { $colNum = 5; $type = "count"; $divide = '\*'; }
	case /^j$/i     { $colNum = 7; $type = "count"; $divide = '\*'; }
	case /^d$/i     { $colNum = 6; $type = "count"; $divide = '\*'; }
	case /^c$/i     { $colNum = 8; $type = "count"; $divide = '\*'; }
	case /^div/     { $colNum = 18; }
	case /blast/    { $colNum = 12; }
	case /charge/   { $colNum = 15; $divide = "charge"; }
	case /features/ { if ($keepVal != 3) {
	                           print "Using subset=all to plot features...\n\n";
				   $keepVal = 0;
			  }
			  $colNum = 9; $type = "count"; }
	case /status/   { if ($keepVal != 3) {
	                           print "Using subset=all to plot read statuses...\n\n";
				   $keepVal = 0;
			  }
			  $colNum = 11; $type = "count"; }
	case /size/     { if ($keepVal != 3) {
	                           print "Using subset=unique to plot read sizes...\n\n";
				   $keepVal = 2; 
			  }
			  $colNum = 17; }
	else {
	    print "Sorry, I don't recognize the plot option $stat. Acceptable values are:\n\t";
	    print join("\n\t", ("ntraw",'aaraw','nttrim','aatrim','ntcdr3','aacdr3','v','j','d','c','div','blastdiv','charge','features','status','size'));
	    die"\nPlease see program help for more information.\n\n";
	}
    }

    push @$linesRef, @{&plot($project, $keepVal, $listRef, $colNum, $type, $divide, $subset, $stat)}; 
    #I am passing by reference from the main program so I don't need to return the modified array

}


sub plot {

    my ($project, $subset, $listRef, $column, $type, $divide, $subsetType, $statType) = @_;

    #for results
    my @values;
    my %categories;

    my $prj_name = basename($project);
    open DATA, "$project/output/tables/$prj_name\_all_seq_stats.txt" or die "Can't read $project/output/tables/$prj_name\_all_seq_stats.txt: $!\n\n";
    my $header = <DATA>;
    while (<DATA>) {
	chomp;
	my @arr = split;
	
	#error checking
	if ($subset == 2 && $#arr < 17) {
	    die "Cannot parse unique sequences from $prj_name; please run 1.4-dereplicate_seqeuences.pl and try again\n\n";
	}
	if ($column==18 && $#arr<18) {
	    if ($#arr == 16) {
		#assume this means that 2.1 was run without running 1.4 first and just roll with it
		$column = 16;
	    } else {
		die "Cannot parse cluster sizes from $prj_name; please run 2.1-calulate_id-div.pl and try again\n\n";
	    }
	}

	#don't count things that shouldn't be counted
	next if $arr[$column] eq "NA";

	#check if this sequence matches the subset
	my $value = 0; $value++ if $arr[11] eq "good"; $value++ if $#arr>=16 && $arr[16] =~ /T/;
	next unless ($value >= $subset || -defined($listRef->{$arr[0]}));
	
	
	#special cases
	#features
	if ($column == 9) {
	    if ($arr[9] =~ /T/) {
		if ($arr[10] =~ /T/) {
		    $arr[9] = "both";
		} else {
		    $arr[9] = "in-del";
		}
	    } else {
		if ($arr[10] =~ /T/) {
		    $arr[9] = "stop";
		} else {
		    $arr[9] = "orf";
		}
	    }
	}
	#status
	if ($column == 11 && $#arr>=16 && $arr[16] =~ /T/) { $arr[11] = "unique"; }
	#divergence
	if ($column == 12 || $column == 16 || $column ==18) {
	    ( $arr[$column] ) = $arr[$column] =~ /^(.+?)%?$/;
	}
	
	
	#get data
	if ($divide =~ /\d+/) { 
	    $arr[$column] /= $divide;
	} elsif ($divide eq "charge") {
	    my $pos = () = $arr[$column] =~ /[KR]/gi;
	    my $neg = () = $arr[$column] =~ /[DE]/gi;
	    $arr[$column] = $pos - $neg;
	} elsif ($divide ne "") {
	    my @spl = split(/$divide/, $arr[$column]);
	    $arr[$column] = $spl[0];
	}
	
	#put it in the right place
	if ($type eq "count") {
	    $categories{$arr[$column]}++;
	} else {
	    push @values, $arr[$column];
	}
	
    }
    close DATA;
    
    if ($type eq "hist") { 

	#get range
	my $min = floor(min(@values)); 
	my $max = ceil( max(@values));
	
	#choose bins
	if ($max > 1000) {
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
    my $groupName = "$prj_name-$subsetType-$statType";
    for my $cat ( sort mysort keys %categories ) {
	#line format: group x y
	#for status and features, group and x are reversed compared to genes
	#if ($column == 9 || $column == 11) {
	#    push @lines, "$cat\t$groupName\t$categories{$cat}";
	#} else {
	push @lines, "$groupName\t$cat\t$categories{$cat}";
	#}
    }
    
    return \@lines;

}


sub mysort {
    no warnings qw(numeric uninitialized);
    my %specificOrder = ( 'orf'=>0, 'wrong_length'=>1, 'noV'=>2, 'noJ'=>3, 'noCDR3'=>4, 'indel'=>5, 'stop'=>6, 'good'=>7, 'unique'=>8, 'both'=>9 );
    $specificOrder{$a}<=>$specificOrder{$b} || $a<=>$b || $a cmp $b;
}
