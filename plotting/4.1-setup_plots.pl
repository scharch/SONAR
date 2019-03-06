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
                        * ntcdr3    [*junction* length in nucleotides (ie IMGT CDR3 +6)]
                        * aacdr3    [*junction* length in amino acids (ie IMGT CDR3 +2)]
                        * v         [V gene usage]
                        * j         [J gene usage]
                        * d         [D gene usage, ignores unassigned reads]
                        * c         [subtype usage, ignores unassigned reads]
                        * div       [germline divergence (nt %) from muscle, requires running 1.4 first]
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

     --plotName  => Specify custom names for each statistic/dataset being plotted.
                       Use multiple times, in the same order as the arguments to 
                       `compare`. Default is project-subset-statistic.

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
                       * showlegend (F)
                       * legendpos ("right")
                       * h [figure height] (chosen by system)
                       * w [figure width] (chosen by system)
                       * dpi (150)

     --list    => If --subset manual is specificed, this opition should point to a
                     text file with one sequence ID per line and statistics will be
                     calculated only for the listed sequences.

 Created by Chaim A. Schramm 2015-09-01.
 Modified by CAS 2017-06-09 to add plotName option and to fix compatability with
                          docopt usage by 4.2-plot_histograms.R
 Modified to plot from new AIRR-format rearrangements.tsv file by CAS 2018-10-19.

 Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
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
my @plotNames   = ( );
my ($type, $subset, $outFile, $listFile, $help) = ( "ntraw", "orf", "", "", 0 );
GetOptions("statistic=s" => \$type,
	   "subset=s"    => \$subset,
	   "compare=s"   => \@compare,
	   "output=s"    => \$outFile,
	   "plotOpt=s"   => \@plotOptions,
	   "plotNames=s" => \@plotNames,
	   "list=s"      => \$listFile,
           "help!"       => \$help
    );

if ($help) { pod2usage(1); }

&logCmdLine($0,@saveArgs);

#change format of options array to be compatible with DocOpt
# doing it this way instead of asking user to input --option value to avoid ambiguities parsing options for this script.
for my $opt (@plotOptions) {
    my @aa= split /=/,$opt;
    $opt = "--$aa[0] '$aa[1]'";
}
my $commandOpts = join(" ",@plotOptions);

#set up some defaults
#docopt throws an error if options are repeated unexpectedly, so check for user specification
$outFile = "output/plots/plot_of_$type\_for_$subset.png" if $outFile eq "";
$commandOpts .= " --title ".basename( getcwd ) unless $commandOpts =~/title/;
$commandOpts .= " --logx T" if $type =~ /size/ && $commandOpts !~ /logx/;
$commandOpts .= " --logy T" if $type =~ /size/ && $commandOpts !~ /logy/;
$commandOpts .= " --mids T" if $type =~ /(raw|trim)/ && $commandOpts !~ /mids/;
$commandOpts .= " --bars stack" if $type =~ /(features|status)/ && $commandOpts !~ /bars/;

#kluge
my $logX=0;
$logX=1 if $commandOpts =~ /logx=T/;


#array to store output
my @lines = ();

#basic plot:
&chooseSubset( $subset, $listFile, $type, \@lines, getcwd, shift(@plotNames) );

#add overlays
for my $over (@compare) {
    if (-d $over) { &chooseSubset( $subset, $listFile, $type, \@lines, $over, shift(@plotNames) );
    } elsif ($over =~ /(all|orf|unique|manual)/i) { &chooseSubset( $over, $listFile, $type, \@lines, getcwd );
    } else { &chooseSubset( $subset, $listFile, $over, \@lines, getcwd ); }
}    


open OUT, ">work/internal/dataForPlotting.txt" or die "Can't write to work/internal/dataForPlotting.txt: $!\n\n";
print OUT join("\n",@lines)."\n";
close OUT;

#call 4.2 with appropriate options
#add single quotes so spaces and parens are handled as expected
system( "4.2-plot_histograms.R work/internal/dataForPlotting.txt $outFile $commandOpts" );

sub chooseSubset {

    my ($subset, $listFile, $stat, $linesRef, $project, $customName) = @_;

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

    &getPlotData($stat, $linesRef, $project, $keepVal, \%lookUp, $subset, $customName )

}


sub getPlotData {

    my ($stat, $linesRef, $project, $keepVal, $listRef, $subset, $customName) = @_;
 
    my $prj_name = basename($project);
    open DATA, "$project/output/tables/$prj_name\_rearrangements.tsv" or die "Can't read $project/output/tables/$prj_name\_rearrangements.tsv.\nIf this project was started with an older version of SONAR, please use utilities/convertToAIRR.py\n   to convert the old all_seq_stats.txt file to the new format and then try again.\n\n";

    my $header = <DATA>;
    chomp $header;
    my @colnames = split(/\t/, $header);
    my %colFinder;
    while (my ($pos, $name) = each @colnames) {
	$colFinder{$name} = $pos;
    }

    my $colNum = -1;
    my $type = "hist"; 
    my $divide = "";

    { no warnings qw(uninitialized);
      
      switch ($stat) {
	  case /ntraw/    { $colNum = $colFinder{"length_raw"}; }
	  case /nttrim/   { $colNum = $colFinder{"length_trimmed"}; } 
	  case /ntcdr3/   { $colNum = $colFinder{"junction_length"}; }
	  case /aaraw/    { $colNum = $colFinder{"length_raw"};      $divide = 3; }
	  case /aatrim/   { $colNum = $colFinder{"length_trimmed"};  $divide = 3; }
	  case /aacdr3/   { $colNum = $colFinder{"junction_length"}; $divide = 3; }
	  case /^v$/i     { $colNum = $colFinder{"v_call"}; $type = "count"; $divide = '\*'; }
	  case /^j$/i     { $colNum = $colFinder{"j_call"}; $type = "count"; $divide = '\*'; }
	  case /^d$/i     { $colNum = $colFinder{"d_call"}; $type = "count"; $divide = '\*'; }
	  case /^c$/i     { $colNum = $colFinder{"c_call"}; $type = "count"; $divide = '\*'; }
	  case /^div/     { $colNum = $colFinder{"v_identity"}; }
	  case /blast/    { $colNum = $colFinder{"blast_identity"}; }
	  case /charge/   { $colNum = $colFinder{"junction_aa"}; $divide = "charge"; }
	  case /features/ { if ($keepVal != 3) {
	                           print "Using subset=all to plot features...\n\n";
				   $keepVal = 0;
			    }
			    $colNum = $colFinder{'indels'}; $type = "count"; }
	  case /status/   { if ($keepVal != 3) {
	                           print "Using subset=all to plot read statuses...\n\n";
				   $keepVal = 0;
			    }
			    $colNum = $colFinder{'status'}; $type = "count"; }
	  case /size/     { if ($keepVal != 3) {
	                           print "Using subset=unique to plot read sizes...\n\n";
				   $keepVal = 2; 
			    }
			    $colNum = $colFinder{'cluster_count'}; }
	  else {
	      print "Sorry, I don't recognize the plot option $stat. Acceptable values are:\n\t";
	      print join("\n\t", ("ntraw",'aaraw','nttrim','aatrim','ntcdr3','aacdr3','v','j','d','c','div','blastdiv','charge','features','status','size'));
	      die"\nPlease see program help for more information.\n\n";
	  }
      }

      if ( ! -defined($colNum) || $colNum<0 ) {
	  die("ERROR: could not find data for $stat in $project/output/tables/$prj_name\_rearrangements.tsv\n\n");
      }
    }


    #now parse actual data

    #for results
    my @values;
    my %categories;

    while (<DATA>) {
	chomp;
	my @arr = split(/\t/, $_);
	
	#don't count things that shouldn't be counted
	next if ! -defined($arr[$colNum]) || $arr[$colNum] eq "" || $arr[$colNum] eq "NA";

	#check if this sequence matches the subset
	if ($keepVal == 3) {
	    next unless -defined($listRef->{$arr[0]});
	} elsif ($keepVal == 2 ) {
	    next unless $arr[$colFinder{'status'}] eq "unique";
	} elsif ($keepVal == 1 ) {
	    next unless $arr[$colFinder{'status'}] =~ /(good|unique)/;
	}
	
	
	#special cases
	if ($stat eq "features") {
	    if ($arr[$colFinder{'indels'}] =~ /T/) {
		if ($arr[$colFinder{'stop_codon'}] =~ /T/) {
		    $arr[$colFinder{'indels'}] = "both";
		} else {
		    $arr[$colFinder{'indels'}] = "in-del";
		}
	    } else {
		if ($arr[$colFinder{'stop_codon'}] =~ /T/) {
		    $arr[$colFinder{'indels'}] = "stop";
		} else {
		    $arr[$colFinder{'indels'}] = "orf";
		}
	    }
	} elsif ($stat =~ /(blast|div)/) {
	    #convert fractional identity to percent divergence
	    $arr[$colNum] = 100 * ( 1 - $arr[$colNum] )
	}
	
	
	#get data
	if ($divide =~ /\d+/) { 
	    $arr[$colNum] /= $divide;
	} elsif ($divide eq "charge") {
	    my $pos = () = $arr[$colNum] =~ /[KR]/gi;
	    my $neg = () = $arr[$colNum] =~ /[DE]/gi;
	    $arr[$colNum] = $pos - $neg;
	} elsif ($divide ne "") {
	    my @spl = split(/$divide/, $arr[$colNum]);
	    $arr[$colNum] = $spl[0];
	}
	
	#put it in the right place
	if ($type eq "count") {
	    $categories{$arr[$colNum]}++;
	} else {
	    push @values, $arr[$colNum];
	}
	
    }
    close DATA;
    
    if ($type eq "hist") { 

	#get range
	my $min = floor(min(@values)); 
	my $max = ceil( max(@values));
	
	#choose bins
	if ($stat =~ /size/ || $logX == 1) {
	    #size plot, use log-scaled x
	    $min = floor(log10($min));
	    $max = ceil( log10($max));
	    for (my $v=$min; $v<$max+0.2; $v+=0.2) {
		$categories{sprintf("%d",10**$v)} = 0;
	    }
	} elsif ($stat =~ /(ntraw|nttrim|aaraw)/) {
	    #nt raw, nt trim, or aa raw read lengths
	    $min = 25 * floor( $min/25 );
	    for (my $v=$min; $v<$max+25; $v+=25) {
		$categories{$v} = 0;
	    }
	} elsif ($stat =~ /aatrim/) {
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
	    for my $bin (@bins) {
		if ($bin <= $obs) {
		    $categories{$bin}++;
		    last;
		}
	    }
	}

    }
    
    $customName = "$prj_name-$subset-$stat" unless -defined($customName);
    for my $cat ( sort mysort keys %categories ) {
	#line format: group x y
	#for status and features, group and x are reversed compared to genes
	#if ($column == 9 || $column == 11) {
	#    push \@linesRef, "$cat\t$customName\t$categories{$cat}";
	#} else {
	push @$linesRef, "$customName\t$cat\t$categories{$cat}";
	#}
    }

}


sub mysort {
    no warnings qw(numeric uninitialized);
    my %specificOrder = ( 'orf'=>0, 'wrong_length'=>1, 'noV'=>2, 'noJ'=>3, 'noCDR3'=>4, 'nonproductive'=>5, 'indel'=>6, 'stop'=>7, 'good'=>8, 'unique'=>9, 'both'=>10 );
    $specificOrder{$a}<=>$specificOrder{$b} || $a<=>$b || $a cmp $b;
}
