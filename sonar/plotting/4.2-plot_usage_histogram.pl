#!/usr/bin/env perl -w

=head 1

 Usage: plotUsageHist.pl <genes.list> <germ_stat_01> [<germ_stat_02> <germ_stat_03> ...]

 plots a (multi-)bargraph for germline gene usage
 collapses alleles into genes and recalculates percentage

 also spits out the reads table for good luck
=cut

use strict;
use File::Basename;
use GD::Graph::bars;
use Pod::Usage;
use POSIX;

if ($#ARGV<0) { pod2usage(1); }

my @files=@ARGV;
#my @mycolors = $files[1]=~/traditional/ ? qw(lred blue green) : qw(blue lred green);
     #got rid of resorting, now everything stays in input order, so don't need case checking above!
my @mycolors=qw(purple green cyan blue lred);# green);

my %genes;
my @orderedFams;
my %profiles;
my %totals;

my %reads;
my $ymax=0;

#first get the gene names
my $list = shift(@files);
open LIST, $list;
while (<LIST>) {
    chomp;
    #(my $family) = $_ =~ /V(\d+D?)[\-\/]/;
    #if (! -defined $genes{$family}) {push @orderedFams, $family;}
    #$genes{$family}{$_} = 0;

    push @orderedFams, $_; #just use order in list file (including spaces) rather than trying for auto sort
    $genes{$_} = 0;        # which means we can forget about needing to break out the families
}

#next, read in all the data
my @orderedSamples;
for my $file (@files) {

    if (! -e $file) {
	print "Warning: Cannot find file $file, skipping...\n";
	next;
    }

    my %unused;
    my ($name,$path,$suffix) = fileparse($file,".txt");
    my @fields = split (/_.?germ_stat/,$name);
    #my $title = join ", ", @fields[0,1];
    my $title =$fields[0];
    push @orderedSamples, $title;

    open STAT, $file or die "Can't read $file: $!\n";
    while (<STAT>) {
	next if $_ =~ /^subject/;

	my @arr = split;
	#if ($arr[0] =~ /(.*V(\d+D?)[\-\/].*)\*/) {
	$arr[0] =~ s/\*.*//;
	$totals{$title} += $arr[1];
	if (-defined $genes{$arr[0]}) {
	    #$genes{$arr[0]}++;
	    $profiles{$title}{$arr[0]} += $arr[1];
	} else { $unused{$arr[0]}++; }
	#} else { print "Warning: Skipped unrecognized gene $arr[0] in $file...\n"; }

    }
    close STAT;
    if (scalar(keys %unused) > 0) { 
	my $list = join(", ", sort keys %unused);
	#warn("Found assigned genes in $file that were not in the list used:\n\t$list\n\n");
    }
    #my $found = join(", ", keys %{$profiles{$title}});
    #print "Found genes: $found\n\n";

    #for my $g (sort keys %{$profiles{$title}}) { print "$g\t$profiles{$title}{$g}\n"; }

    #get read data
    $reads{$title}{'mapped'} = $totals{$title};
    my $dir = $path; $dir =~ s/analysis\/data/logs/;
    open DATA, "$dir/1-split.log" or print "Can't find read data for $title from $dir/1-split.log: $!\n";
    my $line = <DATA>;
    if ($line =~ /total: (\d+); good: (\d+);/) {
	$reads{$title}{'raw'} = $1;
	$reads{$title}{'length'} = $2;
    }
    close DATA;

}


#now format it for plotting
#my @orderedFams = sort {$a <=> $b} keys %genes;
my @data = ([@orderedFams]);
#for my $fam (@orderedFams) {
#    push @{$data[0]}, (sort keys %{$genes{$fam}});
#    push @{$data[0]}, ""; #separate families with a blank column
#}
print "Sample\tRaw\tAppropriate Length\tMapped\n";
for my $sample (@orderedSamples) { #sort keys %profiles) {
    push @data, [];
    print "$sample\t$reads{$sample}{'raw'}\t$reads{$sample}{'length'}\t$reads{$sample}{'mapped'}\n";
    for my $gene (@orderedFams) {
	#for my $gene (sort keys %{$genes{$fam}}) {
	    if (! defined $profiles{$sample}{$gene}) { $profiles{$sample}{$gene} = 0; }
	    my $percent = $profiles{$sample}{$gene} * 100 /$totals{$sample};
	    push @{$data[$#data]}, $percent;
	    #print "$sample, $gene: $percent\n";
	    if ($percent > $ymax) { $ymax = $percent; }	    
	#}
	#push @{$data[$#data]}, 0; #this is the empty between families
    }
}


my $graph = GD::Graph::bars->new(1536, 384);
#my $graph = GD::Graph::bars->new(512, 384);
$graph->set(
    transparent       => 0,
    bgclr             => 'white',
    x_ticks           => 0,
    x_label           => 'Germline Gene',
    x_label_position  => 0.5,
    x_labels_vertical => 1,
    y_label           => 'Percent of all reads',
    y_number_format   => "%d%%",
#    bar_width         => 6,
#    bar_spacing       => 2,
    bargroup_spacing  => 10,
    y_long_ticks      => 1,
    y_max_value       => 5 * ceil($ymax/5),
    y_tick_number     => ceil($ymax/5),
#    legend_placement  => 'RT',
    dclrs             => \@mycolors,
    title             => 'Germline Gene Usage'
    ) or warn $graph->error;

$graph->set_title_font('/usr/share/fonts/truetype/msttcorefonts/arial.ttf', 60);
$graph->set_x_axis_font('/usr/share/fonts/truetype/msttcorefonts/arial.ttf', 38);
$graph->set_y_axis_font('/usr/share/fonts/truetype/msttcorefonts/arial.ttf', 38);
$graph->set_x_label_font('/usr/share/fonts/truetype/msttcorefonts/arial.ttf', 38);
$graph->set_y_label_font('/usr/share/fonts/truetype/msttcorefonts/arial.ttf', 48);
$graph->set_legend_font('/usr/share/fonts/truetype/msttcorefonts/arial.ttf', 48);
$graph->set_legend(sort keys %profiles);

my $image = $graph->plot(\@data) or die $graph->error;

open(IMG, '>germline_usage.png') or die $!;
print IMG $image->png;
close IMG;
