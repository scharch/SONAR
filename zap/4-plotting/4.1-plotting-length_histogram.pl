#!/usr/bin/perl -w

=head 1

 Usage: plotLengthHist.pl <original01.fa> [<original02.fa> <original03.fa> ... output.png -s]

 plots a (multi-)bargraph for 454 read length

 Parameters:
     fasta files with the 454 reads
     output.png - Optional name for output. Must be *.png to be recognized
     -s - Optional flag to get length from sequence (default is from header)
=cut

use strict;
use File::Basename;
use GD::Graph::bars;
use Pod::Usage;
use Switch;
use POSIX;
use Bio::SeqIO;
use List::Util qw/sum max/;

if ($#ARGV<0) { pod2usage(1); }

my @files; my $manual = 0; my $output="output.png";
for my $arg (@ARGV) {
    switch ($arg) {
	case /\.png$/ { $output = $arg; }
	case "-s"     { $manual = 1; }
	else          { push @files, $arg; }
    }
}

my @mycolors = qw(lred dred pink lblue dblue cyan lgreen dgreen dyellow);


my @data = (['1-50','51-100','101-150','151-200','201-250','251-300','301-350','351-400','401-450','451-500','501-550','551-600','601-650','651-700','701-750','751-800','801-850','851-900','901-950','951-1000']);
my @profiles; my $ymax = 0;

for my $file (@files) {

    if (! -e $file) {
	print "Warning: Cannot find file $file, skipping...\n";
	next;
    }

    my ($stem,$path,$suffix) = fileparse($file,qr/\.[^.]*/);
    push @profiles, $stem;

    my @bins = (0) x 20;

    if ($manual) {
	my $in = Bio::SeqIO->new(-file=>$file);
	while (my $seq = $in->next_seq) {
	    my $len = int($seq->length / 50);
	    next if $len >= 20;
	    $bins[$len]++;
	}
    } else {
	open FA, $file or die "Can't read $file: $!\n";
	while (<FA>) {
	    if ($_ =~ /length=(\d+)/) {
		my $len = int($1 / 50);
		next if $len >= 20;
		$bins[$len]++;
	    }
	}
	close FA;
    }
    
    my $totalReads = sum(@bins);
    #print join("+",@bins). "  = ". sum(@bins) . "\n";
    map { $_ *= (100/$totalReads)} @bins;
    #print join("+",@bins). "  = ". sum(@bins) . "\n";
    
    if (max(@bins) > $ymax) { $ymax = max(@bins); }

    push @data, \@bins;

}


my $graph = GD::Graph::bars->new(1536, 384);
$graph->set(
    transparent       => 0,
    bgclr             => 'white',
    x_ticks           => 0,
    x_label           => 'Length of V gene alignement',#454 Read Length',
    x_label_position  => 0.5,
    x_labels_vertical => 1,
    y_label           => 'Percent of Reads',
    y_number_format   => "%d%%",
#    bar_width         => 6,
#    bar_spacing       => 2,
    bargroup_spacing  => 10,
    y_long_ticks      => 1,
    y_max_value       => 5 * ceil($ymax/5),
    y_tick_number     => ceil($ymax/5),
#    legend_placement  => 'RT',
    dclrs             => \@mycolors,
    title             => 'V gene Length Distribution'#454 Read Length Distribution'
    ) or warn $graph->error;

$graph->set_title_font('/usr/share/fonts/truetype/msttcorefonts/arial.ttf', 60);
$graph->set_x_axis_font('/usr/share/fonts/truetype/msttcorefonts/arial.ttf', 38);
$graph->set_y_axis_font('/usr/share/fonts/truetype/msttcorefonts/arial.ttf', 38);
$graph->set_x_label_font('/usr/share/fonts/truetype/msttcorefonts/arial.ttf', 38);
$graph->set_y_label_font('/usr/share/fonts/truetype/msttcorefonts/arial.ttf', 48);
$graph->set_legend_font('/usr/share/fonts/truetype/msttcorefonts/arial.ttf', 48);
$graph->set_legend(@profiles);

my $image = $graph->plot(\@data) or die $graph->error;

open(IMG, ">$output") or die $!;
print IMG $image->png;
close IMG;
