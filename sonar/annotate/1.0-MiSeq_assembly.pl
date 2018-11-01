#!/usr/bin/perl
#    barcode file should be in the format: name\tsequence\n
use strict;
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);
use version qw/logCmdLine/;

#########checking parameters#######
if(@ARGV%2>0||!@ARGV){die "Usage: illumina_filtering.pl 
	-f forward_read_seq file, either compressed .gz or fastq file 
	-r reverse_strand, either compressed .gz or fastq file
        -q quality encoding, default 33 (Sanger, Illumina1.8+)
	-ut usearch_truncation_quality_threshold_before_mapping
	-maxee maximum number of expected errors in the paired reads, default:1.0 If the quality of your sequencing run is low, you may adjust to a higher cutoff
	-trimf number of nucleotides trimed off before assembly for forward mate, default:1
	-trimr number of nucleotides trimed off before assmebly for reverse mate, default:1
	-o output folder
	-maxdiff the maximum number of mismatches allowed in the overlapping region between forward and reverse reads. Default: 10
	-maxdiffpct the maximum percent of mismatches allowed in the overlapping region between forward and reverse reads. Default: 10
	-threads default:20
	-minl minimal length for the merged reads default:300
	-maxl maximal length for the merged reads default:600
	-phix (yes or no) [CURRENTLY DEPRECATED] whether to remove phix genome in the dataset, require usearchv9 and above. Default is no.
Example:
1.0-MiSeq_assembly.pl -usearch usearch -f forward_read.fastq.gz -r reverse_read.fastq.gz -o ./

Created by Zizhang Sheng 2015-12-02.
Modified to use vsearch and clean out dead code by Chaim A Schramm, 2018-07-30.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

	";}

&logCmdLine($0,@ARGV);

my %para=@ARGV;
if(!$para{'-f'}){die "no input f sequence file\n";}
if(!$para{'-r'}){die "no input r sequence file\n";}
if(!$para{'-o'}){$para{'-o'}='.';}
if(!$para{'-ut'}){$para{'-ut'}=3;}
if(!$para{'-q'}){$para{'-q'}=33;}
if(!$para{'-maxee'}){$para{'-maxee'}=1;}
if(!$para{'-trimf'}){$para{'-trimf'}=1;}
if(!$para{'-trimr'}){$para{'-trimr'}=1;}
if(!$para{'-maxdiff'}){$para{'-maxdiff'}=10;}
if(!$para{'-maxdiffpct'}){$para{'-maxdiffpct'}=$para{'-maxdiff'};}
if(!$para{'-threads'}){$para{'-threads'}=20;}
if(!$para{'-minl'}){$para{'-minl'}=300;}
if(!$para{'-maxl'}){$para{'-maxl'}=600;}
if(!$para{'-phix'}){$para{'-phix'}='no';}

if(-e "./merge_report.txt"){unlink "./merge_report.txt";}
my $outputfolder='preprocessed';
#my $outputfolder='0-original';

my $programDir = ppath()

##########processing reads##############
open SAT,">statistics.txt";

#uncompress files
if($para{'-f'}=~/.gz$/){
`gunzip -c $para{'-f'} >forward.fastq`;
}
else{
	`cp $para{'-f'} forward.fastq`;
}
if($para{'-r'}=~/.gz$/){
`gunzip -c $para{'-r'} >revers.fastq`;
}
else{
	`cp $para{'-r'} revers.fastq`;
}

#plot quality
print "Analyzing sequencing quality\n";
`$programDir/fastx_quality_stats -i forward.fastq -o forward_quality.txt -Q $para{'-q'}`;
`$programDir/fastx_quality_stats -i revers.fastq -o revers_quality.txt -Q $para{'-q'}`;
`$programDir/fastq_quality_boxplot_graph.sh -i forward_quality.txt -o forward_quality.png`;
`$programDir/fastq_quality_boxplot_graph.sh -i revers_quality.txt -o revers_quality.png`;

my $lines=`wc -l forward.fastq`;
$lines=~/(\d+)/;
print SAT "Total Raw: ",int($1/4),"\n";

print "Trimming low quality segments\n";
system("$programDir/fastx_trimmer -t $para{'-trimf'} -i forward.fastq -o forward1.fastq");
system("$programDir/fastx_trimmer -t $para{'-trimr'}  -i revers.fastq -o revers1.fastq");
system("mv forward1.fastq forward.fastq");
system("mv revers1.fastq revers.fastq");

my $lines=`wc -l forward.fastq`;
$lines=~/(\d+)/;
print SAT "Total pre-merging quality control: ",int($1/4),"\n";

#Pairing with vsearch
print "Pairing\n";

system("$programDir/vsearch -fastq_mergepairs forward.fastq -reverse revers.fastq -fastq_maxdiffpct $para{'-maxdiffpct'} -fastq_maxdiffs $para{'-maxdiff'} -fastq_truncqual $para{'-ut'} -fastaout good.fna -fastq_eeout -fastq_maxee $para{'-maxee'} -threads $para{'-threads'} -fastq_minmergelen $para{'-minl'}  -fastq_maxmergelen $para{'-maxl'} -fasta_width 0 2>&1 | tee merge_report.txt");

if($para{'-phix'}=~/y/i){
    warn("PhiX filtering is not currently available in VSearch\n");
    #system("$para{'-usearch'} -filter_phix good.fna -output filtered_reads.fna -threads $para{'-threads'} &>phix.txt");
    #system("mv filtered_reads.fna good.fna");
}

#$lines=`wc -l merged.fastq`;
#$lines=~/(\d+)/;
#print SAT "Total merged: ",$1/4,"\n";

$lines=`grep -c '>' good.fna`;
print SAT "Total good: ",$lines,"\n";
	
if(!-e $para{'-o'}){system("mkdir $para{'-o'}");}
system("mkdir ./$para{'-o'}/$outputfolder");
system("mv statistics.txt forward_quality.txt revers_quality.txt good.fna forward_quality.png phix.txt revers_quality.png merge_report.txt ./$para{'-o'}/$outputfolder");
system("mv merged.fastq ./$para{'-o'}/$outputfolder");
close SAT;
unlink 'forward.fastq','revers.fastq';

