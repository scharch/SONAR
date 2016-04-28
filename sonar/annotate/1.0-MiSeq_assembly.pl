#!/usr/bin/perl
#    barcode file should be in the format: name\tsequence\n
use strict;
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);

#########checking parameters#######
if(@ARGV%2>0||!@ARGV){die "Usage: illumina_filtering.pl 
	-usearch path to usearch program
	-f forward_read_seq file, either compressed .gz or fastq file 
	-r reverse_strand, either compressed .gz or fastq file  
	-ut usearch_truncation_quality_threshold_before_mapping
	-maxee maximum number of expected errors in the paired reads, default:1.0 If the quality of your sequencing run is low, you may adjust to a higher cutoff
	-trimf number of nucleotides trimed off before assembly for forward mate, default:1
	-trimr number of nucleotides trimed off before assmebly for reverse mate, default:1
	-o output folder

Example:
1.0-MiSeq_assembly.pl -usearch usearch -f forward_read.fastq.gz -r reverse_read.fastq.gz -o ./

Created by Zizhang Sheng 2015-12-02.

Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

	";}
my %para=@ARGV;
if(!$para{'-f'}){die "no input f sequence file\n";}
if(!$para{'-r'}){die "no input r sequence file\n";}
if(!$para{'-o'}){$para{'-o'}='.';}
if(!$para{'-ut'}){$para{'-ut'}=3;}
if(!$para{'-maxee'}){$para{'-maxee'}=1;}
if(!$para{'-usearch'}){$para{'-usearch'}=ppath('usearch');}
if(!$para{'-fastx_quality_stats'}){$para{'-fastx_quality_stats'}=ppath('fastx') . '/fastx_quality_stats';}
if(!$para{'-fastq_quality_boxplot_graph.sh'}){$para{'-fastq_quality_boxplot_graph.sh'}=ppath('fastx')."/fastq_quality_boxplot_graph.sh";}
if(!$para{'-fastx_trimmer'}){$para{'-fastx_trimmer'}=ppath('fastx').'/fastx_trimmer';}
if(!$para{'-trimf'}){$para{'-trimf'}=1;}
if(!$para{'-trimr'}){$para{'-trimr'}=1;}

my $score="\!\"\#\$\%\&\047\(\)\*\+\,\-\.\/0123456789\:\;\<\=\>\?\@ABCDEFGHI";#Miseq quality scores
	  my %score=();
	  my $i=40;
	  while($score){
	  	$score{chop $score}=$i;
	  	$i--;
	  }
##########processing reads##############
open SAT,">statistics.txt";

#uncompress files
if($para{'-f'}=~/.gz$/){
`gunzip -c $para{'-f'} >forward.fastq`;
}
else{
	`mv $para{'-f'} forward.fastq`;
}
if($para{'-r'}=~/.gz$/){
`gunzip -c $para{'-r'} >revers.fastq`;
}
else{
	`mv $para{'-r'} revers.fastq`;
}

#plot quality
print "Analyzing sequencing quality\n";
`$para{'-fastx_quality_stats'} -i forward.fastq -o forward_quality.txt`;
`$para{'-fastx_quality_stats'} -i revers.fastq -o revers_quality.txt`;
`$para{'-fastq_quality_boxplot_graph.sh'} -i forward_quality.txt -o forward_quality.png`;
`$para{'-fastq_quality_boxplot_graph.sh'} -i revers_quality.txt -o revers_quality.png`;

my $lines=`wc -l forward.fastq`;
my @lines=split/[ \t]+/,$lines;
print SAT "Total Raw: ",int($lines[1]/4),"\n";
#&trim();
#trimming low quality segments
print "Trimming low quality segments\n";
system("$para{'-fastx_trimmer'}  -t $para{'-trimf'} -i forward.fastq -o forward1.fastq");
system("$para{'-fastx_trimmer'} -t $para{'-trimr'}  -i revers.fastq -o revers1.fastq");
system("mv forward1.fastq forward.fastq");
system("mv revers1.fastq revers.fastq");
&find_match_pair('forward.fastq','revers.fastq');
my $lines=`wc -l forward.fastq`;
my @lines=split/[ \t]+/,$lines;
print SAT "Total pre-merging quality control: ",int($lines[1]/4),"\n";

#Pairing with usearch
print "Pairing\n";
 	system("$para{'-usearch'} -fastq_mergepairs forward.fastq -reverse revers.fastq -fastq_truncqual $para{'-ut'} -fastqout merged.fastq");
 	system("$para{'-usearch'} -fastq_filter merged.fastq -fastaout $para{'-o'}_good.fna -fastq_maxee $para{'-maxee'} -eeout ");	
  $lines=`wc -l merged.fastq`;
  @lines=split/[ \t]+/,$lines;
	print SAT "Total merged: ",$lines[1]/4,"\n";
  $lines=`grep -c '>' $para{'-o'}_good.fna`;
  
	print SAT "Total good: ",$lines,"\n";
	if(!-e $para{'-o'}){system("mkdir $para{'-o'}");}
	system("mkdir ./$para{'-o'}/preprocessed");
	system("mv statistics.txt forward_quality.txt revers_quality.txt $para{'-o'}_good.fna forward_quality.png revers_quality.png ./$para{'-o'}/preprocessed");
  system("mv merged.fastq ./$para{'-o'}/preprocessed");
  close SAT;
unlink <*fastq>;
##################
sub rm_polymer{
    my $file=shift;
    my $i=0;
    my $seq='';
    my $count_lowcom=0;
    my $qual='';
    my $seq='';
    my $name='';
    my $count_trim=0;
    open YY,">trimed.fastq";
    open HH,"$file";	
    while(<HH>){
       if($_=~/^\@M/){
        $i=4;
        chomp;
        if($seq){
        	($seq,$qual)=&trim_lowquality($seq,$qual,$para{'-ws'},$para{'-lc'});
        	if($seq){
             print YY "$name\n$seq\n\+\n$qual\n";	
           }
          else{
           $count_trim++;	
           }
           ($name,$seq,$qual)=();
        }
        else{
        	if($name){
        	   $count_lowcom++;	
        	}
        }
        $name=$_;
      }
      elsif($i==3){
      	chomp;
      	$seq=$_;
      }	
      elsif($i==1){
      	chomp;
        $qual=$_;
      }
      $i--;
    }
  print "Total low complex for $file:$count_lowcom\nTotal trimming for $file:$count_trim\n";
	system("mv trimed.fastq $file");
}
sub trim{
	`echo ''>adaptor.fa`;
  system("fastq-mcf -q 0 -w 1 -x 20 -H -X -k 0 -o forward1.fastq  -o reverse1.fastq adaptor.fa forward.fastq revers.fastq");	
	system("mv forward1.fastq forward.fastq");
	system("mv reverse1.fastq revers.fastq");
	unlink 'adaptor.fa';
}

###################
sub trim_lowquality{
	  my ($seq,$qual,$widsize,$cutoff)=@_;
	  
    for(my $j=0;$j<length($qual);$j++){
       	my $str=substr($qual,$j,$widsize);
       	my $meanqual=0;
       	my $length=length($str);      	
       	while($str){
       		my $mark=chop $str;
       		 $meanqual+=$score{"$mark"};
       	}
       	#if($length==0){last;}
       	$i=$j;
       	if($meanqual/$length<$cutoff){
       		 last;
       	}
       	
    }
    
    $seq=substr($seq,0,$i);
    $seq=~s/(([N]{15,})|([T]{15,})|([C]{15,})|([G]{15,})|([A]{15,}))$//;
    $qual=substr($qual,0,length($seq));
    if(length($seq)>15){
       	return $seq,$qual;
    }
    else{
       return 0,0;
    }
}
###################
sub find_match_pair{
	  my ($f,$r)=@_;
	  my %f_reads=();
	  my %r_reads=();
    my %r_name=();
	  open RR,"$r" or die "no filtred reverse reads\n";
	  open FF,"$f" or die "no filtred forward reads\n";
	  my $id='';
    while(<RR>){
    	if($_=~/^\@M/){
    		my @line=split/[\t \/]/,$_;
    		$id=$line[0];
    		$r_name{$id}=$_;
    	}
    	else{
    		 $r_reads{$id}.=$_;
    	}
    }
    $r=~s/\.fastq//g;
	  $f=~s/\.fastq//g;	
    open RO,">$r\_matched.fastq";
	  open FO,">$f\_matched.fastq";
    	my $mark=0;
     while(<FF>){
       if($_=~/^\@M/){
       	my @line=split/[ \t\/]/,$_;
       	$id=$line[0];
       	 if($r_reads{$id}){
       	  print RO "$r_name{$id}$r_reads{$id}";
       	  print FO "$_";
       	  $mark=1;	
       	  }
       	  else{
       	     $mark=0;	
       	  }
       }
       elsif($mark==1){
       print FO "$_";	
      }	
    }
    system("mv $r\_matched.fastq $r\.fastq");
    system("mv $f\_matched.fastq $f\.fastq");
    #return "$r\_matched.fastq","$f\_matched.fastq";
}
