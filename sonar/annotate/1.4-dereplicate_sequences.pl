#!/usr/bin/perl
# performing usearch for sequences
use strict;
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);
use version qw/logCmdLine/;


#########checking parameters#######
my $usage="
Usage:
This script performs two steps of clustering to remove sequences potentially containing sequencing errors. The first step finds duplicates, remove reads with coverage lower then cutoff which may contain sequencing errors, calculate read coverage for each cluster. The second step further cluster the filtered sequences using lower sequence identity cutoff. At the meantime, the second clustering will use reads with high coverage as centroids by assuming biological reads are coming from cDNA with many identical copies while reads containing sequencing errors has very low coverage. please install usearch v7 or higher verion. 

options:
	-id	 percent sequence identity used for the second step of clustering, default:0.99
	-min1	minimun sequencing coverage of a read to be kept in the first step of clustering, default:0
	-min2	minimun sequencing coverage of a read to be kept in the seconde step of clustering, default:3
	-f	sequence file in fasta format
	-t	number of threads to run the script. Default:1
	-s	Optional (Default: 1). A value of 2 is used when clustering a sequence file of interest, unrelated to the pipeline.
Example:
1.4-dereplicate_sequences.pl -pu usearch -min1 2 -min2 3 -f ./test.fa -t 5

Created by Zizhang Sheng.
Edited to switch to VSearch by Chaim A Schramm 2018-07-30.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
 ";
 
foreach(@ARGV){if($_=~/^[\-]{1,2}(h|help)/){die "$usage";}}
if(@ARGV%2>0){die "Number of parameters are not right\n$usage";
 }

&logCmdLine($0,@ARGV);

my %para=@ARGV;
if(!$para{'-min1'}){$para{'-min1'}=0;}
if(!$para{'-min2'}){$para{'-min2'}=3;}
$para{'-pu'}=ppath('usearch');
if(!$para{'-pu'}){die "please give the correct path to usearch program in the PPvars.pm file\n";}
if(!$para{'-f'}){
	my @files=<./output/sequences/nucleotide/*goodVJ.fa>;
	if(-e "$files[0]"){$para{'-f'}=$files[0];}
	else{die "no input seq file\n";}
}
if(!$para{'-s'}){$para{'-s'}=1;}
if(!$para{'-t'}){$para{'-t'}=1;}
if(!$para{'-id'}){$para{'-id'}=0.99;}
#########do calculation##########

my $output=&usearch($para{'-f'},$para{'-p'},$para{'-id'});
my $ids_unique=&changename($para{'-f'},$output);#
my $file=$para{'-f'};
$file=~s/goodVJ/goodCDR3/;
if(-e "$file"&& $para{'-s'}==1){
    print "Finding nucleotide sequences of CDR3s of unique sequences......\n";	
    &read_fasta($file,$ids_unique);
}
$file=~s/nucleotide/amino_acid/;
if(-e "$file"&& $para{'-s'}==1){
    print "Finding amino acid sequences of CDR3s of unique sequences......\n";	
    &read_fasta($file,$ids_unique);
}

$file=~s/goodCDR3/goodVJ/;
if(-e "$file"&& $para{'-s'}==1){
    print "Finding amino acid sequences of unique sequences......\n";	
    &read_fasta($file,$ids_unique);
}


########subrutines#############
sub usearch{#do the two steps of clustering
    my ($file,$id)=@_;	
    my $file_out=$file;
    $file_out=~s/\.fa.*//;	  	
    my %derep=();
    my %final_good=();

    #first step on higher identity
    system("$para{'-pu'} -derep_fulllength $file -threads $para{'-t'} -output $file_out\_dedup.fa -sizein -sizeout -uc $file_out.cluster -minuniquesize $para{'-min1'}");
    if(-z "$file_out\_dedup.fa"){die "No duplicate sequence found in your input sample.\n";}
    
    #second clustering step
    #vsearch cluster_size seems to be the equivalent of the --sortby size option in usearch
    system("$para{'-pu'} -cluster_size $file_out\_dedup.fa -id $para{'-id'} -sizein -sizeout -uc $file_out.cluster -centroids $file_out\_nonredundant.fa");
    system("$para{'-pu'} -sortbysize $file_out\_nonredundant.fa -output $file_out\_unique.fa -minsize $para{'-min2'} -fasta_width 0");
    if(-z "$file_out\_unique.fa"){die "No cluster found for your sequences.\n";}
    system("rm $file_out\_dedup.fa");
    system("rm $file_out\_nonredundant.fa");
    
    unlink "usearchlog.txt";
    return "$file_out\_unique.fa";
    
}


#####################
sub write_raw{
	 my ($good,$output)=@_;
	 open OT,">$output";
	 foreach(sort keys %{$good}){
	 	my $id=$_;
	 	foreach(@{$good->{$id}}){
	 	  print YY "$_\n";
	  }
	}
	close OT;
}

######################
sub changename{#change sequence names back to the input sequence name
    my ($input,$seive)=@_;
    open SEI,"$seive" or die "usearch file $seive not found\n";
    my %state=();
    my $unique_column=0;
    my @stat_file=<./output/tables/*all_seq_stats.txt>;
    my %id=();
    my $id='';
    my $mark=0;
    my %size=0;
    while(<SEI>){
	if($_=~/>([^\t \;]+)/){
	    chomp;
	    $id=$1;
	    $id=~s/centroid\=//;
	    $id{$id}=1;
	    if($_=~/size\=([0-9]+)/){
		$size{$id}=$1;
	    }
	}
    }
    close SEI;
    
    if(-e "$stat_file[0]"&& $para{'-s'}==1){#write statistic info to ./output/table/project_all_seq_stats.txt
	&rm_r($stat_file[0]);
	open STi,"$stat_file[0]";	
	my $title=<STi>;
	my @categories = split/\t+/,$title;
	if($title!~/Unique\tsize/){	  
	    $unique_column = scalar(@categories);
	    chomp $title;
	    $title.="\tUnique\tsize\n";
	} else {
	    #replace existing data based on new run / changed clustering criteria
	    for my $c (0 .. $#categories) {
		if ($categories[$c] =~ /Unique/) { 
		    $unique_column = $c;
		    last;
		}
	    }
	}
	open STo,">stats.txt";
	print STo "$title";

	while(<STi>){
	    if($_=~/^ID\t/||$_!~/[\d\w]/){next;}
	    chomp;
	    ~s/[\r\n]//g;
	    my @l=split/\t+/,$_;
	    $state{$l[0]}=$_;
	}

    }
    open YY,"$input" or die "Original seq file $input not found\n";
    open ZZ,">tempusearch.fa";
    
    while(<YY>){
	if($_=~/>([^\t \;\n\r]+)/){
	    chomp;
	    my $k=$1;
	    $k=~s/centroid\=//;
	    if($id{$k}==1){
		$mark=1;	
		print ZZ "$_ size=$size{$k}\n";
		
	    }
	    else{
		$mark=0;	
	    }
	}
	elsif($mark==1){
	    print ZZ "$_";	
	}
    }	
    close YY;
    close ZZ;

    if($unique_column > 0) {
	foreach(sort {$a<=>$b} keys(%state)){ 
	    my $uniq = "NA";
	    my $size = "NA";
	    if ($id{$_}) {
		$uniq = "T";
		$size = $size{$_};
	    }
      print STo "$state{$_}\t$uniq\t$size\n";
	}		
    
	close STi;
	close STo;
	system("mv stats.txt $stat_file[0]");
    }

    if(-d "./output/sequences/nucleotide"&&$seive!~/output\/sequences\/nucleotide/&&$para{'-s'}==1){#move output files to standard pipeline folders
	system("mv tempusearch.fa ./output/sequences/nucleotide/$seive");
    }
    else{
	system("mv tempusearch.fa $seive");
    }
    return \%id;
}

sub rm_r{#remove \r at line end
      my $file=shift;
   open HH,"$file" or die "rm_r didn't find the file $file\n";
   open YY,">rmtem.txt";
   while(<HH>){
      ~s/\r/\n/g;
      print YY "$_";
   }
    close HH;
    close YY;
  system("mv rmtem.txt $file");	
}

sub read_fasta{
    my ($file,$seqid)=@_;
    open HH,"$file";
    $file=~s/\.fa$/\_unique.fa/;
    open YY,">$file";
    my $mark=0;
    while(<HH>){
    	  if($_=~/>([^\t \;]+)/){
    	  	if($seqid->{$1}){
    	  	   print YY $_;
    	  	   $mark=1;	
    	  	}
    	  	else{
    	  	  $mark=0;	
    	  	}
    	  }
    	  elsif($mark==1){
    	  	print YY $_;
    	  }
    }	
	close HH;
	close YY;
}
