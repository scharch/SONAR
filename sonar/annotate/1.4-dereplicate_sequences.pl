#!/usr/bin/perl
# performing usearch for sequences
use strict;
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);

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
Example:
1.4-dereplicate_sequences.pl -pu usearch -min1 2 -min2 3 -f ./test.fa -t 5

Created by Zizhang Sheng.

Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
 ";
 
foreach(@ARGV){if($_=~/[\-]{1,2}(h|help)/){die "$usage";}}
if(@ARGV%2>0){die "Number of parameters are not right\n$usage";
 }
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
if(!$para{'-t'}){$para{'-t'}=1;}
if(!$para{'-id'}){$para{'-id'}=0.99;}
#########do calculation##########

my $output=&usearch($para{'-f'},$para{'-p'},$para{'-id'});
my $ids_unique=&changename($para{'-f'},$output);
my $file=$para{'-f'};
$file=~s/goodVJ/goodCDR3/;
if(-e "$file"){
    print "Finding nucleotide sequences of CDR3s of unique sequences......\n";	
    &read_fasta($file,$ids_unique);
}
$file=~s/nucleotide/amino_acid/;
if(-e "$file"){
    print "Finding amino acid sequences of CDR3s of unique sequences......\n";	
    &read_fasta($file,$ids_unique);
}

$file=~s/goodCDR3/goodVJ/;
if(-e "$file"){
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
	  	system("$para{'-pu'} -derep_fulllength $file -threads $para{'-t'} -fastaout $file_out\_unique.fa -sizeout -uc $file_out.cluster ");#first step on higher identity
	  	if(-z "$file_out\_unique.fa"){die "No duplicate sequence found in your input sample.\n";}
	  	system("$para{'-pu'} -sortbysize $file_out\_unique.fa -minsize $para{'-min1'} -fastaout $file_out.nonredundant.fa");
	  	if(-z "$file_out.nonredundant.fa"){die "No duplicate sequence found in your input sample.\n";}
	  	system("$para{'-pu'} -cluster_fast $file_out.nonredundant.fa -sort size -id $para{'-id'} -sizein -sizeout -uc $file_out.cluster -centroids $file_out\_unique.fa ");#second step on higher identity
	  	if(-z "$file_out\_unique.fa"){die "No cluster found for your sequences.\n";}
	  	system("$para{'-pu'} -sortbysize $file_out\_unique.fa -minsize $para{'-min2'} -fastaout $file_out.nonredundant.fa");	  	
	  	system("mv $file_out.nonredundant.fa $file_out\_unique.fa");
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
    
    if(-e "$stat_file[0]"){#write statistic info to ./output/table/project_all_seq_stats.txt
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
	    $state{$l[0]}=[@l];
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
	    my @line = @{$state{$_}};
	    $line[$unique_column] = $uniq;
	    $line[$unique_column + 1] = $size;

	    print STo join("\t",@line) . "\n";
	}		
    
	close STi;
	close STo;
	system("mv stats.txt $stat_file[0]");
    }

    if(-d "./output/sequences/nucleotide"&&$seive!~/output\/sequences\/nucleotide/){#move output files to standard pipeline folders
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
