#!/usr/bin/perl
# performing usearch for sequences
use strict;

#########checking parameters#######
my $usage="
Usage:
This script performs two steps of clustering to remove sequences potentially containing sequencing errors. please install usearch v7 or higher verion. 
 options:
 \t-pu\tpath of the usearch program
 \t-min1\tminimun coverage of a read to be kept in the first step of clustering, default:2
 \t-min2\tminimun coverage of a read to be kept in the seconde step of clustering, default:3
 \t-f\tsequence file in fasta format
 \t-t\tnumber of threads to run the script. Default:1
Example:
1.4-dereplicate_sequences.pl -pu usearch -min1 2 -min2 3 -f ./test.fa -t 5

Created by Zizhang Sheng.

Copyright (c) 2011-2015 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
 ";
foreach(@ARGV){if($_=~/[\-]{1,2}(h|help)/){die "$usage";}}
if(@ARGV<1||@ARGV%2>0){die "Number of parameters are not right\n";
 }
my %para=@ARGV;
if(!$para{'-min1'}){$para{'-min1'}=2;}
if(!$para{'-min2'}){$para{'-min2'}=3;}
if(!$para{'-pu'}){die "please give the correct path to usearch program\n";}
if(!$para{'-f'}){die "no input seq file\n";}
if(!$para{'-t'}){$para{'-t'}=1;}
#########do calculation##########
my $output=&usearch($para{'-f'},$para{'-p'},$para{'-i'},$para{'-c'});
&changename($para{'-f'},$output);

########subrutines#############
sub usearch{#do the two steps of clustering
    my ($file,$program,$id,$co)=@_;	
    my $file_out=$file;
    $file_out=~s/\.fa.*//;	  	
    my %derep=();
	  	my %final_good=();
	  	system("$para{'-pu'} -derep_fulllength $file -threads $para{'-t'} -fastaout $file_out.nonredundant.fa -sizeout -uc $file_out.cluster ");#first step on higher identity
	  	system("$para{'-pu'} -sortbysize $file_out.nonredundant.fa -minsize $para{'-min1'} -fastaout $file_out.nonredundant.fa");
	  	system("$para{'-pu'} -cluster_smallmem $file_out.nonredundant.fa -sortedby size -id $id -sizein -sizeout -uc $file_out.cluster -centroids $file_out\_unique.fa ");#second step on higher identity
	  	%final_good=&derep("$file_out.cluster",$para{'-min2'},\%derep);#remove low coverage clusters
	  	system("rm $file_out.nonredundant.fa");
  	  unlink "usearchlog.txt";
	  return "$file_out\_unique.fa";
}


#####################
sub derep{#This is used to remove low abundant clusters after the second step of clustering
    my ($file,$cutoff,$ids)=@_;
    open HH,"$file" or die "Usearch clustering file $file not found\n";
    my %ids=();
    if($ids){
    	while(<HH>){
    		my @line=split/[ \t]+/,$_;
    		$line[8]=~s/\;.+//;
    		$line[9]=~s/\;.+//;
        if($_=~/^C/){ 
        	 if($line[2]<$cutoff){
        	 	foreach(@{$ids{$line[8]}}){
        	   delete $ids->{$_};	
        	 }
        	}
        }
        elsif($_=~/^S/){
	    	  push @{$ids{$line[8]}},$line[8];
	      }
	      else{
	    	  push @{$ids{$line[9]}},$line[8];
	      }
      }
      return %{$ids};
    }
    else{
    	while(<HH>){
    		my @line=split/[ \t]+/,$_;
    		$line[8]=~s/\;.+//;
    		$line[9]=~s/\;.+//;
      if($_=~/^C/){       	
	     	 if($line[2]<$cutoff){
	     	     delete $ids{$line[8]};	     
	     	}
	    }	
	    elsif($_=~/^S/){
	    	push @{$ids{$line[8]}},$line[8];
	    }
	    else{
	    	push @{$ids{$line[9]}},$line[8];
	    }
	   	
     }
     return %ids;
    }	
	
}
#####################
sub write_raw{
	 my ($good,$output)=@_;
	 open YY,">$output";
	 foreach(sort keys %{$good}){
	 	my $id=$_;
	 	foreach(@{$good->{$id}}){
	 	  print YY "$_\n";
	  }
	}
}

######################
sub changename{#change sequence names back to the input sequence name
  my ($input,$seive)=@_;
  open HH,"$seive" or die "usearch file $seive not found\n";
	open YY,"$input" or die "Original seq file $input not found\n";
	open ZZ,">tempusearch.fa";
	my %id=();
	my $id='';
	my $mark=0;
	my %size=0;
	while(<HH>){
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
	while(<YY>){
		if($_=~/>([^\t \;\n\r]+)/){
			chomp;
			my $k=$1;
			$k=~s/centroid\=//;
		   if($id{$k}==1){
		   	  $mark=1;	
		   	  print ZZ "$_\n";
		   }
		 	 else{
		 	    $mark=0;	
		 	}
		}
		elsif($mark==1){
		  print ZZ "$_";	
		}
	}
	if(-d "./output/sequences/nucleotide"){#move output files to standard pipeline folders
		system("mv tempusearch.fa ./output/sequences/nucleotide/$seive");
	}
	else{
	  system("mv tempusearch.fa $seive");
  }
}

