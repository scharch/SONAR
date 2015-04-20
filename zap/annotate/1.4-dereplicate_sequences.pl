#!/usr/bin/perl
# performing usearch for sequences
use strict;

if(@ARGV<1){die "Usage: usearch.pl\n This script performs two steps of clustering to remove sequences potentially containing sequencing errors. please install usearch v7 or higher verion. options\n\t-pp path of the usearch program\n\t-min1 minimun coverage of a read to be kept in the first step of clustering, default:2\n\t-min2 minimun coverage of a read to be kept in the seconde step of clustering, default:3\n\t-f sequence file\n\t-threads number of cpu cores to use. Default:1\n";}

my %para=@ARGV;
if(!$para{'-min1'}){$para{'-min1'}=2;}
if(!$para{'-min2'}){$para{'-min2'}=3;}
if(!$para{'-pp'}){die "please give the correct path to usearch program\n";}
if(!$para{'-f'}){die "no input seq file\n";}
if(!$para{'-threads'}){$para{'-threads'}=1;}

my $output=&usearch($para{'-f'},$para{'-p'},$para{'-i'},$para{'-c'});
&changename($para{'-f'},$output);
#####################
sub usearch{
    my ($file,$program,$id,$co)=@_;	
    my $file_out=$file;
    $file_out=~s/\.fa//;	  	
    my %derep=();
	  	my %final_good=();
	  	system("$para{'-pp'} -derep_fulllength $file -threads $para{'-threads'} -fastaout $file_out.nonredundant.fa -sizeout -uc $file_out.cluster ");#> usearchlog.txt
	  	system("$para{'-pp'} -sortbysize $file_out.nonredundant.fa -minsize $para{'-min1'} -fastaout $file_out.nonredundant.fa");
	  	system("$para{'-pp'} -cluster_smallmem $file_out.nonredundant.fa -sortedby size -id $id -sizein -sizeout -uc $file_out.cluster -centroids $file_out.unique.fa ");#-query_cov $co -target_cov $co> usearchlog.txt
	  	%final_good=&derep("$file_out.cluster",$para{'-min2'},\%derep);
	  	system("rm $file_out.nonredundant.fa");

	  unlink "usearchlog.txt";
	  return "$file_out.unique.fa";
}


#####################
sub derep{
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
	 #open HH,"$file" or die "input raw seq file $file not found\n";
	 open YY,">$output";
	 #my $mark=0;
	 foreach(sort keys %{$good}){
	 	my $id=$_;
	 	foreach(@{$good->{$id}}){
	 	  print YY "$_\n";
	  }
	}
}

######################
sub changename{
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
	system("mv tempusearch.fa $seive");
}

