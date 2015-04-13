#!/usr/bin/perl
# performing usearch for sequences
use strict;

if(@ARGV<1){die "Usage: usearch.pl please install usearch v7 or higher verion. options\n\t-pp path of the usearch program\n\t-p [1,default cluster_fast,use longest seq as seed or center of cluster;2,use sequences with high coverage as seed or center of clusters;3, remove identical sequences in the dataset]\n\t-f sequence file\n\t-i sequence identity cutoff, default:0.975\n\t-low remove sequences with coverage lower than cutoff\n\t-threads number of cpu cores to use. Default:1\n";}

my %para=@ARGV;
if(!$para{'-p'}){$para{'-p'}=1;}
if(!$para{'-pp'}){die "please give the correct path to usearch program\n";}
if(!$para{'-f'}){die "no input seq file\n";}

if(!$para{'-i'}){$para{'-i'}=0.975;}
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
	  if($program==1){
	  	system("$para{'-pp'} -cluster_fast $file -threads $para{'-threads'}  -id $id  -uc $file_out.cluster -sizeout -centroids $file_out.unique.fa > usearchlog.txt");#-clusters tep-query_cov $co -target_cov $co
	  }
	  elsif($program==2){#good sequences potentially don't have errors
	  	system("$para{'-pp'} -derep_fulllength $file -threads $para{'-threads'} -fastaout $file_out.nonredundant.fa -sizeout -uc $file_out.cluster ");#> usearchlog.txt
	  	system("$para{'-pp'} -sortbysize $file_out.nonredundant.fa -minsize 1 -fastaout $file_out.nonredundant.fa");
	  	system("$para{'-pp'} -cluster_smallmem $file_out.nonredundant.fa -sortedby size -id $id -sizein -sizeout -uc $file_out.cluster -centroids $file_out.unique.fa ");#-query_cov $co -target_cov $co> usearchlog.txt
	  	system("rm $file_out.nonredundant.fa");
	  }
	  elsif($program==3){
	  	system("$para{'-pp'} -derep_fulllength $file -threads $para{'-threads'} -fastaout $file_out.unique.fa -uc $file_out.cluster -sizeout > usearchlog.txt");
	  	
	  }
	  unlink "usearchlog.txt";
	  if($para{'-low'}){print "Removing low quality reads\n";
	  &remove_lowquality("$file_out.unique.fa","$file_out.cluster",$para{'-low'});}
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

#####################
sub remove_lowquality{
     my($seq,$cluster,$cutoff)=@_;
     open HH,"$seq" or die "sequence file not right\n";
     open YY,"$cluster"or die "cluster file not right\n";
     open ZZ,">temp_cutoff.fa";
     my %mark=();	
     my $count=0;
	   while(<YY>){
	     if($_=~/^C/){
	       my @line=split/[ \t]+/,$_;	
	     	 if($line[2]>=$cutoff){
	     	 $line[8]=~s/\;.+//;
	     	 $count++;
	     	   $mark{$line[8]}=1;	
	     	}
	    }	
	  }
	  print "$count sequences passed coverage cutoff\n";
	  my $mark=0;
	  while(<HH>){
	  	if($_=~/>(.+)/){
	  		my @line=split/[ \t\;]/,$1;
	  		if($mark{$line[0]}){
	  			print ZZ "$_";
	  			$mark=1;
	  		}
	  		else{
	  			$mark=0;
	  		}
	  	}
	  	elsif($mark==1){
	  		print ZZ $_;
	  	}
	  }
	  system("mv temp_cutoff.fa $seq");
}