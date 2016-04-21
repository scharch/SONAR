#!/usr/bin/env perl
use strict;
use threads;
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);

my $usage="Usage: 
This script is used to calculate sequence identity between germline V, antibody gene and reads or between antibody CDR3 and read CDR3.
Options:
	-f\tsequence file, if reads come from different germline V genes, the name of each read should 
	  \tcontain a field saying it's germline V gene for somatic hypermutation calculation. Otherwise, no 
	  \tsomatic hypermutation is calculated. fasta file output from pipeline is ok. For example, 
	  \t'>00000089 V_gene=IGHV1-2*02'. Please seperate '00000089 and 'V_gene' with space or tab.
	-g\tgermline V gene file, optional.Default: IgHV.fa in the program folder.
	-a\tfasta file has the sequence of interested antibody. optional.
	-t\tthreads, default:5
	-npt\tnumber of sequences per thread. default:1000
	-p\tprotein or DNA sequence. default: DNA
	-ap\t name of the program for sequence alignment. muscle or clustalo or mafft. required. Based on our 
	   \texperience, muscle is ~2 fold faster than clustalo. Clustalo version of 1.2.0 or higher is required.
	-pu\tpath to usearch program to remove duplicates in the read file. Optional.
  -CDR3\tWhether calculating sequence identity between CDR3s
  
Example:
2.1-calculate_id-div.pl -f test.fa -g germline.fa -a antibody.fa -t 5 -npt 1000 -p DNA -ap muscle -pu usearch

Created by Zizhang Sheng.

Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.	
 ";
foreach(@ARGV){if($_=~/[\-]{1,2}(h|help)/){die "$usage";}}
if(@ARGV%2>0){die "$usage"; }
my %para=@ARGV;
if(!$para{'-t'}){$para{'-t'}=5;}
if(!$para{'-npt'}){$para{'-npt'}=1000;}
$para{'-ap'}=ppath($para{'-ap'});
$para{'-pu'}=ppath($para{'-pu'});
if(!$para{'-ap'}){
	 if(ppath('muscle')){
	   $para{'-ap'}=ppath('muscle');	
	}else{
	 die "please select a program for sequence alignment\n";
	}
	}
if(!$para{'-g'}){
	$para{'-g'}=ppath('sonar')."/germDB/IgHV.fa";
	if(-e $para{'-g'}){
	  	print "using $para{'-g'} as germline file\n";
	 }
	else{
	   warn "No calculation for hypermutation\n";
	 }
	}
if(!$para{'-p'}){$para{'-p'}='DNA';}
if($para{'-a'}) {
    if ( ! -e $para{'-a'}){warn "file $para{'-a'} doesn't exist.\n";}
} else { warn "No native sequence file specified; only divergence will be calculated...\n"; }
if(!$para{'-f'}|| ! -e $para{'-f'}){
	my @files=<./output/sequences/nucleotide/*goodVJ_unique.fa>;
	if(-e "$files[0]"){
	   $para{'-f'}=$files[0];	
	   print "Using sequence file $files[0]\n";
	 }
	else{
	    @files=<./output/sequences/nucleotide/*goodVJ.fa>;
	    if(-e "$files[0]"){
		$para{'-f'}=$files[0];	
		print "Using sequence file $files[0]\n";
	    }else{
		die "Sequence file $para{'-f'} doesn't exist.\n";
	    }
	}
}
#if($para{'-g'}&& ! -e $para{'-g'}){die "file $para{'-g'} doesn't exist.\n";}
my %germ_db=();

############Reading seqs####################
&rm_r($para{'-g'});
if($para{'-a'}){&rm_r($para{'-a'});}
my ($germV,$germg)=&readfasta($para{'-g'});
my ($anti,$antigerm)= $para{'-a'} ? &readfasta($para{'-a'}) : ({},{});

###############Processing###################
print "processing $para{'-f'}...\n";
my $changefile=$para{'-f'};
$changefile=~s/\.fa.*//;
my @filepath=split/\//,$changefile;
$changefile=pop @filepath;
my $output_id=$changefile;
my $output_cov='';
if($para{'-CDR3'}){
  	$output_id.="_CDR3-id.tab";
}
else{
    $output_id.="_id-div.tab";
}
if(-d "./output/tables/"){#detect output folders
  open YY,">./output/tables/$output_id";
  open ZZ,">./output/tables/$changefile\_coverage.tab";  
  $output_id="./output/tables/$output_id";
  $output_cov="./output/tables/$changefile\_coverage.tab";
}
else{
	open YY,">$output_id";
  open ZZ,">$changefile\_coverage.tab";
  $output_cov="$changefile\_coverage.tab";
}

if($para{'-g'}){print YY "ID\tgerm_div";}
else{print YY "ID";}
foreach(sort keys %{$anti}){
    print YY "\t$_";
}
print YY "\n";
my $file_calculation=$para{'-f'};
if($para{'-pu'}){#dereplicate
  system("$para{'-pu'} -derep_fulllength $para{'-f'} -threads $para{'-t'} -fastaout $changefile\_unique.fa -uc $changefile.cluster -sizeout > usearchlog.txt");
  $file_calculation="$changefile\_unique.fa";
}
open READs,"$file_calculation"or die "$file_calculation not found\n";#read sequences and assign to threads
    my $i=1;
    my $j=0;
    my $count_thread=0;
    my %name=();
    my %seq=();
    my $id='';
    my %readgerm=();
    while(<READs>){
        chomp;
        if($_=~/\>/){
               if($i==1+$para{'-npt'}){
               	 while(threads->list(threads::running)>=$para{'-t'}){#threading
               	 	foreach(threads->list(threads::joinable)){
                     my @result=$_->join();
                     print  YY $result[0];
                     print ZZ $result[1];
                      }
    		            sleep(1);
    	           } 
                threads->create({'context' => 'list'},\&calculation,\%seq,\%readgerm,$germV,$anti);
             	  $i=0; 
             	  print "processing thread $count_thread\n";
             	  $count_thread++;
            }
            if($i==0){
                %seq=();
             	  %readgerm=();	
            }
                $i++;     
             my @line=split/[ \t,]+/,$_;
             $id=substr($line[0],1,);
             #$line[1]=~s/\,.+//g;
             #$readgerm{$id}=$line[1];
	       for my $field (@line) { if ($field=~/V_gene=/) { $readgerm{$id}=$field; } }
            $readgerm{$id}=~s/\-//g;
            $readgerm{$id}=~s/V\_gene\=//g;
            
        }
        else{
            $seq{$id}.=$_;
        }
    }
    
 if(%seq){

  threads->create({'context' => 'list'},\&calculation,\%seq,\%readgerm,$germV,$anti);
 	
}   

while(threads->list()){#waiting for all threads finish
    foreach(threads->list(threads::joinable)){
        my @result=$_->join();
        print  YY $result[0];
        print ZZ $result[1];
    }
    sleep(1);
}

if($para{'-pu'}){#add redundant sequences back
  &recover($changefile,$output_id,$output_cov);
  unlink "$changefile\_unique.fa","$changefile.cluster","usearchlog.txt";
}
&add_column_to_statistic($output_id);
##################################
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
###################################
sub calculation{#main subrutine in each thread to do calculation
    my ($seq,$germ,$germV,$anti)=@_;
    my $identity='';
    my $coverage='';
  my %seq1=%{$seq};
  my %germ1=%{$germ};
  my %germV1=%{$germV};
  my $num_germ=()=keys %germV1;
    foreach(sort keys %{$seq}){#print "$_\n";
    	my $germ_seq='';
    	if($num_germ==1){my @k=values %germV1;$germ_seq=$k[0];}
    	else{$germ_seq=$germV1{$germ1{$_}};}
        my ($iden,$cov)=&paired_identity($_,$seq1{$_},$germ_seq,$anti);
        $identity.=$iden."\n";
        $coverage.=$cov."\n";
    }
    #chomp($identity);
    return $identity,$coverage;
}
#################################
sub paired_identity{#calculate sequence identity
    my($id,$seq,$germ,$anti)=@_;
    my @identity=();
    my @coverage=();
    push @identity,$id;
    push @coverage,$id;
    if($germ ){#calculate germline divergence and coverage if germline sequence is given
    	my ($seqg,$seqr)=&aln($germ,$seq,$para{'-p'});
      my $div=sprintf("%.2f",100-&identity($seqg,$seqr));
      my $cov=sprintf("%.2f",&coverage($seqg,$seqr));
      push @coverage,$cov; 
      push @identity,$div;
    }
    elsif($para{'-g'}&&!$germ){
    	push @coverage,'NA'; 
      push @identity,'NA';
    }
    
    if($anti){#calculate identity and coverage to given antibody sequence
    foreach(sort keys %{$anti}){
    	my ($seqg,$seqr)=&aln($anti->{$_},$seq,$para{'-p'});
        push @identity,&identity($seqg,$seqr);
      my $cov=sprintf("%.2f",&coverage($seqg,$seqr));
        push @coverage,$cov; 
     }
    }
    my $identity=join "\t",@identity;
    my $coverage=join "\t",@coverage;
    return $identity,$coverage;
}
#################################
sub readfasta{# read in sequences and germline assign info from fasta file
    my $file=shift;
    my %seq=();
    my %seqgerm=();
    my $id='';
    open HH,"$file" or warn "Sequence file $file does not exist\n";
    while(<HH>){
        chomp;
        if($_=~/>(.+)/){
            my @ll=split/[\t ]+/,$1;
            $id=$ll[0];
            $ll[1]=~s/\-//g;
            $id=~s/\-//g;
            $seqgerm{$id}=$ll[1];
        }
        else{
            $seq{$id}.=$_;
        }
    }
    return \%seq,\%seqgerm;
}

#####################################
sub aln{#pairwise sequence alignment using input program
    my @seq=@_;
    if($para{'-ap'} =~/muscle/i){
      return &muscle(@seq);	
    }
    elsif($para{'-ap'} =~/clustalo/i){
    	return &clustalo(@seq);
    }
    elsif($para{'-ap'} =~/mafft/i){
      	return &mafft(@seq);
    }
    else{
      die "$para{'-ap'} is not a recognized program\n";	
    }
	
}
#################################
sub mafft{#sequence alignment using mafft 
    my ($seq1,$seq2,$type)=@_;    	
	  my %seq=();
	  if($seq1=~/\*/||$seq2=~/\*/){print "$seq1\n$seq2\n";}
    if($seq1!~/[A-Za-z]/||$seq2!~/[A-Za-z]/){warn "$seq1\n$seq2\n no sequence for alignment\n";return '','';last;}
    my @aln=`echo ">seq1\n$seq1\n>seq2\n$seq2\n"| $para{'-ap'} --quiet -`;
    while(!$aln[3]){@aln=`echo ">seq1\n$seq1\n>seq2\n$seq2\n"| $para{'-ap'} --quiet - >@aln`;}
    foreach(@aln){
        chomp;
        if($_=~/>(.+)/){
            $id=$1;
        }
        else{
            $seq{$id}.=$_;
        }
    }
        $seq{'seq1'}=~s/([a-z]+)/\U$1/g;
    $seq{'seq2'}=~s/([a-z]+)/\U$1/g;
    if($seq2&&$seq1&&(!$seq{'seq1'}||!$seq{'seq2'})){die "Mafft aln error: >seq1\n$seq1>seq2\n$seq2\n";}
    return $seq{'seq1'},$seq{'seq2'};
}
###################
sub clustalo{#sequence alignment using clustalo
    my ($seq1,$seq2,$type)=@_;    	
	  my %seq=();
	  if($seq1=~/\*/||$seq2=~/\*/){print "$seq1\n$seq2\n";}
    if($seq1!~/[A-Za-z]/||$seq2!~/[A-Za-z]/){warn "$seq1\n$seq2\n no sequence for alignment\n";return '','';last;}
    my @aln=`echo ">seq1\n$seq1\n>seq2\n$seq2\n"| $para{'-ap'} -i - --wrap 999999`;
    my $id='';
    while(!$aln[3]){@aln=`echo ">seq1\n$seq1\n>seq2\n$seq2\n"| $para{'-ap'} -i - --wrap 999999`;}
    foreach(@aln){
        chomp;
        if($_=~/>(.+)/){
            $id=$1;
        }
        else{
            $seq{$id}.=$_;
        }
    }

    if($seq2&&$seq1&&(!$seq{'seq1'}||!$seq{'seq2'})){die "Clustalo aln error: >seq1\n$seq1>seq2\n$seq2\n";}
    return $seq{'seq1'},$seq{'seq2'};
}
#################################
sub muscle{#sequence alignment using muscle
    my ($seq1,$seq2,$type)=@_;
    my %seq=();
    if($seq1=~/\*/||$seq2=~/\*/){print "$seq1\n$seq2\n";}
    if($seq1!~/[A-Za-z]/||$seq2!~/[A-Za-z]/){warn "$seq1\n$seq2\n no sequence for alignment\n";return '','';last;}
    my @aln=`echo ">seq1\n$seq1\n>seq2\n$seq2\n"| $para{'-ap'} -seqtype $type -quiet -gapopen -1000`;
    my $id='';
    while(!$aln[3]){@aln=`echo ">seq1\n$seq1\n>seq2\n$seq2\n"| $para{'-ap'} -seqtype $type -quiet -gapopen -1000`;}
    foreach(@aln){
        chomp;
        if($_=~/>(.+)/){
            $id=$1;
        }
        else{
            $seq{$id}.=$_;
        }
    }
    if($seq2&&$seq1&&(!$seq{'seq1'}||!$seq{'seq2'})){die "muscle aln error: >seq1\n$seq1>seq2\n$seq2\n";}
    return $seq{'seq1'},$seq{'seq2'};
}
#################################
sub identity{#include gaps
    my ($seq1,$seq2,$type)=@_;
    if(!$seq1||!$seq2||$seq2!~/[A-Z]/||$seq1!~/[A-Z]/){warn "Nothing in the sequence for identity calculation:\n>seq1\n$seq1\n>seq2\n$seq2#\n";return 0;last;}
    my $score=0;
    my $leng=0;
    my $misc='';
    if($type eq 'DNA'){$misc='N';}
    else{$misc='X';}
    my $start=0;#remove terminal gaps
    my $end=length($seq1);
    if($seq1=~/^(\-+)/){
        $start=length($1);
    }
    elsif($seq2=~/^(\-+)/){
        $start=length($1);
    }
    
    if($seq2=~/(\-+)$/){$end=$end-length($1);}
    elsif($seq1=~/(\-+)$/){
        $end=$end-length($1);
    }
    $seq1=substr($seq1,$start,$end-$start);
    $seq2=substr($seq2,$start,$end-$start);
    for(my $j=0;$j<length($seq1);$j++){
        if(substr($seq1,$j,1) ne $misc &&substr($seq2,$j,1) ne $misc){
            $leng++;
            if(substr($seq1,$j,1) eq substr($seq2,$j,1)){
                $score++;
            }
        }
    }
    $score=sprintf("%.2f",100*$score/$leng);
    if($score !~/[0-9]/){$score='NA';}
    return $score;
}

################################
sub coverage{#calculate sequence coverage between pair of seqeunces
    my ($seq1,$seq2)=@_;
        if(!$seq1||!$seq2){warn "Nothing in the sequence for coverage calculation\n";}
    my $start=0;
    my $end=length($seq1);
    if($seq1=~/^[A-Z]/){}
    elsif($seq1=~/^(\-+)/){
        $start=length($1);
    }
    if($seq1=~/[A-Z]$/){}
    else{$seq1=~/(\-+)$/;
        $end=$end-length($1);
    }
    my $str1=substr($seq1,$start,$end-$start);
    my $str2=substr($seq2,$start,$end-$start);
    $str1=~s/[\-\.]//g;
    $str2=~s/[\-\.]//g;
    my $c=length($str2)/length($str1);
    if($c>1){$c=1/$c;}
    my $d=sprintf ("%3.2f",100*$c);
    if($d!~/[0-9]/){$d='NA';}
 return $d;       
}
################################
sub recover{#add back all the identity and coverage calculation for redundant reads to the final output
   my ($filename,$output_id,$output_cov)=@_;
   my %ids=();
   &rm_r("$filename.cluster");
  open HH,"$filename.cluster" or die "Usearch clustering file not exist";
  
   while(<HH>){
    if($_=~/^H/){
    	chomp;
      my @line=split/\t/,$_;
      push @{$ids{$line[9]}},$line[8];
     }	
   }  	
	open CD,"$output_id" or die "$output_id not exist\n";#add back identity for redundant reads
	open ZZ,">tempchange.txt";
	foreach(<CD>){
	  my @line=split/\t/,$_;
	  my $line=$_;
	  foreach(@{$ids{$line[0]}}){
	  	  my $line1=$line;
	  	  $line1=~s/$line[0]/$_/;
	  	  print ZZ "$line1";
	  }	
	}
	  close ZZ;
		close CD;
		system("cat tempchange.txt >>$output_id");
   unlink "tempchange.txt";
   	open CD,"$output_cov" or die "$output_cov not exist\n";#add back coverage for redundant reads
	open ZZ,">tempchange.txt";
	foreach(<CD>){
	  my @line=split/\t/,$_;
	  my $line=$_;
	  foreach(@{$ids{$line[0]}}){
	  	  my $line1=$line;
	  	  $line1=~s/$line[0]/$_/;
	  	  print ZZ "$line1";
	  }	
	}
	  close ZZ;
		close CD;
		system("cat tempchange.txt >>$output_cov");
   unlink "tempchange.txt";
}

sub add_column_to_statistic{
    my ($identity_file)=@_;
    open HH,"$identity_file" or die "Sequence identity file $identity_file not found\n";
    my $l=<HH>;
    my @l=split/\t/,$l;
    my %identity=();
    if($l[1]!~/germ_div/){
      	last;
    }	
    else{
	while(<HH>){
	    chomp;
	    @l=split/\t/,$_;
	    $identity{$l[0]}=$l[1];	
	}
    }
    close HH;
    my @stats=<./output/tables/*all_seq_stats.txt>;
    if(-e "$stats[0]"){
	print "$stats[0]##\n";
	&rm_r($stats[0]);
	open ST,"$stats[0]";
	
	my $line=<ST>;
	if($line!~/V\_div/){
	    chomp $line;		 	
	    open STo,">temstat.txt";
	    print STo "$line\tV_div\n";
	    while(<ST>){
		chomp;
		if($_=~/^[\d\w]/){
		    my @li=split/\t/,$_;
		    if($identity{$li[0]})	{
			print STo "$_\t$identity{$li[0]}\%\n";#print "$_ $identity{$li[0]}\n";
		    }
		    else{
			print STo "$_\tNA\n";
		    }
		}
	    } 	
	    close ST;
	    close STo;
	    system("mv temstat.txt $stats[0]"); 	
	}
    } else { print "\n\nCould not find master table in output/tables/ \n\t-cannot add germline divergence to the all_seq_stats table for use with 4.1_setup_plots.pl\n\n"; }
    
}

