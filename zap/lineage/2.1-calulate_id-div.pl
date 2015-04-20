#!/usr/bin/perl
use strict;
use threads;
# This script is used to calculate sequence identity between germline V, antibody gene and reads
if(@ARGV%2>0||@ARGV==0){die "Usage: 2Didentity.pl 
	-s readfile, if reads come from different germline V genes, the name of each read should contain a field saying it's germline V gene for somatic hypermutation calculation. Otherwise, no somatic hypermutation is calculated. fasta file output from pipeline is ok. For example, '>00000089 V_gene=IGHV1-2*02'. Please seperate '00000089 and 'V_gene' with space or tab.
	-g germline V gene file 
	-a fasta file has the sequence of interested antibody
	-t threads, default:5
	-npt number of sequences per thread. default:1000
	-p protein or DNA sequence. default: DNA
	-ap absolute path to the alignment program. muscle or clustalo or mafft. required. Based on our experience, muscle is ~2 fold faster than clustalo.
	-usearch path to usearch program to remove duplicates in the read file. Optional.
	\n";}
my %para=@ARGV;
if(!$para{'-t'}){$para{'-t'}=5;}
if(!$para{'-npt'}){$para{'-npt'}=1000;}
if(!$para{'-ap'}){die "please select a program for sequence alignment\n";}
if(!$para{'-g'}){warn "No calculation for hypermutation\n";}
if(!$para{'-p'}){$para{'-p'}='DNA';}
if($para{'-a'}&& ! -e $para{'-a'}){die "file $para{'-a'} doesn't exist. ):\n";}
if($para{'-s'}&& ! -e $para{'-s'}){die "file $para{'-s'} doesn't exist. ):\n";}
if($para{'-g'}&& ! -e $para{'-g'}){die "file $para{'-g'} doesn't exist. ):\n";}
my %germ_db=();
############Reading seqs####################
&rm_r($para{'-g'});
&rm_r($para{'-a'});
my ($germV,$germg)=&readfasta($para{'-g'});
my ($anti,$antigerm)=&readfasta($para{'-a'});

###############Processing###################
print "processing $para{'-s'}...\n";
my $changefile=$para{'-s'};
$changefile=~s/\.fa//;
open YY,">$changefile\_identity.txt";
open ZZ,">$changefile\_coverage.txt";
if($para{'-g'}){print YY "ID\tgerm_div";}
else{print YY "ID";}
foreach(sort keys %{$anti}){
    print YY "\t$_";
}
print YY "\n";
my $file_calculation=$para{'-s'};
if($para{'-usearch'}){
system("$para{'-usearch'} -derep_fulllength $para{'-s'} -threads $para{'-t'} -fastaout $changefile.unique.fa -uc $changefile.cluster -sizeout > usearchlog.txt");
$file_calculation="$changefile.unique.";
}
open READs,"$file_calculation";
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
               	 while(threads->list(threads::running)>=$para{'-t'}){
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
             my @line=split/[ \t]+/,$_;
             $id=substr($line[0],1,);
             $line[1]=~s/\,.+//g;
             $readgerm{$id}=$line[1];
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

while(threads->list()){
    foreach(threads->list(threads::joinable)){
        my @result=$_->join();
        print  YY $result[0];
        print ZZ $result[1];
    }
    sleep(1);
}

if($para{'-usearch'}){
  &recover($changefile);
  unlink "$changefile.unique.fa","$changefile.cluster";
}
##################################
sub rm_r{
      my $file=shift;
   open HH,"$file" or die "rm_r didn't find the file $file\n";
   open YY,">rmtem.txt";
   while(<HH>){
      ~s/\r/\n/g;
      print YY "$_";
   }
    close HH;
  system("mv rmtem.txt $file");	
}
###################
sub calculation{
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
###################
sub paired_identity{
    my($id,$seq,$germ,$anti)=@_;
    my @identity=();
    my @coverage=();
    push @identity,$id;
    push @coverage,$id;
    if($germ){
    	my ($seqg,$seqr)=&aln($germ,$seq,$para{'-p'});
      my $div=sprintf("%.2f",100-&identity($seqg,$seqr));
      my $cov=sprintf("%.2f",&coverage($seqg,$seqr));
      push @coverage,$cov; 
    push @identity,$div;
    }
    if($anti){
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
###################
sub readfasta{
    my $file=shift;
    my %seq=();
    my %seqgerm=();
    my $id='';
    open HH,"$file" or die "file $file not exist\n";
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
#################################
sub separte_seq{
    my ($seq,$section)=@_;
    my @keys=keys %{$seq};
    my %seq_w=();
    my $subset=sprintf("%d",@keys/$section);
    if(@keys%$section>0){$section++;}
    foreach(1..$section){
        my $start=($_-1)*$subset;
        my $end=$_*$subset-1;
        if($_*$subset>$#keys){
            $end=$#keys;
        }
        @{$seq_w{$_}}=@keys[$start..$end];
    }
    return %seq_w;
}
###################
sub aln{
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
sub mafft{
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
sub clustalo{
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
sub muscle{
    my ($seq1,$seq2,$type)=@_;
    my %seq=();
    #if($type eq 'DNA'){$type='nucleo';}
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
    return $score;
}
#################################
sub identity1{#exclude gaps
    my ($seq1,$seq2,$type)=@_;
    if(!$seq1||!$seq2||$seq2!~/[A-Z]/||$seq1!~/[A-Z]/){warn "Nothing in the sequence for identity calculation:\n$seq1#\n$seq2#\n";return 0;last;}
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
        if(substr($seq1,$j,1) ne $misc &&substr($seq2,$j,1) ne $misc&&substr($seq1,$j,1) ne '-' &&substr($seq2,$j,1) ne '-'){
            $leng++;
            if(substr($seq1,$j,1) eq substr($seq2,$j,1)){
                $score++;
            }
        }
    }
    $score=sprintf("%.2f",100*$score/$leng);
    return $score;
}
################################
sub coverage{
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
return sprintf ("%3.2f",100*$c);    
    
}
################################
sub recover{#add back all the redundant reads
   my $filename=shift;
   my %ids=();
   system("rm_r.pl $filename.cluster");
  open HH,"$filename.cluster" or die "Usearch clustering file not exist";
  
   while(<HH>){
    if($_=~/^H/){
    	chomp;
      my @line=split/\t/,$_;
      push @{$ids{$line[9]}},$line[8];
     }	
   }  	
	open CD,"$filename\_identity.txt" or die "$filename\_identity.txt not exist\n";
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
		system("cat tempchange.txt >>$filename\_identity.txt");
   unlink "tempchange.txt";
   	open CD,"$filename\_coverage.txt" or die "$filename\_identity.txt not exist\n";
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
		system("cat tempchange.txt >>$filename\_coverage.txt");
   unlink "tempchange.txt";

}