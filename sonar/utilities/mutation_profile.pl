#!/usr/bin/perl
use strict;

if(@ARGV<1||@ARGV%2>0){die "
	usage: script 
	-f seq_file 
	-r reference_Vgene name
	-rdb reference germline gene database
	-repair whether to find sequencing errors in the sequences, default :yes
	";}
my %para=@ARGV;
if(!$para{'-f'}){die "sequence file not exist\n";}

open HH,"$para{'-f'}" or die;#read sequence file
my %codon_to_aa=('ATT','I','ATC','I','ATA','I','CTT','L','CTC','L','CTA','L','CTG','L','TTA','L','TTG','L','GTT','V','GTC','V','GTA','V','GTG','V','TTT','F','TTC','F','ATG','M','TGT','C','TGC','C','GCT','A','GCC','A','GCA','A','GCG','A','GGT','G','GGC','G','GGA','G','GGG','G','CCT','P','CCC','P','CCA','P','CCG','P','ACT','T','ACC','T','ACA','T','ACG','T','TCT','S','TCC','S','TCA','S','TCG','S','AGT','S','AGC','S','TAT','Y','TAC','Y','TGG','W','CAA','Q','CAG','Q','AAT','N','AAC','N','CAT','H','CAC','H','GAA','E','GAG','E','GAT','D','GAC','D','AAA','K','AAG','K','CGT','R','CGC','R','CGA','R','CGG','R','AGA','R','AGG','R','TAA','*','TAG','*','TGA','*');
my $germ='IGHV1-46*01';
if($para{'-r'}){$germ=$para{'-r'}}
if(!$para{'-repair'}){$para{'-repair'}='yes';}
if(!$para{'-rdb'}){$para{'-rdb'}='/Users/sheng/work/HIV/scripts/db/germline/human_HL_short_novel_allels.fasta';}
my %mutation=();
my %mutation_total=();
my %seq_germ=();
my $id='';
my %seq=();
my %germ=();
my %total=();
my %germ_group=();
my $total=0;
my $i=1;
my %V_length=();
my %mutation_CM=();
my %germ_seqs=();
my %germ_seqsaa=();
my %germ_allelaa=();
while(<HH>){
	chomp;
    if($_=~/>([^ \t]+)/){
       $id=$i;	
       $i++;
       if($_=~/V_gene\=([^ \t]+)/){      	
       		$germ=$1; 
       		my @germs=split/\,/,$1;
       		my %mark=();
       		my $t='';
       		foreach(@germs){  
       		#push @{$germ_group{$_}},$id;
       		$t=$_;
       		$t=~s/\*.+//g;
       		#if(!$mark{$_}){$total{$_}++;$mark{$_}=1;}
       				$mark{$t}++; 
       		} 
       		my $germs_count=()=keys %mark;
       		if($germs_count>1){$id=''; }
       		else{$total++;@{$germ{$id}}=@germs; }
       				
       }
    }	
    elsif($id){
      $seq{$id}.=$_;	
    }
	
}

my $total1=$total;
#open G,"/Users/sheng/work/HIV/scripts/db/germline/human_HL_short_novel_allels.fasta";#read germline database
open G,"$para{'-rdb'}";#read germline database
while(<G>){
	chomp;
    if($_=~/>([^ \t]+)/){
       $id=$1;	
    }	
    else{
      $seq_germ{$id}.=$_;	
    }
}

foreach(sort keys %seq_germ){
	  $id=$_;
	  $id=~s/\*.+//;
	  $germ_seqs{$id}.=">$_\n$seq_germ{$_}\n";
	  my $aa=&translation($seq_germ{$_});
	  $germ_seqsaa{$id}.=">$_\n$aa\n";
	  $germ_allelaa{$_}=$aa;
}

############calculation###############
$para{'-f'}=~s/\.fa.*//;
open MP,">$para{'-f'}_mutations.txt";
my $ignored=0;
if($para{'-repair'} =~/yes/){
			open GS,">$para{'-f'}\_goodseq.fa";
		}
		print "Processing $para{'-f'}....\n";
foreach(sort {$a<=>$b} keys %seq){
	  my $id=$_;
	  my $aa_germ='';
	  my $germid='';
	  my $aa_seq=&translation($seq{$_});
	  my $aa_seq=$seq{$_};
	      $total1--;
    if($total1%100 == 0){print "$total1 sequences left...\n";}
	  foreach(@{$germ{$id}}){
      $germid=$_;
      if($seq_germ{$germid}){
      last;
      }
	  }
	  if(!$seq_germ{$germid}){warn "$id $seq{$id} germ $germid not found\n";next;}
	  
	  my ($testa,$testg)=&clustalw2($para{'-f'},">germ\n$germ_allelaa{$germid}\n",$aa_seq,'protein');
	  my $iden2=&identity($testa,$testg,'protein'); 
	  if($iden2==100){$ignored++; print "###$id\n"; next;}
	  if($para{'-repair'} =~/yes/){
	  	my @marks=&contain_frameshift($seq_germ{$germid},$seq{$id},$germ_allelaa{$germid});
	  	if($marks[0]>0){
	  		#print ES ">$id V_gene=@{$germ{$id}} $marks[0]\n$seq{$id}\n";
	  		#if($marks[2]){print ES ">$id\_repaired\n$marks[2]\n";} 
	  		
	  		next;}
	  	else{print GS ">$id V_gene=",join(',',@{$germ{$id}}),"\n$marks[2]\n";}	  
	  	if(!$marks[4]){die "no aa sequence\n";}
	    $marks[4]=~s/\-//g;
	    $aa_seq=$marks[4];
	  }
	  $germid=~s/\*.+//; 
	  $aa_germ=$germ_seqsaa{$germid};	  
	  if(!$aa_germ||!$aa_seq){warn "$_ $germid wrong\n$aa_germ\n$aa_seq\n";next;}	  
	  $total{$germid}++;

	  my @aln=&clustalw2($para{'-f'},$aa_germ,$aa_seq,'PROTEIN');
	  my $mutationp=&mutation_freq($germid,@aln);
    print MP "$id\t$germid\t$mutationp\n";
}
#foreach(sort keys %germ_group){
#    my $allele=$_;
#    my $aa_germ=&translation($seq_germ{$_});
#    if(!$aa_germ){next;}
#    foreach(@{$germ_group{$allele}}){print "$_\n";
#        my $aa_seq=&translation($seq{$_});
#        #if(!$aa_146||!$aa_germ||!$aa_seq){die "$_ wrong\n$aa_146\n$aa_germ\n$aa_seq\n"}	
#    	 #my @aln=&clustalo($aa_146,$aa_germ,$aa_seq,'PROTEIN');
#    	 if(!$aa_germ||!$aa_seq){warn "$_ $allele wrong\n$aa_germ\n$aa_seq\n";next;}
#    	 my @aln=&clustalo($aa_germ,$aa_seq,'PROTEIN');print "$aln[0]\n$aln[1]\n";
#    	 &mutation_freq(@aln,$allele);
#    	 $total1--;
#    	 print "$total1 sequences left...\n";
##   }	
#}
print "$ignored sequences ignored because of low SHM\n";
open HH,">$para{'-f'}_mutation_freq.txt";
open YY,">$para{'-f'}_mutation_profile.txt";
open ZZ,">$para{'-f'}_profile_forlogo.txt";
open XX,">$para{'-f'}_sequences.txt";
open CM,">$para{'-f'}_composition.txt";

print YY "Vgene\tpos";
print CM "Vgene\tpos";
foreach('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'){
      print YY "\t$_";	
	    print CM "\t$_";
}
print YY "\n";
print CM "\n";

foreach(sort keys %mutation){
	my $v=$_;
	print ZZ "$v\n";
	print XX "$v\t$total{$v}\n";
	foreach(1..$V_length{$v}){
	   my $pos=$_;	   
	   if(!$mutation_total{$v}{$pos}){	      
	      	print ZZ "($_) startstack\n";
	      	print ZZ "endstack\n"; 
	      	print YY "$v\t$pos";
	      	printf HH "$v\t$pos\t%.2f\n",0;
	   foreach('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'){
	   	    print YY "\t0";
	   }
	   print YY "\n";
	      	next;
	  }
	  	   
	   print ZZ "($pos) startstack\n"; 
	   printf HH "$v\t$pos\t%.2f\n",100*$mutation_total{$v}{$pos}/$total{$v};
	   my $weight=$mutation_total{$v}{$pos}/$total{$v};
	   	print YY "$v\t$pos";
	   foreach('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'){
	   	  if($mutation{$v}{$pos}{$_}){
	   	  	  printf YY "\t%.2f",100*$mutation{$v}{$pos}{$_}/$mutation_total{$v}{$pos};
	   	  }
	   	  else{
	   	  	print YY "\t0";	   	  	
	   	  }
	  }	
	 
	  foreach(sort {$mutation{$v}{$pos}{$a}<=>$mutation{$v}{$pos}{$b}} keys %{$mutation{$v}{$pos}}){
	  	   if($mutation{$v}{$pos}{$_}){	   	  	  
	   	  	  print ZZ " ($pos) ",$weight*100*$mutation{$v}{$pos}{$_}/$mutation_total{$v}{$pos}," ($_) numchar\n";
	   	    }
	  }
	  print ZZ "endstack\n";
	  print YY "\n";
	}
}

foreach(sort keys %mutation_CM){
     my $v=$_;
     foreach(1..$V_length{$v}){
     	my $pos=$_;
     	print CM "$v\t$pos";
     		foreach('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'){
     	     my $aa=$_;
     	     if($mutation_CM{$v}{$pos}{$aa}){
     	          printf CM "\t%.2f",100*$mutation_CM{$v}{$pos}{$_}/$total{$v};	
     	     }
     	     else{
     	          print CM "\t0";	
     	    }
     		}	
     		print CM "\n";     
     }	
}

unlink "$para{'-f'}\_tmp.fa","$para{'-f'}\_tmp.aln","$para{'-f'}\_tmp.dnd";
#######################
sub mutation_freq{
	  my ($V,@seq)=@_;
	  my $target=$seq[$#seq];
	  my $ref=$seq[0];
	  $ref=~s/\-//g;
	  $V=~s/\*.+//;
	  my $mutationpos=',';
	  my $leng=length($ref);
	  $V_length{$V}=$leng;
	  my $j=0;
	  my $i=0;
	  while($j<=length($seq[0])){
	  	$j++;
	  	my $mark=0;
	  	my $taa=substr($target,$j-1,1);
	  	my $gaa=substr($seq[0],$j-1,1);
	  	if($gaa eq '-'){next;}
	  	
	  	foreach(0..($#seq-1)){
	  		 $gaa=substr($seq[$_],$j-1,1);
	  		 if($gaa eq $taa ){$mark=1;last;}
	  	}
	  	$i++;
	  	$mutation_CM{$V}{$i}{$taa}++;
	  	if($mark==0 && $taa!~/[X\-\*]/ && $gaa!~/[X\-\*]/){
	  		$mutation{$V}{$i}{$taa}++;
	  		$mutation_total{$V}{$i}++;
	  		$mutationpos.="$gaa$i$taa\,";
	  	}
	  }
	return $mutationpos;
}

#######################
sub mutation_freq1{
	  my @seq=@_;
	  my $ref=$seq[0];
	  $ref=~s/\-//g;
	  my $leng=length($ref);
	  my $j=0;
	  my $V=$seq[2];
	  $V=~s/\*.+//;
	  my $i=0;
	  while($j<=length($seq[0])){
	  	$j++;
	  	my $aa=substr($seq[0],$j-1,1);
	  	if($aa eq '-'){next;}
	  	$leng--;
	  	$i++;
	  	if((substr($seq[1],$j-1,1) ne substr($seq[0],$j-1,1)) &&(substr($seq[1],$j-1,1)!~/[X\-]/)&&(substr($seq[0],$j-1,1)!~/[X\-]/)){
	  		$mutation{$V}{$i}{substr($seq[1],$j-1,1)}++;if($i ==92){print substr($seq[0],$j-1,1),"\t",substr($seq[1],$j-1,1),"\n";}
	  		$mutation_total{$V}{$i}++;
	  	}
	  }
	
}
##############
sub clustalw2{
    my($fi,$seq1,$seq2,$type)=@_;
    my %alnseq=();
    open CW,">$fi\_tmp.fa" or die "aln file can not be opened\n";
    print CW "$seq1\n>seq\n$seq2\n"; 
    close CW;
    system("clustalw2 -INFILE=$fi\_tmp.fa -OUTFILE=$fi\_tmp.aln -OUTPUT=fasta -QUIET -TYPE=$type");
    open WA,"$fi\_tmp.aln" or die "Clustalw2 aln file not found\n";
    while(<WA>){
        chomp;
        if($_=~/>(.+)/){
            $id=$1;
        }
        else{
            $alnseq{$id}.=$_;
        }

    }
    close WA;
    #unlink "$fi\_tmp.aln";
    $seq2=$alnseq{'seq'};
    delete $alnseq{'seq'};
    #print "cls\t",values %alnseq,"\n",$seq2,"\n";
    return values %alnseq,$seq2;
    
}

###################
sub clustalo{
    #my ($seq1,$seq2,$seq3,$type)=@_; 
    my ($seq1,$seq2,$type)=@_;    	
	  my %seq=();
	  if($seq2=~/\*/){print "$seq1\n$seq2\n";}
    if($seq1!~/[A-Za-z]/||$seq2!~/[A-Za-z]/){warn "$seq1\n$seq2\n no sequence for alignment\n";return '','';last;}
    my @aln=`echo "$seq1>seq\n$seq2"| clustalo -i - --wrap 999999`;
    my $id='';
   
    while(!$aln[3]){@aln=`echo "$seq1>seq\n$seq2"| clustalo -i - --wrap 999999`;}
    foreach(@aln){
        chomp;
        if($_=~/>(.+)/){
            $id=$1;
        }
        else{
            $seq{$id}.=$_;
        }
    }
    $seq2=$seq{'seq'};
    delete $seq{'seq'};
    #if($seq2&&$seq1&&(!$seq{'seq1'}||!$seq{'seq2'})){die "mafft aln error: >seq1\n$seq1>seq2\n$seq2\n";}
    return values %seq,$seq2;
}

################################
sub translation{
    my $seq=shift;
    my $aa1='';
    my @aa=();
    foreach(0){
    for(my $i=$_;$i<length($seq);$i+=3){
        if(length(substr($seq,$i,3))<3){last;}
        if($codon_to_aa{substr($seq,$i,3)}){
            $aa[$_].=$codon_to_aa{substr($seq,$i,3)};
        }
        else{
            $aa[$_].='X';
        }
 
    }
    }
    return $aa[0];
}
############################
sub contain_frameshift{
    my ($germDNA,$readDNA,$germaa)=@_;
    #remove less than two single nucleitide insertions in the read V region sequence
    my $count_insert=0;
    my $mark=0;
    ($germDNA,$readDNA)=&clustalw2($para{'-f'},">germ\n$germDNA\n",$readDNA,'DNA');
    my $start=0;#remove terminal gaps
    my $end=length($germDNA);
    if($germDNA=~/^(\-+)/){
        $start=length($1);
    }
    elsif($readDNA=~/^(\-+)/){
    	  
        $start=length($1)+(3-$1%3);
    }
    
    if($readDNA=~/(\-+)$/){$end=$end-length($1);}
    elsif($germDNA=~/(\-+)$/){
        $end=$end-length($1);
    }
    $germDNA=substr($germDNA,$start,$end-$start);
    $readDNA=substr($readDNA,$start,$end-$start);
    
    while($germDNA=~m/([A-Z]{2})(\-{1,2})([A-Z]{2})/g){#repair single insertion
        my $f=$1;
        my $b=$3;
        my $in=length($2);
        my $pos=$-[0];
            $count_insert++;
            if($count_insert>5){last;}
            if(substr($readDNA,$pos,4+$in)=~/$f([A-Z]{$in})$b/){
                substr($readDNA,$pos+2,$in,'');
                substr($germDNA,$pos+2,$in,'');
                
            }
    }
    while($readDNA=~m/([A-Z]{2})(\-{1,2})([A-Z]{2})/g){#repair single deletion
        my $f=$1;
        my $b=$3;
        my $pos=$-[0];
        my $in=length($2);
        $count_insert++;
        if($count_insert>5){last;}
        if(substr($germDNA,$pos,4+$in)=~/$f([A-Z]{$in})$b/){
            substr($readDNA,$pos+2,$in,$1);
        }
    }
    ($germDNA,$readDNA)=&clustalw2($para{'-f'},">germ\n$germDNA\n",$readDNA,'DNA');
    if($count_insert>5){$mark=1;return 1;last;}

    #check stop codon
    my $readseq=$readDNA;
    $readseq=~s/\-//g;
    #if(length($readseq)%3>0){return 5;last;}
    $readseq=&translation($readseq);
    if($readseq=~/(\*)/){
        $readseq=$readDNA;
        $readseq=~s/\-//g;
        my $shift=length($readseq)%3;
        $readseq=substr($readseq,$shift,);
        $readseq=&translation($readseq);
        if($readseq=~/(\*)/){
            return 2;last;
        }
        else{
            $readDNA=~s/\-//g;
            $readDNA=substr($readDNA,$shift,);
            ($germDNA,$readDNA)=&clustalw2($para{'-f'},">germ\n$germDNA\n",$readDNA,'DNA');
        }
    }
    
    if($germDNA=~/([A-Z]{3})\-\-([A-Z]{3})/){return 4;last;}
    #check aa/nt>0.7
    $germDNA=~s/\-//g;
    $readDNA=~s/\-//g;
    my ($aant,$refaa_aln,$targetaa_aln)=&aa_nt($germDNA,$readDNA,$germaa,$readseq);
    if($aant>0){$mark=3;}
    return $mark,$readseq,$readDNA,$refaa_aln,$targetaa_aln;
}

############################
sub aa_nt{
    my($refdna,$targetdna,$refaa,$targetaa)=@_; 
    my $window=25;
    ($refaa,$targetaa)=&clustalw2($para{'-f'},">germ\n$refaa\n",$targetaa,'protein');
    my $pos=0;
    my $i=0;
    my $prepos=0;
    my $preref=0;
    my $pretarget=0;
    my $mark=0;
    my $dna_identity=0;
    while($i<length($refaa)){
        if(substr($refaa,$i,1)=~/[A-Z]/){
            $pos++;
        }
        if($pos>=$window &&($pos%$window==0||$pos==length($refaa))){
            my $refseg=substr($refaa,$prepos,$window);
            my $targetseg=substr($targetaa,$prepos,$window);
            if($refseg!~/[A-Z]/||$targetseg!~/[A-Z]/){$i++;next;}
            my $identity_aa=&identity($refseg,$targetseg,'protein');
            $refseg=~s/\-//g;
            $targetseg=~s/\-//g;
            if(!$targetseg||!$refseg){last;}
            my $refsegdna=substr($refdna,$preref,3*length($refseg));
            my $targetsegdna=substr($targetdna,$pretarget,3*length($targetseg));
            if(!$refsegdna||!$targetsegdna){$i++;next;}
            my $identity_dna=&identity(&clustalw2($para{'-f'},">germ\n$refsegdna\n",$targetsegdna,'DNA'),'DNA');
            if($identity_dna==0&&$identity_aa>=0){$mark=1;last;}
            if($identity_aa/$identity_dna<0.7){$mark=1;last;}
            $preref+=3*length($refseg);
            $pretarget+=3*length($targetseg);
            $prepos=$i+1;
        }
        $i++;
    }
    return $mark,$refaa,$targetaa;
}
#############################
sub DNA_aln{
    my @seq=@_;
    my $seq1_DNA='';
    my $pos=0;
    #if(length($seq[0])*3 != length($seq[1])){die "DNA and protein length does not match!\n";}
    for(my $i=0;$i<length($seq[0]);$i++){
        if(substr($seq[0],$i,1) ne '-'){
            $seq1_DNA.=substr($seq[1],$pos*3,3);
            $pos++;
        }
        else{
            $seq1_DNA.='---';
        }
    }
    return $seq1_DNA;
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
    if($leng){
       $score=sprintf("%.2f",100*$score/$leng);
      return $score;
    }
    else{
      return 0;	
    }
}