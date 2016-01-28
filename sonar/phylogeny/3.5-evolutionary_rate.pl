#!/usr/bin/env perl
use strict;
use FindBin;
use lib "$FindBin::Bin/../PPvars" qw(ppath);

################################################
#This is a mater script that allow user to generate xml configuration files and perform MCMC simulation to estimate evolutionary rate. User has to install beast package first. The script accepts fasta format sequence file as input. Be aware that only temporal data with at least from two time points, are allowed to calculate evolutionary rate. All sequences in the dataset must be non-redundant. Please remove duplicate sequences before running the script. If there are duplicate sequences in the dataset, it's better keep the one with earlier time point. Currently, only GTR substitution model was incorporated but our tests showed the rate calculation is robust to substitution models. it's users' responsibility to check if the MCMC simulation converged using Tracer in Beast2 Package.
################################################

my $usage="Usage:
Evolutionary rate calculation on an antibody lineage from at least two time points. Only nucleotide sequences accepted at this time
Options:
	-f\tinput antibody nucleotide VDJ sequence file with fasta format, must be aligned sequences
	-CDR\tsequences of the three CDRs
	-FW\tsequences of the four FWs
	-CDRb\tInstead of providing sequences for CDRs and FWs, user can also provide boundaries of the three CDR regions in the DNA sequences 
	   \tconnected with :. For example, 76:99:151:174:289:399
	-codon_pos\twhether calculate rate for the first,second and third codon positions seperately. Currently only support seperating 1+2 
	   \tfrom 3 codon position. The calculation is for full length sequence only. Default is 0 (not seperate them,otherwise give it 1).
	-spliter\twhich marker used to seperate seq name, default:_
	-nc\twhich column is time info after parsing using -spliter in the sequence name, start from 1. The time 
	   \tinfo must be numerical. default:1
	-pop\tmodel of coalecent population size, constant or bayesian_skyline, default is later
	-o\toutput file prefix
	-n\thow many sequences to randomly choose from each time point. default: 25
	-chainLength\tnumber of steps to run, default:10000000
	-storeEvery\tstore mcmc state info every n steps, default:5000
	-log\tstore a tree every n step, default:1000
	-beast\tpath to beast2 program
	-lower_limit\tlower limit of when the lineage start evolve, for example, all sequences week 11 and 50 post infection, the suspected 
	   \ttime of lineage activation is week 1, then use -lower_limit 1. Default:1. if the date includes year and month,please use this format yyyy-mm-dd.
	-upper_limit\tupper limit of when the lineage end evolve, Default:Infinity. if the date includes year and month,please use this 
	   \tformat yyyy-mm-dd.
	-burnin\t percentage of steps to be thrown away before rate estimation (burnin in beast). Default:10 

Example:
3.5-evolutionary_rate.pl -f test.fa -CDR test_CDR.fa -FW test_FW.fa -o test -beast beast -burnin 10

Created by Zizhang Sheng.

Copyright (c) 2011-2015 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
	";
foreach(@ARGV){if($_=~/[\-]{1,2}(h|help)/){die "$usage";}}
if(@ARGV%2>0||@ARGV==0){die "$usage"; }	
my %para=@ARGV;
if(!$para{'-f'}){die "No input sequence file\n";}
$para{'-beast'}=ppath('beast');
if(!$para{'-beast'}){die "Please set up absolute path to beast program\n";}
if(!$para{'-n'}){$para{'-n'}=25;}
if(!$para{'-spliter'}){$para{'-spliter'}='_';}
if(!$para{'-nc'}){$para{'-nc'}=1;}
if(!$para{'-o'}){$para{'-o'}=$para{'-f'};$para{'-o'}=~s/\..*//;}
if(!$para{'-chainLength'}){$para{'-chainLength'}=10000000;}
if(!$para{'-pop'}){$para{'-pop'}='bayesian_skyline';}
if(!$para{'-log'}){$para{'-log'}=1000;}
if(!$para{'-storeEvery'}){$para{'-storeEvery'}=5000;}
if(!$para{'-burnin'}){$para{'-burnin'}=10;}
if(!$para{'-lower_limit'}){$para{'-lower_limit'}=1;}
if(!$para{'-upper_limit'}){$para{'-upper_limit'}='Infinity';}
if($para{'-lower_limit'}=~/[0-9]{4}\-[0-9]{2}\-[0-9]{2}/){$para{'-lower_limit'}=&year_to_digital($para{'-lower_limit'});}
if($para{'-upper_limit'}=~/[0-9]{4}\-[0-9]{2}\-[0-9]{2}/){$para{'-upper_limit'}=&year_to_digital($para{'-upper_limit'});}
if(!$para{'-codon_pos'}){$para{'-codon_pos'}=0;}


######Processing######
my $seq_pos12='';
my $seq_pos3='';
my $result='';
my %rate_partitions=();
my @partition_names=();
my ($seq,$seq_name_selected)=&read_fasta($para{'-f'},$para{'-n'});#read sequences and select n sequences for each time point
if($para{'-CDR'}){$rate_partitions{'CDR'}=&read_fasta($para{'-CDR'});}
if($para{'-FW'}){$rate_partitions{'FW'}=&read_fasta($para{'-FW'});}

if(!$para{'-CDR'}&&!$para{'-FW'}&&$para{'-CDRb'}){#find CDR seqences
	  my @regions=split/\:/,$para{'-CDRb'};
	  if(@regions!=6&&@regions!=4){die "Provided CDR boundaries are not correct\n";}
	  ($rate_partitions{'CDR'},$rate_partitions{'FW'})=&extract_CDR(@regions);
	}
	
if($para{'-codon_pos'}){($rate_partitions{'pos12'},$rate_partitions{'pos3'})=&extract_codon_pos();}#find sequences for codon positions

&write_xml();
system("$para{'-beast'} $para{'-o'}.xml");
&parsing_results("$para{'-o'}.log");
if(-d "./output/rate/"){
  system("mv $para{'-o'}.xml $para{'-o'}.xml.state $para{'-o'}.log $para{'-o'}.trees ./output/rate/");	
}

print $result;
############################
sub extract_codon_pos{
	  my %pos12=();
	  my %pos3=();
	  foreach(@$seq_name_selected){
	  	my $name=$_;
	  	chomp $seq->{$name};
	  	for(my $i=0;$i<length($seq->{$name});$i+=3){
	  		$pos3{$name}.=substr($seq->{$name},$i+2,1);
	  		$pos12{$name}.=substr($seq->{$name},$i,2);
	  		
	  	}
	  }
	return \%pos12,\%pos3;
}
############################
sub extract_CDR{
    my (@region)=@_;
    my %seqCDR=();
    my %seqFW=();
    foreach(@$seq_name_selected){
        $seqCDR{$_}=substr($seq->{$_},$region[0]-1,$region[1]-$region[0]+1).substr($seq->{$_},$region[2]-1,$region[3]-$region[2]+1).substr($seq->{$_},$region[4]-1,$region[5]-$region[4]+1);
        $seqFW{$_}=substr($seq->{$_},0,$region[0]-1).substr($seq->{$_},$region[1],$region[2]-1-$region[1]).substr($seq->{$_},$region[3],$region[4]-1-$region[3]).substr($seq->{$_},$region[5],);	
    }	
	 return \%seqCDR,\%seqFW;
}
############################
sub parsing_results{#calculate mean rate and 95% highest probability density range (or confidence interval)
    my ($file,$codon)=@_;
    open HH,"$file" or die "resulte file $file not found\n";#read in log file
    my %rate=();
    my %rate_colomn=();
	  while(<HH>){
	  	if($_=~/^Sample/){
	  		 my @line=split/	/,$_;
	  		 my $i=0;
	  		 foreach(@line){
	  		   if($_=~/rate.mean/){#find colomn for mean rate 
	  		       	$rate_colomn{'overall'}=$i;
	  		    }
	  		   elsif($_=~/rate\.\_(.+)\.mean/){#find colomn for framework region rate 
	  		       	$rate_colomn{$1}=$i;		
	  		    }
	  		  $i++;
	  		}
	  	}
	  elsif($_=~/^[0-9]/){
	     my @line1=split/\t/,$_;
	     foreach(keys %rate_colomn){
	       push @{$rate{$_}},$line1[$rate_colomn{$_}];	
	    }
	   }
	  	
	  }
	  
	  #output rates
	  foreach(sort keys %rate){
	  	my @rate=@{$rate{$_}};
	    for(my $j=0;$j<@rate/$para{'-burnin'};$j++){
	     shift @rate;	
	    }
	    my @rate_s=&statistic(@rate);
	    $result.="Mean rate $_: $rate_s[0]\n95% HPD interval: [$rate_s[3],$rate_s[4]]\n\n";
	  }

}
#############################
sub year_to_digital{
    my ($date)=@_;
    my %month=('01',31,'02',28,'03',31,'04',30,'05',31,'06',30,'07',31,'08',31,'09',30,'10',31,'11',30,'12',31);	
	  my @line=split/\-/,$date;
	  my $dates=0;
	  foreach(sort {$a<=>$b} keys %month){
	      if($_<$line[1]){
	      	 $dates+=$month{$_};
	      }	
	  }
	  $dates+=$line[2];
	  $dates=$dates/365;
	  if($dates>=1){$dates=0.99999;}
	  return $line[0].substr($dates,1,);
}
#########################
sub statistic{# estimate mean, median and standard deviations of rate
    my @array=@_;
    @array=sort {$a <=>$b} @array;
    my $mean=0;
    my $sum=0;
    my $median=0;
    my $stddev=0;
    my $variance=0;
    my $stderr=0;
    foreach(@array){
    	$sum+=$_;
    }	
    $mean=$sum/@array;
    $median=$array[int(@array/2)];
    foreach(@array){
      $variance+=($_-$mean)**2;	
    }
    if(@array/1==1){$stddev=0;$variance=0;}
    else{
      $variance=$variance/(@array-1);
      $stddev=sqrt($variance);
      $stderr=$stddev/sqrt(@array/1);
    }
    return $median,$mean,$stderr,&HPD(@array);
}
#########################
sub HPD{#find 95% high probability density region
    my @array=@_;	
	  my $interval=9999999;
	  my $low=0;
	  my $high=0;
	  if(@array<2){}
	  else{
	   for(my $i=0;($i+1)/@array<=0.05;$i++){
	  	my $up=sprintf("%d",@array*0.95)+$i;
	  	if($up==@array){$up--;}
	  	if($up>@array){print "something wrong\n";last;}	  	    
	  	 if($array[$up]-$array[$i]<$interval){
	  	    $low=$array[$i];
	  	    $high=$array[$up];
	  	    $interval=$high-$low;
	  	}
	   }
	  }
	return $low,$high;
}
############################
sub write_xml{#write cofiguration file for beast
my $type='';
open YY,">$para{'-o'}.xml";
my %distribution=();
my $key='';
print YY '<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate=\'Standard\' beautistatus=\'\' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" version="2.0">',"\n\n\n    <data\n";#print title
print YY "id=\"$para{'-o'}\"\nname=\"alignment\">\n";
my $taxas=@$seq_name_selected;
foreach(@$seq_name_selected){
	print YY "                        <sequence id=\"seq_$_\" taxon=\"$_\" totalcount=\"4\" value=\"$seq->{$_}\"/>\n";	
	}
print YY "                    </data>\n";
my $seqtaxon=1;
foreach(sort keys %rate_partitions){
	$key=$_;
	print YY "<data\nid=\"$para{'-o'}\_$key\">\n";
	foreach(@$seq_name_selected){
		if(!$rate_partitions{$key}->{$_}){die "Please make sure you have the same sequencs in full length and $key files\n";}
	print YY "                        <sequence id=\"seq_$_",$seqtaxon,"\" taxon=\"$_\" totalcount=\"4\" value=\"$rate_partitions{$key}->{$_}\"/>\n";	
	}
	print YY "                    </data>\n";
	$seqtaxon++;
	}


print YY '
<map name="Beta">beast.math.distributions.Beta</map>
<map name="Exponential">beast.math.distributions.Exponential</map>
<map name="InverseGamma">beast.math.distributions.InverseGamma</map>
<map name="LogNormal">beast.math.distributions.LogNormalDistributionModel</map>
<map name="Gamma">beast.math.distributions.Gamma</map>
<map name="Uniform">beast.math.distributions.Uniform</map>
<map name="prior">beast.math.distributions.Prior</map>
<map name="LaplaceDistribution">beast.math.distributions.LaplaceDistribution</map>
<map name="OneOnX">beast.math.distributions.OneOnX</map>
<map name="Normal">beast.math.distributions.Normal</map>';

print YY "

<run chainLength=\"$para{'-chainLength'}\" id=\"mcmc\" spec=\"MCMC\">
    <state id=\"state\" storeEvery=\"$para{'-storeEvery'}\">
        <tree id=\"Tree.t:$para{'-o'}\" name=\"stateNode\">
            <trait id=\"dateTrait.t:$para{'-o'}\" spec=\"beast.evolution.tree.TraitSet\" traitname=\"date\">
                ";
my $taxon='';
foreach(@$seq_name_selected){
	my @line=split/[$para{'-spliter'}]/,$_;
	$line[$para{'-nc'}-1]=~s/^0+//;
	$taxon.="$_=".&year_to_digital($line[$para{'-nc'}-1]).",\n"
	}
	chomp $taxon;
	chop $taxon;
	print YY $taxon;
print YY "                <taxa id=\"TaxonSet.$para{'-o'}\" spec=\"TaxonSet\">
                    <data
idref=\"$para{'-o'}\"
name=\"alignment\"/>
                </taxa>
            </trait>
            <taxonset idref=\"TaxonSet.$para{'-o'}\"/>
        </tree>
        <parameter dimension=\"4\" id=\"freqParameter.s:$para{'-o'}\" lower=\"0.0\" name=\"stateNode\" upper=\"1.0\">0.25</parameter>
        <parameter id=\"rateAC.s:$para{'-o'}\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>
        <parameter id=\"rateAG.s:$para{'-o'}\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>
        <parameter id=\"rateAT.s:$para{'-o'}\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>
        <parameter id=\"rateCG.s:$para{'-o'}\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>
        <parameter id=\"rateGT.s:$para{'-o'}\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>
        <parameter id=\"gammaShape.s:$para{'-o'}\" name=\"stateNode\">1.0</parameter>
        <parameter id=\"ucldMean.c:$para{'-o'}\" name=\"stateNode\">1.0</parameter>
        <parameter id=\"ucldStdev.c:$para{'-o'}\" lower=\"0.0\" name=\"stateNode\" upper=\"5.0\">0.5</parameter>
        <stateNode dimension=\"",2*$taxas-2,"\" id=\"rateCategories.c:$para{'-o'}\" spec=\"parameter.IntegerParameter\">1</stateNode>
";	
foreach(sort keys %rate_partitions){
	$key=$_;
print YY "        <parameter id=\"ucldMean.c:$para{'-o'}_$key\" name=\"stateNode\">1.0</parameter>
        <parameter id=\"ucldStdev.c:$para{'-o'}_$key\" lower=\"0.0\" name=\"stateNode\" upper=\"5.0\">0.5</parameter>
        <stateNode dimension=\"$taxas\" id=\"rateCategories.c:$para{'-o'}_$key\" spec=\"parameter.IntegerParameter\">1</stateNode>\n";

}#rate setting up

if($para{'-pop'} =~/constant/){
	  print YY "<parameter id=\"popSize.t:$para{'-o'}\" name=\"stateNode\">0.3</parameter>
    </state>

    <init estimate=\"false\" id=\"RandomTree.t:$para{'-o'}\" initial=\"\@Tree.t:$para{'-o'}\" spec=\"beast.evolution.tree.RandomTree\" taxa=\"\@$para{'-o'}\">
        <populationModel id=\"ConstantPopulation0.t:$para{'-o'}\" spec=\"ConstantPopulation\">
            <parameter id=\"randomPopSize.t:$para{'-o'}\" name=\"popSize\">1.0</parameter>
        </populationModel>
    </init>

    <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">
        <distribution id=\"prior\" spec=\"util.CompoundDistribution\">
            <distribution id=\"CoalescentConstant.t:$para{'-o'}\" spec=\"Coalescent\">
                <populationModel id=\"ConstantPopulation.t:$para{'-o'}\" popSize=\"\@popSize.t:$para{'-o'}\" spec=\"ConstantPopulation\"/>
                <treeIntervals id=\"TreeIntervals.t:$para{'-o'}\" spec=\"TreeIntervals\" tree=\"\@Tree.t:$para{'-o'}\"/>
            </distribution>
            <prior id=\"GammaShapePrior.s:$para{'-o'}\" name=\"distribution\" x=\"\@gammaShape.s:$para{'-o'}\">
                <Exponential id=\"Exponential.0\" name=\"distr\">
                    <parameter id=\"RealParameter.0\" lower=\"0.0\" name=\"mean\" upper=\"0.0\">1.0</parameter>
                </Exponential>
            </prior>
            <prior id=\"PopSizePrior.t:$para{'-o'}\" name=\"distribution\" x=\"\@popSize.t:$para{'-o'}\">
                <OneOnX id=\"OneOnX.0\" name=\"distr\"/>
            </prior>\n";

}
else{
	print YY "        <parameter dimension=\"5\" id=\"bPopSizes.t:$para{'-o'}\" lower=\"0.0\" name=\"stateNode\" upper=\"380000.0\">380.0</parameter>
        <stateNode dimension=\"5\" id=\"bGroupSizes.t:$para{'-o'}\" spec=\"parameter.IntegerParameter\">1</stateNode>
    </state>

    <init estimate=\"false\" id=\"RandomTree.t:$para{'-o'}\" initial=\"\@Tree.t:$para{'-o'}\" spec=\"beast.evolution.tree.RandomTree\" taxa=\"\@$para{'-o'}\">
        <populationModel id=\"ConstantPopulation0.t:$para{'-o'}\" spec=\"ConstantPopulation\">
            <parameter id=\"randomPopSize.t:$para{'-o'}\" name=\"popSize\">1.0</parameter>
        </populationModel>
    </init>

    <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">
        <distribution id=\"prior\" spec=\"util.CompoundDistribution\">
            <distribution groupSizes=\"\@bGroupSizes.t:$para{'-o'}\" id=\"BayesianSkyline.t:$para{'-o'}\" popSizes=\"\@bPopSizes.t:$para{'-o'}\" spec=\"BayesianSkyline\">
                <treeIntervals id=\"BSPTreeIntervals.t:$para{'-o'}\" spec=\"TreeIntervals\" tree=\"\@Tree.t:$para{'-o'}\"/>
            </distribution>
            <distribution id=\"MarkovChainedPopSizes.t:$para{'-o'}\" jeffreys=\"true\" parameter=\"\@bPopSizes.t:$para{'-o'}\" spec=\"beast.math.distributions.MarkovChainDistribution\"/>
            <prior id=\"GammaShapePrior.s:$para{'-o'}\" name=\"distribution\" x=\"\@gammaShape.s:$para{'-o'}\">
                <Exponential id=\"Exponential.0\" name=\"distr\">
                    <parameter id=\"RealParameter.0\" lower=\"0.0\" name=\"mean\" upper=\"0.0\">1.0</parameter>
                </Exponential>
            </prior>\n";
	
}

print YY "            <prior id=\"RateACPrior.s:$para{'-o'}\" name=\"distribution\" x=\"\@rateAC.s:$para{'-o'}\">
                <Gamma id=\"Gamma.0\" name=\"distr\">
                    <parameter id=\"RealParameter.01\" lower=\"0.0\" name=\"alpha\" upper=\"0.0\">0.05</parameter>
                    <parameter id=\"RealParameter.02\" lower=\"0.0\" name=\"beta\" upper=\"0.0\">10.0</parameter>
                </Gamma>
            </prior>
            <prior id=\"RateAGPrior.s:$para{'-o'}\" name=\"distribution\" x=\"\@rateAG.s:$para{'-o'}\">
                <Gamma id=\"Gamma.01\" name=\"distr\">
                    <parameter id=\"RealParameter.03\" lower=\"0.0\" name=\"alpha\" upper=\"0.0\">0.05</parameter>
                    <parameter id=\"RealParameter.04\" lower=\"0.0\" name=\"beta\" upper=\"0.0\">20.0</parameter>
                </Gamma>
            </prior>
            <prior id=\"RateATPrior.s:$para{'-o'}\" name=\"distribution\" x=\"\@rateAT.s:$para{'-o'}\">
                <Gamma id=\"Gamma.02\" name=\"distr\">
                    <parameter id=\"RealParameter.05\" lower=\"0.0\" name=\"alpha\" upper=\"0.0\">0.05</parameter>
                    <parameter id=\"RealParameter.06\" lower=\"0.0\" name=\"beta\" upper=\"0.0\">10.0</parameter>
                </Gamma>
            </prior>
            <prior id=\"RateCGPrior.s:$para{'-o'}\" name=\"distribution\" x=\"\@rateCG.s:$para{'-o'}\">
                <Gamma id=\"Gamma.03\" name=\"distr\">
                    <parameter id=\"RealParameter.07\" lower=\"0.0\" name=\"alpha\" upper=\"0.0\">0.05</parameter>
                    <parameter id=\"RealParameter.08\" lower=\"0.0\" name=\"beta\" upper=\"0.0\">10.0</parameter>
                </Gamma>
            </prior>
            <prior id=\"RateGTPrior.s:$para{'-o'}\" name=\"distribution\" x=\"\@rateGT.s:$para{'-o'}\">
                <Gamma id=\"Gamma.04\" name=\"distr\">
                    <parameter id=\"RealParameter.09\" lower=\"0.0\" name=\"alpha\" upper=\"0.0\">0.05</parameter>
                    <parameter id=\"RealParameter.010\" lower=\"0.0\" name=\"beta\" upper=\"0.0\">10.0</parameter>
                </Gamma>
            </prior>
            <prior id=\"MeanRatePrior.c:$para{'-o'}\" name=\"distribution\" x=\"\@ucldMean.c:$para{'-o'}\">
                <Uniform id=\"Uniform.0\" name=\"distr\" upper=\"Infinity\"/>
            </prior>\n";
            $distribution{'uniform'}=1;
            foreach(sort keys %rate_partitions){
	             $key=$_;	             
               print YY "            <prior id=\"MeanRatePrior.c:$para{'-o'}\_$key\" name=\"distribution\" x=\"\@ucldMean.c:$para{'-o'}\_$key\">
                <Uniform id=\"Uniform.0$distribution{'uniform'}\" name=\"distr\" upper=\"Infinity\"/>
            </prior>\n";
            $distribution{'uniform'}++;
            }
            
            
     
print YY "            <prior id=\"ucldStdevPrior.c:$para{'-o'}\" name=\"distribution\" x=\"\@ucldStdev.c:$para{'-o'}\">
                <Exponential id=\"Exponential.01\" name=\"distr\">
                    <parameter estimate=\"false\" id=\"RealParameter.011\" name=\"mean\">0.3333</parameter>
                </Exponential>
            </prior>\n";
$distribution{'exp'}=2; 
$distribution{'realp'}=12;            
           foreach(sort keys %rate_partitions){
	             $key=$_;	             
	             print YY "            <prior id=\"ucldStdevPrior.c:$para{'-o'}\_$key\" name=\"distribution\" x=\"\@ucldStdev.c:$para{'-o'}\_$key\">
                <Exponential id=\"Exponential.0$distribution{'exp'}\" name=\"distr\">
                    <parameter estimate=\"false\" id=\"RealParameter.0$distribution{'realp'}\" name=\"mean\">0.3333</parameter>
                </Exponential>
            </prior>\n";
            $distribution{'exp'}++;
            $distribution{'realp'}++;
	         }
     
print YY "            <distribution id=\"all.prior\" monophyletic=\"true\" spec=\"beast.math.distributions.MRCAPrior\" tree=\"\@Tree.t:$para{'-o'}\">
                <taxonset id=\"all\" spec=\"TaxonSet\">\n";
foreach(@$seq_name_selected){
   print YY "                    <taxon id=\"$_\" spec=\"Taxon\"/>\n";	
}

print YY "               </taxonset>
                <Uniform id=\"Uniform.0$distribution{'uniform'}\" lower=\"$para{'-lower_limit'}\" name=\"distr\" upper=\"$para{'-upper_limit'}\"/>
            </distribution>
        </distribution>
        <distribution id=\"likelihood\" spec=\"util.CompoundDistribution\">
            <distribution data=\"\@$para{'-o'}\" id=\"treeLikelihood.$para{'-o'}\" spec=\"TreeLikelihood\" tree=\"\@Tree.t:$para{'-o'}\">
                <siteModel gammaCategoryCount=\"4\" id=\"SiteModel.s:$para{'-o'}\" shape=\"\@gammaShape.s:$para{'-o'}\" spec=\"SiteModel\">
                    <parameter estimate=\"false\" id=\"mutationRate.s:$para{'-o'}\" name=\"mutationRate\">1.0</parameter>
                    <parameter estimate=\"false\" id=\"proportionInvariant.s:$para{'-o'}\" lower=\"0.0\" name=\"proportionInvariant\" upper=\"1.0\">0.0</parameter>
                    <substModel id=\"gtr.s:$para{'-o'}\" rateAC=\"\@rateAC.s:$para{'-o'}\" rateAG=\"\@rateAG.s:$para{'-o'}\" rateAT=\"\@rateAT.s:$para{'-o'}\" rateCG=\"\@rateCG.s:$para{'-o'}\" rateGT=\"\@rateGT.s:$para{'-o'}\" spec=\"GTR\">
                        <parameter estimate=\"false\" id=\"rateCT.s:$para{'-o'}\" lower=\"0.0\" name=\"rateCT\">1.0</parameter>
                        <frequencies frequencies=\"\@freqParameter.s:$para{'-o'}\" id=\"estimatedFreqs.s:$para{'-o'}\" spec=\"Frequencies\"/>
                    </substModel>
                </siteModel>
                <branchRateModel clock.rate=\"\@ucldMean.c:$para{'-o'}\" id=\"RelaxedClock.c:$para{'-o'}\" rateCategories=\"\@rateCategories.c:$para{'-o'}\" spec=\"beast.evolution.branchratemodel.UCRelaxedClockModel\" tree=\"\@Tree.t:$para{'-o'}\">
                    <LogNormal S=\"\@ucldStdev.c:$para{'-o'}\" id=\"LogNormalDistributionModel.c:$para{'-o'}\" meanInRealSpace=\"true\" name=\"distr\">
                        <parameter estimate=\"false\" id=\"RealParameter.0$distribution{'realp'}\" lower=\"0.0\" name=\"M\" upper=\"1.0\">1.0</parameter>
                    </LogNormal>
                </branchRateModel>
            </distribution>\n";
            $distribution{'realp'}++;
            $distribution{'uniform'}++;
            
foreach(sort keys %rate_partitions){
	             $key=$_;	  
   print YY "<distribution id=\"treeLikelihood.$para{'-o'}\_$key\" siteModel=\"\@SiteModel.s:$para{'-o'}\" spec=\"TreeLikelihood\" tree=\"\@Tree.t:$para{'-o'}\">
                <data
idref=\"$para{'-o'}\_$key\"/>
                <branchRateModel clock.rate=\"\@ucldMean.c:$para{'-o'}\_$key\" id=\"RelaxedClock.c:$para{'-o'}\_$key\" rateCategories=\"\@rateCategories.c:$para{'-o'}\_$key\" spec=\"beast.evolution.branchratemodel.UCRelaxedClockModel\" tree=\"\@Tree.t:$para{'-o'}\">
                    <LogNormal S=\"\@ucldStdev.c:$para{'-o'}\_$key\" id=\"LogNormalDistributionModel.c:$para{'-o'}\_$key\" meanInRealSpace=\"true\" name=\"distr\">
                        <parameter estimate=\"false\" id=\"RealParameter.0$distribution{'realp'}\" lower=\"0.0\" name=\"M\" upper=\"1.0\">1.0</parameter>
                    </LogNormal>
                </branchRateModel>
            </distribution>\n";
            $distribution{'realp'}++;
            }
           
print YY "
        </distribution>
    </distribution>
";
print YY "
    <operator id=\"treeScaler.t:$para{'-o'}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" tree=\"\@Tree.t:$para{'-o'}\" weight=\"3.0\"/>

    <operator id=\"treeRootScaler.t:$para{'-o'}\" rootOnly=\"true\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" tree=\"\@Tree.t:$para{'-o'}\" weight=\"3.0\"/>

    <operator id=\"UniformOperator.t:$para{'-o'}\" spec=\"Uniform\" tree=\"\@Tree.t:$para{'-o'}\" weight=\"30.0\"/>

    <operator id=\"SubtreeSlide.t:$para{'-o'}\" spec=\"SubtreeSlide\" tree=\"\@Tree.t:$para{'-o'}\" weight=\"15.0\"/>

    <operator id=\"narrow.t:$para{'-o'}\" spec=\"Exchange\" tree=\"\@Tree.t:$para{'-o'}\" weight=\"15.0\"/>

    <operator id=\"wide.t:$para{'-o'}\" isNarrow=\"false\" spec=\"Exchange\" tree=\"\@Tree.t:$para{'-o'}\" weight=\"3.0\"/>

    <operator id=\"WilsonBalding.t:$para{'-o'}\" spec=\"WilsonBalding\" tree=\"\@Tree.t:$para{'-o'}\" weight=\"3.0\"/>

    <operator delta=\"0.01\" id=\"FrequenciesExchanger.s:$para{'-o'}\" spec=\"DeltaExchangeOperator\" weight=\"0.1\">
        <parameter idref=\"freqParameter.s:$para{'-o'}\"/>
    </operator>

    <operator id=\"RateACScaler.s:$para{'-o'}\" parameter=\"\@rateAC.s:$para{'-o'}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>

    <operator id=\"RateAGScaler.s:$para{'-o'}\" parameter=\"\@rateAG.s:$para{'-o'}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>

    <operator id=\"RateATScaler.s:$para{'-o'}\" parameter=\"\@rateAT.s:$para{'-o'}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>

    <operator id=\"RateCGScaler.s:$para{'-o'}\" parameter=\"\@rateCG.s:$para{'-o'}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>

    <operator id=\"RateGTScaler.s:$para{'-o'}\" parameter=\"\@rateGT.s:$para{'-o'}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>

    <operator id=\"gammaShapeScaler.s:$para{'-o'}\" parameter=\"\@gammaShape.s:$para{'-o'}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>

    <operator id=\"ucldMeanScaler.c:$para{'-o'}\" parameter=\"\@ucldMean.c:$para{'-o'}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"1.0\"/>

    <operator id=\"ucldStdevScaler.c:$para{'-o'}\" parameter=\"\@ucldStdev.c:$para{'-o'}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"3.0\"/>

    <operator id=\"CategoriesRandomWalk.c:$para{'-o'}\" parameter=\"\@rateCategories.c:$para{'-o'}\" spec=\"IntRandomWalkOperator\" weight=\"10.0\" windowSize=\"1\"/>

    <operator id=\"CategoriesSwapOperator.c:$para{'-o'}\" intparameter=\"\@rateCategories.c:$para{'-o'}\" spec=\"SwapOperator\" weight=\"10.0\"/>

    <operator id=\"CategoriesUniform.c:$para{'-o'}\" parameter=\"\@rateCategories.c:$para{'-o'}\" spec=\"UniformOperator\" weight=\"10.0\"/>

    <operator id=\"relaxedUpDownOperator.c:$para{'-o'}\" scaleFactor=\"0.75\" spec=\"UpDownOperator\" weight=\"3.0\">
        <parameter idref=\"ucldMean.c:$para{'-o'}\" name=\"up\"/>
        <tree idref=\"Tree.t:$para{'-o'}\" name=\"down\"/>
    </operator>\n";
    
foreach(sort keys %rate_partitions){
	             $key=$_;	
print YY"    <operator id=\"ucldMeanScaler.c:$para{'-o'}\_$key\" parameter=\"\@ucldMean.c:$para{'-o'}\_$key\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"1.0\"/>

    <operator id=\"ucldStdevScaler.c:$para{'-o'}\_$key\" parameter=\"\@ucldStdev.c:$para{'-o'}\_$key\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"3.0\"/>

    <operator id=\"CategoriesRandomWalk.c:$para{'-o'}\_$key\" parameter=\"\@rateCategories.c:$para{'-o'}\_$key\" spec=\"IntRandomWalkOperator\" weight=\"10.0\" windowSize=\"1\"/>

    <operator id=\"CategoriesSwapOperator.c:$para{'-o'}\_$key\" intparameter=\"\@rateCategories.c:$para{'-o'}\_$key\" spec=\"SwapOperator\" weight=\"10.0\"/>

    <operator id=\"CategoriesUniform.c:$para{'-o'}\_$key\" parameter=\"\@rateCategories.c:$para{'-o'}\_$key\" spec=\"UniformOperator\" weight=\"10.0\"/>

    <operator id=\"relaxedUpDownOperator.c:$para{'-o'}\_$key\" scaleFactor=\"0.75\" spec=\"UpDownOperator\" weight=\"3.0\">
        <parameter idref=\"ucldMean.c:$para{'-o'}\_$key\" name=\"up\"/>
        <tree idref=\"Tree.t:$para{'-o'}\" name=\"down\"/>
    </operator>\n";
}


if($para{'-pop'} =~/constant/){
  print YY "    <operator id=\"PopSizeScaler.t:$para{'-o'}\" parameter=\"\@popSize.t:$para{'-o'}\" scaleFactor=\"0.75\" spec=\"ScaleOperator\" weight=\"3.0\"/>

    <logger fileName=\"$para{'-o'}.log\" id=\"tracelog\" logEvery=\"$para{'-log'}\" model=\"\@posterior\" sanitiseHeaders=\"true\" sort=\"smart\">
        <log idref=\"posterior\"/>
        <log idref=\"likelihood\"/>
        <log idref=\"prior\"/>
        <log idref=\"treeLikelihood.$para{'-o'}\"/>
        <log id=\"TreeHeight.t:$para{'-o'}\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"\@Tree.t:$para{'-o'}\"/>
        <parameter idref=\"freqParameter.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"rateAC.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"rateAG.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"rateAT.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"rateCG.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"rateGT.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"gammaShape.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"ucldMean.c:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"ucldStdev.c:$para{'-o'}\" name=\"log\"/>
        <log branchratemodel=\"\@RelaxedClock.c:$para{'-o'}\" id=\"rate.c:$para{'-o'}\" spec=\"beast.evolution.branchratemodel.RateStatistic\" tree=\"\@Tree.t:$para{'-o'}\"/>\n";
foreach(sort keys %rate_partitions){
	             $key=$_;	
    print YY "        <parameter idref=\"ucldMean.c:$para{'-o'}\_$key\" name=\"log\"/>
        <parameter idref=\"ucldStdev.c:$para{'-o'}\_$key\" name=\"log\"/>
        <log branchratemodel=\"\@RelaxedClock.c:$para{'-o'}\_$key\" id=\"rate.c:$para{'-o'}\_$key\" spec=\"beast.evolution.branchratemodel.RateStatistic\" tree=\"\@Tree.t:$para{'-o'}\"/>\n";
}

print YY "
        <log idref=\"all.prior\"/>
        <parameter idref=\"popSize.t:$para{'-o'}\" name=\"log\"/>
        <log idref=\"CoalescentConstant.t:$para{'-o'}\"/>
    </logger>

    <logger id=\"screenlog\" logEvery=\"$para{'-log'}\">
        <log idref=\"posterior\"/>
        <log arg=\"\@posterior\" id=\"ESS.0\" spec=\"util.ESS\"/>
        <log idref=\"likelihood\"/>
        <log idref=\"prior\"/>
    </logger>

    <logger fileName=\"$para{'-o'}.trees\" id=\"treelog.t:$para{'-o'}\" logEvery=\"$para{'-log'}\" mode=\"tree\">
        <log branchratemodel=\"\@RelaxedClock.c:$para{'-o'}\" id=\"TreeWithMetaDataLogger.t:$para{'-o'}\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"\@Tree.t:$para{'-o'}\"/>
    </logger>

</run>

</beast>
";

}
else{
print YY "
    <operator id=\"popSizesScaler.t:$para{'-o'}\" parameter=\"\@bPopSizes.t:$para{'-o'}\" scaleFactor=\"0.75\" spec=\"ScaleOperator\" weight=\"15.0\"/>

    <operator id=\"groupSizesDelta.t:$para{'-o'}\" integer=\"true\" spec=\"DeltaExchangeOperator\" weight=\"6.0\">
        <intparameter idref=\"bGroupSizes.t:$para{'-o'}\"/>
    </operator>

    <logger fileName=\"$para{'-o'}.log\" id=\"tracelog\" logEvery=\"$para{'-log'}\" model=\"\@posterior\" sanitiseHeaders=\"true\" sort=\"smart\">
        <log idref=\"posterior\"/>
        <log idref=\"likelihood\"/>
        <log idref=\"prior\"/>
        <log idref=\"treeLikelihood.$para{'-o'}\"/>
        <log id=\"TreeHeight.t:$para{'-o'}\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"\@Tree.t:$para{'-o'}\"/>
        <parameter idref=\"freqParameter.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"rateAC.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"rateAG.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"rateAT.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"rateCG.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"rateGT.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"gammaShape.s:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"ucldMean.c:$para{'-o'}\" name=\"log\"/>
        <parameter idref=\"ucldStdev.c:$para{'-o'}\" name=\"log\"/>
        <log branchratemodel=\"\@RelaxedClock.c:$para{'-o'}\" id=\"rate.c:$para{'-o'}\" spec=\"beast.evolution.branchratemodel.RateStatistic\" tree=\"\@Tree.t:$para{'-o'}\"/>\n";
foreach(sort keys %rate_partitions){
	             $key=$_;	
    print YY "        <parameter idref=\"ucldMean.c:$para{'-o'}\_$key\" name=\"log\"/>
        <parameter idref=\"ucldStdev.c:$para{'-o'}\_$key\" name=\"log\"/>
        <log branchratemodel=\"\@RelaxedClock.c:$para{'-o'}\_$key\" id=\"rate.c:$para{'-o'}\_$key\" spec=\"beast.evolution.branchratemodel.RateStatistic\" tree=\"\@Tree.t:$para{'-o'}\"/>\n";
}     

print YY "
        <log idref=\"BayesianSkyline.t:$para{'-o'}\"/>
        <parameter idref=\"bPopSizes.t:$para{'-o'}\" name=\"log\"/>
        <log idref=\"bGroupSizes.t:$para{'-o'}\"/>
        <log idref=\"all.prior\"/>
    </logger>

    <logger id=\"screenlog\" logEvery=\"$para{'-log'}\">
        <log idref=\"posterior\"/>
        <log arg=\"\@posterior\" id=\"ESS.0\" spec=\"util.ESS\"/>
        <log idref=\"likelihood\"/>
        <log idref=\"prior\"/>
    </logger>

    <logger fileName=\"$para{'-o'}.trees\" id=\"treelog.t:$para{'-o'}\" logEvery=\"$para{'-log'}\" mode=\"tree\">
        <log branchratemodel=\"\@RelaxedClock.c:$para{'-o'}\" id=\"TreeWithMetaDataLogger.t:$para{'-o'}\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"\@Tree.t:$para{'-o'}\"/>
    </logger>

</run>

</beast>";
}
}
######################

######################
sub read_fasta{ #read in fasta files
     my ($file,$n)=@_;	
     my %seq=();
     my %timepoint=();
     my $id='';
     my @selected=();
	open HH,"$file" or die "seq file $file not found :(\n";
  while(<HH>){
	  chomp;
  	if($_=~/>([^\t ]+)/){
  		my @id=split/$para{'-spliter'}/,$1; 		
  		$id=$1;
  		my $timepoint=$id[$para{'-nc'}-1];
  		if($timepoint=~/[^\d\-]/){die "Time info contain no digit characters\n";}
  		if($timepoint=~/[0-9]{4}\-[0-9]{2}\-[0-9]{2}/){$timepoint=&year_to_digital($timepoint);}
  		 push @{$timepoint{$timepoint}},$id;
  	}
	 else{
	  $seq{$id}.=$_;	
	  }
   }
  if($n>0){
    foreach(sort {$a<=>$b} keys %timepoint){
       my @seq_name=&random_sele($n,@{$timepoint{$_}});
       print "Time point $_: ",@seq_name/1," selected\n";
     	 push @selected, @seq_name;
    }
  
	  return \%seq,\@selected;	
	}
	else{
	 return \%seq;
	}
}
###############
sub random_sele{	#randomly select sequences
   my ($num_seq,@ids)=@_;	
   if($num_seq>@ids){return @ids;}
   else{
   	   &fisher_yates_shuffle(\@ids);
	     return @ids[0..($num_seq-1)];
	}
	
}
#####################
# randomly permutate @array in place
sub fisher_yates_shuffle
{
    my $array = shift;
    my $jj=3;
    while($jj>0){
    	$jj--;
    my $i = @$array;
    while ( --$i )
    {
        my $j = int rand( $i+1 );
        @$array[$i,$j] = @$array[$j,$i];
    }
  }
}
