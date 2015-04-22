#!/usr/bin/perl
use strict;
#use lib ("/Users/sheng/work/HIV/scripts/github/zap/zap/");
use PPvars qw(ppath);
################################################
#This is a mater script that allow user to generate xml configuration files and perform MCMC simulation to estimate evolutionary rate. User has to install beast package first. The script accepts fasta format sequence file as input. Be aware that only temporal data with at least from two time points, are allowed to calculate evolutionary rate. All sequences in the dataset must be non-redundant. Please remove duplicate sequences before running the script. If there are duplicate sequences in the dataset, it's better keep the one with earlier time point. Currently, only GTR substitution model was incorporated but our tests showed the rate calculation is robust to substitution models. it's users' responsibility to check if the MCMC simulation converged using Tracer in Beast2 Package.
################################################

my $usage="Usage:
Evolutionary rate calculation on a lineage from at least two time points. Only nucleotide sequences accepted at this time
Options:
	-f\tinput antibody nucleotide VDJ sequence file with fasta format
	-CDR\tsequences of the three CDRs
	-FW\tsequences of the four FWs
	-spliter\twhich marker used to seperate seq name, default:_
	-nc\twhich column is time info after parsing using -spliter in the sequence name, start from 1. The time 
	   \tinfo must be numerical. default:1
	-pop\tmodel of coalecent population size, constant or bayesian_skyline, default is later
	-o\toutput file prefix
	-n \thow many sequences to randomly choose from each time point. default: 25
	-chainLength\tnumber of steps to run, default:10000000
	-storeEvery\tstore log info every n steps, default:1000
	-log\tstore a tree every n step, default:1000
	-beast\tpath to beast2 program
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
if(!$para{'-storeEvery'}){$para{'-storeEvery'}=1000;}
if(!$para{'-burnin'}){$para{'-burnin'}=10;}
######Processing######
my $seq_CDR='';
my $seq_FW='';
my ($seq,$seq_name_selected)=&read_fasta($para{'-f'},$para{'-n'});#read sequences and select n sequences for each time point
if($para{'-CDR'}){$seq_CDR=&read_fasta($para{'-CDR'});}
if($para{'-FW'}){$seq_FW=&read_fasta($para{'-FW'});}
&write_xml();
system("$para{'-beast'} $para{'-o'}.xml");
&parsing_results("$para{'-o'}.log");
if(-d "./output/rate/"){
  system("mv $para{'-o'}.xml $para{'-o'}.xml.state $para{'-o'}.log $para{'-o'}.trees ./output/rate/");	
}

############################
sub parsing_results{#calculate mean rate and 95% highest probability density range (or confidence interval)
    my $file=shift;
    open HH,"$file" or die "resulte file $file not found\n";#read in log file
    my $mean_c=0;
    my $mean_CDR=0;
    my $mean_FW=0;	
    my @rate=();
    my @rate_CDR=();
    my @rate_FW=();
	  while(<HH>){
	  	if($_=~/^Sample/){
	  		 my @line=split/	/,$_;
	  		 my $i=0;
	  		 foreach(@line){
	  		   if($_=~/rate.mean/){#find colomn for mean rate 
	  		       	$mean_c=$i;
	  		    }
	  		   elsif($_=~/rate.+FW.mean/){#find colomn for framework region rate 
	  		       	$mean_FW=$i;		
	  		    }
	  		   elsif($_=~/rate.+CDR.mean/){#find colomn for CDRs region rate
	  		   	   $mean_CDR=$i;
	  		  }
	  		  $i++;
	  		}
	  	}
	  elsif($_=~/^[0-9]/){
	     my @line1=split/\t/,$_;
	     push @rate,$line1[$mean_c];
	     if($mean_CDR>0){
	        push @rate_CDR,$line1[$mean_CDR];	
	     }	
	     if($mean_FW>0){
	        push @rate_FW,$line1[$mean_FW];	
	     }	
	   }
	  	
	  }
	  
	  #output rates
	  for(my $j=0;$j<@rate/$para{'-burnin'};$j++){
	     shift @rate;	
	  }
	  
	  my @rate_s=&statistic(@rate);
	  print "Mean rate: $rate_s[0]\n95% HPD interval: [$rate_s[3],$rate_s[4]]\n";
	  if(@rate_CDR){
	  	for(my $j=0;$j<@rate/$para{'-burnin'};$j++){
	     shift @rate_CDR;	#10% steps of burn in
	    }
	    @rate_s=&statistic(@rate_CDR);
	    print "Mean rate for CDRs: $rate_s[0]\n95% HPD interval: [$rate_s[3],$rate_s[4]]\n";
	  }
	  if(@rate_FW){
	  	for(my $j=0;$j<@rate/$para{'-burnin'};$j++){
	     shift @rate_FW;	#10% steps of burn in
	    }
	    @rate_s=&statistic(@rate_FW);
	    print "Mean rate for FWs: $rate_s[0]\n95% HPD interval: [$rate_s[3],$rate_s[4]]\n";
	  }
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
open YY,">$para{'-o'}.xml";

print YY '<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate=\'Standard\' beautistatus=\'\' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" version="2.0">',"\n\n\n    <data\n";#print title
print YY "id=\"$para{'-o'}\"\nname=\"alignment\">\n";
my $taxas=@$seq_name_selected;
foreach(@$seq_name_selected){
	print YY "                        <sequence id=\"seq_$_\" taxon=\"$_\" totalcount=\"4\" value=\"$seq->{$_}\"/>\n";	
	}
print YY "                    </data>\n";
if($para{'-FW'}){
	print YY "<data\nid=\"$para{'-o'}\_FW\">\n";
	foreach(@$seq_name_selected){
		if(!$seq_FW->{$_}){die "Please make sure you have the same sequencs in full length and FW files\n";}
	print YY "                        <sequence id=\"seq_$_",1,"\" taxon=\"$_\" totalcount=\"4\" value=\"$seq_FW->{$_}\"/>\n";	
	}
	print YY "                    </data>\n";
}
if($para{'-CDR'}){
  print YY "<data\nid=\"$para{'-o'}\_CDR\">\n"; 
	foreach(@$seq_name_selected){
		if(!$seq_CDR->{$_}){die "Please make sure you have the same sequencs in full length and CDR files\n";}
	print YY "                        <sequence id=\"seq_$_",2,"\" taxon=\"$_\" totalcount=\"4\" value=\"$seq_CDR->{$_}\"/>\n";	
	}
	print YY "                    </data>\n";
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
	$taxon.="$_=$line[$para{'-nc'}-1],\n"
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
if($para{'-FW'}){
print YY "        <parameter id=\"ucldMean.c:$para{'-o'}_FW\" name=\"stateNode\">1.0</parameter>
        <parameter id=\"ucldStdev.c:$para{'-o'}_FW\" lower=\"0.0\" name=\"stateNode\" upper=\"5.0\">0.5</parameter>
        <stateNode dimension=\"$taxas\" id=\"rateCategories.c:$para{'-o'}_FW\" spec=\"parameter.IntegerParameter\">1</stateNode>\n";

}
if($para{'-CDR'}){
print YY "        <parameter id=\"ucldMean.c:$para{'-o'}_CDR\" name=\"stateNode\">1.0</parameter>
        <parameter id=\"ucldStdev.c:$para{'-o'}_CDR\" lower=\"0.0\" name=\"stateNode\" upper=\"5.0\">0.5</parameter>
        <stateNode dimension=\"$taxas\" id=\"rateCategories.c:$para{'-o'}_CDR\" spec=\"parameter.IntegerParameter\">1</stateNode>\n";

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
      if($para{'-FW'}){print YY "            <prior id=\"MeanRatePrior.c:$para{'-o'}\_FW\" name=\"distribution\" x=\"\@ucldMean.c:$para{'-o'}\_FW\">
                <Uniform id=\"Uniform.01\" name=\"distr\" upper=\"Infinity\"/>
            </prior>\n";}
      if($para{'-CDR'}){print YY "            <prior id=\"MeanRatePrior.c:$para{'-o'}\_CDR\" name=\"distribution\" x=\"\@ucldMean.c:$para{'-o'}\_CDR\">
                <Uniform id=\"Uniform.02\" name=\"distr\" upper=\"Infinity\"/>
            </prior>\n";}
print YY "            <prior id=\"ucldStdevPrior.c:$para{'-o'}\" name=\"distribution\" x=\"\@ucldStdev.c:$para{'-o'}\">
                <Exponential id=\"Exponential.01\" name=\"distr\">
                    <parameter estimate=\"false\" id=\"RealParameter.011\" name=\"mean\">0.3333</parameter>
                </Exponential>
            </prior>\n";
      if($para{'-FW'}){print YY "            <prior id=\"ucldStdevPrior.c:$para{'-o'}\_FW\" name=\"distribution\" x=\"\@ucldStdev.c:$para{'-o'}\_FW\">
                <Exponential id=\"Exponential.02\" name=\"distr\">
                    <parameter estimate=\"false\" id=\"RealParameter.012\" name=\"mean\">0.3333</parameter>
                </Exponential>
            </prior>\n";}
      if($para{'-CDR'}){print YY "            <prior id=\"ucldStdevPrior.c:$para{'-o'}\_CDR\" name=\"distribution\" x=\"\@ucldStdev.c:$para{'-o'}\_CDR\">
                <Exponential id=\"Exponential.03\" name=\"distr\">
                    <parameter estimate=\"false\" id=\"RealParameter.013\" name=\"mean\">0.3333</parameter>
                </Exponential>
            </prior>\n";}      
print YY "            <distribution id=\"all.prior\" monophyletic=\"true\" spec=\"beast.math.distributions.MRCAPrior\" tree=\"\@Tree.t:$para{'-o'}\">
                <taxonset id=\"all\" spec=\"TaxonSet\">\n";
foreach(@$seq_name_selected){
   print YY "                    <taxon id=\"$_\" spec=\"Taxon\"/>\n";	
}

print YY "               </taxonset>
                <Uniform id=\"Uniform.03\" lower=\"1.0\" name=\"distr\" upper=\"Infinity\"/>
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
                        <parameter estimate=\"false\" id=\"RealParameter.014\" lower=\"0.0\" name=\"M\" upper=\"1.0\">1.0</parameter>
                    </LogNormal>
                </branchRateModel>
            </distribution>\n";
 if($para{'-FW'}){
 
    print YY "<distribution id=\"treeLikelihood.$para{'-o'}\_FW\" siteModel=\"\@SiteModel.s:$para{'-o'}\" spec=\"TreeLikelihood\" tree=\"\@Tree.t:$para{'-o'}\">
                <data
idref=\"$para{'-o'}\_FW\"/>
                <branchRateModel clock.rate=\"\@ucldMean.c:$para{'-o'}\_FW\" id=\"RelaxedClock.c:$para{'-o'}\_FW\" rateCategories=\"\@rateCategories.c:$para{'-o'}\_FW\" spec=\"beast.evolution.branchratemodel.UCRelaxedClockModel\" tree=\"\@Tree.t:$para{'-o'}\">
                    <LogNormal S=\"\@ucldStdev.c:$para{'-o'}\_FW\" id=\"LogNormalDistributionModel.c:$para{'-o'}\_FW\" meanInRealSpace=\"true\" name=\"distr\">
                        <parameter estimate=\"false\" id=\"RealParameter.015\" lower=\"0.0\" name=\"M\" upper=\"1.0\">1.0</parameter>
                    </LogNormal>
                </branchRateModel>
            </distribution>\n";
 
}
 if($para{'-CDR'}){
 
    print YY "<distribution id=\"treeLikelihood.$para{'-o'}\_CDR\" siteModel=\"\@SiteModel.s:$para{'-o'}\" spec=\"TreeLikelihood\" tree=\"\@Tree.t:$para{'-o'}\">
                <data
idref=\"$para{'-o'}\_CDR\"/>
                <branchRateModel clock.rate=\"\@ucldMean.c:$para{'-o'}\_CDR\" id=\"RelaxedClock.c:$para{'-o'}\_CDR\" rateCategories=\"\@rateCategories.c:$para{'-o'}\_CDR\" spec=\"beast.evolution.branchratemodel.UCRelaxedClockModel\" tree=\"\@Tree.t:$para{'-o'}\">
                    <LogNormal S=\"\@ucldStdev.c:$para{'-o'}\_CDR\" id=\"LogNormalDistributionModel.c:$para{'-o'}\_CDR\" meanInRealSpace=\"true\" name=\"distr\">
                        <parameter estimate=\"false\" id=\"RealParameter.016\" lower=\"0.0\" name=\"M\" upper=\"1.0\">1.0</parameter>
                    </LogNormal>
                </branchRateModel>
            </distribution>\n";
 
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
    
if($para{'-FW'}){
print YY"    <operator id=\"ucldMeanScaler.c:$para{'-o'}\_FW\" parameter=\"\@ucldMean.c:$para{'-o'}\_FW\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"1.0\"/>

    <operator id=\"ucldStdevScaler.c:$para{'-o'}\_FW\" parameter=\"\@ucldStdev.c:$para{'-o'}\_FW\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"3.0\"/>

    <operator id=\"CategoriesRandomWalk.c:$para{'-o'}\_FW\" parameter=\"\@rateCategories.c:$para{'-o'}\_FW\" spec=\"IntRandomWalkOperator\" weight=\"10.0\" windowSize=\"1\"/>

    <operator id=\"CategoriesSwapOperator.c:$para{'-o'}\_FW\" intparameter=\"\@rateCategories.c:$para{'-o'}\_FW\" spec=\"SwapOperator\" weight=\"10.0\"/>

    <operator id=\"CategoriesUniform.c:$para{'-o'}\_FW\" parameter=\"\@rateCategories.c:$para{'-o'}\_FW\" spec=\"UniformOperator\" weight=\"10.0\"/>

    <operator id=\"relaxedUpDownOperator.c:$para{'-o'}\_FW\" scaleFactor=\"0.75\" spec=\"UpDownOperator\" weight=\"3.0\">
        <parameter idref=\"ucldMean.c:$para{'-o'}\_FW\" name=\"up\"/>
        <tree idref=\"Tree.t:$para{'-o'}\" name=\"down\"/>
    </operator>\n";
}
if($para{'-CDR'}){
print YY "     <operator id=\"ucldMeanScaler.c:$para{'-o'}\_CDR\" parameter=\"\@ucldMean.c:$para{'-o'}\_CDR\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"1.0\"/>

    <operator id=\"ucldStdevScaler.c:$para{'-o'}\_CDR\" parameter=\"\@ucldStdev.c:$para{'-o'}\_CDR\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"3.0\"/>

    <operator id=\"CategoriesRandomWalk.c:$para{'-o'}\_CDR\" parameter=\"\@rateCategories.c:$para{'-o'}\_CDR\" spec=\"IntRandomWalkOperator\" weight=\"10.0\" windowSize=\"1\"/>

    <operator id=\"CategoriesSwapOperator.c:$para{'-o'}\_CDR\" intparameter=\"\@rateCategories.c:$para{'-o'}\_CDR\" spec=\"SwapOperator\" weight=\"10.0\"/>

    <operator id=\"CategoriesUniform.c:$para{'-o'}\_CDR\" parameter=\"\@rateCategories.c:$para{'-o'}\_CDR\" spec=\"UniformOperator\" weight=\"10.0\"/>

    <operator id=\"relaxedUpDownOperator.c:$para{'-o'}\_CDR\" scaleFactor=\"0.75\" spec=\"UpDownOperator\" weight=\"3.0\">
        <parameter idref=\"ucldMean.c:$para{'-o'}\_CDR\" name=\"up\"/>
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
if($para{"-FW"}){
    print YY "        <parameter idref=\"ucldMean.c:$para{'-o'}\_FW\" name=\"log\"/>
        <parameter idref=\"ucldStdev.c:$para{'-o'}\_FW\" name=\"log\"/>
        <log branchratemodel=\"\@RelaxedClock.c:$para{'-o'}\_FW\" id=\"rate.c:$para{'-o'}\_FW\" spec=\"beast.evolution.branchratemodel.RateStatistic\" tree=\"\@Tree.t:$para{'-o'}\"/>\n";
}
if($para{"-CDR"}){
    print YY "        <parameter idref=\"ucldMean.c:$para{'-o'}\_CDR\" name=\"log\"/>
        <parameter idref=\"ucldStdev.c:$para{'-o'}\_CDR\" name=\"log\"/>
        <log branchratemodel=\"\@RelaxedClock.c:$para{'-o'}\_CDR\" id=\"rate.c:$para{'-o'}\_CDR\" spec=\"beast.evolution.branchratemodel.RateStatistic\" tree=\"\@Tree.t:$para{'-o'}\"/>\n";
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
if($para{'-FW'}){
    print YY "        <parameter idref=\"ucldMean.c:$para{'-o'}\_FW\" name=\"log\"/>
        <parameter idref=\"ucldStdev.c:$para{'-o'}\_FW\" name=\"log\"/>
        <log branchratemodel=\"\@RelaxedClock.c:$para{'-o'}\_FW\" id=\"rate.c:$para{'-o'}\_FW\" spec=\"beast.evolution.branchratemodel.RateStatistic\" tree=\"\@Tree.t:$para{'-o'}\"/>\n";
}
if($para{"-CDR"}){
    print YY "        <parameter idref=\"ucldMean.c:$para{'-o'}\_CDR\" name=\"log\"/>
        <parameter idref=\"ucldStdev.c:$para{'-o'}\_CDR\" name=\"log\"/>
        <log branchratemodel=\"\@RelaxedClock.c:$para{'-o'}\_CDR\" id=\"rate.c:$para{'-o'}\_CDR\" spec=\"beast.evolution.branchratemodel.RateStatistic\" tree=\"\@Tree.t:$para{'-o'}\"/>\n";
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
  	if($_=~/>(.+)/){
  		my @id=split/$para{'-spliter'}/,$1; 		
  		$id=$1;
  		 push @{$timepoint{$id[$para{'-nc'}-1]}},$id;
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