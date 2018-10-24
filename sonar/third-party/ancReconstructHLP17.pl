#!/usr/bin/perl
#Kenneth Hoehn
#23/8/2016
#Simulate sequences under HLP16 model
#Parameters read in through config.txt

use strict;
use warnings;
use hotSpotUtilsPDL;

#Read in parameters from congif file
my $config = $ARGV[0];
open(C,$config) or die("Couldn't open config file ($config)");

my $nsim;
my $kappa;
my @omegas;
my @motifs;
my $partfile;
my $treefile;
my @hs;
my $freqs;
my $length;
my $outdir;
my $igphyml;
my $seqfile;
my $context=0;
my $rooted=0;
my $rootid;
my $ancstate=0;
my $statsfile;
my $stem;
my $statslhood;
my $ambigfile;
my $aamodel;
my @aaints;
my $aaslope="N";
my $aatree = "N";
my $alrt=0;
my $hint = "N";
my $hslope = "N";
my $rootonly = 0;
my $rootpis = 1;
my $bstats;

while(<C>){
	my $line = $_;
	chomp($line);
	if($line =~ /nsim\s+(\S+)/){
		$nsim = $1;
	}
	if($line =~ /omegas\s+(\S+)/){
		@omegas = split(",",$1);
	}
	if($line =~ /kappa\s+(\S+)/){
		$kappa = $1;
	}
	if($line =~ /motifs\s+(\S+)/){
		@motifs = split(",",$1);
	}
	if($line =~ /hs\s+(\S+)/){
		@hs = split(",",$1);
	}
	if($line =~ /freqs\s+(\S+)/){
		$freqs = $1;
	}
	if($line =~ /^tree\s+(\S+)$/){
		$treefile = $1;
	}
	if($line =~ /fullcontext\s+(\S+)/){
		$context=$1;
	}
	if($line =~ /outdir\s+(\S+)/){
		$outdir=$1;
	}
	if($line =~ /rooted\s+(\S+)/){
		$rooted=$1;
	}
	if($line =~ /length\s+(\S+)/){
		$length=$1;
	}
	if($line =~ /rootid\s+(\S+)/){
		$rootid=$1;
	}
	if($line =~ /part\s+(\S+)/){
		$partfile=$1;
	}
	if($line =~ /igphyml\s+(\S+)/){
		$igphyml=$1;
	}
	if($line =~ /seqfile\s+(\S+)/){
		$seqfile=$1;
	}
	if($line =~ /ancstate\s+(\S+)/){
		$ancstate=$1;
	}
	if($line =~ /stats\s+(\S+)/){
		print "match $&\n";
		$statsfile=$1;
		print "statsfile $1 $statsfile\n";
	}
	if($line =~ /stem\s+(\S+)/){
		$stem=$1;
	}
	if($line =~ /ambigfile\s+(\S+)/){
		$ambigfile=$1;
	}
	if($line =~ /aatree\s+(\S+)/){
		$aatree=$1;
	}
	if($line =~ /alrt\s+(\S+)/){
		$alrt=$1;
	}
	if($line =~ /rootonly\s+(\S+)/){
		$rootonly=$1;
	}
	if($line =~ /rootpis\s+(\S+)/){
		$rootpis=$1;
	}
}

#check to see if stats file was specified in command line
for(my $i = 1; $i < scalar(@ARGV); $i++){
	my $line = $ARGV[$i];
	if($line =~ /-stats/){
		$statsfile=$ARGV[$i+1];
	}
}

#Read in igphyml stats file, if specified
print "stats: $statsfile\n";
if(defined $statsfile && $statsfile ne "N"){
	open(STATS,$statsfile)or die("Couldn't open $statsfile\n");
	my @stats = <STATS>;
	#print "@stats\n";
	@motifs = (0)x0;
	@omegas = (0)x0;
	@hs = (0)x0;
	$freqs = "stats";
	foreach my $l (@stats){
	  chomp($l);
	  #print "$l\n";
	  if($l =~ /Motif:\s+(\S+)#\d+\s+\d\s+\d\s+(\S+)/){
	  	push(@motifs,$1);
	  	push(@hs,$2);
	  	print "Read motif h $1 = $2 from $statsfile\n";
	  }elsif($l =~ /Motif:\s+(\S+)\s+\d\s+\d\s+(\S+)/){
	  	push(@motifs,$1);
	  	push(@hs,$2);
	  	print "Read motif h $1 = $2 from $statsfile\n";
	  }
	  if($l =~ /\. Omega\s+(\d+)\s+\S+:\s+(\S+)/){
	  	$omegas[$1]=$2;
	  	print "Read omega $1 = $2 from $statsfile\n";
	  }
	  if($l =~ /. Nonsynonymous\/synonymous ratio:\s+(\S+)/){
	  	$omegas[0]=$1;
	  	print "Read old school omega 0 = $1 from $statsfile\n";
	  }
	  if($l =~ /\. Transition\/transversion ratio:\s+(\S+)/){
	  	$kappa=$1;
	  	print "Read kappa $kappa from $statsfile\n";
	  }
	  if($l =~ /\. Log-likelihood:\s+(\S+)/){
	  	$statslhood=$1;
	  	print "Read stats lhood $statslhood from $statsfile\n";
	  }
	  if($l =~ /. Partition file:\s+(\S+)/){
	  	$partfile=$1;
	  	print "Read partfile $partfile from $statsfile\n";
	  }
	  if($l =~ /.\s+(\S+)\s+aa_intercept\s+(\d+)\s+(\S+):\s+(\S+)/){
	  	$aamodel = $1;
	  	$aaints[$2]=$4;
	  	print "Read $aamodel $aaints[$2] from $statsfile\n";
	  }
	  if($l =~ /.\s+(\S+)\s+aa_slope:\s+(\S+)/){
	  	$aaslope = $2;
	  	print "Read $aamodel $aaslope from $statsfile\n";
	  }
	  if($l =~ /Motif:\s+S5F\s+h_intercept:\s+(\S+)/){
	  	$hint = $1;
	  	print "Read motif hint S5F = $1 from $statsfile\n";
	  }
	  if($l =~ /Motif:\s+S5F\s+h_slope:\s+(\S+)/){
	  	$hslope = $1;
	  	print "Read motif hslope S5F = $1 from $statsfile\n";
	  }
	  $bstats .= " $l";
	}
	print "Reading in frequency parameters from IgPhyML stats file\n";
}

#check command line args to see if any should be over-ridden
for(my $i = 1; $i < scalar(@ARGV); $i++){
	my $line = $ARGV[$i];
	if($line =~ /-nsim/){
		$nsim = $ARGV[$i+1];
	}
	if($line =~ /-omegas/){
		@omegas = split(",",$ARGV[$i+1]);
	}
	if($line =~ /-kappa/){
		$kappa = $ARGV[$i+1];
	}
	if($line =~ /-part/){
		$partfile = $ARGV[$i+1];
	}
	if($line =~ /-motifs/){
		@motifs = split(",",$ARGV[$i+1]);
	}
	if($line =~ /-hs/){
		@hs = split(",",$ARGV[$i+1]);
	}
	if($line =~ /-freqs/){
		$freqs = $ARGV[$i+1];
	}
	if($line =~ /-tree/){
		$treefile = $ARGV[$i+1];
	}
	if($line =~ /-fullcontext/){
		$context=$ARGV[$i+1];
	}
	if($line =~ /-outdir/){
		$outdir=$ARGV[$i+1];
	}
	if($line =~ /-rooted/){
		$rooted=$ARGV[$i+1];
	}
	if($line =~ /-length/){
		$length=$ARGV[$i+1];
	}
	if($line =~ /-rootid/){
		$rootid=$ARGV[$i+1];
	}
	if($line =~ /-igphyml/){
		$igphyml=$ARGV[$i+1];
	}
	if($line =~ /-seqfile/){
		$seqfile=$ARGV[$i+1];
	}
	if($line =~ /-ancstate/){
		$ancstate=$ARGV[$i+1];
	}
	if($line =~ /-stats/){
		$statsfile=$ARGV[$i+1];
	}
	if($line =~ /-stem/){
		$stem=$ARGV[$i+1];
	}
	if($line =~ /-ambigfile/){
		$ambigfile=$ARGV[$i+1];
	}
	if($line =~ /-aatree/){
		$aatree=$ARGV[$i+1];
	}
	if($line =~ /-alrt/){
		$alrt=$ARGV[$i+1];
	}
	if($line =~ /-rootonly/){
		$rootonly=$ARGV[$i+1];
	}
	if($line =~ /-rootpis/){
		$rootpis=$ARGV[$i+1];
	}
}

#check that all necessary parameters are specified
if(!defined $kappa){die("kappa needs to be specified")}
#if(scalar(@omegas)==0){die("omegas needs to be specified")}
#if(scalar(@motifs)==0){die("motifs needs to be specified")}
if(scalar(@hs)==0 && $hint eq "N"){die("hs needs to be specified")}
if(!defined $freqs){die("freqs needs to be specified")}
if(!defined $outdir){die("outdir needs to be specified")}
if(!defined $length){die("length needs to be specified")}
if(!defined $rootid){die("rootid needs to be specified")}
if(!defined $igphyml){die("igphyml needs to be specified")}
if(!defined $stem){die("stem needs to be specified")}
if(!defined $seqfile){die("seqfile needs to be specified");}
if(!defined $statsfile){$statsfile="N";}
if(!defined $partfile){$partfile="N";}
if(!defined $ambigfile){$ambigfile="N";}

print "\nReconstruction Settings\n";
print "kappa\t$kappa\n";
print "omegas\t@omegas\n";
print "motifs\t@motifs\n";
print "hs\t@hs\n";
print "freqs\t$freqs\n";
print "outdir\t$outdir\n";
print "length\t$length\n";
print "rootid\t$rootid\n";
print "seqfile\t$seqfile\n";
print "stats\t$statsfile\n";
print "ambigfile\t$ambigfile\n";
print "outfile format: $outdir/$stem\_\n";
print "aatree: $aatree\n";
print "alrt: $alrt\n";
print "seqfile: $seqfile\n";
print "rootonly: $rootonly\n";
print "rootpis: $rootpis\n";
print "\n";

if(scalar(@hs) ne scalar(@motifs)){die(scalar(@hs)." h values but ".scalar(@motifs)." motifs!\n")}

#Read in tree
open(TREE,$treefile)or die("Couldn't open $treefile");
my $tree = <TREE>;
chomp($tree);
my %root;
if($rooted==1){
	%root = %{readInRootedNewick($tree,0,$alrt)};
}else{
	%root = %{readInUnrootedNewick($tree,$rootid,0)};
}

my $t = printtreestring(\%root,"",$alrt).";";
print "Tree! $t\n";

#Check to see that the root sequence is in good order
my $seqs;
my @rs;
if(defined $seqfile || $seqfile ne "N"){
	$seqs = getfasta($seqfile);
	if(!exists($seqs->{$rootid})){
		die("$rootid not found in sequence file: $seqfile\n");
	}
	if(length($seqs->{$rootid}) ne $length*3){
		die("Specified root sequence is not specified length: $length ".(length($seqs->{$rootid})/3)."\n");
	}
}else{
	die("Need sequences and root!");
}

#read in ambiguous character file
my %ambig;
if($ambigfile ne "N"){
open(AM,$ambigfile) or die("Can't open ambiguous character file\n");
while(<AM>){
	my $line = $_;
	chomp($line);
	my @in = split(" ",$line);
	if(!exists($ambig{$in[0]})){
		my %new;
		$ambig{$in[0]} = \%new;
	}
	if(!exists($ambig{$in[0]}->{$in[1]})){
		my @new = (0)x61;
		$ambig{$in[0]}->{$in[1]} = \@new;
	}
	$ambig{$in[0]}->{$in[1]}->[$in[2]]=$in[3];
}
}
#Set up partition model from file
my @part = ((-1)x$length);
my @rpart = ((-1)x$length);
my @siteind; #site partition index
my @omegaind; #index in omega vector
my @branchind; #which branch label
my @hbranchind; #which hotspot model
my $nparts=0;
my $nrparts=0;
my $nqmats=1;
if($partfile ne "N"){
	open(P,$partfile) or die("Couldn't open $partfile");
	my $h = <P>;
	my @info = split(" ",$h);
	$nparts = $info[1];
	$nrparts=$info[0];
	my $brmodel=$info[3];
	print "nrparts: $nrparts\n";
	#read in site partitions
	for(0..($nrparts-1)){
		my $line = <P>;
		chomp($line);
		print "line: $line\n";
		my @in1 = split(" ",$line);
		my @in2 = split(",",$in1[2]);
		print "@in2\n";
		for(my $i=0;$i<scalar(@in2);$i++){
			my @in3=split("\\.\\.",$in2[$i]);
			for(my $j=$in3[0];$j<=$in3[1];$j++){
				if($j >= $length){die("Partition $nparts extends beyond specified sequnece length!")}
				$part[$j]=$in1[0];
				if($in1[1] eq "FWR"){
					$rpart[$j]=0;
				}elsif($in1[1] eq "CDR"){
					$rpart[$j]=1;
				}else{
					die($in1[1]);
				}
			}
		}
		#$nrparts++;
	}
	if($brmodel==1){
		$nqmats=<P>; #number of qmatrices
		print "n qmats: $nqmats";
		chomp($nqmats);
		for(my$i=0;$i<$nqmats;$i++){
			my $line=<P>;
			my @in = split(" ",$line);
			push(@siteind,$in[0]);
			push(@omegaind,$in[1]);
			push(@branchind,$in[2]);
			push(@hbranchind,$in[3]);
		}
	}else{
		$nqmats=$nparts;
		for(my$i=0;$i<$nqmats;$i++){
			push(@siteind,$i);
			push(@omegaind,$i);
			push(@branchind,0);
			push(@hbranchind,0);
		}
	}
}else{ #default of a single partition
	@part = ((0)x$length);
	@rpart = ((0)x$length);
	$nparts=1;
	$nrparts=1;
	push(@siteind,0);
	push(@omegaind,0);
	push(@branchind,0);
	push(@hbranchind,0);
}
print "omegaind: @omegaind\n";
print "omegas: @omegas\n";
#if($nparts != scalar(@omegas)){die("$nparts partitions, but ".(scalar(@omegas)." omegas!"))}
print "Omega partition index:\n";
for(my$i=0;$i<$length;$i++){
	if($part[$i] == -1){
		die("Position $i unspecified in partition file.");
	}else{
		print "$part[$i]";
	}
}
print "\n";
print "Region partition index:\n";
for(my$i=0;$i<$length;$i++){
	if($rpart[$i] == -1){
		#die("Position $i unspecified in partition file.");
	}else{
		print "$rpart[$i]";
	}
}
print "\n";


#read in FWR_CDR file
#open(FWRCDR,"$igphyml/src/motifs/HTABLE_CDR_FWR_UNK.single.csv") or die();
my @rootCDR = (0)x61;
my @rootFWR = (0)x61;
my @rootUNK = (0)x61;
#for(my $i=0;$i<183;$i++){
#	my $num = <FWRCDR>;
	#print "$num\n";
#	chomp($num);
#	if($i < 61){
#		 $rootCDR[$i]=$num;
#	}elsif($i < 122){
#	 	$rootFWR[$i-61]=$num;
#	}else{
	 #die($i)
	 #print "$i\t$num\n";
#	 	$rootUNK[$i-122]=$num;
#	}
#}
print "@rootUNK\n";

#Make codon indexes and get frequencies
my @transfreq;
my @codons;
my %codoni;
my $index = 0;
my @chars = ("t","c","a","g");
my %freqs;
my $fsum=0;
foreach my $a (@chars){
	foreach my $b (@chars){
		foreach my $c (@chars){
			if($a.$b.$c ne "tga" && $a.$b.$c ne "taa" && $a.$b.$c ne "tag"){
				push(@codons,$a.$b.$c);
				$codoni{$a.$b.$c}=$index;
				my $match = uc $a.$b.$c;
				if(defined $statsfile && $statsfile ne "N"){
					if($bstats =~ /f\($match\)=(\d+\.*\d*)/){
						$freqs{lc $match} = $1;
						if($freqs{lc $match} == 0){
							print "Zero frequency caught\n";
							$freqs{lc $match}=1e-10;
						}
					}else{die($match)}
				}elsif($freqs eq "uniform"){
					$freqs{lc $match} = 1/61;
				}else{
					die("freqs not set properly\n");
				}
				$fsum += $freqs{$a.$b.$c};
				$index++;
			}
		}
	}
}
foreach my $k (keys %freqs){
	$freqs{$k} = $freqs{$k}/$fsum;
	$transfreq[$codoni{$k}]=$freqs{$k};
}
print "Codon frequencies: @transfreq\n";

#read in aatable
my @aatable;
if($aaslope ne "N"){
	open(HTABLE,"$igphyml/src/motifs/HTABLE_AAtable\_$aamodel") or die("Couldn't open $igphyml/src/motifs/HTABLE_AAtable\_$aamodel\n");
	print "reading $igphyml/src/motifs/HTABLE_AAtable\_$aamodel\n";
	while(<HTABLE>){
		my $l = $_;
		chomp($l);
		push(@aatable,$l);
	}
	close(HTABLE);
}

#Make B matrix
print "Setting up B matrix\n";
my @Bmat = (0)x(61*61);
my $fi;my $ti;my $li;my $ri;

if($hint eq "N"){
for(my $mi = 0; $mi < scalar(@motifs); $mi++){
	my $motif = $motifs[$mi];
	print "Reading in $motif table\n";
	my @htable;
	open(HTABLE,"$igphyml/src/motifs/HTABLE_$motif") or die("Couldn't open $igphyml/src/motifs/HTABLE_$motif\n");
	while(<HTABLE>){
		my $l = $_;
		chomp($l);
		push(@htable,$l);
	}
	close(HTABLE);

	for($fi=0;$fi<61;$fi++){
		for($ti=0;$ti<61;$ti++){
			my @htotals = (0)x(1);
			for($li=0;$li<61;$li++){
				for($ri=0;$ri<61;$ri++){
					$htotals[0] += $transfreq[$li]*$transfreq[$ri]*$htable[$fi*61*61*61+$ti*61*61+$li*61+$ri];
				}
			}
			my $hsum = $hs[$mi]*$htotals[0];
			$Bmat[61*$fi+$ti]+=$hsum;
		}
	}
	if(scalar(@Bmat) != (61*61)){die("@Bmat")}
}
}else{
	print "Reading in S5F table\n";
	my @htable;
	open(HTABLE,"$igphyml/src/motifs/HTABLE_S5F") or die("Couldn't open $igphyml/src/motifs/HTABLE_S5F\n");
	while(<HTABLE>){
		my $l = $_;
		chomp($l);
		push(@htable,$l);
	}
	close(HTABLE);

	for($fi=0;$fi<61;$fi++){
		for($ti=0;$ti<61;$ti++){
			for($li=0;$li<61;$li++){
				for($ri=0;$ri<61;$ri++){
					my $val = $hint + $hslope*$htable[$fi*61*61*61+$ti*61*61+$li*61+$ri];
					if($val < -1){$val=-1;}
					$Bmat[61*$fi+$ti] += $transfreq[$li]*$transfreq[$ri]*$val;
				}
			}
			if($Bmat[61*$fi+$ti] < -1.0000000001){die("$fi $ti $Bmat[61*$fi+$ti]")}
		}
	}
	if(scalar(@Bmat) != (61*61)){die("@Bmat")}
}

#Make Q matrices
my @Qs;
print "omegaind: @omegaind\n";
print "omegas: @omegas\n";
for(my $i=0;$i<$nqmats;$i++){
	if($aaslope eq "N"){
		print "Making Q matrix $i, omega: ".$omegas[$omegaind[$i]]."\n";
		push(@Qs,getQmat_HLP16(\@Bmat,$kappa,$omegas[$omegaind[$i]],\%freqs,\@codons,0))
	}else{
		print "Making Q matrix $i, aaint: $aaints[$i], aaslope: $aaslope\n";
		push(@Qs,getQmat_SE_HLP16(\@Bmat,$kappa,$aaints[$i],\%freqs,\@codons,0,$aaslope,\@aatable));
	}
}

foreach my $q(@Qs){
	print mexp($q)->at(2,3)."\n";
}
#die();
#collect subtaxa of each node
getSubTaxa(\%root);

#assign node parents
assignparents(\%root,"NA");

setUnlabled(\%root);


print "About to do reconstructions\n";



#Do felsenstein's pruning algorithm
my $lhood = Pruning_Lhood(\%root,$seqs,\%codoni,\@Qs,\@part,$nqmats,\%ambig,\@rpart,\@rootFWR,\@rootCDR,\%freqs,\@codons,\@rootUNK,$rootpis,$nrparts);
#Fill in Up's
Fill_upp(\%root,$seqs,\@Qs,\@part,$nparts,$rootonly);

#Do marginal ancestral state reconstructions
Marginal_ASR(\%root,$seqs,\@Qs,\@part,$nparts,$rootonly);

my $mcodons = printML_codon(\%root,"",\@codons,$rootonly);
print "$mcodons\n";

my $codonTableSingle = codonTableSingle();
my $codonTableTriple = codonTable();

my $maa = printML_aa(\%root,"",\@codons,$codonTableSingle,$codonTableTriple,$rootonly);
print "$maa\n";

print "StatsFile: $statslhood, Reconstruction: $lhood\n";


open(OUT,">$outdir/$stem.MLcodons.fa") or die();
print OUT uc "$mcodons\n";
close(OUT);

open(OUT,">$outdir/$stem.MLaas.fa") or die();
print OUT uc "$maa\n";
close(OUT);
