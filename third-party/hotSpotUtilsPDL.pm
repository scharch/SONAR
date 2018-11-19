#Utilities for doing analysis with BCR phylogenetics, especially with regards to 
#hotspot selection models
#include PDL dependent functions
#29/March/2016
#Kenneth Hoehn
#Uses new form of h parameterization 16/6/2016

use strict;
use warnings;
use PDL;
use PDL::LinearAlgebra::Trans;
use hotSpotUtils;

sub Fill_upp{
    my $node = $_[0];
    my $seqs = $_[1];
    my $Qs = $_[2]; #array reference of q matrixes
    my $partition = $_[3]; #array of partition indexes
    my $nparts = $_[4]; #number of unique partitions
    my $rootonly = $_[5]; #only reconstruct root?

    my @keys = keys %$seqs;
    my $length = length($seqs->{$keys[0]})/3;
    my @bigempty;
    for(my $i=0; $i < $length; $i++){
        my @empty = (-200)x61;
        push(@bigempty,\@empty);
    }
    $node->{"uppmat"}=\@bigempty;

   if($node->{"level"}==0){#if at root node
    #don't do anything, because this isn't a real node :)
   }elsif($node->{"level"}==1){ #one below root
    my $other; #other node to either the left or right
    if($node->{"up"}->{"left"} eq $node){$other="right"}
    elsif($node->{"up"}->{"right"} eq $node){$other="left"}
    else{die("something weird happened")}
    my @Pxz;
      for(my $i=0;$i<$nparts;$i++){
        push(@Pxz,mexp($Qs->[$i]*$node->{"up"}->{$other}->{"dist"}));
      }
      for(my $i=0; $i < $length;$i++){
          for(my $j=0;$j<61;$j++){
              my $sumxz = log($Pxz[$partition->[$i]]->at(0,$j))+$node->{"up"}->{$other}->{"mat"}->[$i][0];
              for(my $k=1;$k<61;$k++){
                  if($Pxz[$partition->[$i]]->at($k,$j)==0){print $node->{"dist"}." $k $j\n"}
                  my $pxz = log($Pxz[$partition->[$i]]->at($k,$j)) + $node->{"up"}->{$other}->{"mat"}->[$i][$k];
                  $sumxz = $sumxz + log(1+exp($pxz-$sumxz));
              }
              $node->{"uppmat"}->[$i][$j] = $sumxz;
          }
      }
   }else{#if not at the root node
    if($rootonly==1){return;}
    my $other; #other node to either the left or right
    if($node->{"up"}->{"left"} eq $node){$other="right"}
    elsif($node->{"up"}->{"right"} eq $node){$other="left"}
    else{die("something weird happened")}
    my @Pxy; my @Pyv;
    for(my $i=0;$i<$nparts;$i++){
      push(@Pxy,mexp($Qs->[$i]*$node->{"up"}->{"dist"}));
      push(@Pyv,mexp($Qs->[$i]*$node->{"up"}->{$other}->{"dist"}));
    }
    #pxy
    for(my $i=0; $i < $length;$i++){
        for(my $j=0;$j<61;$j++){
            my $sumxy = log($Pxy[$partition->[$i]]->at($j,0))+$node->{"up"}->{"uppmat"}->[$i][0];
            for(my $k=1;$k<61;$k++){
                if($Pxy[$partition->[$i]]->at($k,$j)==0){print $node->{"up"}->{"dist"}." $k $j\n"}
                my $pxy = log($Pxy[$partition->[$i]]->at($j,$k)) + $node->{"up"}->{"uppmat"}->[$i][$k];
                $sumxy = $sumxy + log(1+exp($pxy-$sumxy));
            }
            $node->{"uppmat"}->[$i][$j] = $sumxy;
        }
    }
    #pyv
    for(my $i=0; $i < $length;$i++){
        for(my $j=0;$j<61;$j++){
            my $sumyv = log($Pyv[$partition->[$i]]->at(0,$j))+$node->{"up"}->{$other}->{"mat"}->[$i][0];
            for(my $k=1;$k<61;$k++){
                if($Pyv[$partition->[$i]]->at($k,$j)==0){print $node->{"up"}->{$other}->{"dist"}." $k $j ".$Pyv[$partition->[$i]]->at($k,$j)."\n"}
                my $pyv = log($Pyv[$partition->[$i]]->at($k,$j)) + $node->{"up"}->{$other}->{"mat"}->[$i][$k];
                $sumyv = $sumyv + log(1+exp($pyv-$sumyv));
            }
            $node->{"uppmat"}->[$i][$j] += $sumyv;
        }
    }
   }
   #tally up stuff
   my $upphood = 0;
   for(my $i=0; $i < $length;$i++){
    my $sumyv;
        for(my $j=0;$j<61;$j++){
            if($j==0){$sumyv=$node->{"uppmat"}->[$i][$j];}
            else{$sumyv = $sumyv + log(1+exp($node->{"uppmat"}->[$i][$j]-$sumyv));}
        }
        $upphood += $sumyv;
    }
    print "Upphood\t".$node->{"subtaxa"}."\t$upphood\n";

   if(exists($node->{"left"})){ #recurse!
    Fill_upp($node->{"right"},$seqs,$Qs,$partition,$nparts,$rootonly);
    Fill_upp($node->{"left"},$seqs,$Qs,$partition,$nparts,$rootonly);
   }  
}

sub Lk_at_each{
  my $node = $_[0];
  my $seqs = $_[1];
  my $Qs = $_[2]; #array reference of q matrixes
  my $partition = $_[3]; #array of partition indexes
  my $nparts = $_[4]; #number of unique partitions

  if($node->{"level"} != 0){
    my @Pyv;
    for(my $i=0;$i<$nparts;$i++){
      push(@Pyv,mexp($Qs->[$i]*$node->{"dist"}));
    }
    my $lhood = 0;
    for(my $i=0; $i < scalar(@{$node->{"mat"}});$i++){
      my $sitelhood;
      for(my $j=0;$j<61;$j++){
            my $sumyv = $node->{"uppmat"}->[$i][$j]+log($Pyv[$partition->[$i]]->at(0,$j))+$node->{"mat"}->[$i][0];
            for(my $k=1;$k<61;$k++){
              my $pyv = $node->{"uppmat"}->[$i][$j]+log($Pyv[$partition->[$i]]->at($k,$j))+$node->{"mat"}->[$i][$k];
              $sumyv = $sumyv + log(1+exp($pyv-$sumyv));
            }
            if($j==0){$sitelhood = $sumyv;}
            else{$sitelhood = $sitelhood + log(1+exp($sumyv-$sitelhood));}
        }
        $lhood += $sitelhood;
    }
    print $node->{"subtaxa"}."\t$lhood\n";
  }else{

  }
  if(exists($node->{"left"})){
    Lk_at_each($node->{"right"},$seqs,$Qs,$partition,$nparts);
    Lk_at_each($node->{"left"},$seqs,$Qs,$partition,$nparts);
  }
}

sub Marginal_ASR{
  my $node = $_[0];
  my $seqs = $_[1];
  my $Qs = $_[2]; #array reference of q matrixes
  my $partition = $_[3]; #array of partition indexes
  my $nparts = $_[4]; #number of unique partitions
  my $rootonly=$_[5];
  my @keys = keys %$seqs;
  my $length = length($seqs->{$keys[0]})/3;

  if($node->{"level"} != 0){
    my @charmat;
    my @Pyv;
    for(my $i=0;$i<$nparts;$i++){
      push(@Pyv,mexp($Qs->[$i]*$node->{"dist"}));
  }
  my $lhood = 0;
  for(my $i=0; $i < $length;$i++){
    my @sitemat;#relative lhoods of each v at the site
    my $sitelhood;
    my $maxchar;
    my $maxlhood;
    for(my $v=0;$v<61;$v++){
        my $lhoodv;
        for(my $y=0;$y<61;$y++){
          my $val = $node->{"uppmat"}->[$i][$y]+log($Pyv[$partition->[$i]]->at($v,$y))+$node->{"mat"}->[$i][$v];
          if($y==0){$lhoodv=$val;}
          else{$lhoodv = $lhoodv + log(1+exp($val-$lhoodv));}
       }
       push(@sitemat,$lhoodv);
       if($v == 0){$sitelhood = $lhoodv;}
       else{$sitelhood=$sitelhood + log(1+exp($lhoodv-$sitelhood));}
    }
    push(@charmat,\@sitemat);
    $lhood += $sitelhood;
  }
  $node->{"Codon_lhoods"}=\@charmat;
  print "Anc recon:\t".$node->{"subtaxa"}."\tLikelihood:\t$lhood\n";
  }
  if($node->{"level"} == 1 && $rootonly==1){return;}
  if(exists($node->{"left"})){
    Marginal_ASR($node->{"right"},$seqs,$Qs,$partition,$nparts,$rootonly);
    Marginal_ASR($node->{"left"},$seqs,$Qs,$partition,$nparts,$rootonly);
  }
}



sub Pruning_Lhood{
    my $node = $_[0];
    my $seqs = $_[1];
    my $codoni = $_[2];
    my $Qs = $_[3]; #array reference of q matrixes
    my $partition = $_[4]; #array of partition indexes
    my $nqmats = $_[5]; #number of unique partitions
    my $ambig_char = $_[6];
    my $rpartIndex = $_[7];
    my $rootFWR = $_[8];
    my $rootCDR = $_[9];
    my $freqs = $_[10];
    my $codons = $_[11];
    my $rootUNK=$_[12];
    my $rootpis = $_[13];
    my $nrparts = $_[14];

    my @keys = keys %$seqs;
    my $length = length($seqs->{$keys[0]})/3;
    my @bigempty;
    for(my $i=0; $i < $length; $i++){
        my @empty = (-200)x61;
        push(@bigempty,\@empty);
    }
    $node->{"mat"}=\@bigempty;
    if(exists($node->{"left"})){ #internal node
        my $r = Pruning_Lhood($node->{"right"},$seqs,$codoni,$Qs,$partition,$nqmats,$ambig_char,$rpartIndex,$rootFWR,$rootCDR,$freqs,$codons,$rootUNK,$rootpis,$nrparts);
        my $l = Pruning_Lhood($node->{"left"},$seqs,$codoni,$Qs,$partition,$nqmats,$ambig_char,$rpartIndex,$rootFWR,$rootCDR,$freqs,$codons,$rootUNK,$rootpis,$nrparts);
        my @Prs; my @Pls;
        for(my $i=0;$i<$nqmats;$i++){
          push(@Prs,mexp($Qs->[$i]*$node->{"right"}->{"dist"}));
          push(@Pls,mexp($Qs->[$i]*$node->{"left"}->{"dist"}));
        }
        for(my $i=0; $i < $length;$i++){
            for(my $j=0;$j<61;$j++){
                my $sumr = log($Prs[$partition->[$i]+$nrparts*$node->{"right"}->{"olab"}]->at(0,$j))+$node->{"right"}->{"mat"}->[$i][0];
                my $suml = log($Pls[$partition->[$i]+$nrparts*$node->{"left"}->{"olab"}]->at(0,$j))+$node->{"left"}->{"mat"}->[$i][0];
                for(my $k=1;$k<61;$k++){
                    if($Prs[$partition->[$i]]->at($k,$j)==0){print $node->{"dist"}." $k $j\n"}
                    my $pr = log($Prs[$partition->[$i]+$nrparts*$node->{"right"}->{"olab"}]->at($k,$j)) + $node->{"right"}->{"mat"}->[$i][$k];
                    my $pl = log($Pls[$partition->[$i]+$nrparts*$node->{"left"}->{"olab"}]->at($k,$j)) + $node->{"left"}->{"mat"}->[$i][$k];
                    $sumr = $sumr + log(1+exp($pr-$sumr));
                    $suml = $suml + log(1+exp($pl-$suml));
                }
                $node->{"mat"}->[$i][$j] = $sumr+$suml;
            }
        }
    }else{ #external node
        my @s = split("",$seqs->{$node->{"id"}});
        my @t = @{transarrayCodon(\@s,$codoni)};
        for(my $i=0; $i < scalar(@{$node->{"mat"}});$i++){
            my $val=log(1);
            if($t[$i] ne "NA"){ #adjust for ambiguous sites
                $node->{"mat"}->[$i][$t[$i]]=$val;
            }else{
                #fill in equilibrium frequencies for ambiguous state
                for(my $j = 0; $j < 61; $j++){
                  my $val;
                  if($node->{"level"} == 1){
                    if($rootpis==1){
                      $val = $freqs->{$codons->[$j]};
                    }elsif($rpartIndex->[$i] == 0){
                      #$val = $rootFWR->[$j];
                      #$val = $freqs->{$codons->[$j]};
                      $val=$rootUNK->[$j];
                    }else{
                      #$val = $rootCDR->[$j];
                      $val=$rootUNK->[$j];
                    }
                  }else{
                    if(!exists($ambig_char->{$node->{"id"}}->{$i}->[$j])){
                      print $node->{"id"}."\t$i\t$j\n";
                      print $ambig_char->{$node->{"id"}}."\n";
                      print $ambig_char->{$node->{"id"}}->{$i}."\n";
                      print $ambig_char->{$node->{"id"}}->{$i}->[$j]."\n";
                      die();
                    }
                    $val = $ambig_char->{$node->{"id"}}->{$i}->[$j];
                  }
                  if($val == 0){$val=-200}
                  else{$val = log($val)}
                  $node->{"mat"}->[$i][$j] =  $val;
                }
            }
        }
    }
   my $lhood=0;
   for(my $i=0; $i < $length;$i++){
       my $slhood = $node->{"mat"}->[$i]->[0];
       for(my $j=1;$j<61;$j++){
           $slhood = $slhood + log(1+exp($node->{"mat"}->[$i]->[$j]-$slhood));
       }
       $lhood += $slhood;
   }
   if($node->{"level"} != 0){
       print $node->{"id"}."\t".$node->{"up"}->{"id"}."\t$lhood\n";
     }else{
       print $node->{"id"}."\t"."NONE"."\t$lhood\n";
    }
    return($lhood);
}


sub Pruning_Lhood_g{
    my $node = $_[0];
    my $seqs = $_[1];
    my $codoni = $_[2];
    my $Qs = $_[3]; #array reference of q matrixes
    my $codons = $_[4];

    my @keys = keys %$seqs;
    my $length = length($seqs->{$keys[0]});
    my $nstates = scalar(@$codons);
    my @bigempty;
    for(my $i=0; $i < $length; $i++){
        my @empty = (-200)x61;
        push(@bigempty,\@empty);
    }
    $node->{"mat"}=\@bigempty;

    #if(!exists($node->{"id"})){print $node->{"subtaxa"}."\n";die("node ID not found :-/ Check that all IDs in fasta file are found in the tree\n")}
    if(exists($node->{"left"})){ #internal node
        my $r = Pruning_Lhood_g($node->{"right"},$seqs,$codoni,$Qs,$codons);
        my $l = Pruning_Lhood_g($node->{"left"},$seqs,$codoni,$Qs,$codons);
        
        my $Prs = mexp($Qs*$node->{"right"}->{"dist"});
        my $Pls = mexp($Qs*$node->{"left"}->{"dist"});
        
        for(my $i=0; $i < $length;$i++){
            for(my $j=0;$j<$nstates;$j++){
                my $sumr = log($Prs->at(0,$j))+$node->{"right"}->{"mat"}->[$i][0];
                my $suml = log($Pls->at(0,$j))+$node->{"left"}->{"mat"}->[$i][0];
                for(my $k=1;$k<$nstates;$k++){
                    if($Prs->at($k,$j)==0){print $node->{"dist"}." $k $j\n"}
                    my $pr = log($Prs->at($k,$j)) + $node->{"right"}->{"mat"}->[$i][$k];
                    my $pl = log($Pls->at($k,$j)) + $node->{"left"}->{"mat"}->[$i][$k];
                    $sumr = $sumr + log(1+exp($pr-$sumr));
                    $suml = $suml + log(1+exp($pl-$suml));
                }
                $node->{"mat"}->[$i][$j] = $sumr+$suml;
            }
        }
    }else{ #external node
        my @s = split("",$seqs->{$node->{"id"}});
        my @t = @{transarrayChar(\@s,$codoni)};
        for(my $i=0; $i < scalar(@{$node->{"mat"}});$i++){
            my $val=log(1);
            if($t[$i] ne "NA"){ #adjust for ambiguous sites
                $node->{"mat"}->[$i][$t[$i]]=$val;
            }else{
                die("ambiguous states not allowed yet!");
            }
        }
    }
   my $lhood=0;
   for(my $i=0; $i < $length;$i++){
       my $slhood = $node->{"mat"}->[$i]->[0];
       for(my $j=1;$j<61;$j++){
           $slhood = $slhood + log(1+exp($node->{"mat"}->[$i]->[$j]-$slhood));
       }
       $lhood += $slhood;
   }
   if($node->{"level"} != 0){
    #   print $node->{"id"}."\t".$node->{"up"}->{"id"}."\t$lhood\n";
     }else{
     #  print $node->{"id"}."\t"."NONE"."\t$lhood\n";
    }
    return($lhood);
}

#Make Q matrix for HLP16
sub getQmat_HLP16{
    my $bij = $_[0];
    my $kappa =$_[1];
    my $omega = $_[2];
    my %freqs = %{$_[3]};
    my @codons = @{$_[4]};
    my $print = $_[5];
    my %tr = %{codonTable()};
    my %q;
    my $small = 1e-200;
    for(my $i=0; $i < scalar(@codons); $i++){
        my $from = $codons[$i];
        for(my $j=0; $j < scalar(@codons); $j++){
            my $to = $codons[$j];
            if($from eq $to){$q{$from.$to}=0;next;}
            my @diff = @{diffPos($from,$to)};
            if(scalar(@diff) > 1){
                $q{$from.$to}=$small;
            }else{
                my @dc = sort {$a cmp $b} (substr($from,$diff[0],1),substr($to,$diff[0],1));
                if("@dc" eq "a g" || "@dc" eq "c t"){
                    if($tr{$from} eq $tr{$to}){
                        $q{$from.$to} = $kappa*$freqs{$to}*(1+$bij->[$i*61+$j]);
                    }else{
                        $q{$from.$to} = $omega*$kappa*$freqs{$to}*(1+$bij->[$i*61+$j]);
                    }
                }
                else{
                    if($tr{$from} eq $tr{$to}){
                        $q{$from.$to} = $freqs{$to}*(1+$bij->[$i*61+$j]);
                    }else{
                        $q{$from.$to} = $omega*$freqs{$to}*(1+$bij->[$i*61+$j]);
                    }
                }
            }
        }
    }
    my $trate = 0;
    foreach my $from (@codons){
        my $fsum = 0;
        foreach my $to (@codons){
            $fsum += $q{$from.$to};
        }
        $q{$from.$from}=-$fsum;
        $trate += $freqs{$from}*$fsum;
    }
    if($print){ print "$trate\n";}
    #check that rows sum to 0
    foreach my $from (@codons){
        my $fsum = 0;
        foreach my $to (@codons){
            $fsum += $q{$from.$to};
        }
        if($fsum > 0.0001){die("fsum > 0.0001f irst $from $fsum\n")}
    }
    #scale to a relative rate of 1
    my $rate = 0;
    foreach my $from (@codons){
        foreach my $to (@codons){
            $q{$from.$to} = $q{$from.$to}/$trate;
        }
        $rate -= $freqs{$from}*$q{$from.$from};
    }
    if($print){print "Mean rate: ".$rate."\n";} 
    #check that rows sum to 0
    foreach my $from (@codons){
        my $fsum = 0;
        foreach my $to (@codons){
            $fsum += $q{$from.$to};
        }
        if($fsum > 0.0001){die("$from $fsum\n")}
    }
 
    my @pdls;
    foreach my $from (@codons){
        my @row;
        foreach my $to (@codons){
            push(@row,$q{$from.$to})
        }
        push(@pdls,pdl([@row]));
    }
    my $Q = $pdls[0];
    for(my $i = 1; $i < 61; $i++){
        $Q = $Q->glue(1,$pdls[$i]);
    }
    return($Q);
}

#Make Q matrix for HLP16
sub getQmat_SE_HLP16{
    my $bij = $_[0];
    my $kappa =$_[1];
    my $aaint = $_[2];
    my %freqs = %{$_[3]};
    my @codons = @{$_[4]};
    my $print = $_[5];
    my %tr = %{codonTable()};
    my $aaslope = $_[6];
    my $aatable = $_[7];
    my %q;
    my $small = 1e-200;
    for(my $i=0; $i < scalar(@codons); $i++){
        my $from = $codons[$i];
        for(my $j=0; $j < scalar(@codons); $j++){
            my $to = $codons[$j];
            if($from eq $to){$q{$from.$to}=0;next;}
            my @diff = @{diffPos($from,$to)};
            if(scalar(@diff) > 1){
                $q{$from.$to}=$small;
            }else{
                my @dc = sort {$a cmp $b} (substr($from,$diff[0],1),substr($to,$diff[0],1));
                if("@dc" eq "a g" || "@dc" eq "c t"){
                    if($tr{$from} eq $tr{$to}){
                        $q{$from.$to} = $kappa*$freqs{$to}*(1+$bij->[$i*61+$j]);
                    }else{
                        if($aatable->[$i*61+$j] < 0){die($aatable->[$i*61+$j])}
                        $q{$from.$to} = ($aaint*exp($aaslope*log($aatable->[$i*61+$j])))*$kappa*$freqs{$to}*(1+$bij->[$i*61+$j]);
                    }
                }
                else{
                    if($tr{$from} eq $tr{$to}){
                        $q{$from.$to} = $freqs{$to}*(1+$bij->[$i*61+$j]);
                    }else{
                        if($aatable->[$i*61+$j] < 0){die($aatable->[$i*61+$j])}
                        $q{$from.$to} = ($aaint*exp($aaslope*log($aatable->[$i*61+$j])))*$freqs{$to}*(1+$bij->[$i*61+$j]);
                    }
                }
            }
        }
    }
    my $trate = 0;
    foreach my $from (@codons){
        my $fsum = 0;
        foreach my $to (@codons){
            $fsum += $q{$from.$to};
        }
        $q{$from.$from}=-$fsum;
        $trate += $freqs{$from}*$fsum;
    }
    if($print){ print "$trate\n";}
    #check that rows sum to 0
    foreach my $from (@codons){
        my $fsum = 0;
        foreach my $to (@codons){
            $fsum += $q{$from.$to};
        }
        if($fsum > 0.0001){die("fsum > 0.0001f irst $from $fsum\n")}
    }
    #scale to a relative rate of 1
    my $rate = 0;
    foreach my $from (@codons){
        foreach my $to (@codons){
            $q{$from.$to} = $q{$from.$to}/$trate;
        }
        $rate -= $freqs{$from}*$q{$from.$from};
    }
    if($print){print "Mean rate: ".$rate."\n";} 
    #check that rows sum to 0
    foreach my $from (@codons){
        my $fsum = 0;
        foreach my $to (@codons){
            $fsum += $q{$from.$to};
        }
        if($fsum > 0.0001){die("$from $fsum\n")}
    }
 
    my @pdls;
    foreach my $from (@codons){
        my @row;
        foreach my $to (@codons){
            push(@row,$q{$from.$to})
        }
        push(@pdls,pdl([@row]));
    }
    my $Q = $pdls[0];
    for(my $i = 1; $i < 61; $i++){
        $Q = $Q->glue(1,$pdls[$i]);
    }
    return($Q);
}

#Mutate ancestral sequence according a Q matrix
sub simulate_HLP16_br{
  my $node = $_[0];
  my $aseq = $_[1];
  my $Qs = $_[2];
  my $part = $_[3];
  my $npart = $_[4];
  my $nqmats = $_[5];
  my @Ps;
  for(my $i=0;$i<$nqmats;$i++){
    push(@Ps,mexp($Qs->[$i]*$node->{"dist"}));
  }
  my @seq;
  my @aseq = @$aseq;
  if($node->{"level"} != 0){ #root node is same as supplied root
    print "olab:\t".$node->{"olab"}."\n";
    for(my $i = 0; $i < scalar(@aseq); $i++){
      print $part->[$i]+$npart*$node->{"olab"};
      my $char = $aseq[$i];
      my $rowsum = 0;
      for(my $ind = 0; $ind < 61; $ind++){
        #print "$part->[$i]\t$i\t$ind\t$char\n";
        $rowsum+=$Ps[$part->[$i]+$npart*$node->{"olab"}]->at($ind,$char);
      }
      if($rowsum < 0.99){die("$rowsum\t$char\t".$node->{"dist"}."\n")}
      my $n = rand($rowsum);
      my $sum = 0;
      for(my $ind = 0; $ind < 61; $ind++){
        $sum += $Ps[$part->[$i]+$npart*$node->{"olab"}]->at($ind,$char);
        if($sum >= $n){
          $seq[$i] = $ind;
          last;
        }
      }
    }
    print "\n";
  }else{
    for(my $i = 0; $i < scalar(@aseq); $i++){
      $seq[$i]=$aseq[$i];
    }
  }
  $node->{"sequence"} = \@seq;
  if(exists($node->{"left"})){
    simulate_HLP16_br($node->{"left"},\@seq,$Qs,$part,$npart,$nqmats);
    simulate_HLP16_br($node->{"right"},\@seq,$Qs,$part,$npart,$nqmats);
  }
}

#Mutate ancestral sequence according a Q matrix
sub simulate_HLP16{
  my $node = $_[0];
  my $aseq = $_[1];
  my $Qs = $_[2];
  my $part = $_[3];
  my $npart = $_[4];
  my $rate= $_[5];
  $rate = defined($rate)?$rate:1;
  my @Ps;
  for(my $i=0;$i<$npart;$i++){
    push(@Ps,mexp($Qs->[$i]*$node->{"dist"}*$rate));
  }
  my @seq;
  my @aseq = @$aseq;
  for(my $i = 0; $i < scalar(@aseq); $i++){
    my $char = $aseq[$i];
    my $rowsum = 0;
    for(my $ind = 0; $ind < 61; $ind++){
      #print "$part->[$i]\t$i\t$ind\t$char\n";
      $rowsum+=$Ps[$part->[$i]]->at($ind,$char);
    }
    if($rowsum < 0.99){die("$rowsum\t$char\t".$node->{"dist"}."\n")}
    my $n = rand($rowsum);
    my $sum = 0;
    for(my $ind = 0; $ind < 61; $ind++){
      $sum += $Ps[$part->[$i]]->at($ind,$char);
      if($sum >= $n){
        $seq[$i] = $ind;
        last;
      }
    }
  }
  $node->{"sequence"} = \@seq;
  if(exists($node->{"left"})){
    simulate_HLP16($node->{"left"},\@seq,$Qs,$part,$npart,$rate);
    simulate_HLP16($node->{"right"},\@seq,$Qs,$part,$npart,$rate);
  }
}

#Mutate ancestral sequence according a Q matrix
sub simulate_1_omega{
  my $node = $_[0];
  my $aseq = $_[1];
  my $Q = $_[2];
  my $P = mexp($Q*$node->{"dist"});
  my @seq;
  my @aseq = @$aseq;
  for(my $i = 0; $i < scalar(@aseq); $i++){
    my $char = $aseq[$i];
    my $rowsum = 0;
    for(my $ind = 0; $ind < 61; $ind++){
      $rowsum+=$P->at($ind,$char);
    }
    if($rowsum < 0.99){die("$rowsum\t$char\t".$node->{"dist"}."\n")}
    my $n = rand($rowsum);
    my $sum = 0;
    for(my $ind = 0; $ind < 61; $ind++){
      $sum += $P->at($ind,$char);
      if($sum >= $n){
        $seq[$i] = $ind;
        last;
      }
    }
  }
  $node->{"sequence"} = \@seq;
  if(exists($node->{"left"})){
    simulate_1_omega($node->{"left"},\@seq,$Q);
    simulate_1_omega($node->{"right"},\@seq,$Q);
  }
}

#Mutate ancestral sequence according a Q matrix
sub simulate_1_omega_full_context{
  my $node = $_[0];
  my $aseq = $_[1];
  my $hs = $_[2];
  my $omega = $_[3];
  my $kappa = $_[4];
  my $freqs = $_[5];
  my $tr = $_[6];
  my $codons = $_[7];
  my $rmotifs = $_[8];
  my $lmotifs = $_[9];
 
  my @seq;
  my @aseq = @$aseq;
  for(my $i = 0; $i < scalar(@aseq); $i++){
    #get Bij values    
    my $modeqs;
    if($i == 0){
        $modeqs = makeBmat_known_left($hs,$rmotifs,$lmotifs,$codons,$freqs,$codons->[$aseq[$i+1]]);
    }elsif($i==(scalar(@aseq)-1)){
        $modeqs = makeBmat_known_right($hs,$rmotifs,$lmotifs,$codons,$freqs,$codons->[$aseq[$i-1]]);
    }else{
        $modeqs = makeBmat_known($hs,$rmotifs,$lmotifs,$codons,$codons->[$aseq[$i-1]],$codons->[$aseq[$i+1]]);
    }
    #make new Q matrix given the nucleotide context
    my $Q = getQmat_HLP16($modeqs,$kappa,$omega,$freqs,$codons,0);
    #exponentiate
    my $P = mexp($Q*$node->{"dist"});

    my $char = $aseq[$i];
    my $rowsum = 0;
    for(my $ind = 0; $ind < 61; $ind++){
      $rowsum+=$P->at($ind,$char);
    }
    if($rowsum < 0.99){die("$rowsum\t$char\t".$node->{"dist"}."\n")}
    my $n = rand($rowsum);
    my $sum = 0;
    for(my $ind = 0; $ind < 61; $ind++){
      $sum += $P->at($ind,$char);
      if($sum >= $n){
        $seq[$i] = $ind;
        last;
      }
    }
  }
  $node->{"sequence"} = \@seq;
 # print "done\n";
  if(exists($node->{"left"})){
    simulate_1_omega_full_context($node->{"left"},\@seq,$hs,$omega,$kappa,$freqs,$tr,$codons,$rmotifs,$lmotifs);
    simulate_1_omega_full_context($node->{"right"},\@seq,$hs,$omega,$kappa,$freqs,$tr,$codons,$rmotifs,$lmotifs);
  }
}

sub makeBmat_freq{
    my $hotness = $_[0];
    my $rhotspot = $_[1];
    my $lhotspot = $_[2];
    my @codons = @{$_[3]};
    my $freqs = $_[4];
 
    #Modify frequencies
    my @res;
    foreach my $from (@codons){
        foreach my $to (@codons){
            #No modification if more than one mutation (will equal 0 in the Q matrix anyway)
            my @diff = @{diffPos($from,$to)};
            if(scalar(@diff) != 1){
                push(@res,0);
                next;
            }
            #what is the probability that this mutation occurred in a hotspot motif?
            my $htotal = 0; 
            foreach my $left (@codons){
                foreach my $right (@codons){
                    my $h = 0;
                    my @diff = @{diffPos($left.$from.$right,$left.$to.$right)};
                    my @c = split("",$left.$from.$right);
                    my $l = $c[$diff[0]].$c[$diff[0]+1].$c[$diff[0]+2];
                    my $r = $c[$diff[0]-2].$c[$diff[0]-1].$c[$diff[0]];         
                    if($rhotspot->{uc $r} == 1 or $lhotspot->{uc $l}==1){$h=1}
                    $htotal += $freqs->{$left}*$freqs->{$right}*$h
                }
            }
            #if($from eq "gta" || $from eq "gtt"){print "$from $to $htotal $freqs->{$to}\n";}
             push(@res,$htotal*$hotness);
        }
    }
    return(\@res);
}

sub makeBmat_known{
    my $hotness = $_[0];
    my $rhotspot = $_[1];
    my $lhotspot = $_[2];
    my @codons = @{$_[3]};
    my $left = $_[4];
    my $right = $_[5];

    #Modify frequencies
    my @res;
    foreach my $from (@codons){
        foreach my $to (@codons){
            #No modification if more than one mutation (will equal 0 in the Q matrix anyway)
            my @diff = @{diffPos($from,$to)};
            if(scalar(@diff) != 1){
                push(@res,0);
                next;
            }
            #what is the probability that this mutation occurred in a hotspot motif?
            my $h = 0;
            @diff = @{diffPos($left.$from.$right,$left.$to.$right)};
            my @c = split("",$left.$from.$right);
            my $l = $c[$diff[0]].$c[$diff[0]+1].$c[$diff[0]+2];
            my $r = $c[$diff[0]-2].$c[$diff[0]-1].$c[$diff[0]];         
            if($rhotspot->{uc $r} == 1 or $lhotspot->{uc $l}==1){$h=1}
            push(@res,$h*$hotness);
        }
    }
    return(\@res);
}

sub makeBmat_known_left{
    my $hotness = $_[0];
    my $rhotspot = $_[1];
    my $lhotspot = $_[2];
    my @codons = @{$_[3]};
    my $freqs = $_[4];
    my $right = $_[5];
 #   my $right = $codons[$rightc];

    #Modify frequencies
    my @res;
    foreach my $from (@codons){
        foreach my $to (@codons){
            #No modification if more than one mutation (will equal 0 in the Q matrix anyway)
            my @diff = @{diffPos($from,$to)};
            if(scalar(@diff) != 1){
                push(@res,0);
                next;
            }
            #what is the probability that this mutation occurred in a hotspot motif?
            my $htotal = 0; 
            foreach my $left (@codons){
             #   foreach my $right (@codons){
                    my $h = 0;
                    my @diff = @{diffPos($left.$from.$right,$left.$to.$right)};
                    my @c = split("",$left.$from.$right);
                    my $l = $c[$diff[0]].$c[$diff[0]+1].$c[$diff[0]+2];
                    my $r = $c[$diff[0]-2].$c[$diff[0]-1].$c[$diff[0]];         
                    if($rhotspot->{uc $r} == 1 or $lhotspot->{uc $l}==1){$h=1}
                    $htotal += $freqs->{$left}*$h
            #  }
            }
             push(@res,$htotal*$hotness);
        }
    }
    return(\@res);
}

sub makeBmat_known_right{
    my $hotness = $_[0];
    my $rhotspot = $_[1];
    my $lhotspot = $_[2];
    my @codons = @{$_[3]};
    my $freqs = $_[4];
    my $left = $_[5];
  #  my $left = $codons[$leftc];

    #Modify frequencies
    my @res;
    foreach my $from (@codons){
        foreach my $to (@codons){
            #No modification if more than one mutation (will equal 0 in the Q matrix anyway)
            my @diff = @{diffPos($from,$to)};
            if(scalar(@diff) != 1){
                push(@res,0);
                next;
            }
            #what is the probability that this mutation occurred in a hotspot motif?
            my $htotal = 0; 
           # foreach my $left (@codons){
                foreach my $right (@codons){
                    my $h = 0;
                    my @diff = @{diffPos($left.$from.$right,$left.$to.$right)};
                    my @c = split("",$left.$from.$right);
                    my $l = $c[$diff[0]].$c[$diff[0]+1].$c[$diff[0]+2];
                    my $r = $c[$diff[0]-2].$c[$diff[0]-1].$c[$diff[0]];         
                    if($rhotspot->{uc $r} == 1 or $lhotspot->{uc $l}==1){$h=1}
                    $htotal += $freqs->{$right}*$h
              }
            #}
             push(@res,$htotal*$hotness);
        }
    }
    return(\@res);

}

1