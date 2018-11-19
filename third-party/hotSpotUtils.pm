#Utilities for doing analysis with BCR phylogenetics, especially with regards to 
#hotspot selection models
#20/Feb/2016
#Kenneth Hoehn

use strict;
use warnings;

sub translate{
  my $s = $_[0];
  my $ct3 = codonTable();
  my $ct1 = codonTableSingle();
  my $aa = "";
  #if(length($s)%3==0){
    for(my$i=0;$i<length($s);$i+=3){
      my $ss = substr($s,$i,3);
      if(exists($ct3->{lc $ss})){
        if(!exists($ct1->{uc $ct3->{lc $ss}})){
          $aa .= '*';
        }else{
          my $val = $ct1->{uc $ct3->{lc $ss}};
          if($ss =~ /[a-z]/){
            $val = lc $val;
          }
          $aa .= $val;
        }
      }elsif($ss =~ /---/){
        $aa .="-";
      }elsif($ss =~ /[A-Z]--/){
        my $c = "-";my $counter=$i+2;#move forward and fix the codon
        my $gaps = 1;
        while($c eq "-"){
          $counter++;$gaps++;
          $c=substr($s,$counter,1);
        }
        if($counter == length($s)){last;}
        $s = substr($s,0,$i+1).substr($s,$counter,2).substr($s,$i+1,$gaps).substr($s,$counter+2);
        $i-=3;
      }elsif($ss =~ /[A-Z][A-Z]\-/){
        my $c = "-";my $counter=$i+2;#move forward and fix the codon
        my $gaps = 0;
        while($c eq "-"){
          $counter++;$gaps++;
          $c=substr($s,$counter,1);
        }
        if($counter == length($s)){last;}
        #print "$i\t$counter\t".length($s)."\n";
        $s = substr($s,0,$i+2).substr($s,$counter,1).substr($s,$i+2,$gaps).substr($s,$counter+1);
        $i-=3;
      }elsif($ss =~ /[A-Z]-[A-Z]/){#fix split codon
        $s = substr($s,0,$i+1).substr($s,$i+2,2)."-".substr($s,$i+4);
        $i-=3;
      }
      elsif($ss =~ /-[A-Z][A-Z]/){
        $i-=2;
        $aa .= "-";
      }elsif($ss =~ /--[A-Z]/){
        $i-=1;
        $aa .= "-";
      }
      else{
        $aa .= "X";
      }
    }
  return $aa;
}


#Print ML amino acid sequence at all nodes
sub aaEachSite{
  my $node = $_[0];
  my $string = $_[1];
  my $codons = $_[2];
  my $aatable1 = $_[3];
  my $aatable3 = $_[4];
  if(exists($node->{"left"})){
    aaEachSite($node->{"left"}, $string,$codons,$aatable1,$aatable3);
    aaEachSite($node->{"right"},$string,$codons,$aatable1,$aatable3);
  }

  if($node->{"level"}!=0){
  my @sequence;
    for(my $i=0; $i < scalar(@{$node->{"Codon_lhoods"}});$i++){
      my %aas;
      for(my $j=0; $j < 61; $j++){ #tally up relative likelihoods of amino acids
        my $n = $node->{"Codon_lhoods"}->[$i][$j];
        my $aa =  $aatable3->{$codons->[$j]}; #triple letter amino acid
        if(exists($aas{$aa})){
          $aas{$aa} = $aas{$aa} + log(1+exp($n - $aas{$aa}))
        }else{
          $aas{$aa} = $n
        }
      }
      my $max = -inf;
      my $maxchar = -1;
      foreach my $key (keys %aas){
        my $n = $aas{$key};
        if($n > $max){
          $max = $n; 
          $maxchar=$key;
        }
      }
      $node->{"AA_$i"} = $maxchar;
    }
    }
}

#print out for figtree
sub printtreestring_AAeach_figtree{
  my $node = $_[0];
  my $seqsf = $_[1];
  my $file = $_[2];
  my $slength = $_[3];
  my $alrt = $_[4];
  if(!defined($alrt)){$alrt=0;}

  my $seqs = getfasta($seqsf);
  my @keys = keys %$seqs;
  open(FIG,">$file") or die($file);
  print FIG "#NEXUS\nBegin taxa;\n\tDimensions ntax=".scalar(@keys).";\n\tTaxlabels\n";
  foreach my $k (@keys){
    print FIG "\t\t$k\n";
  }
  print FIG "\t\t;\nEnd;\nBegin trees;Tree TREE1 = [\&R] ";
  my $t="";
  my $tree = printtreestring_AAeach($node,$t,$slength,$alrt);
  print FIG "$tree;\nEnd;";
}

#Print out the tree in Newick format to a string
sub printtreestring_AAeach{
  my $node = $_[0];
  my $string = $_[1];
  my $slength = $_[2];
  my $alrt = $_[3];
  if(exists($node->{"left"})){
    $string = $string."(";
    $string=printtreestring_AAeach($node->{"left"},$string,$slength,$alrt);
    $string = $string.",";
    $string=printtreestring_AAeach($node->{"right"},$string,$slength,$alrt);
    $string=$string.")";
  }
  if($node->{"level"} != 0){
    my $a = "";
    if($alrt == 1){
      if($node->{"alrt"} ne ""){
        $string = $string.$node->{"id"}."[\&aLRT=".$node->{"alrt"}.",AA_0\=".$node->{"AA_0"};
      }else{
        $string = $string.$node->{"id"}."[\&AA_0\=".$node->{"AA_0"};
      }
    }else{
      $string = $string.$node->{"id"}."[\&AA_0\=".$node->{"AA_0"};
    }
    for(my $j = 1; $j < $slength; $j++){
      $string = $string.",AA_$j\=".$node->{"AA_$j"};
    }
    $string = $string."]:".$node->{"dist"};
  }
  return($string);
}

#Read in a fasta file
sub getfasta{
  my $in = $_[0];
  chomp($in);
  open(FASTAFILEIN, $in) or die("Couldn't open $in.\n");  # Create a new file
  my %seqs;
  my $id;
  my $s = '';
  while(<FASTAFILEIN>){
    my $line=$_;
    chomp($line);
    if($line =~/^\s*$/) {next;}     
    elsif($line =~/^\s*#/) {next;}  
    elsif($line =~ /^\>\s*(\S+)\s*$/) {
      my $temp = $1;
      $temp =~ s/\s//sg;
      if(length($s) > 0) {
        $s =~ s/\s//sg;
        $seqs{$id} = $s;
        $s = '';
      }
      $id = $temp;
    }else{
     $s .= $line;
   }
  }
  $s =~ s/\s//sg;
  $seqs{$id} = $s;
  close(FASTAFILEIN);
  return \%seqs;
}

#raw distance between sequences of the same length, not counting gaps on either
sub distngap{ 
  my ($one,$two) = @_;
  my $dist = 0;
  $one = uc $one;
  $two = uc $two;
  if(length($one) != length($two)){die("Seqs not of same length! 1st:$one\n2nd:$two\n")}
  for(my $i = 0; $i < length($one); $i++ ){
    if(substr($one,$i,1) ne substr($two,$i,1) && substr($one,$i,1) ne '-' && substr($two,$i,1) ne '-' && substr($one,$i,1) ne 'N' && substr($two,$i,1) ne 'N'){
      $dist++;
    }
  }
  return($dist);
}

sub distngap_av{ 
  my ($one,$two) = @_;
  $one = uc $one;
  $two = uc $two;
  my $dist = 0;
  my $sites=0;
  if(length($one) != length($two)){die("Seqs not of same length! 1st:$one\n2nd:$two\n")}
  for(my $i = 0; $i < length($one); $i++ ){
    if(substr($one,$i,1) ne '-' && substr($two,$i,1) ne '-' && substr($one,$i,1) ne 'N' && substr($two,$i,1) ne 'N'){
      if(substr($one,$i,1) ne substr($two,$i,1)){
        $dist++;
      }
      $sites++;
    }
  }
  my @ret=($dist,$sites);
  return \@ret;
}

#raw distance between sequences of the same length, not counting gaps on either
sub distngap_aa{ 
  my ($one,$two) = @_;
  my $dist = 0;
  if(length($one) != length($two)){die("Seqs not of same length! 1st:$one\n2nd:$two\n")}
  for(my $i = 0; $i < length($one); $i++ ){
    if(substr($one,$i,1) ne substr($two,$i,1) && substr($one,$i,1) ne '-' && substr($two,$i,1) ne '-'){
      $dist++;
    }
  }
  return($dist);
}

#Position of differences between two strings
sub diffPos{
  my $s1 = $_[0];
  my $s2 = $_[1];
  if(length($s1) != length($s2)){die("$s1 $s2 no same length")}
  my @diffpos;
  for(my $i = 0; $i < length($s1); $i++){
    if(substr($s1,$i,1) ne substr($s2,$i,1)){
      push(@diffpos,$i);
    }
  }
  return(\@diffpos);
}

#raw distance between sequences of the same length, not counting gaps
sub ardistngap{ 
  my ($first,$second) = @_;
  my @one = @{$first};
  my @two = @{$second};
  my $dist = 0;
  if(scalar(@one) != scalar(@two)){die("Not of equal length:\n@one\n\n@two")}
  for(my $i = 0; $i < scalar(@one); $i++ ){
     if($one[$i] ne "-" && $two[$i] ne "-" && $one[$i] ne "N" && $two[$i] ne "N"){
      if($one[$i] ne $two[$i]){
        $dist++;
      }
    }
  }
  return($dist);
}


#read in an unrooted newick tree
#and convert to tree rooted at leftmost taxa
sub readInUnrootedNewick{
  my $in = $_[0];
  my $rootname = $_[1];
  my $printtree = $_[2];
  if($in =~ /\($rootname\:(\d+\.*\d*[Ee]*-*\d*),/){
    #extract leftmost taxa as the root
    my $rdist = $1;
    my $next = "(".$';#remove root taxon
    my $oin = $next;

    my %root; 
    my %rright;
    my %rleft;
  
    #read in the rest of the tree as rooted
    $rright{"dist"}=0;
    rootedNewick($next,\%rright,0);
    getdists(\%rright); #parse distance

    my $t = printtreestring(\%rright,"").";";
    if($t ne $oin){print "Tree read in incorrectly!\nOriginal Tree:\n$oin\nRead tree:\n$t\n";die();}
    else{print "Tree read in correctly\n";}
  
    #set up root node
    $root{"dist"} = 0;
    $root{"id"}="";
    $root{"level"} = -1;
    $rright{"dist"} = $rdist;
    $root{"right"} = \%rright;

    $rleft{"dist"}= 0.0000000000000000000001;
    $rleft{"level"}= 0;
    $rleft{"id"}=$rootname;
    $root{"left"} = \%rleft;

    relevel(\%root,1); #increase levels by 1

    getdivergence(\%root,0); #get divergences

    #optionally print out trees to make sure everything is good
    if($printtree){
      $t = printtreestring(\%root,"").";";
      print("$in\n$oin\n$t\n");
    }
    return(\%root);
  }else{die("Input tree not formatted correctly - root is not leftmost taxa.\n");}
}

sub relevel{
  my $node = $_[0];
  my $increase = $_[1];
  if(exists($node->{"right"})){
    relevel($node->{"right"},$increase);
    relevel($node->{"left"},$increase);
  }
  $node->{"level"}=$node->{"level"}+$increase;
}

#Read in a rooted Newick tree
sub rootedNewick {
  my $in = $_[0];
  my $node = $_[1];
  my $level = $_[2];
  my $first; my $id;
  if($in =~ /(,|\(|\)|;)/){
    $first = $&;
    $node->{"id"} = $`;
    $node->{"level"} = $level;
    $in = $';
    if($first eq ","){$in = $first.$in}
  }else{
    die($in);
  }
  if($first eq "("){#left
    my %n;
    $node->{"left"} = \%n;
    $in = rootedNewick($in,\%n,$level+1);
  }
  elsif($first  eq ","){#up
    return($in);
  }
  elsif($first  eq ")"){#up
    return($in);
  }
  elsif($first  eq ";"){#up
    return($in);
  }
  my $second;
  my $pre;
  if($in =~ /(,|\(|\)|;)/){
    $second = $&;
    $in = $';
    $pre = $`;
  }else{
    die($in);
  }
  if($second eq ","){#right
    my %n;
    $node->{"right"} = \%n;
    $in = rootedNewick($in,\%n,$level+1);
  }elsif($second  eq ")"){#up
  }
  if($in =~ /(,|\(|\)|;)/){
    $node->{"id"} = $`;
    $in = $';
    if($& eq ","){$in = $&.$'}
  }else{
    die($in);
  }
  return($in);
}


#read in a rooted newick tree
sub readInRootedNewick{
  my $in = $_[0];
  my $printtree = $_[1];
  my $alrt = $_[2];
  if(!defined($alrt)){$alrt=0}

   my %root; #set up root node
   my $oin = $in;
   $root{"dist"}=0;

   my $nin = "";
   if($alrt == 2){
    while($in =~ /(\[?\&?AA\_\d+\=[A-Z][a-z]+)([,|\]])/){
      $nin.=$`.$1;
      if($2 eq "]"){$nin.="]"}
      else{$nin.="%";}
      $in=$';
    }
    $nin.=$in;
    $in=$nin;
    $nin="";
    while($in =~ /(\[?\&?aLRT\=-?\d+\.\d+)([,|\]])/){
      $nin.=$`.$1;
      if($2 eq "]"){$nin.="]"}
      else{$nin.="%";}
      $in=$';
    }
    $nin.=$in;
    $in=$nin;
    $nin="";
    while($in =~ /(\[?\&?CO\_\d+\=[A-Z]+)([,|\]])/){
      $nin.=$`.$1;
      if($2 eq "]"){$nin.="]"}
      else{$nin.="%";}
      $in=$';
    }
    $nin.=$in;
    $in=$nin;
   }
   if($alrt == 3){
    while($in =~ /(\[\&\!name\=\"1\"\])/){
      $nin.=$`.$1;
      $in=$';
    }
    $nin.=$in;
    $in=$nin;
    $nin="";
   }
   if($alrt == 4){
    while($in =~ /(#\d+)/){
      $nin.=$`.$1;
      $in=$';
    }
    $nin.=$in;
    $in=$nin;
    $nin="";
   }
   if($alrt == 5){
    #print "trying\n";
    while($in =~ /(\[\&pZero\=\d?.?\d+),cZero\=\d+,cOne\=\d+,oLab=\d+,num=\d+]/){
      #print "$1\n";
      $nin.=$`.$1."]";
      $in=$';
    }
    $nin.=$in;
    $in=$nin;
    $nin="";
   }
    if($alrt == 6){
    #print "trying\n";
    while($in =~ /\[\&CDRdiff=-?\d*.?\d+,LR=-?\d*.?\d+,Num=-?\d+,altid=\"\S{1,11}\"\]/){
      #print "$1\n";
      my $m=$&;
      my $post=$';
      my$pre=$`;
      $m=~s/,/|/g;
      $nin.=$pre.$m;
      $in=$post;
    }
    $nin.=$in;
    $in=$nin;
    $nin="";
   }
   if($alrt == 7){
    #print "trying\n";
    while($in =~ /\[\&Num\=(\d+\_\d+)\]/){
      #print "$1\n";
      my $m=$&;
      my $post=$';
      my$pre=$`;
      $m=~s/,/|/g;
      $nin.=$pre.$m;
      $in=$post;
    }
    $nin.=$in;
    $in=$nin;
    $nin="";
   }
   rootedNewick($in,\%root,0);
   getdists(\%root,$alrt); #parse distance
   getdivergence(\%root,0); #get divergences

   my $t = printtreestring(\%root,"",$alrt).";";
   #print "$t\n";
   #if($alrt eq 0){
   # if($t ne $oin){print "Tree read in incorrectly!\nOriginal Tree:\n$oin\nRead tree:\n$t\n";die();}
   # else{print "Tree read in correctly\n";}
   #}
  
   if($printtree){
    print "processed tree:\n";
     $t = printtreestring(\%root,"").";";
     print("$oin\n$t\n");
   }
   return(\%root);
}

#Once tree is read in, need to get the branch lengths
sub getdists{
  my $node = $_[0];
  my $alrt = $_[1];
  if(!defined($alrt)){$alrt=0;}
  if(exists($node->{"left"})){
    getdists($node->{"left"},$alrt);
    getdists($node->{"right"},$alrt);
  }
  if(!exists($node->{"id"})){
    die("Node ID doens't exist!");
  }else{
    if($node->{"id"}=~/\:/){
      $node->{"dist"} = $';
      my $pre = $`;
      if($alrt==1){
        if($pre =~ /^(-?\d+\.?\d+)#(\d+)$/){
          $node->{"alrt"}=$1;
          $node->{"label"}=$2;
          #print "$1\t$2\n";
          $pre="";
        }elsif($pre =~ /^(\S*)#(\d+)$/){
          $pre=$1;
          $node->{"label"}=$2;
          $node->{"alrt"}="";
        }elsif($pre =~ /^(-?\d+\.?\d+)$/){
          $node->{"alrt"}=$1;
          $pre="";
        }else{
          $node->{"alrt"}="";
        }
      }elsif($alrt==2){
        if($pre =~ /(\S*)(\[.*\])/){
          $node->{"aas"}=$2;
          $pre="$1";
        }else{
          $node->{"aas"}="";
        }
      }
      elsif($alrt==3){
        if($pre =~ /(\S*)\[\&\!name\=\"(\d+)\"\]/){
          $node->{"olab"}=$2;
          $pre="$1";
        }else{
          $node->{"olab"}="";
        }
      }
      elsif($alrt==4){
        if($pre =~ /(\S*)#(\d+)/){
          $node->{"olab"}=$2;
          $pre="$1";
        }else{
          $node->{"olab"}="";
        }
      }
      elsif($alrt==5){
        #print "Pre: $pre\n";
        #die();
        if($pre =~ /(\S*)\[\&pZero\=(-?\d?.?\d+)]/){
          $node->{"olab"}=$2;
          $pre="$1";
        }else{
          $node->{"olab"}="";
        }
      }
       elsif($alrt==6){
        #print "Pre: $pre\n";
        #die();
        if($pre =~ /(\S*)(\[\&CDRdiff=-?\d*.?\d+\|LR=-?\d*.?\d+\|Num=-?\d+\|altid=\"\S{1,11}\"\])/){
          $node->{"olab"}=$2;
          $pre="$1";
        }else{
          $node->{"olab"}="";
        }
      }
      elsif($alrt==7){
        if($pre =~ /(\S*)\[\&Num=(\d+\_\d+)\]/){
          $node->{"enum"}=$2;
          $pre="$1";
        }else{
          $node->{"enum"}="";
        }
      }
      $node->{"id"} = $pre;
      #print $node->{"label"}."\n";
      #print "$pre\n";
      if($node->{"dist"} == 0){
        $node->{"dist"} = 0.0000000000000000000001;#print "zero length branch length\n";
      }
    }else{
      if($node->{"level"} != 0){
        die($node->{"id"}." level ".$node->{"level"}." is formatted incorrectly!");
      }
    }
  }
}

sub setUnlabled{
  my $node=$_[0];
  if($node->{"level"} != 0){
    if(!exists($node->{"olab"})){
      $node->{"olab"}=0;
    }elsif($node->{"olab"} eq ""){
        $node->{"olab"}=0;
      }
  }
  if(exists($node->{"left"})){
    setUnlabled($node->{"left"});
    setUnlabled($node->{"right"});
  }
}

#re-root the tree at a given node
sub reRoot{
  my $oldroot= $_[0]; #current root
  my $new = $_[1]; #id of new root

  my $old = $oldroot->{"left"};
  #get reference to new root node
  my $root = getnode($oldroot,$new);

  #add blank zero length node to new root
  my %blank;
  $blank{"right"}=$root->{"up"}; #tree is right of new
  $blank{"left"}=$root;  #root is left of new
  $root->{"up"}=\%blank;  #blank is up of root
  #re-assign blank to either left or right of tree, depending on where root was
  if(exists($blank{"right"}->{"left"})){
    print $blank{"right"}->{"left"}."\t$root\n";
    if($blank{"right"}->{"left"} eq $root){ 
      print "new root on left\n";
      $blank{"right"}->{"left"} = \%blank;
    }
  }
  if(exists($blank{"right"}->{"right"})){
    if($blank{"right"}->{"right"} eq $root){
      print "new root on right\n";
      $blank{"right"}->{"right"} = \%blank;
    }
  }
  my $dist = $root->{"dist"};
  $blank{"dist"}=$dist;
  $root->{"dist"}=0.0000000000000000000001;

  $root->{"level"}=1;
  $blank{"level"}=0;
  $blank{"id"}="blank";
  $blank{"subtaxa"}="NONE";
  
  #remove zero-length node from old root
  $old->{"up"}->{"right"}->{"up"}=$old;
  $old->{"dist"}=$oldroot->{"right"}->{"dist"};

  #change, left, right, and anc branches for each node
  reRoot_recurse($blank{"right"},\%blank,0,$dist);

  #re-calculate divergence and subtaxa
  getdivergence(\%blank,0); #get divergences

  return(\%blank);
}

sub reRoot_recurse{
  my $node = $_[0];
  my $anc = $_[1];
  my $level = $_[2];
  my $dist = $_[3];

  print $node->{"subtaxa"}."\t".$anc->{"subtaxa"}."\t".$node->{"dist"}."\n";

  my $oanc = $node->{"up"};
  $node->{"level"}=$level;
  $level++;
  my $ndistl; my $ndistr;
  if(exists($node->{"left"})){#if internal node
    my $oanc = $node->{"up"};
    my $ol = $node->{"left"};
    my $or = $node->{"right"};
    if($node->{"left"} eq $anc){
      $node->{"right"}=$oanc;
      $node->{"left"}=$or;
      $ndistr=$node->{"dist"};
      $ndistl=$node->{"left"}->{"dist"};
    }elsif($node->{"right"} eq $anc){
      $node->{"left"}=$oanc;
      $node->{"right"}=$ol;
      $ndistl=$node->{"dist"};
      $ndistr=$node->{"right"}->{"dist"};
    }elsif($node->{"up"} eq $anc){
      $node->{"left"}=$ol;
      $node->{"right"}=$or;
      $ndistl=$node->{"left"}->{"dist"};
      $ndistr=$node->{"right"}->{"dist"};
    }else{die("something weird happened")}
    reRoot_recurse($node->{"left"},$node,$level,$ndistl);
    reRoot_recurse($node->{"right"},$node,$level,$ndistr);
  }
  $node->{"dist"}=$dist;
  $node->{"up"}=$anc;
}


#re-root the tree at a given node
sub reRoot_internal{
  my $oldroot= $_[0]; #current root
  my $root = $_[1]; #hash reference of new root

  my $old = $oldroot->{"left"};

  #add blank zero length node to new root
  my %blank;
  $blank{"right"}=$root->{"up"}; #tree is right of new
  $blank{"left"}=$root;  #root is left of new
  $root->{"up"}=\%blank;  #blank is up of root
  #re-assign blank to either left or right of tree, depending on where root was
  if(exists($blank{"right"}->{"left"})){
    print $blank{"right"}->{"left"}."\t$root\n";
    if($blank{"right"}->{"left"} eq $root){ 
      print "new root on left\n";
      $blank{"right"}->{"left"} = \%blank;
    }
  }
  if(exists($blank{"right"}->{"right"})){
    if($blank{"right"}->{"right"} eq $root){
      print "new root on right\n";
      $blank{"right"}->{"right"} = \%blank;
    }
  }
  my $dist = $root->{"dist"};
  $blank{"dist"}=$dist;
  $root->{"dist"}=0.0000000000000000000001;

  $root->{"level"}=1;
  $blank{"level"}=0;
  $blank{"id"}="blank";

  #remove zero-length node from old root
  $old->{"up"}->{"right"}->{"up"}=$old;
  $old->{"dist"}=$oldroot->{"right"}->{"dist"};
  $blank{"subtaxa"}="NONE";

  #change, left, right, and anc branches for each node
  reRoot_recurse($blank{"right"},\%blank,0,$dist);
  reRoot_recurse($blank{"left"},\%blank,0,0.0000000000000000000001);

  #re-calculate divergence and subtaxa
  return(\%blank);
}



#Get the divergences for each node in the tree
sub getdivergence{
  my $node = $_[0];
  my $div = $_[1];
  $node->{"divergence"} = $div + $node->{"dist"};
  if(exists($node->{"left"})){
    getdivergence($node->{"left"},$node->{"divergence"});
    getdivergence($node->{"right"},$node->{"divergence"});
  }
}

#Print out the tree in Newick format to a string
sub printtreestring{
  my $node = $_[0];
  my $string = $_[1];
  my $alrt = $_[2];
  if(!defined($alrt)){$alrt=0}
  if(exists($node->{"left"})){
    $string = $string."(";
    $string=printtreestring($node->{"left"},$string,$alrt);
    $string = $string.",";
    $string=printtreestring($node->{"right"},$string,$alrt);
    $string=$string.")";
  }
  if($node->{"level"} != 0){
    if($alrt == 0){
      $string = $string.$node->{"id"}.":".$node->{"dist"};
    }elsif($alrt == 1){
      $string = $string.$node->{"id"}.$node->{"alrt"}.":".$node->{"dist"};
    }else{
        $string = $string.$node->{"id"}.":".$node->{"dist"};
    }
  }
  return($string);
}

sub printtreestring_label{
  my $node = $_[0];
  my $string = $_[1];
  my $alrt = $_[2];
  if(!defined($alrt)){$alrt=0}
  if(exists($node->{"left"})){
    $string = $string."(";
    $string=printtreestring_label($node->{"left"},$string,$alrt);
    $string = $string.",";
    $string=printtreestring_label($node->{"right"},$string,$alrt);
    $string=$string.")";
  }
  if($node->{"id"} ne ""){$node->{"id"} .= "_".$node->{"sequence"}->[0];}
  if($node->{"level"} != 0){
    if($alrt == 0){
      $string = $string.$node->{"id"}.":".$node->{"dist"};
    }elsif($alrt == 1){
      #$string = $string.$node->{"id"}."_".$node->{"sequence"}->[0].$node->{"alrt"}.":".$node->{"dist"};
      die();
    }else{
        $string = $string.$node->{"id"}.":".$node->{"dist"};
    }
  }
  return($string);
}

#Print out the tree in Newick format to a string
sub printtreestring_freq{
  my $node = $_[0];
  my $string = $_[1];
  if(exists($node->{"left"})){
    $string = $string."(";
    $string=printtreestring_freq($node->{"left"},$string);
    $string = $string.",";
    $string=printtreestring_freq($node->{"right"},$string);
    $string=$string.")";
  }
  if($node->{"level"} != 0){
    $string = $string.$node->{"id"}."[&WRC=".$node->{"WRC_c"}.",WRC_count=".$node->{"WRC_count"}.",GYW=".$node->{"GYW_c"}.",GYW_count=".$node->{"GYW_count"}.",WA=".$node->{"WA_c"}.",WA_count=".$node->{"WA_count"}.",TW=".$node->{"TW_c"}.",TW_count=".$node->{"TW_count"}.",SYC=".$node->{"SYC_c"}.",SYC_count=".$node->{"SYC_count"}.",GRS=".$node->{"GRS_c"}.",GRS_count=".$node->{"GRS_count"}.",WRCGYW=".($node->{"WRC_c"}+$node->{"GYW_c"}).",WATW=".($node->{"WA_c"}+$node->{"TW_c"}).",SYCGRS=".($node->{"SYC_c"}+$node->{"GRS_c"})."]:".$node->{"dist"};
  }
  return($string);
}


#print out for figtree
sub printtreestring_figtree{
  my $node = $_[0];
  my $seqsf = $_[1];
  my $file = $_[2];

  my $seqs = getfasta($seqsf);
  my @keys = keys %$seqs;
  open(FIG,">$file") or die($file);
  print FIG "#NEXUS\nBegin taxa;\n\tDimensions ntax=".scalar(@keys).";\n\tTaxlabels\n";
  foreach my $k (@keys){
    print FIG "\t\t$k\n";
  }
  print FIG "\t\t;\nEnd;\nBegin trees;Tree TREE1 = [\&R] ";
  my $t="";
  my $tree = printtreestring_freq($node,$t);
  print FIG "$tree;\nEnd;";

}

#Calculate hotspot frequency at each node
sub getHotSpotFreq{
  my $node = $_[0];
  my $codons = $_[1];
  my $motif = $_[2]; #Motif in IUPAC format
  
  #http://www.bioinformatics.org/sms/iupac.html
  my %iupac = (
    "A" => "A",
    "C" => "C",
    "G" => "G",
    "T" => "T",
    "R" => "[A|G]",
    "Y" => "[C|T]",
    "S" => "[G|C]",
    "W" => "[A|T]",
    "K" => "[G|T]",
    "M" => "[A|C]",
    "B" => "[C|G|T]",
    "D" => "[A|G|T]",
    "H" => "[A|C|T]",
    "V" => "[A|C|G]",
    "N" => "[A|C|G|T]"
  );
  
  my $regex = $iupac{substr($motif,0,1)};
  for(my $i = 1; $i < length($motif); $i++){
    $regex = $regex.$iupac{substr($motif,$i,1)}
  }
  print "$regex\n";

  getHotSpotFreq_recur($node,$codons,$regex,$motif);

}
sub getHotSpotFreq_recur{
  my $node = $_[0];
  my $codons = $_[1];
  my $regex = $_[2];
  my $motif = $_[3];#
  my $length = length($motif);

  if(exists($node->{"left"})){
    getHotSpotFreq_recur($node->{"left"},$codons,$regex,$motif,$length);
    getHotSpotFreq_recur($node->{"right"},$codons,$regex,$motif,$length);
  }

  if($node->{"level"}!=0){
  if(!exists($node->{"sequence"})){die($node->{"subtaxa"})}

  my @c = @$codons;
  my @seqtemp = @{untransarrayCodon($node->{"sequence"},$codons)};
  my $sequence = "";
  foreach my $s (@seqtemp){$sequence = $sequence.$s;}
  my @seq = split("",uc $sequence);
  my $count = 0;
  for(my $i=0;$i<(scalar(@seq)-$length+1);$i++){
    my $mer = "";
    for(my $j = 0; $j < $length; $j++){
      $mer .= $seq[$i+$j];
    }
    if($mer =~ /$regex/){
      $count++;
    }
  }
 # print "@seq\t$count\n";
  $node->{"$motif\_c"}=$count/length($sequence);
  $node->{"$motif\_count"}=$count;
  }
}


#translate codons to indexes in Q matrix
sub transarrayCodon{
  my @in = @{$_[0]};
  my %codoni = %{$_[1]};
  my @trans;
  for(my $i = 0; $i < scalar(@in); $i+=3){
    if(!exists($codoni{lc $in[$i].$in[$i+1].$in[$i+2]})){
      push(@trans,"NA");
    }else{
      push(@trans,$codoni{lc $in[$i].$in[$i+1].$in[$i+2]});
    }
  }
  return(\@trans);
}

#re-translate indexes to codons
sub untransarrayCodon{
  my @in = @{$_[0]};
  my @codons = @{$_[1]};
  my @trans;
  for(my $i = 0; $i < scalar(@in); $i++){
    if($in[$i] ne "NA"){
      push(@trans,$codons[$in[$i]]);
    }else{
      print "NA found $i\n";
      push(@trans,"NNN");
    }
  }
  return(\@trans);
}

#translate codons to indexes in Q matrix
sub transarrayChar{
  my @in = @{$_[0]};
  my %codoni = %{$_[1]};
  my @trans;
  for(my $i = 0; $i < scalar(@in); $i++){
    if(!exists($codoni{lc $in[$i]})){
      push(@trans,"NA");
    }else{
      push(@trans,$codoni{lc $in[$i]});
    }
  }
  return(\@trans);
}

#re-translate indexes to codons
sub untransarrayChar{
  my @in = @{$_[0]};
  my @codons = @{$_[1]};
  my @trans;
  for(my $i = 0; $i < scalar(@in); $i++){
    if($in[$i] ne "NA"){
      push(@trans,$codons[$in[$i]]);
    }else{
      print "NA found $i\n";
      push(@trans,"NNN");
    }
  }
  return(\@trans);
}

#translate sequence back to nucleotide
sub printseqsCodon{
  my $node = $_[0];
  my $string = $_[1];
  my $codons = $_[2];
  if(exists($node->{"left"})){
    $string = printseqsCodon($node->{"left"},$string,$codons);
    $string = printseqsCodon($node->{"right"},$string,$codons);
  }
  elsif($node->{"id"} ne ""){
    if(!exists($node->{"sequence"})){die($node->{"id"})}
    my @seq = @{untransarrayCodon($node->{"sequence"},$codons)};
    my $sequence = "";
    foreach my $s (@seq){$sequence = $sequence.$s;}
    $string = $string.">".$node->{"id"}."\n$sequence\n";
  }
  return($string);
}

#translate sequence back to nucleotide
sub printseqsGen{
  my $node = $_[0];
  my $string = $_[1];
  my $codons = $_[2];
  if(exists($node->{"left"})){
    $string = printseqsCodon($node->{"left"},$string,$codons);
    $string = printseqsCodon($node->{"right"},$string,$codons);
  }
  elsif($node->{"id"} ne ""){
    if(!exists($node->{"sequence"})){die($node->{"id"})}
    my @seq = @{untransarrayCodon($node->{"sequence"},$codons)};
    my $sequence = "";
    foreach my $s (@seq){$sequence = $sequence.$s;}
    $string = $string.">".$node->{"id"}."\n$sequence\n";
  }
  return($string);
}

#Print sequences at all nodes
sub printAncSeqsCodon{
  my $node = $_[0];
  my $string = $_[1];
  my $codons = $_[2];
  if(exists($node->{"left"})){
    $string = printAncSeqsCodon($node->{"left"},$string,$codons);
    $string = printAncSeqsCodon($node->{"right"},$string,$codons);
  }

  if(!exists($node->{"sequence"})){die($node->{"subtaxa"})}
  my @seq = @{untransarrayCodon($node->{"sequence"},$codons)};
  my $sequence = "";
  foreach my $s (@seq){$sequence = $sequence.$s;}
  $string = $string.">".$node->{"level"}.";".$node->{"subtaxa"}.";".$node->{"divergence"}."\n$sequence\n";
  return($string);
}


#Print ML codon sequence at all nodes
sub printML_codon{
  my $node = $_[0];
  my $string = $_[1];
  my $codons = $_[2];
  my $rootonly = $_[3];
  if($rootonly == 1 && $node->{"level"} > 1){return $string;}

  if(exists($node->{"left"})){
    $string = printML_codon($node->{"left"},$string,$codons,$rootonly);
    $string = printML_codon($node->{"right"},$string,$codons,$rootonly);
  }

  if($node->{"level"}!=0){
  #print ">".$node->{"enum"}."\n";
  print ">".$node->{"level"}.";".$node->{"subtaxa"}.";".$node->{"divergence"}."\n";
  my @sequence;
  for(my $i=0; $i < scalar(@{$node->{"Codon_lhoods"}});$i++){
    my $sum;
    my $max = -inf;
    my $maxchar = -1;
    for(my $j=0; $j < 61; $j++){
      my $n = $node->{"Codon_lhoods"}->[$i][$j];
      if($n > $max){
        $max = $n; 
        $maxchar=$j;
      }
      if($j==0){$sum=$n}
      else{$sum+=log(1+exp($n-$sum));}
    }
    push(@sequence,$maxchar);
    my @test=($maxchar);
    my @seq = @{untransarrayCodon(\@test,$codons)};
    print "$i\t".($seq[0])."\t".(exp($max-$sum))."\t$max\t$sum\n";
  }
  $node->{"sequence"} = \@sequence;
  
  my @seq = @{untransarrayCodon($node->{"sequence"},$codons)};
  my $sequence = "";
  foreach my $s (@seq){$sequence = $sequence.$s;}
  $string = $string.">".$node->{"level"}.";".$node->{"subtaxa"}.";".$node->{"divergence"}."\n$sequence\n";
  }

  return($string);
}

#Print ML codon sequence at all nodes
sub printMAR_codon{
  my $node = $_[0];
  my $string = $_[1];
  my $codons = $_[2];
  my $rootonly = $_[3];
  if($rootonly == 1 && $node->{"level"} > 1){return $string;}

  if(exists($node->{"left"})){
    $string = printMAR_codon($node->{"left"},$string,$codons,$rootonly);
    $string = printMAR_codon($node->{"right"},$string,$codons,$rootonly);
  }

  if($node->{"level"}!=0){
  my @sequence;
  for(my $i=0; $i < scalar(@{$node->{"Codon_lhoods"}});$i++){
    for(my $j=0; $j < 61; $j++){
      my $n = $node->{"Codon_lhoods"}->[$i][$j];
      $string .= $node->{"id"}."\t$i\t$j\t$n\n";
    }
  }
  }
  return($string);
}


#Print ML amino acid sequence at all nodes
sub printML_aa{
  my $node = $_[0];
  my $string = $_[1];
  my $codons = $_[2];
  my $aatable1 = $_[3];
  my $aatable3 = $_[4];
  my $rootonly = $_[5];

  if($rootonly == 1 && $node->{"level"} > 1){return $string;}
  if(exists($node->{"left"})){
    $string = printML_aa($node->{"left"}, $string,$codons,$aatable1,$aatable3,$rootonly);
    $string = printML_aa($node->{"right"},$string,$codons,$aatable1,$aatable3,$rootonly);
  }

  if($node->{"level"}!=0){
  my @sequence;
    for(my $i=0; $i < scalar(@{$node->{"Codon_lhoods"}});$i++){
      my %aas;
      for(my $j=0; $j < 61; $j++){ #tally up relative likelihoods of amino acids
        my $n = $node->{"Codon_lhoods"}->[$i][$j];
        my $aa = $aatable1->{uc $aatable3->{$codons->[$j]}}; #single letter amino acid
        if(exists($aas{$aa})){
          $aas{$aa} = $aas{$aa} + log(1+exp($n - $aas{$aa}))
        }else{
          $aas{$aa} = $n
        }
      }

      my $max = -inf;
      my $maxchar = -1;
      foreach my $key (keys %aas){
        my $n = $aas{$key};
        if($n > $max){
          $max = $n; 
          $maxchar=$key;
        }
      }
      push(@sequence,$maxchar);
    }
    #$node->{"sequence"} = \@sequence;
    my $sequence = "";
    foreach my $s (@sequence){$sequence = $sequence.$s;}
    $string = $string.">".$node->{"level"}.";".$node->{"subtaxa"}.";".$node->{"divergence"}."\n$sequence\n";
    }
  return($string);
}


sub getTotalLength {
  my $node = $_[0];
  my $dist = $_[1];
# print $node->{"level"}."\t".$node->{"dist"}."\t$dist\n";
  if(exists($node->{"left"})){
    $dist = getTotalLength($node->{"left"},$dist);
    $dist = getTotalLength($node->{"right"},$dist);
  }
  $dist += $node->{"dist"};
  return($dist);  
}

sub getInternalLength {
  my $node = $_[0];
  my $dist = $_[1];
# print $node->{"level"}."\t".$node->{"dist"}."\t$dist\n";
  if(exists($node->{"left"})){
    $dist += $node->{"dist"};
    $dist = getInternalLength($node->{"left"},$dist);
    $dist = getInternalLength($node->{"right"},$dist);
  }
  return($dist);  
}


#make subtaxa labels
sub getSubTaxa{
  my $node = $_[0];
  #if a tip
  if(!exists($node->{"right"})){
    $node->{"subtaxa"}=$node->{"id"};
  }else{
    #if internal node
    getSubTaxa($node->{"right"});
    getSubTaxa($node->{"left"});
    my @l = split(",",$node->{"left"}->{"subtaxa"});
    my @r = split(",",$node->{"right"}->{"subtaxa"});
    my @total = (@l,@r);
    @total = sort @total;
    $node->{"subtaxa"}=join(",",@total);
  }
}

#get most likely ancestral sequence at each node
sub getMLAnc{
  my $node = $_[0];
  if(exists($node->{"right"})){
    getMLAnc($node->{"right"});
    getMLAnc($node->{"left"});
  }
  my @sequence;
  if($node->{"level"}!=0){
  for(my $i=0; $i < scalar(@{$node->{"mat"}});$i++){
    my $max = -inf;
    my $maxchar = -1;
    for(my $j=0; $j < 61; $j++){
      my $n = ($node->{"Codon_lhoods"}->[$i][$j]);
      if($n > $max){
        $max = $n; 
        $maxchar=$j;
      }
    }
    push(@sequence,$maxchar);
  }
  $node->{"sequence"} = \@sequence;
  }
}

sub translateSeqs{
  my $seqs = $_[0];
  my %tr = %{codonTable()};
  my %new;
  foreach my $k (keys %$seqs){
    my @s = split("",$seqs->{$k});
    my @n;
    for(my $i=0;$i < scalar(@s);$i+=3){
      push(@n,$tr{lc $s[$i].$s[$i+1].$s[$i+2]});
    }
    $new{$k}=join(",",@n);
  } 
  return \%new;
}

sub translateSeqsSingle{
  my $seqs = $_[0];
  my %tr = %{codonTable()};
  my %tr1 = %{codonTableSingle()};
  my %new;
  foreach my $k (keys %$seqs){
    my @s = split("",$seqs->{$k});
    my @n;
    for(my $i=0;$i < scalar(@s);$i+=3){
      push(@n,$tr1{uc $tr{lc $s[$i].$s[$i+1].$s[$i+2]}});
    }
    $new{$k}=join("",@n);
  } 
  return \%new;
}

#change a specified edge length
sub changeedge{
  my $node = $_[0];
  my $id = $_[1];
  my $length = $_[2];
  if(exists($node->{"left"})){
    changeedge($node->{"left"},$id,$length);
    changeedge($node->{"right"},$id,$length);
  }
  if(!exists($node->{"id"})){
    die("Node ID doens't exist!");
  }elsif($node->{"subtaxa"} eq $id){
    $node->{"dist"}=$length;
  }
}

#assign parents to each node
sub assignparents{
  my $node = $_[0];
  my $parent = $_[1];
  if($node->{"level"} != 0){
    $node->{"up"}=$parent;
  }
  if(exists($node->{"left"})){
    assignparents($node->{"left"},$node);
    assignparents($node->{"right"},$node);
  }
}

#get a reference to a particular node
sub getnode{
  my $node = $_[0];
  my $id=$_[1];
  my $nid = "";
  if(exists($node->{"left"})){
    my $nid1 = getnode($node->{"left"},$id);
    my $nid2 = getnode($node->{"right"},$id);
    if($nid1 ne ""){$nid=$nid1}
    if($nid2 ne ""){$nid=$nid2}
  }
 if($node->{"id"} eq $id){
  return $node;
 }else{
  return $nid;
 }
}

#get a reference to a particular node
sub getnode_subtaxa{
  my $node = $_[0];
  my $id=$_[1];
  my $nid = "";
  if(exists($node->{"left"})){
    my $nid1 = getnode_subtaxa($node->{"left"},$id);
    my $nid2 = getnode_subtaxa($node->{"right"},$id);
    if($nid1 ne ""){$nid=$nid1}
    if($nid2 ne ""){$nid=$nid2}
  }
 # print $node->{"subtaxa"}."\n";
 if($node->{"subtaxa"} eq $id){
  return $node;
 }else{
  return $nid;
 }
}


#return codon translation table
sub codonTable{
  my %codons = (
"ttt" =>  "Phe",
"ttc" =>  "Phe",
"tta" =>  "Leu",
"ttg" =>  "Leu",
"ctt" =>  "Leu",
"ctc" =>  "Leu",
"cta" =>  "Leu",
"ctg" =>  "Leu",
"att" =>  "Ile",
"atc" =>  "Ile",
"ata" =>  "Ile",
"atg" =>  "Met",
"gtt" =>  "Val",
"gtc" =>  "Val",
"gta" =>  "Val",
"gtg" =>  "Val",
"tct" =>  "Ser",
"tcc" =>  "Ser",
"tca" =>  "Ser",
"tcg" =>  "Ser",
"cct" =>  "Pro",
"ccc" =>  "Pro",
"cca" =>  "Pro",
"ccg" =>  "Pro",
"act" =>  "Thr",
"acc" =>  "Thr",
"aca" =>  "Thr",
"acg" =>  "Thr",
"gct" =>  "Ala",
"gcc" =>  "Ala",
"gca" =>  "Ala",
"gcg" =>  "Ala",
"tat" =>  "Tyr",
"tac" =>  "Tyr",
"taa" =>  "STOP",
"tag" =>  "STOP",
"cat" =>  "His",
"cac" =>  "His",
"caa" =>  "Gln",
"cag" =>  "Gln",
"aat" =>  "Asn",
"aac" =>  "Asn",
"aaa" =>  "Lys",
"aag" =>  "Lys",
"gat" =>  "Asp",
"gac" =>  "Asp",
"gaa" =>  "Glu",
"gag" =>  "Glu",
"tgt" =>  "Cys",
"tgc" =>  "Cys",
"tga" =>  "STOP",
"tgg" =>  "Trp",
"cgt" =>  "Arg",
"cgc" =>  "Arg",
"cga" =>  "Arg",
"cgg" =>  "Arg",
"agt" =>  "Ser",
"agc" =>  "Ser",
"aga" =>  "Arg",
"agg" =>  "Arg",
"ggt" =>  "Gly",
"ggc" =>  "Gly",
"gga" =>  "Gly",
"ggg" =>  "Gly"
);
  return \%codons;
}

#return codon translation table
sub codonTableSingle{
  my %codons = (
"ALA" => "A",
"CYS" => "C",
"ASP" => "D",
"GLU" => "E",
"PHE" => "F",
"GLY" => "G",
"HIS" => "H",
"ILE" => "I",
"LYS" => "K",
"LEU" => "L",
"MET" => "M",
"ASN" => "N",
"PRO" => "P",
"GLN" => "Q",
"ARG" => "R",
"SER" => "S",
"THR" => "T",
"VAL" => "V",
"TRP" => "W",
"TYR" => "Y");
  return \%codons;
}

1