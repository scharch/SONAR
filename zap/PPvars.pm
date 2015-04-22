#!/usr/bin/perl

package PPvars;
use strict;

use vars '@ISA', '@EXPORT', '$NAME', '$VERSION', '$DATE', '$AUTHOR';
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(ppath);



my %ppath=();
$ppath{'usearch'}="usearch";# absolute path to usearch
$ppath{'clustalo'}="clustalo";# absolute path to clustalo
$ppath{'clustalw'}="clustalw2";# absolute path to clustalw
$ppath{'mafft'}="mafft";# absolute path to mafft
$ppath{'muscle'}="muscle";# absolute path to neighbor
$ppath{'beast'}="/Applications/BEAST/bin/beast";# absolute path to beast

sub ppath{
    my $p=shift;	
	if($ppath{$p}){
	  return $ppath{$p};	
	}
	else{
	  return '';	
	}
}




1;

