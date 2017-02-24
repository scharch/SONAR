#!/usr/bin/env perl

=head 1 SUMMARY

 quickConsensus.pl
 
 This script generates a consensus sequences from a (phylip-style) MSA

 
 Usage: quickConsensus.pl <alignment.phy>

     Invoke with -h or --help to print this documentation. 

 Parameters:
     alignment.phy - the MSA

 Created by Chaim A. Schramm 2013-03-29
 Added to SONAR by CA Schramm 2017-02-24
 Copyright (c) 2011-2017 Columbia University and Vaccine Research Center, National 
                          Institutes of Health, USA. All rights reserved.
=cut

use warnings;
use strict;
use diagnostics;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);
use version qw/logCmdLine/;
use Bio::AlignIO;

if ($#ARGV<0 || $ARGV[0] =~ /-h/) { pod2usage(1); }

&logCmdLine($0,@saveArgs);

my $in = Bio::AlignIO->new( -file=>$ARGV[0], -format=>"phylip");
my $aln = $in->next_aln;
my $name = $ARGV[0]; $name =~ s/\.phy//;
print ">$name\n" . $aln->consensus_string() . "\n";

