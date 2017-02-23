#!/usr/bin/perl -w

=head 1 SUMMARY

 quickConsensus.pl
 Usage: quickConsensus.pl <alignment.phy>
 
 This script generates a consensus sequences from a (phylip-style) MSA
 
 Parameters:
     alignment.phy - the MSA

 Created by Chaim A. Schramm 2013-03-29
 Copyright (c) 2013 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
=cut

use strict;
use diagnostics;
use Pod::Usage;
use Bio::AlignIO;

if ($#ARGV<0 || $ARGV[0] =~ /-h/) { pod2usage(1); }

my $in = Bio::AlignIO->new( -file=>$ARGV[0], -format=>"phylip");
my $aln = $in->next_aln;
my $name = $ARGV[0]; $name =~ s/\.phy//;
print ">$name\n" . $aln->consensus_string() . "\n";

