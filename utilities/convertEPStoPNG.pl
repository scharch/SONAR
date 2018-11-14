#!/usr/bin/env perl

=head 1 SUMMARY
 
 convertEPStoPNG.pl

 This program uses LYNX (text broswer) to script the conversion of GPPS plots in
       EPS format to high-resolution PNG files using the psviewer.org website.
       This is necessary because the standard tools have trouble with the way
       these files are constructed (adapted from WebLogo (weblogo.berkeley.edu/logo.cgi)
       for reasons that I don't really understand.

 Usage: convertEPStoPNG.pl file1.eps [ file2.eps file3.eps ... ]

     Invoke with -h or --help to print this documentation. 
 
 Required Parameters:
     file.eps  =>  File(s) to convert. Output will be saved as file.png in the
                         same directory.
 
 Created by Chaim A. Schramm 2017-02-24.

 Copyright (c) 2011-2017 Columbia University and Vaccine Research Center, National
                          Institutes of Health, USA. All rights reserved.
=cut


use warnings;
use strict;
use diagnostics;
use Pod::Usage;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);
use version qw/logCmdLine/;


if ($#ARGV<0 || $ARGV[0] =~ /-h/) { pod2usage(1); }
&logCmdLine($0,@ARGV);


for my $file (@ARGV) {
    my $outDir    = dirname($file);
    my $splitFile = join("\nkey ", split(//, $file));
    my $splitDir  = $outDir eq "" ? "" : "\nkey Home\nkey " . join("\nkey ", split(//, $outDir)) . "\nkey /";

    open SCRIPT, ">tmp.lx" or die "Can't create temporary script tmp.lx: $!\n";
    print SCRIPT <<FIN;
key g
key p
key s
key v
key i
key e
key w
key e
key r
key .
key o
key r
key g
key /
key c
key o
key n
key v
key e
key r
key t
key p
key s
key t
key o
key p
key n
key g
key .
key a
key s
key p
key x
key ^J
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key $splitFile
key Down Arrow
key Down Arrow
key Down Arrow
key Right Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Right Arrow
key y
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key <delete>
key <delete>
key <delete>
key 3
key 0
key 0
key 3
key Down Arrow
key <delete>
key <delete>
key <delete>
key 1
key 0
key 0
key 0
key Down Arrow
key Right Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Down Arrow
key Right Arrow
key D
key Up Arrow
key Right Arrow$splitDir
key ^J
key q
key y
FIN

    close SCRIPT;
    system("lynx -cmd_script=tmp.lx");

} #next file
