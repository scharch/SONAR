#!/usr/bin/perl
#    barcode file should be in the format: name\tsequence\n
use strict;
use FindBin;
use lib "$FindBin::Bin/../";
use PPvars qw(ppath);

#########checking parameters#######
my $usage="

Usage: SONAR_master.pl command parameters....
	Example:
	SONAR_master.pl dereplicate_sequences -pu usearch -min1 2 -min2 3 -f ./test.fa -t 5
	SONAR_master.pl [dereplicate_sequences] -h for help information.
Command options:
	MiSeq_assembly	Merge Illumina MiSeq paired reads to a single read
	blast_V	\tAssign germline V genes to each read. It can also assign V(D)J with one step (see help). 
	blast_J	\tAssignemnt of D and J genes
	finalize_assignments	finalize the assignemnt of V(D)J genes
	dereplicate_sequences	Remove duplicate sequences 
	calculate_id-div	Calculate somatic hypermutation level and sesquence divergence to target antibodies
	get_island	\tExtract sequences within an interested island on 2D plot.
	intradonor_analysis	Find lineage related sequences of known antibodies on phylogenetic tree
	cluster_into_groups	Cluster sequences into clones or lineages
	merge_timepoints	Recluster selected sequences from multiple timepoints to determine birthdays and persistence times
	runDNAML	\tBiuld phylogenetic tree using DNAML in the PHYLIP package.
	pick_intermediates	Select developmental intermediates for the specified native mAb 
	cluster_tree	\tCluster sequences in a phylgentic tree so that all clusters are monophyletic clades
	evolutionary_rate	Generate xml configuration files for BEAST2 evolutionary rate estimation.
	setup_plots	\tParse SONAR output tables and generates files for easy plotting with 4.2-generate_plots.R
	plot_histogram	\tGenerate nice histograms to demonstrate the features of the repertore.
	plot_identity_divergence	Generate identity-divergence plot
	plot_tree	\tGenerate phylogenetic tree color coded longitudinally
	checkClusterBlast	Monitor the SGE array jobs submitted by scripts 1.1 and 1.2 and resubmit failed jobs.
	checkForFrameshift	Detect and discard sequences with frameshift errors
	fastq_quality_boxplot_graph	Plot the sequencing quality of a dataset
	getFastaFromList	extracting a subset of sequences from a large fasta file
	getListFromFasta	getting a list of the sequence identifiers from a large fasta file

Created by Zizhang Sheng.

Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
 ";
my %scripts=(
							'assign_VDJ','1',
							'MiSeq_assembly','1.0-MiSeq_assembly.pl',
							'blast_V','1.1-blast_V.py',
							'blast_J','1.2-blast_J.py',
							'finalize_assignments','1.3-finalize_assignments.py',
							'dereplicate_sequences','1.4-dereplicate_sequences.pl',
							'calculate_id-div','2.1-calculate_id-div.pl',
							'get_island','2.2-get_island.py',
							'intradonor_analysis','2.3-intradonor_analysis.py',
							'cluster_into_groups','2.4-cluster_into_groups.py',
							'merge_timepoints','3.1-merge_timepoints.pl',
							'run_DNAML','3.2-run_DNAML.py',
							'pick_intermediates','3.3-pick_intermediates.pl',
							'cluster_tree','3.4-cluster_tree.pl',
							'evolutionary_rate','3.5-evolutionary_rate.pl',
							'setup_plots','4.1-setup_plots.pl',
							'plot_histograms','4.2-plot_histograms.R',
							'plot_identity_divergence','4.3-plot_identity_divergence.R',
							'display_tree','4.4-display_tree.py',
							'checkClusterBlast','checkClusterBlast.py',
							'getFastaFromList','getFastaFromList.py',
							'checkForFrameshift','checkForFrameshift.py',
							'getListFromFasta','getListFromFasta.py',
							'fastq_quality_boxplot_graph','fastq_quality_boxplot_graph.sh'
							);


my $command=shift @ARGV;
if($command=~/\-{1,2}(h|help)$/){die "$usage";}

foreach(@ARGV){if($_=~/\-{1,2}(h|help)$/){&help($command);die "\n";}}
if(!$scripts{$command}){die "No $command found. Check spelling\n";}


#####################################
sub help{
    my ($script)=@_;print "$scripts{$script}\n";
    system("python ./utilities/$scripts{$script} -h"); 
  }
  
sub tasks{
	  my ($script)=@_;
}



