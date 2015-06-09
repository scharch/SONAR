
#~/bin/bash

#This script sets up paths for the SOAnAR (Antibodyomics) pipeline


#print useful information
echo -e "
General Prerequisites: 
\tPython with BioPython
\tPerl with BioPerl
\tUSearch v8
\tMuscle v5

Module 1 (annotation) additional prerequisites:
\t**legacy** blast (will be updated to BLAST+ in final release)

Module 2 (lineage) additional prerequisites:
\t(none)

Module 3 (phylogeny) additional prerequisites:
\tDNAML
\tBEAST

Module 4 (plotting) additional prerequisites:
\tR with ggplot2, MASS, and grid
\tete2 Python module

"



#setup local directories
while [ "$verified" != "Y" ] && [ "$verified" != "y" ]; do

#these are the default values
    homeDir=$(cd ~; pwd)
    pipeDir=$(pwd)
    clustal="clustalo"
    usearch="usearch"
    muscle="muscle"
    blast="blastall"
    dnaml="dnaml"
    beast="/Applications/BEAST/bin/beast";

#Ask user for input
    echo ""
    read -e -p "Please enter the default root directory: " -i $homeDir homeDir
    read -e -p "Please enter the directory where SOAnAR is installed: " -i $pipeDir pipeDir
    read -e -p "Please enter the path to ClustalO: " -i $clustal clustal
    read -e -p "Please enter the path to Muscle: " -i $muscle muscle
    read -e -p "Please enter the path to USearch: " -i $usearch usearch
    read -e -p "Please enter the path to blastall: " -i $blast blastall
    read -e -p "Please enter the path to DNAML: " -i $dnaml dnaml
    read -e -p "Please enter the path to BEAST: " -i $beast beast
    
    #verify inputs
    echo -e "\n\nYou have entered the following values:
\tDefault root directory: $homeDir
\tDirectory where SOAnAR is installed: $pipeDir
\tPath to ClustalO: $clustal
\tPath to Muscle: $muscle
\tPath to USearch: $usearch
\tPath to blastall: $blast
\tPath to DNAML: $dnaml
\tPath to BEAST:  $beast\n"

    read -e -p "Is this correct [y/n]? " verified

done



#Check for cluster and get cluster paths
read -e -p "Will you be using an HPC cluster [y/n]?" cluster
if [ "$cluster" == "Y" ] || [ "$cluster" == "y" ]; then
    
    while [ "$checked" != "Y" ] && [ "$checked" != "y" ]; do
	#defaults
	qsub="qsub"
	clustMuscle="muscle"
	clustBlast="blastall"

	#Ask user for input
	echo ""
	read -e -p "Please enter the path to qsub: " -i $qsub homeDir
	read -e -p "Please enter the path to Muscle on the cluster: " -i $clustMuscle clustMuscle
	read -e -p "Please enter the path to blastall on the cluster: " -i $clustBlast clustBlast

	#and check
	echo -e "\n\nYou have entered the following values:
\tPath to qsub:  $qsub
\tPath to Muscle on the cluster: $clustMuscle
\tPath to blastall on the cluster: $clustBlast\n"

	read -e -p "Is this correct [y/n]? " checked

    done

fi #with cluster specific questions



#add SOAnAR to path?
read -e -p "Would you like to permanently add SOAnAR to the path (by adding a line to .bashrc) [y/n]?" path
if [ "$path" == "Y" ] || [ "$path" == "y" ]; then
    echo "export \$PATH=$pipeDIR/annotate:$pipeDIR/lineage:$pipeDIR/phylogeny:$pipeDIR/plotting:$pipeDIR/utilities:\$PATH" >> ~/.bashrc
    pypath = $(cd $pipeDir; cd ..; pwd)
    echo "export \$PYTHONPATH=$pypath:\$PYTHONPATH" >> ~/.bashrc
fi




#Now create output:

#create a file for the Python portion of the pipeline
echo "

HOME_FOLDER    = \"$homeDir\"
SCRIPT_FOLDER     = \"$pipeDir\"
clustal	       = \"$clustal\"
blastall_cmd   = \"$blast\"
usearch        = \"$usearch\"

clusterExists  = True
qsub           = \"$qsub\"
cluster_muscle = \"$clustMuscle\"
cluster_blast  = \"$clustBlast\"

" > paths.py



#create a file for the Perl portion of the pipeline
echo "
#!/usr/bin/perl

package PPvars;
use strict;

use vars '@ISA', '@EXPORT', '\$NAME', '\$VERSION', '\$DATE', '\$AUTHOR';
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(ppath);



my %ppath=();
\$ppath{'usearch'}=\"$usearch\";# absolute path to usearch
\$ppath{'clustalo'}=\"$clustalo\";# absolute path to clustalo
\$ppath{'muscle'}=\"$muscle\";# absolute path to neighbor
\$ppath{'beast'}=\"$beast\";# absolute path to beast

sub ppath{
    my \$p=shift;	
	if(\$ppath{\$p}){
	  return \$ppath{\$p};	
	}
	else{
	  return '';	
	}
}




1;
" > PPVars.pm
