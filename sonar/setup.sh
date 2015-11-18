
#!/bin/bash

#This script sets up paths for the SONAR (Antibodyomics) pipeline


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
    blast="blastn"
    dnaml="dnaml"
    beast="/Applications/BEAST/bin/beast";

#Ask user for input
    echo ""
    if [[ ${BASH_VERSION[0]} < 4 ]]
    then
	read -e -p "Please enter the default root directory (default $homeDir): " var; homeDir=${var:-$homeDir}
	read -e -p "Please enter the directory where SONAR is installed (default $pipeDir): " var; pipeDir=${var:-$pipeDir}
	read -e -p "Please enter the path to ClustalO (default $clustal): " var; clustal=${var:-$clustal}
	read -e -p "Please enter the path to Muscle (default $muscle): " var; muscle=${var:-$muscle}
	read -e -p "Please enter the path to USearch (default $usearch): " var; usearch=${var:-$usearch}
	read -e -p "Please enter the path to blastn (default $blast): " var; blast=${var:-$blast}
	read -e -p "Please enter the path to DNAML (default $dnaml): " var; dnaml=${var:-$dnaml}
	read -e -p "Please enter the path to BEAST (default $beast):" var; beast=${var:-$beast}
    else
	read -e -p "Please enter the default root directory: " -i $homeDir homeDir
	read -e -p "Please enter the directory where SONAR is installed: " -i $pipeDir pipeDir
	read -e -p "Please enter the path to ClustalO: " -i $clustal clustal
	read -e -p "Please enter the path to Muscle: " -i $muscle muscle
	read -e -p "Please enter the path to USearch: " -i $usearch usearch
	read -e -p "Please enter the path to blastn: " -i $blast blast
	read -e -p "Please enter the path to DNAML: " -i $dnaml dnaml
	read -e -p "Please enter the path to BEAST: " -i $beast beast
    fi

    #verify inputs
    echo -e "\n\nYou have entered the following values:
\tDefault root directory: $homeDir
\tDirectory where SONAR is installed: $pipeDir
\tPath to ClustalO: $clustal
\tPath to Muscle: $muscle
\tPath to USearch: $usearch
\tPath to blastn: $blast
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
	clustBlast="blastn"

	#Ask user for input
	echo ""
	if [[ ${BASH_VERSION[0]} < 4 ]]
	then
	    read -e -p "Please enter the path to qsub (default $qsub): " var; qsub=${var:-$qsub}
	    read -e -p "Please enter the path to Muscle on the cluster (default $clustMuscle): "  var; clustMuscle=${var:-$clustMuscle}
	    read -e -p "Please enter the path to blastn on the cluster (default $clustBlast): " var; clustBlast=${var:-$clustBlast}
	else
	    read -e -p "Please enter the path to qsub: " -i $qsub qsub
	    read -e -p "Please enter the path to Muscle on the cluster: " -i $clustMuscle clustMuscle
	    read -e -p "Please enter the path to blastn on the cluster: " -i $clustBlast clustBlast
	fi

	#and check
	echo -e "\n\nYou have entered the following values:
\tPath to qsub:  $qsub
\tPath to Muscle on the cluster: $clustMuscle
\tPath to blastn on the cluster: $clustBlast\n"

	read -e -p "Is this correct [y/n]? " checked

    done

fi #with cluster specific questions



#add SONAR to path?
read -e -p "Would you like to permanently add SONAR to the path (by adding a line to .bashrc) [y/n]?" path
if [ "$path" == "Y" ] || [ "$path" == "y" ]; then
    echo "export PATH=$pipeDir/annotate:$pipeDir/lineage:$pipeDir/phylogeny:$pipeDir/plotting:$pipeDir/utilities:\$PATH" >> ~/.bashrc
    echo "export PERL5LIB=$pipeDir:\$PERL5LIB" >> ~/.bashrc
    pypath=$(cd $pipeDir; cd ..; pwd)
    echo "export PYTHONPATH=$pypath:\$PYTHONPATH" >> ~/.bashrc
fi




#Now create output:

#create a file for the Python portion of the pipeline
echo "

HOME_FOLDER    = \"$homeDir\"
SCRIPT_FOLDER  = \"$pipeDir\"
clustal	       = \"$clustal\"
blast_cmd      = \"$blast\"
usearch        = \"$usearch\"

clusterExists  = True
qsub           = \"$qsub\"
cluster_muscle = \"$clustMuscle\"
cluster_blast  = \"$clustBlast\"

" > paths.py



#create a file for the Perl portion of the pipeline
echo "
#!/usr/bin/env perl

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
\$ppath{'sonar'}=\"$pipeDir\";# absolute path to SONAR
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
