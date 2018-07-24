===============================================================
SONAR - Software for Ontogenic aNalysis of Antibody Repertoires
===============================================================

Introduction
============

Antibodyomics 2.0 pipeline developed collaboratively by the Vaccine Research Center, NIH and the Shapiro lab at Columbia University.
Please cite Schramm et al Frontiers Immunol 2016
For GSSPs, please cite Sheng et al Frontiers Immunol 2017

Prerequisites
=============

** SONAR has now been update to Python3 **

#General Prerequisites:
* Python3 with Biopython and docopt
* Perl5 with BioPerl, Statistics::Basic, List::Util, and Algorithm::Combinatorics
* R with docopt, ggplot2, MASS, and grid
* USearch v8 (newer versions will break parsing of output)
* Muscle v5
* NCBI BLAST+

#Optional Prerequisites:
*For preprocessing of raw sequencing data: FastX-Toolkit
*For evolutionary rate calculations: BEAST2
*For displaying trees: ete3, PyQT4, and PyQt4.QtOpenGL python packages
*For comparing GSSPs: pandas python package



More information coming. In the meantime, please see:

Doria-Rose et al Nature 2014 doi: 10.1038/nature13036

Wu et al Cell 2015 doi: 10.1016/j.cell.2015.03.004

