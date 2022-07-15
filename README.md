SONAR - Software for Ontogenic aNalysis of Antibody Repertoires [![Build Status](https://travis-ci.com/scharch/SONAR.svg?branch=master)](https://travis-ci.com/scharch/SONAR)
=====

<a href="https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html"><img src="https://github.com/airr-community/software-wg/blob/master/swtools_badge/sw_tools_badge_v1.0.png" alt="AIRR SW Tools v1.0 valid"/></a>

Introduction
-----

SONAR is the Antibodyomics2 pipeline developed collaboratively by the Structural Bioinformatics Section of the Vaccine Research Center, NIAID, NIH and the Shapiro lab at Columbia University. It is designed to process longitudinal NGS samples by extracting sequences related to a known antibody or antibodies of interest and reconstructing the ontogeny of the lineage. For more details, please see [below](#Using-SONAR). For examples of papers using SONAR, please see [here](#Papers).

If you use SONAR, please cite:

Schramm et al "SONAR: A High-Throughput Pipeline for Inferring Antibody Ontogenies from Longitudinal Sequencing of B Cell Transcripts." **Frontiers Immunol**, 2016. [PMCID: PMC5030719](https://www.frontiersin.org/articles/10.3389/fimmu.2016.00372/full)

For GSSPs, please cite:

Sheng et al "Gene-Specific Substitution Profiles Describe the Types and Frequencies of Amino Acid Changes during Antibody Somatic Hypermutation." **Frontiers Immunol**, 2017. [PMCID: PMC5424261](https://www.frontiersin.org/articles/10.3389/fimmu.2017.00537/full).

Installation
-----

**SONAR requires Python >=3.6**

### Docker
SONAR is available as an automatically updated Docker image. To use Docker:
```
$> docker pull scharch/sonar
$> docker run -it -e DISPLAY=$DISPLAY -v ~:/host/home scharch/sonar
$root@1abcde234> cd /host/home
$root@1abcde234> /SONAR/sample_data/runVignette.sh
<< *OR* >>
$root@1abcde234> cd /host/home/*path*/*to*/*data*/
$root@1abcde234> sonar 1.1
.
.
.
```

### Installing locally

#### General Prerequisites:
* Python3 with Biopython, airr, and docopt
* Perl5 with BioPerl, Statistics::Basic, List::Util, and Algorithm::Combinatorics
* R with docopt, ggplot2, MASS, grid, and ptinpoly

#### Optional Prerequisites:
* For using the master script: fuzzywuzzy python package
* For single cell (paired heavy/light) clonality: networkx python package
* For inferring ancestor sequences with IgPhyML: PDL::LinearAlgebra::Trans perl package
* For displaying trees: ete3, PyQT4, and PyQt4.QtOpenGL python packages
* For comparing GSSPs: pandas python package

For details on how to install the prerequisites, follow the recipe used in the [Dockerfile](Dockerfile).

Then clone the github repo and run the setup utility:
```
$> git clone https://github.com/scharch/SONAR.git
$> cd SONAR
$> ./setup.py
$> cp sonar ~/bin
```

If you wish, you may verify/test the installation by running
```
$> python3 tests/run_tests.py
```

Using SONAR
-----
To see a summary of all SONAR scripts and what they do, simply run `sonar -h`. Alternatively, take advantange of the fuzzy matching to find the scripts in a particular module, eg `sonar annotate`. All sonar scripts will print detailed options and usage when passed the `-h` flag. For a detailed summary, please see [the vignette](vignette.pdf).

Support
-----
I am more than happy to assist with basic SONAR usage. Please file all bugs reports and requests for help as GitHub [issues](https://github.com/scharch/SONAR/issues). I will typically respond within a day or two, though it may take me up to a month to push out bug fixes, depending on the criticallity and complexity of the bug, as well as other obligations.

Change Log
-----
### New in version 4.3
* `1.3-finalize_assignments.py` now produces a `complete_vdj` column in the AIRR TSV.
* `2.1-calculate_id-div.py` and `2.4-cluster_into_groups.py` now support the `--species` option.
* `2.4-cluster_into_groups.py` now only uses AIRR TSV input/output.
* In `2.4-cluster_into_groups.py`, addition of "native" mAbs has been generalized to the inclusion of multiple rearrangement TSVs from any source. This also partially replaces `3.1-merge_timepoints.pl`. The `--names` parameter has been added to distinguish the source of each input rearrangement set.
* `2.4-cluster_into_groups.py` now has a `--geneClusters` options to group V genes into closely related/hard to distinguish cluster, allowing more sensitivity in clonal assigment, especially for amplicons with only partial V genes. This option also turns off J gene matching. Clusters have been generated for the default human/rhesus gene sets (see `sample_data/functionalClusters`); for custom gene sets, clusters can be specificed with `--cutomClusters`.
* `2.4-cluster_into_groups.py` now supports a `--singlecell` mode, which assigns clonality based on a joint analysis of CDRH3 and CDRL3 clusterings. In this mode, the default sequence identity threshold is 80%, providing increased sensitivity, while the joint analysis keeps the specificity reasonable.
* The filtering function has been removed from `getFastaFromAIRR.py` and put in a new standalone script, `filterAIRR.py`. In addition, the syntax for filters has been rationalized and hopefully simplified.

[_OLDER CHANGES_](Changes.md)

Papers
-----
* Doria-Rose et al "Developmental pathway for potent V1V2-directed HIV-neutralizing antibodies." **Nature**, 2014. [PMCID: PMC4395007](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4395007/)
* Wu et al "Maturation and Diversity of the VRC01-Antibody Lineage over 15 Years of Chronic HIV-1 Infection." **Cell** 2015. [PMCID: PMC4706178](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4706178/)
* Bhiman et al "Viral variants that initiate and drive maturation of V1V2-directed HIV-1 broadly neutralizing antibodies." **Nat Med**, 2015. [PMCID: PMC4637988](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4637988/)
* Doria-Rose et al "New Member of the V1V2-Directed CAP256-VRC26 Lineage That Shows Increased Breadth and Exceptional Potency." **J Virol**, 2016. [PMCID: PMC4702551](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4702551/)
* Sheng et al "Effects of Darwinian Selection and Mutability on Rate of Broadly Neutralizing Antibody Evolution during HIV-1 Infection." **PLoS Comput Biol**, 2016. [PMCID: PMC4871536](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4871536/)
* Huang et al "Identification of a CD4-Binding-Site Antibody to HIV that Evolved Near-Pan Neutralization Breadth." **Immunity**, 2016. [PMCID: PMC5770152](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5770152/) 
* Sheng et al "Gene-Specific Substitution Profiles Describe the Types and Frequencies of Amino Acid Changes during Antibody Somatic Hypermutation." **Frontiers Immunol**, 2017. [PMCID: PMC5424261](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5424261/)
* Krebs et al "Longitudinal Analysis Reveals Early Development of Three MPER-Directed Neutralizing Antibody Lineages from an HIV-1-Infected Individual." **Immunity**, 2019. [PMCID: PMC6555550](https://www.ncbi.nlm.nih.gov/pubmed/30876875)
