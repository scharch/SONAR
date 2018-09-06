Suggested workflow for generating GSSPs using SONAR
===

1. Run Module 1 (annotation) on NGS data. In this case, we only want to dereplicate sequences without clustering, as we will do error control _within_ lineages.
   ```
   $> 1.1-blast_V.py -fasta ngs.fq
   $> 1.2-blast_J.py
   $> 1.3-finalize_assignments.py
   $> 1.4-dereplicate_sequences.pl -id 1 -min1 0 -min2 1
   ```
1. Now that we have assigned V and J, we can cluster sequences into pseudolineages.
   ```
   $> 2.4-cluster_into_groups.py
   ```
1. Repick the representative sequence from each lineage -- use the read which is closest to the consenus of the lineage. In the paper, we find that this step isn't needed for robust results. If you skip it, however, we recommend manually removing lineages with only one member sequence before proceeding.
   ```
   $> 5.1-repick_lineage_representative.py
   ```
1. Now check for frameshifts compared to germline (especially important for 454 or PacBio data) and translate to amino acids. Remove sequences with no amino acid changes from germline.
   ```
   $> 5.2-filter_sequences.py output/sequences/nucleotide/<project>_consensusLineageRepresentatives.fa \
                              output/sequences/amino_acid/<project>_finalForProfiles.fa
   ```
1. Now we are ready to construct the GSSPs. We recommend using at least 300 sequences per GSSP for robust results, but sometimes use 100 for rarer V genes. Here we calculate profiles just using the standard IMGT database, however we strongly recommend using a personalized germline allele database constructed using [TIgGER](http://tigger.readthedocs.io/en/0.2.8/), [partis](https://github.com/psathyrella/partis), or [IgDiscover](https://bitbucket.org/igdiscover/igdiscover). Incomplete sequences should be masked using the -m argument.
   ```
   $> 5.3-make_profiles.py output/sequences/amino_acid/<project>_finalForProfiles.fa -o <project>_profiles.txt -n 300 -p 0 -g /<path-to-sonar>/germDB/IgHKLV_cysTruncated.AA.fa -a
   ```
1. We can generate logo plots of the GSSPs. (Many programs seem to have problems opening the resulting EPS files; [This site](http://convertepstojpg.com/) can be used to convert them to standard image formats.)
   ```
   $> 4.5-create_GSSP_logo.pl <project>_profiles.txt output/plots/<project>-GSSPs
   ```
1. Finally, we can compare GSSPs generated from multiple datasets. This generates a matrix of weighted average Jensen-Shannon divergences between each GSSP, which can be plotted using `cmdscale` in R. It also produces tables of the weighted average entropy of each GSSP and the rarity of every possible mutation in each GSSP.
   ```
   $> 5.4-compare_profiles.py GSSP-comparisons <project1>_profiles.txt <project2>_profiles.txt <project3>_profiles.txt ...
   ```
1. Note that GSSPs can be constructed from repertoires of nonproductive rearrangements, as well. To do so, run 1.4 on all sequences with an assigned J gene, and use splitFunctionalAndNonfunctional.py in the sonar/utilities folder to get only those lineages made up entirely of nonproductive sequences.
   ```
   $> 1.4-dereplicate_sequences.pl -f output/sequences/nucleotide/<project>_allJ.fa -id 1 -min1 0 -min2 1 -s 2
   $> getListFromFasta.py -f output/sequences/nucleotide/<project>_allJ_unique.fa | \
      getFastaFromList.py -f output/sequences/nucleotide/<project>_allCDR3.fa -o output/sequences/nucleotide/<project>_allCDR3_unique.fa
   $> 2.4-cluster_into_groups.py -full output/sequences/nucleotide/<project>_allJ_unique.fa -cdr3 output/sequences/nucleotide/<project>_allCDR3_unique.fa
   $> splitFunctionalAndNonfunctional.py
   ```
At this point, the rest of the mGSSP pipeline (5.1-5.3) can be run separately for the functional and nonfunctional repertoires. Note that the --keep flag should be passed to 5.2-filter_sequences.py when applied to nonfunctional repertoires:
   ```
   $> 5.2-filter_sequences.py output/sequences/nucleotide/<project>_nonFunctionalRepresentatives.fa \
                              output/sequences/amino_acid/<project>_nonFunctional_finalForProfiles.fa --keep
   ```
