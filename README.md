# finemap-scripts
Scripts used for the formatting of files in order to run FINEMAP to calculate causal variant probabilities and downstream analyses
The scripts and README files shown here detail the creation of files required to run FINEMAP for multiple regions of genetic association.

# Running FINEMAP
FINEMAP requires three input files for each associated region, a pairwise variant correlation matrix (Pearson’s r), a file giving the z-score for each variant and a file containing priors for the number of causal variants in the region.

**The correlation (r) matrix (README_3_correlation.md)**  
Correlation matrices for the full IBD dataset were generated using the ‘corrcoeff’ function in the numPy package (python-2.7.8), which outputs a square matrix of r values from allelic dosages of the alternate allele in the IBD data to match the allele used to calculate the Z-score for the association. The order of the variants in this file is kept the same as that of the z-score file. The allele used to calculate the r correlation matrix must be the same as the allele used to calculate the Z-score. 
The ‘–r square’ option in PLINKv9 was used for the 1KGP and UK10K data while forcing the correct effect allele to be used for the correlation calculation.
The use of alternative reference panels for the r correlation matrix presents a few issues regarding the matching of variant id’s with the IBD dataset to ensure the order and number of variants included in the analysis remains constant and that variants are not lost due to naming differences. This especially affects insertion/deletions, which can be named using a number of different conventions from rsids to ‘chr:position:I/D’ format. 

In our analysis, the scripts create_r2_matrix.py and run_create_r2.sh were used to create the correlation matrix for the full IBD dataset.
**See the README_3_correlation.md file for details on how these files were created. 

**The Z-score file (README_4_inputfiles.md)**

The is a linear file with no header and one row per variant and two columns: snp id and directional Z-score. 
The calculation of the Z-scores and the method to ensure that the z-score and correlation information are in the same direction see README_2_alleles.md

```
$ head HD71_1.z
rs791336 2.37497791971
rs791337 2.32597283332
rs811008 2.13098958761
rs791338 2.11735181982
```

**The priors file**

Create a region.k file for each region containing one line with the prior probabilities for the region having 1, 2, 3, n signals. The number of priors dictates the maximum number of causal variants (independent signals) considered by FINEMAP.
We calculated the priors for the number of causal variants to reflect the proportion of regions with one, two and three independent signals (or causal variants) in the IBD regions included in our analysis. For example, for a region where we will consider a maximum of 3 signals the priors might be:
```
0.816 0.15 0.033
```

**Create a parameter file containing specifics for the FINEMAP analysis **
The **make_inputmaster_c3.sh** script, provided here, creates a parameter file to tell FINEMAP where the files are and what they are called for each input and output files. 
The file can contain multiple lines with details for each region that requires running FINEMAP to determine likely causal variants. The 'z', 'ld' and 'k' files are the Z-score, pairwise correlation and priors file mentioned above, respectively. The 'snp', 'config' and 'log' files are all outputs from FINEMAP. The 'n-ind' is the total number of samlples in the dataset.
```
z;ld;snp;config;k;log;n-ind
/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior/HD110_1.z;/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior/HD110_1.ld;/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior_test/prior1/HD110_1.snp;/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior_test/prior1/HD110_1.config;/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior_test/prior1/HD110_1.k;/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior_test/prior1/HD110_1.log;68428
/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior/HD150_1.z;/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior/HD150_1.ld;/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior_test/prior1/HD150_1.snp;/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior_test/prior1/HD150_1.config;/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior_test/prior1/HD150_1.k;/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior_test/prior1/HD150_1.log;68428
```  
 
**Commands to run finemap across multiple regions using a parameter file (--in-files option)**
Can be run over multiple directories using a shell wrapper as in run_finemap_actual_67k_master.sh to run different parameters in the same script. 

```
./finemap --sss --in-files ${datadir}/input_master_c3alt_cln_multi_IBD --log --prior-k --n-configs-top 100 --corr-threshold 0.98
```
* **in-files** this is the directory and name of the parameter file.
* **prior-k** tells finemap to use the region.k files for the priors.
* **n-configs-top 100** only return the top 100 by probability configurations of causal variants.
