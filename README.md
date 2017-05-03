# finemap-scripts
Scripts used for the formatting of files and running of FINEMAP and downstream analysis

# Running FINEMAP
First create a region.k file for each region with the priors for the region having 1, 2, 3, n signals. The number of priors dictates the maximum number of causal variants (independent signals) considered by FINEMAP.

**make_inputmaster.sh**
This creates a parameter file to tell FINEMAP where the files are and what they are called for each input and output files.

**run_finemap_actual_67k_master.sh** 
```
./finemap --sss --in-files ${datadir}/input_master_c3alt_cln_multi_${i} --log --prior-k --n-configs-top 100 --corr-threshold 0.98
```
* **in-files** this is the directory and name of the parameter file.
* **prior-k** tells finemap to use the region.k files for the priors.
* **n-configs-top 100** only return the top 100 by probability configurations of causal variants.
