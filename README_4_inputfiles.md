**These scripts use the output files from README 2 and 3 to create the input files to run FINEMAP: z-score and .k files**
* **get_zscores2.sh. Required the .key file for each region and the snporder list.**
```
#!/bin/sh
#for i in `cat cred_set_regionnless50_116.list` ;  do
for i in `cat /psorsfm/IBD/caviarbf/phase2_signals.list` ;  do
key=/psorsfm/IBD/caviarbf/${i}_dosage_results_in1KGP6.key
bimorder=/psorsfm/IBD/caviarbf/1KGP/snporder/HD${i}.snporder_IBDformat6.list
#for regions with signals that belong to different traits this will have to be run for each trait to get the correct pvalues and use the correct key file
#snp id in $2 as rsid - need to match this to the psoriasis results

#the files are the same and do not really need to be duplicated for different LD dataset
outdir=/psorsfm/IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior
awk -vkey="${key}" '
BEGIN{
OFS=" "
FS=" "
while(getline<key){
 zscoreA1[$1]=$9 #based on A1 - this is OR in the new files
 zscoreA2[$1]=$10  #based on A2
 zscoreR[$1]=$11 #based on risk
 variance[$1]=$13
 snp[$1]=$1
 snp1kgp[$1]=$12
 }
}
{
rsid=$1
IBDsnp=snp[rsid]
print IBDsnp, zscoreA2[rsid]
}' ${bimorder} > ${outdir}/HD${i}.z
done
```

* **Make priors files for all regions - .k files**
**run_makek.sh - this basically takes the prior probability for the number of independent signals in the region and creates a separate file for each region. If you want to run the analysis allowing up to 3 independent signals, then there will be three priors - one for: 1 signal, 2 signals and 3 signals.**
```
0.816 0.15 0.033
```
Scripts: 
```
#!/bin/sh
#export PATH=/software/python-2.7.8/bin:$PATH
#export LD_LIBRARY_PATH=/software/python-2.7.8/lib:$LD_LIBRARY_PATH
#source ~/virtualenvs/myenv/bin/activate

files=/IBD/finemap_Nov14/paintor/finemap
workdir=/psorsfm/IBD/FINEMAPv2/actual_Feb16
#for i in `cat ${workdir}/r_files_67k/signallist_final_280616` ; do
for i in `cat /psorsfm/IBD/caviarbf/phase2_signals.list` ; do
    
region=${i}

#thie following takes the required priors from the file provided by Christian Benner
#file=${files}/prior_k.file
#the following file is created from the proportion of signals seen in the actual IBD data included in the analysis
file=/psorsfm/IBD/FINEMAPv2/actual_Feb16/r_files_67k/actual_signals_priors
#for f in  r_files_379 r_files_2k ; do
for f in r_files_67k ; do

awk '{print $1, $2, $3}' ${file} > ${workdir}/${f}/c3_altprior/HD${i}.k
awk '{print $1, $2}' ${file} > ${workdir}/${f}/c2_new/HD${i}.k
awk '{print $1}' ${file} > ${workdir}/${f}/c1_new/HD${i}.k

done
   done
```
* ** Create the files to run FINEMAP- make_inputmaster.sh**
This has been uploaded as a separate file. The output is a file containing the parameters and file locations of all the input files for all regions being run. This saves you having to place all the files for a specific run in the same directory and duplicating files.
You should now have all the files to run FINEMAP. 
