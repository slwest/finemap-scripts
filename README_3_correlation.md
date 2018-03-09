**the scripts in this file require that those in README_2 have already been run. 
**The following scripts use the .key file produced for each region to produce a SNP order file so that the correlation matrix for the 1KGP reference panel and the matrix for the actual IBD association data contain the same SNPs in the same order. This allows us to run and compare the reference panels for the anlaysis.
** Additional scripts included here are used to create the actual correlation matrices. 

* **1. Create the plink files in 1KGP dataset containing only the SNPs for the regions - force the allele to A2 from the key file.
* getLD_plink.sh
```
#!/bin/sh
module load hgi/plink/1.90b2f
#using a file with regions separated by chromosomes 
for i in 1 3 4 5 6 7 8 10 16 18 22 ; do
    #calculate LD
for f in `cat /psorsfm/IBD/caviarbf/phase2/chr${i}.phase2_signals.txt` ; do

key=/psorsfm/IBD/caviarbf/${f}_dosage_results_in1KGP6.key
home=/psorsfm/IBD/caviarbf/1KGP/LD

files=/IBD/finemap_Nov14/paintor/1KGP_data/chr${i}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR
snps=/psorsfm/IBD/caviarbf/info_files/HD${f}.snp_1KGPformat.list

#run make_snp.list.sh to determine the SNPs that need to be included in the analysis using the key file and the results from the conditional and standard assoc analysis
#the alleles are being forced to be the same as the key file which matches the direction of effect reflected in the z-score
plink --bfile ${files} --extract ${snps} --a1-allele ${key} 3 12 --allow-no-sex --make-bed --out ${home}/${f}_1KGP_snpsA2
done
done
```
**A2 in the key should now match A1 in the plink .bim file as this allele corresponds to the IBD dosage files used to generate the correlation matrix for the actual IBD data.
* **2. Calculate the pairwise r correlation in each region for 1KGP  
This version means that r was calculated from the files that you are using for the SNP order.
```
plink --bfile ${home}/${f}_1KGP_snpsA2 --a1-allele ${key} 3 12 --allow-no-sex --r square --out ${home}/${f}_1KGP_snpsldrA2
```
* **3. Change the field separator of the matrix to space separated - this outputs the region.ld file in the correct format to run FINEMAP.
**run_diagonal.sh calls diagonal_spacesep.py (provided as a separate file)  
```
#create a list of the LD files to change
ls *.ld | awk '{split($1,part,"."); print part[1]}' > file.list

export PATH=/software/python-2.7.8/bin:$PATH
export LD_LIBRARY_PATH=/software/python-2.7.8/lib:$LD_LIBRARY_PATH

homedir=/psorsfm/IBD/caviarbf
outdir=/psorsfm/IBD/caviarbf/FINEMAP/1KGP/c3

source ~/virtualenvs/myenv/bin/activate
#for i in `cat ${homedir}/cred_set_regionnless50_116.list` ; do
for i in `cat ${homedir}/phase2_signals.list` ; do

python diagonal_spacesep.py ${homedir}/1KGP/LD/${i}_1KGP_snpsldrA2.ld ${outdir}/HD${i}.ld
done
deactivate
```

* **4. Prepare files to create the correlation r files for the actual IBD dataset
**create the SNP order list. Requires the .key file for the region and the plink .bim file from 1KGP created in 2.
```
#!/bin/sh
#this script returns the SNPs in the same order as they are in the plink bim files in the 1KGP directory
#this is required to produce the LD and z-score files in the same order
#for i in `cat cred_set_regionnless50_116.list` ; do
for i in `cat /psorsfm/IBD/phase2_signals.list` ; do
key=${i}_dosage_results_in1KGP6.key
#snps=HD${i}.snp_IBDformat.list
#snp1KGP=HD${i}.snp_1KGPformat.list
#these are the bim files for the region created from the snp lists
bim=/psorsfm/IBD/caviarbf/1KGP/LD/${i}_1KGP_snpsA2.bim
awk -vkey="${key}" '
BEGIN{
FS=" " 
OFS="\t"
while(getline<key){
  IBDsnp[$12]=$1
  KGPsnp[$12]=$12 
  }
}
{
snp=$2
print IBDsnp[snp]
}' ${bim} > 1KGP/snporder/HD${i}.snporder_IBDformat6.list

done
```

* **5. Create the correlation files from the actual IBD data dosage files in the order of the snporder list in 4.
**run_create_r2.sh, which calls create_r2_matrix.py (provided as separate file)
The output file is in the correct format to run as input for FINEMAP
```
export PATH=/software/python-2.7.8/bin:$PATH
export LD_LIBRARY_PATH=/software/python-2.7.8/lib:$LD_LIBRARY_PATH

source ~/virtualenvs/myenv/bin/activate

#for i in `cat /psorsfm/IBD/FINEMAPv2/actual_Feb16/2plus_final_loci.list` ; do 
for i in `cat /psorsfm/IBD/caviarbf/phase2_signals.list` ; do
   region=(${i//_?/ })
    filedir=/lustre/scratch114/projects/psorsfm/IBD/Dosage_files
    cd ${filedir}/unpack
    tar -zxvf ${filedir}/HD${region}.dosage.tar.gz

cd /psorsfm/IBD/FINEMAPv2/actual_Feb16
for f in 67k ; do
       python create_r2_matrix.py ${i} _${f}

done
done
rm ${filedir}/unpack/ibd_rel5_EU_*

deactivate
```

 
