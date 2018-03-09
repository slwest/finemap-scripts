* Create a key file to ensure alleles match between Z-score and correlation files and to hold the z-scores
**Scripts to generate the key files containing alleles and Z-scores generate the following files:**
* **1. phase2_signals.list** 
* **2. ${i}.phase2_signals.txt (where ${i} is IBD, CD or UC)**
* **3. chr${i}.phase2_signals.txt (where ${i} is chr 1-22)**
* **4. List of variants for each region: ibd_rel5_EU_clean_HD${region}_HG19.core.filtered.dosage.snps**
* **5. List of z-scores by allele and variant ID key: {region}_dosage_results_in1KGP6.key**

* **1. Create the list of regions to be analysed**
e.g.
```
:/IBD/caviarbf$ cat phase2_signals.list 
1_1
45_1
56_1
66_1
67_1
71_1
78_1
83_1
84_1
100_1
109_1
110_1
112_1
150_1
162_1
180_1
181_1
63_1
```

* **2. Then get the traits for the signals and separate by CD, UC or IBD phenotype:**
This is only required as for this dataset, the analysis was split by the different phenotypes and so the z-scores need to be calculated from the correct analysis for the association signal (region).
**get_region_traits.sh**
```
#!/bin/sh

#signal_list=cred_set_regionnless50_116.list
signal_list=phase2_signals.list

for i in CD IBD UC; do
    awk -v trait="${i}" -v regions="${signal_list}" '
    BEGIN{
    while(getline<regions){
        region[$1]=$1
    }
}
($3==trait){
hdr=$2
if(region[hdr]!=""){
    print $2
}
}' HDR.traits > ${i}.phase2_signals.txt 
     done
```

The HDR.traits input file lists each region and the phenotype it is associated with:    
```
/IBD/caviarbf$ head HDR.traits
     65 100_1 CD
      1 102_1 IBD
      2 10_2 CD
      4 104_1 IBD
```

* **3. Get the Chromosomes for the each associated region**
The 1KGP data is in PLINK format and is separated by chromosome - this step is to enable queries of the correct 1KGP files for the LD data and SNP order.
**get_region_chrs.sh**
```
#!/bin/sh

#signal_list=cred_set_regionnless50_116.list
signal_list=phase2_signals.list

for i in `seq 1 22` ; do
    awk -v chr="${i}" -v regions="${signal_list}" '
    BEGIN{
    while(getline<regions){
        region[$1]=$1
    }
}
($3==chr){
hdr=$2
if(region[hdr]!=""){    #this line ensures that the signals are included in the list of signals to be analysed
    print $2
}
}' HDR.chr > chr${i}.phase2_signals.txt 
    
     done
```
The top of the input file HDR.chr looks like this; column 3 is the chromosome:
```
     65 100_1 8
      2 10_2 1
      1 102_1 9
      4 104_1 9
     29 104_2 9
```

* **4. Create the SNP lists.**
The **get_dosage_header.sh** script unpacks the dosage files for each associated region in phase2_signals.list and prints the list of variants contained in the file. The order of this list determines the order of the correlation file and thus the order of the Z-score file must be the same too. 
```
#!/bin/sh

#for i in `cat cred_set_regionnless50_116.list` ; do
for i in `cat phase2_signals.list` ; do

    region=(${i//_?/ })

filedir=/IBD/Dosage_files
cd ${filedir}/unpack
tar -zxvf ${filedir}/HD${region}.dosage.tar.gz

awk '
(NR==1){
c=2
while(c<=NF){
        snp[c]=$c
        print snp[c]
        c=c+1
        }
}' ${filedir}/unpack/ibd_rel5_EU_clean_HD${region}_HG19.core.filtered.dosage > ${filedir}/ibd_rel5_EU_clean_HD${region}_HG19.core.filtered.dosage.snps

done

rm ${filedir}/unpack/ibd_rel5_EU_*
```

* **5. Create the .key (ref) file of Z-scores and variant IDs**
**run_create_keyfile.sh** - This script is a shell wrapper for a python script that outputs a .key file that contains the variant IDs as seen in the IBD study data, the variant IDs as they appear in 1000 genomes reference panel (1KGP) and the Z-scores for each variant.
```
export PATH=/software/python-2.7.8/bin:$PATH
export LD_LIBRARY_PATH=/software/python-2.7.8/lib:$LD_LIBRARY_PATH

homedir=/IBD/caviarbf/
resdir=/users/sw20/IBD/assoc_files/

source ~/virtualenvs/myenv/bin/activate

regions=/IBD/caviarbf/cred_set_regionnless50_116.list

for i in `cat phase2_signals.list` ; do #this is the sys arg used in the python script

python create_zscore_caviarin6.py ${i}
done

deactivate
```

The **create_zscore_caviarin6.py** script called by the above wrapper outputs the actual key file and contains calculations for calculating the z-score. This script has been added to the repository separately. 
The OR in the association results file is with reference to allele A1. 
The input file for this python script is the association results file from the dosage files for the dataset and are of the following format:
**ibd_rel5_EU_clean_HD177_HG19.core.UC.assoc.dosage**
```
         SNP  A1  A2     FRQ    INFO      OR      SE       P
 rs146525874  TA   T  1.0000  0.9527      NA      NA      NA
 rs182073174   A   G  0.9996  0.2697      NA      NA      NA
 rs113758151   A   G  0.9999  0.8624      NA      NA      NA
   rs1888484   C   G  0.9997  0.8278      NA      NA      NA
 rs138761116   T   C  1.0000  0.5107      NA      NA      NA
 chr21:40402489   C CTG  0.8127  1.0033  0.9667  0.0183  0.0639
```

**Z-score calculation example:**
```
if OR <1 ##of A1
    zscore_A1 = st.norm.ppf(pvalue/2)
    zscore_A2 = -st.norm.ppf(pvalue/2)
else:
    zscore_A1 = -st.norm.ppf(pvalue/2)
    zscore_A2 = st.norm.ppf(pvalue/2)
```
The output *.key file is space separated and has the following format (where A1 is the risk allele in line 1 and A2 is the risk allele in line 2):
```
SNP A1 A2 FRQ INFO OR SE P chr Zscore_A2 Zscore_risk newSNPid variance_A2
rs76393761   T   C  0.9081  0.8548  1.0715  0.0209 0.0009597 22 -3.3020820284 3.3020820284 rs76393761 0.16690878
rs56083406   G   A  0.9088  0.7700  0.9939  0.0221  0.7829 22 0.275541821359 0.275541821359 rs56083406 0.16576512
```

**scripts to create the correlation files follow in README_3_correlation.md 
