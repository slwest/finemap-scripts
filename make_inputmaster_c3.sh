##script to create a parameter file containing all the regions to be run using finemap providing all the locations of the files and the file header: z;ld;snp;config;k;log;n-ind

pheno=CD
number=52541
#for i in 107 123 165 176 60 73 77 ; do 
#pheno=IBD
#number=68428
#for i in 23 ; do
#pheno=UC
#number=48245
#for i in 39 4 91 ; do

#head -1 input_master_All > input_master_${pheno}
#head -1 input_master_All > input_master_c3alt_actual_${pheno}
#head -1 input_master_All > input_master_c3alt_actual_multi_${pheno}
head -1 input_master_All > input_master_c3alt_actualp2_${pheno}

#for i in 200 1000 3652 5000 8000 10000 15000 20000 30000 ; do
#for i in 67k 2k 379 ; do
for i in 67k ; do
  
    #mkdir IBD/FINEMAPv2/actual_Feb16/r_files_${i}/c3
awk  -v size="${i}" -v phen="${pheno}" -v num="${number}" '
BEGIN{
FS=";"
OFS=";"
}
{
region=$1
#the multi signal files require the "_1" to be added to the name
zfile="IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior/HD"region".z"
#zfile="IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior/HD"region"_1.z"
ld="IBD/FINEMAPv2/actual_Feb16/r_files_" size "/c3_altprior/HD" region ".ld"
#ld="IBD/FINEMAPv2/actual_Feb16/r_files_" size "/c3_altprior/HD" region "_1.ld"
snp="IBD/FINEMAPv2/actual_Feb16/r_files_" size "/c3_altprior/HD"region".snp"
#snp="IBD/FINEMAPv2/actual_Feb16/r_files_" size "/c3_altprior/HD"region"_1.snp"

config="IBD/FINEMAPv2/actual_Feb16/r_files_" size "/c3_altprior/HD"region".config"
#config="IBD/FINEMAPv2/actual_Feb16/r_files_" size "/c3_altprior/HD"region"_1.config"
kfile="IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior/HD"region".k"
#kfile="IBD/FINEMAPv2/actual_Feb16/r_files_67k/c3_altprior/HD"region"_1.k"

lag="IBD/FINEMAPv2/actual_Feb16/r_files_" size "/c3_altprior/HD"region".log"
#lag="IBD/FINEMAPv2/actual_Feb16/r_files_" size "/c3_altprior/HD"region"_1.log"
print zfile, ld, snp, config, kfile, lag, num 
}' IBD/caviarbf/phase2/${pheno}.phase2_signals.txt  >> input_master_c3alt_actualp2_${pheno}
# final_loci${pheno}.list >> input_master_c3alt_actual_multi_${pheno}
#final_locisgl${pheno}.list >> input_master_c3alt_actual_${pheno}

#final_loci${pheno}.list >> input_master_c3_${pheno}
#final_loci${pheno}.list >> input_master_${pheno}
#final_loci${pheno}.list >> input_master_c2_${pheno}
#the input file is only used to provide the correct number of lines for the output
done

