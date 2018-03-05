import pandas as pd
import numpy as np
import sys

#this script produces as subset of the dosage file dependent on the sample list and creates the r matrix in the same order as the snplist
filedir="IBD/Dosage_files/unpack/"
outdir="IBD/FINEMAPv2/actual_Feb16/r_files_67k/"
region=sys.argv[1] #this should be HD23_1 format
outfile="HD" + region + ".ld"
##split the region to return just the beginning
HDR=region.split("_")[0]
#ensure the order of the SNPs is the same as in the z-score file. The snporder list has 1 line per snp and 1 column of snp ids.
snplist="IBD/caviarbf/1KGP/snporder/" + "HD" + region + ".snporder_IBDformat6.list"
snpfile = [line.strip() for line in open(snplist, 'r')]

dosage=pd.read_table(filedir + "ibd_rel5_EU_clean_HD" + HDR + "_HG19.core.filtered.dosage", sep="\t")
#the dosage file is sample_id	snp_1	snp_2	snp_n
#one column per snp showing the dosage of the allele 
myarray=dosage.as_matrix(columns=snpfile)

rarray=np.corrcoef(myarray, rowvar=0)
rarray[rarray>1]=1
rarray[rarray<-1]=-1
##rowvar 0= variables are columns (snps) and obervations as rows (dosages) - the output is a square array of r2 values. 
#save the results into a space delimited file
np.savetxt(outdir + outfile , rarray, fmt='%.16g', delimiter=" ", comments='') 
