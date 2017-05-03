#!/bin/sh

datadir=FINEMAP_files #location of the master parameter files accessed by "--in-files" argument below


export LD_LIBRARY_PATH=<path to gcc-4.9.1/lib64>; <path to linux2.6-glibc2.3-x86_64/lib>
#can reduce the files based on the correlation threshold using the following
for i in IBD CD UC; do
 
#this is the 11th July 2016 version that now allows users to set the input file names individually
program=programs/finemap_v1.1_x86_64
cd ${program}

# --prior-k0 default is 0

#for c=2
./finemap --sss --in-files ${datadir}/input_master_c2alt_cln_multi_${i} --log --prior-k --n-configs-top 100 --corr-threshold 0.98
./finemap --sss --in-files ${datadir}/input_master_c2alt_cln_sgl_${i} --log --prior-k --n-configs-top 100 --corr-threshold 0.98
./finemap --sss --in-files ${datadir}/input_master_c2alt_cln_phase2_${i} --log --prior-k --n-configs-top 100 --corr-threshold 0.98

done
