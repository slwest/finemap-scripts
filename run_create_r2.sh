export PATH=/software/python-2.7.8/bin:$PATH
export LD_LIBRARY_PATH=/software/python-2.7.8/lib:$LD_LIBRARY_PATH

source ~/virtualenvs/myenv/bin/activate

for i in `cat IBD/caviarbf/cred_set_regionnless50_116.list` ; do
   region=(${i//_?/ })
    filedir=/IBD/Dosage_files
    cd ${filedir}/unpack
    tar -zxvf ${filedir}/HD${region}.dosage.tar.gz

cd IBD/FINEMAPv2/actual_Feb16
for f in 67k ; do
  python create_r2_matrix.py ${i} _${f}

done
done
rm ${filedir}/unpack/ibd_rel5_EU_*

deactivate
