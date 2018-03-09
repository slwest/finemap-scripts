import numpy as np
import pandas as pd
import scipy.stats as st
import sys
import re
import datetime
#this script produces the zscore files for each hdsignal and produces a file with one column of zscores  
    
# precompile the regex to save time
RE_MARKER=re.compile(r"chr\d{1,2}:\d+",re.IGNORECASE)

def get_trait(region,trait_dict, resdir):
    for line in trait_dict:
    #the line is the entry in the dictionary
        if line[0]==region:
            mytrait = line[1]
        #    if mytrait=='IBD':
        #        mytrait = 'CD' #there are only P value for CD or UC not IBD
            myregion = region.split("_")[0] #the signal numbers are not included
            filename = resdir + "HD" + myregion + "/ibd_rel5_EU_clean_HD" + myregion + "_HG19.core." + mytrait + ".assoc.dosage"
    return filename

def get_chr(region,trait_dict):
    for line in trait_dict:
        if line[0]==region:
            mychr = str(line[2])
    return mychr

def get_pos(snp,dbsnp_dict, myChr):
    if RE_MARKER.match(snp):
        return snp
    if snp in dbsnp_dict:
        _chr=dbsnp_dict[snp][0]
        _pos=dbsnp_dict[snp][1] 
        return "chr{}:{}".format(_chr, _pos)
#        return "chr" + _chr + ":" + str(_pos)
    return snp #if the 1st 2 fail return the snp as is
 
def main():
    workdir = "/lustre/scratch114/projects/psorsfm/IBD/caviarbf/"
    resdir = "/lustre/scratch113/teams/barrett/users/sw20/IBD/assoc_files/"
    ldfiles = "/lustre/scratch113/teams/barrett/users/sw20/IBD/finemap_Nov14/paintor/1KGP_data/"
    infile = "credible_set_march30_hdr.txt"
    dbsnp = "/lustre/scratch113/teams/barrett/users/sw20/IBD/finemap_Nov14/paintor/dbsnp/chr"
    region = sys.argv[1]
    
    #get the data for the SNPs and add a column called hdsignal 
    credsnpdata = pd.read_table(workdir + infile, sep="\t")
    credsnpdata["hdsignal"] = credsnpdata['HD'].map(str) + "_" + credsnpdata['signal'].astype(str)
    print 'cred data', datetime.datetime.now()
    
    #retrieve the regions and the traits and remove duplicates
    trait_dict=credsnpdata[['hdsignal','trait_reassigned','Chr']].drop_duplicates().to_dict('split')['data']
    myChr = get_chr(region, trait_dict)
    mysnpsfile = open(ldfiles + "chr" + myChr + ".phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.bim")
    print 'get chr', datetime.datetime.now()
    
    dbsnp_db = pd.read_table(dbsnp + myChr + "_dbsnp_142.chrpos.list.gz", sep="\t", compression="gzip")
    #dbsnp_dict=dbsnp_db[['CHROM','POS', 'ID']].to_dict('split')['data'] 
    dbsnp_dict = {}
    for index, row in dbsnp_db.iterrows():
        dbsnp_dict[row.ID] = [row.CHROM, row.POS]
    print 'create dbsnp dict', datetime.datetime.now()
    
   #def get_pos(snp,dbsnp_dict, myChr):
    #    chrpos = None
    #    #add a query in here to check the format of the snp id - start chr or rs - to determine what to match to    
    #    for line in dbsnp_dict:
    #        presnp = snp.split(":")[0]   #print the first part of the ":" delimited text
    #        dbsnp_chrpos = "chr" + str(line[0]) + ":" + str(line[1])
    #        if line[2]==snp:
    #            chrpos = str(line[0]) + ":" + str(line[1])
    #        if presnp== "chr" + myChr:
    #            if dbsnp_chrpos==snp:  
    #                chrpos = snp  
    #    if not chrpos:    #this means if chrpos is empty
    #	chrpos = snp
    #    return chrpos
    
    #myfilename=get_trait(region, trait_dict)
    #get the assoc results for each region for the correct trait and create zscore
    myassoc = open(get_trait(region, trait_dict, resdir)) #these are the assoc results
    lines=myassoc.readlines()
    mysnps={}
    for variant in mysnpsfile:
        rsid=variant.split("\t")[1].strip()
        mysnps[rsid]=True
    
    #mysnps = map(str.strip,mysnpsfile.readlines()) #this is the list of SNPs that are in 1KGP
    #print "mysnps"
    #print mysnps
    mysnpsfile.close()
    header = lines[0].split()
    legend={}
    for index,name in enumerate(header):
        legend[name.strip()]=index #creates a mapping between name and index
    signals={}    
    for line in lines[1:]:
    #    print "line"
    #    print line
        columns= line.split() #the assoc file is multiple white space
        snp = columns[legend['SNP']]
        snp2 = snp
        chrompos = get_pos(snp,dbsnp_dict, myChr) 
    #    if snp not in mysnps: #this bit checks whether the snps are present in the mysnps list
    #        continue
        hds = region
        freq = float(columns[legend['FRQ']])
        if columns[legend['OR']] != 'NA':
            OR = float(columns[legend['OR']])
            pvalue = (float(columns[legend['P']]))
            if OR <1:
                zscore_A1 = st.norm.ppf(pvalue/2)
                zscore_A2 = -st.norm.ppf(pvalue/2)
            else:
                zscore_A1 = -st.norm.ppf(pvalue/2)
                zscore_A2 = st.norm.ppf(pvalue/2)

            zscore_risk = -st.norm.ppf(pvalue/2)
            variance_A2 = 2*(1-freq)*(1-(1-freq))
            if freq > 0.5:
                minor_allele = columns[legend['A2']]
                maj_allele = columns[legend['A1']]
                if len(minor_allele) > 1:
                    snp2 = str(chrompos) + ":I"
                if len(maj_allele) > 1:
                    snp2 = str(chrompos) + ":D" #with reference to the minor allele
                OR_minor = 1/OR
                maf = 1-freq
                variance_minor = 2*maf*(1-maf)
            if freq <= 0.5:
                minor_allele = columns[legend['A1']]
                maj_allele = columns[legend['A2']]
                if len(minor_allele) > 1:
                    snp2 = str(chrompos) + ":I"
                if len(maj_allele) > 1:
                    snp2 = str(chrompos) + ":D" #with reference to the minor allele
                OR_minor = OR
                variance_minor = 2*freq*(1-freq)
            if OR_minor <1:
                zscore_min =  st.norm.ppf(pvalue/2)
            else:
                zscore_min = -st.norm.ppf(pvalue/2)
        if snp2 not in mysnps:  #removes SNPs that are not present in the 1KGP dataset with the new ids
            continue
        if columns[legend['OR']] != 'NA': #removes any variant with NA as OR, which is all MAF<0.01 SNPS
            if hds in signals:
                signals[hds].append('{} {} {} {} {} {}\n'.format(line.strip(), myChr, zscore_A2, zscore_risk,snp2, variance_A2)) 
            else:
        #if hds key is not in the dict then add it in
                signals[hds]=['{} {} {} {} {} {}\n'.format(line.strip(), myChr, zscore_A2, zscore_risk, snp2,  variance_A2)]    
       
#        print 'finito', datetime.datetime.now()
    
    for hds, details in signals.iteritems():
    #    print hds
        f=open(workdir + hds + "_dosage_results_in1KGP6.key",'w')
        #f=open(workdir + hds + "_dosage_results.key",'w')
        header.extend(["chr", "Zscore_A2","Zscore_risk", "newSNPid",  "variance_A2"])
        f.write(" ".join(header)+"\n")
        for detail in details:
            f.write(detail)
        f.close()
    
    myassoc.close()
    #mysnps.close()


if __name__ == "__main__":
    main()
