#!/usr/bin/env python
# script to parse variant annotations

import pandas as pd
import argparse
import json

def parse_arguments():
     parser = argparse.ArgumentParser(description='A script for making annotations from a tabulated VEP file')
     parser.add_argument('--chr',  dest='chr', required=True, help="chromosome")
     args = parser.parse_args()
     return args

def read_file(chrom):
     annots=files["VEP_ANNOTATIONS"].split("*")
     data = pd.read_csv(annots[0] + str(chrom) + annots[1], header=0, sep= "\t", compression='gzip', na_values=' ')
     return data

# read the annotations for the requested chromosome
args=parse_arguments()
files=json.load(open("file_locations"))
chrom = str(args.chr)
data=read_file(chrom)

#make files
f=open(files["SNPS_ANNOTATED"] + str(chrom), 'w')
for gene, gene_data in data.groupby("Feature"):
    gene_counter = 1 
    all_snps_tmp=""
    
    gene_data["BP"]= gene_data["Location"].str.split(":").str[1] # split out the location of the SNP
    gene_tmp = str(gene_data["SYMBOL"].iloc[0]) + "_"  + str(gene_data["Feature"].iloc[0]) # create string of gene name and transcript
    
    # create dataframe of consequences we're interested in
    gene_data_tmp = gene_data.loc[gene_data['Consequence'].str.contains("missense_variant", case=False,na=False) | gene_data['Consequence'].str.contains("synonymous_variant", case=False,na=False) |  gene_data['Consequence'].str.contains("inframe_insertion", case=False,na=False) | gene_data['Consequence'].str.contains("inframe_deletion", case=False,na=False) | gene_data['LoF'].str.contains("HC", case=False,na=False) | gene_data['LoF'].str.contains("LC", case=False,na=False)]
    if len(gene_data_tmp) ==0:     # skip if there are no variants left for this transcript
         continue
    first_snp = str(min(gene_data["BP"])).split("-")[0]
     
    # relabel missense variants based on CADD
    gene_data_tmp.loc[gene_data_tmp['Consequence'].str.contains("missense_variant", case=False,na=False) & (gene_data_tmp["CADD_PHRED"] >= 30), "Consequence"] = "missense_cadd30plus_variant"
    gene_data_tmp.loc[gene_data_tmp['Consequence'].str.contains("missense_variant", case=False,na=False) & (gene_data_tmp["CADD_PHRED"] >= 25), "Consequence"] = "missense_cadd25plus_variant"
    gene_data_tmp.loc[gene_data_tmp['Consequence'].str.contains("missense_variant", case=False,na=False), "Consequence"] = "missense_cadd24less_variant"
    
    # relabel inframe insertions/deletions based on CADD score
    gene_data_tmp.loc[(gene_data_tmp['Consequence'].str.contains("inframe_deletion", case=False,na=False) | gene_data_tmp['Consequence'].str.contains("inframe_insertion", case=False,na=False)) & (gene_data_tmp["CADD_PHRED"] >= 30), "Consequence"] = "inframe_indel_cadd30plus_variant"
    gene_data_tmp.loc[(gene_data_tmp['Consequence'].str.contains("inframe_deletion", case=False,na=False) | gene_data_tmp['Consequence'].str.contains("inframe_insertion", case=False,na=False))  & (gene_data_tmp["CADD_PHRED"] >= 25), "Consequence"] = "inframe_indel_cadd25plus_variant"
    gene_data_tmp.loc[(gene_data_tmp['Consequence'].str.contains("inframe_deletion", case=False,na=False) | gene_data_tmp['Consequence'].str.contains("inframe_insertion", case=False,na=False)), "Consequence"] = "inframe_indel_cadd24less_variant"
    
    # relabel synonymous
    gene_data_tmp.loc[gene_data_tmp['Consequence'].str.contains("synonymous_variant", case=False,na=False), "Consequence"] = "synonymous"
    
    # relabel LoF variants
    gene_data_tmp.loc[gene_data_tmp['LoF_info'].str.contains("50_BP_RULE:FAIL",  case=False,na=False), "LoF_info"] = "FAIL"
    gene_data_tmp.loc[gene_data_tmp['LoF_info'].str.contains("50_BP_RULE:PASS",  case=False,na=False), "LoF_info"] = "PASS"
    gene_data_tmp.loc[gene_data_tmp['LoF'].str.contains("HC", case=False,na=False) & gene_data_tmp['VARIANT_CLASS'].str.contains("SNV", case=False,na=False), "Consequence"] = "LoF_SNV"
    gene_data_tmp.loc[gene_data_tmp['LoF'].str.contains("HC", case=False,na=False) & gene_data_tmp['VARIANT_CLASS'].str.contains("SNV", case=False,na=False), "LoF"] = "LoF_SNV"
    gene_data_tmp.loc[gene_data_tmp['LoF'].str.contains("HC", case=False,na=False), "Consequence"] = "LoF_HC"
    gene_data_tmp.loc[gene_data_tmp['LoF'].str.contains("LC", case=False,na=False), "Consequence"] = "LoF_LC"
    
    first_bp_rule = float(str(gene_data_tmp["BP"][gene_data_tmp["LoF_info"] == "FAIL"].min()).split("-")[0])    
    
    for index, row in gene_data_tmp.iterrows():
        snp_tmp = str(row["Uploaded_variation"])
        consq_tmp = str(row["Consequence"])
        bp_rule = str(row["LoF_info"])
        if float(row["STRAND"])==1:
              if ("LoF" in consq_tmp) & ((bp_rule=="FAIL") | (float(str(row["BP"]).split("-")[0]) >= first_bp_rule)):
                   consq_tmp+="_50bp_fail"
              elif ("LoF" in consq_tmp) & ((bp_rule=="PASS") | (float(str(row["BP"]).split("-")[0]) < first_bp_rule)):
                   consq_tmp+="_50bp_pass"
        elif float(row["STRAND"])==-1:
              if ("LoF" in consq_tmp) & ((bp_rule=="FAIL") | (float(str(row["BP"]).split("-")[0]) <= first_bp_rule)):
                   consq_tmp+="_50bp_fail"
              elif ("LoF" in consq_tmp) & ((bp_rule=="PASS") | (float(str(row["BP"]).split("-")[0]) > first_bp_rule)):
                   consq_tmp+="_50bp_pass"
        
        if gene_counter == len(gene_data_tmp):
            all_snps_tmp += snp_tmp 
        else:
            all_snps_tmp += snp_tmp + ","
        f.write(str(snp_tmp)+"\t"+str(gene_tmp)+"\t"+str(consq_tmp)+"\n")
        
        gene_counter += 1

f.close()