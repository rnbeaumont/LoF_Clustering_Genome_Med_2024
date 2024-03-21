#!/usr/bin/env python
# script to format clinvar VCF file
import pandas as pd
import numpy as np
import json

files=json.load(open("file_locations"))
clinvar_vcf_file = files["CLINVAR_VCF"]
clinvar_formatted_file = files["CLINVAR_FORMATTED"]

# read in the vcf file
clinvar=pd.read_csv(clinvar_vcf_file,sep="\t",skiprows=27)
# split out the info field
INFO=clinvar.INFO.str.split('[;]',expand=True,regex=True)
# create array to track pathogenicity status
clinvar["path"]=False
clinvar["benign"]=False
for col in INFO:
    arr=INFO[col].str.split('[=|,]',expand=True,regex=True)
    clinvar.loc[((arr[0]=="CLNSIG") & 
        ((arr[1]=="Pathogenic") | (arr[1]=="Pathogenic/Likely_pathogenic") | 
        (arr[1]=="Pathogenic/Likely_risk_allele") | (arr[1]=="Likely_pathogenic") | 
        (arr[1]=="Likely_risk_allele"))),"path"]=True
    clinvar.loc[((arr[0]=="CLNSIG") & 
        ((arr[1]=="Benign") | (arr[1]=="Benign/Likely_benign") | 
        (arr[1]=="Conflicting_interpretations_of_pathogenicity") | 
        (arr[1]=="Likely_benign"))),"benign"]=True

clinvar["variantID"]=clinvar["#CHROM"].astype(str)+":"+clinvar["POS"].astype(str)+":"+clinvar["REF"]+":"+clinvar["ALT"]
pathogenic=clinvar["variantID"].loc[clinvar["path"]==True]
benign=clinvar["variantID"].loc[clinvar["benign"]==True]
pathogenic=pathogenic.loc[~pathogenic.isin(clinvar["variantID"].loc[clinvar["benign"]==True])]
pathogenic.to_csv(clinvar_formatted_file,sep="\t")