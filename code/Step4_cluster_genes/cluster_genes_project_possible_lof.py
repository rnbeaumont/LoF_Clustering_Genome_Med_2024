#!/usr/bin/env python
# python script for clustering variant positions within genes
import pandas as pd
import numpy as np
import math
from sklearn.cluster import KMeans  # https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_digits.html#sphx-glr-auto-examples-cluster-plot-kmeans-digits-py
from sklearn.cluster import AgglomerativeClustering # https://scikit-learn.org/stable/auto_examples/cluster/plot_digits_linkage.html#sphx-glr-auto-examples-cluster-plot-digits-linkage-py
from sklearn import decomposition   # for PCA
from sklearn.mixture import GaussianMixture # https://scikit-learn.org/stable/modules/generated/sklearn.mixture.BayesianGaussianMixture.html
import plotly.figure_factory as ff
import plotly.express as px
import json

files=json.load(open("file_locations"))
# set temporary testing defaults. TODO change these to read in the argumenst
variants_file = files["VARIANTS_FILE_PREFIX"]+"_all_variants"
variants_file_pseudotranscript = variants_file+"_single_pseudo_transcript_all_variants"
output_prefix = files["FIGURES_PREFIX"]
gencc_file = files["GENCC"]
gene_names = files["GENE_POS"]
mane_transcripts = files["MANE_TRANSCRIPTS"]
clinvar_file = files["VARIANTS_FILE_PREFIX"]+"_clinvar_variants"
clinvar_file_pseudotranscript = files["VARIANTS_FILE_PREFIX"]+"_single_pseudo_transcript_clinvar_variants"
possible_lof_file = files["VARIANTS_FILE_PREFIX"]+"_possible_lof_positions.txt"
simulated_genes_file_prefix = files["VARIANTS_FILE_PREFIX"]+"_simulated_genes_n_quantile_"
output_prefix="plots/test_bayesian/test_output_clinvar_possible_lof"
n_break_min=4
n_break_max=7   # number of quantiles to break the gene into
n_clust_min=4
n_clust_max=10   # maximum number of clusters to consider

run_mane=1
run_pseudo=0    # NOTE: pseudo not implemented yet TODO
group_lof=1
indel_thresh=1
run_sims=0

def read_gene_names(gene_names,col_list):    # Read gene names from file
    ensembl=pd.read_csv(gene_names,sep="\t",usecols=col_list)
    gene_dict=ensembl.set_index(col_list[0]).to_dict()["Gene name"]
    return gene_dict

genes=read_gene_names(gene_names,["Transcript stable ID","Gene name"])

def read_variants(variants_file,n_break):   # function to read in data frame and add annotation super-groups
    super_groups={
        "LoF_HC": "lof",
        "LoF_HC_50bp_fail": "lof",
        "LoF_HC_50bp_pass": "lof",
        "LoF_LC": "lof",
        "LoF_LC_50bp_fail": "lof",
        "LoF_LC_50bp_pass": "lof",
        "LoF_SNV": "lof_snv",
        "LoF_SNV_50bp_fail": "lof_snv",
        "LoF_SNV_50bp_pass": "lof_snv",
        "missense_cadd24less_variant": "missense",
        "missense_cadd25plus_variant": "missense",
        "missense_cadd30plus_variant": "missense",
        "inframe_indel_cadd24less_variant": "indel",
        "inframe_indel_cadd25plus_variant": "indel",
        "inframe_indel_cadd30plus_variant": "indel"
    }
    variants=pd.read_csv(variants_file,sep=" ")
    variants["annotation_grouped"] = variants.annotation
    variants=variants.replace({"annotation_grouped":super_groups})
    variants["quantile"]=pd.cut(variants.position,[item/float(n_break) for item in range(0,n_break+1)],include_lowest=True,labels=range(1,n_break+1))   # split positions into n_break groups
    variants_indels=variants[variants["annotation_grouped"].isin(["indel"])]
    variants=variants[-(variants["annotation_grouped"].isin(["indel"]))]
    return [variants,variants_indels]

def read_clinvar_variants(clinvar_file,n_break): # function to read in clinvar variants and make them compatible with other variants
    clinvar_vars=pd.read_csv(clinvar_file,sep=" ")
    clinvar_vars["annotation_grouped"]="clinvar_pathogenic"
    clinvar_vars["quantile"]=pd.cut(clinvar_vars.position,[item/float(n_break) for item in range(0,n_break+1)],include_lowest=True,labels=range(1,n_break+1))
    # work out which variants are indels so we can run with and without them
    clinvar_vars["a1"]=clinvar_vars.variantID.str.split(":",expand=True)[2]
    clinvar_vars["a2"]=clinvar_vars.variantID.str.split(":",expand=True)[3]
    clinvar_vars["indel"]=np.where(((clinvar_vars.a1.str.len()>1) | (clinvar_vars.a2.str.len()>1)),1,0)
    clinvar_vars=clinvar_vars.drop(columns=["a1","a2"]) # drop temporary columns
    return clinvar_vars

def read_possible_lof_variants(possible_lof_file,n_break): # function to read in clinvar variants and make them compatible with other variants
    lof_vars=pd.read_csv(possible_lof_file,sep=" ")
    lof_vars["annotation_grouped"]="possible_lof_position"
    lof_vars["quantile"]=pd.cut(lof_vars.position,[item/float(n_break) for item in range(0,n_break+1)],include_lowest=True,labels=range(1,n_break+1))
    # work out which variants are indels so we can run with and without them
    lof_vars["variantID"]="gene"+lof_vars.position.astype(str)
    lof_vars["indel"]=0
    return lof_vars

def read_simulated_genes(simulated_genes_file_prefix,n_break):
    n_var_vec=list(range(5,21))
    n_var_vec.append(30)
    n_var_vec.append(50)
    n_var_vec.append(100)
    # generate vector to use as header
    header=["gene"]
    for quantile in range(1,n_break+1):
        header.append("quantile"+str(quantile))
    genes_df=pd.DataFrame(columns=header)
    for n_var in n_var_vec:
        genes=pd.read_csv(simulated_genes_file_prefix+"_n_var_"+str(n_var),sep=" ")
        genes["gene"]="gene_"+genes["gene"].astype(str)+"_nvar_"+str(n_var)
        genes_df=genes_df.append(genes)
    genes=pd.read_csv(simulated_genes_file_prefix+"_ukb_profile",sep=" ")
    genes_df=genes_df.append(genes)
    return genes_df

def filter_mane(variants,mane_transcripts,col_list):
    # read in ensembl transcripts
    ensembl=pd.read_csv(mane_transcripts,sep="\t",usecols=col_list)
    ensembl=ensembl.rename(columns={col_list[0]: "transcript",col_list[1]: "MANE_status"})
    ensembl.transcript=ensembl.transcript.str.split(".",expand=True)[0]
    ensembl=ensembl[ensembl.MANE_status=="MANE Select"]
    variants=variants[(variants.transcript.isin(ensembl.transcript))]
    return variants

def filter_indel_length(variants_indels):
    variants_indels["a1"]=variants_indels.variantID.str.split(":",expand=True)[2]
    variants_indels["a2"]=variants_indels.variantID.str.split(":",expand=True)[3]
    variants_indels=variants_indels[(variants_indels.a1.str.len()<5) & (variants_indels.a2.str.len()<5)]
    variants_indels=variants_indels.drop(columns=["a1","a2"])
    return variants_indels

def filter_lof_indel_length(variants):
    variants["a1"]=variants.variantID.str.split(":",expand=True)[2]
    variants["a2"]=variants.variantID.str.split(":",expand=True)[3]
    variants=variants[(variants["annotation_grouped"]=="lof") & (variants.a1.str.len()==1) & (variants.a2.str.len()==1)]
    return variants

def count_variants(variants,group_to_count,quantile_groups,pivot_index,pivot_column,mac=False):
    # count total number of variants per transcipt and annotation
    if mac:
        n_variants=variants.groupby(group_to_count)["mac"].sum().reset_index().rename(columns={"mac": "n_variants"})
    else:
        n_variants=variants.groupby(group_to_count).size().reset_index().rename(columns={0: "n_variants"})
    # merge with quantile_groups
    quantile_groups=pd.merge(quantile_groups,n_variants,on=group_to_count,how="outer")
    quantile_groups.counts=quantile_groups.counts/quantile_groups.n_variants
    quantile_groups.counts=quantile_groups.counts.fillna(0)
    quantile_groups.n_variants=quantile_groups.n_variants.fillna(0)
    # convert to table for clustering
    quantile_groups_table=quantile_groups.pivot_table(values="counts",index=pivot_index,columns=pivot_column,aggfunc=np.sum)
    quantile_groups_table=quantile_groups_table.fillna(0)   # convert NaN to 0
    return quantile_groups_table

def pca_decomposition(quantile_groups,quantile_groups_clinvar,quantile_groups_clinvar_no_indel,quantile_groups_possible_lofs,quantile_groups_indels,quantile_groups_snv,quantile_groups_simulated_genes,n_comp=5):
    pca=decomposition.PCA(n_components=n_comp)
    pca.fit(quantile_groups)
    pc_decomposition=pca.transform(quantile_groups)
    pc_decomposition=pd.DataFrame(pc_decomposition)
    pc_decomposition.columns=[col for col in range(1,n_comp+1)]
    ## project into other datasets
    # clinvar
    pc_decomposition_clinvar=pca.transform(quantile_groups_clinvar)
    pc_decomposition_clinvar=pd.DataFrame(pc_decomposition_clinvar)
    pc_decomposition_clinvar.columns=[col for col in range(1,n_comp+1)]
    # clinvar no indels
    pc_decomposition_clinvar_no_indel=pca.transform(quantile_groups_clinvar_no_indel)
    pc_decomposition_clinvar_no_indel=pd.DataFrame(pc_decomposition_clinvar_no_indel)
    pc_decomposition_clinvar_no_indel.columns=[col for col in range(1,n_comp+1)]
    # possible LoF
    pc_decomposition_possible_lofs=pca.transform(quantile_groups_possible_lofs)
    pc_decomposition_possible_lofs=pd.DataFrame(pc_decomposition_possible_lofs)
    pc_decomposition_possible_lofs.columns=[col for col in range(1,n_comp+1)]
    # indels
    pc_decomposition_indels=pca.transform(quantile_groups_indels)
    pc_decomposition_indels=pd.DataFrame(pc_decomposition_indels)
    pc_decomposition_indels.columns=[col for col in range(1,n_comp+1)]
    # lof and indel grouped
    pc_decomposition_snv=pca.transform(quantile_groups_snv)
    pc_decomposition_snv=pd.DataFrame(pc_decomposition_snv)
    pc_decomposition_snv.columns=[col for col in range(1,n_comp+1)]
    # simulated genes
    pc_decomposition_simulated_genes=pca.transform(quantile_groups_simulated_genes)
    pc_decomposition_simulated_genes=pd.DataFrame(pc_decomposition_simulated_genes)
    pc_decomposition_simulated_genes.columns=[col for col in range(1,n_comp+1)]
    # sort headers and pull out only pcs
    quantile_groups_tmp=quantile_groups.copy()
    quantile_groups_tmp_clinvar=quantile_groups_clinvar.copy()
    quantile_groups_tmp_clinvar_no_indel=quantile_groups_clinvar_no_indel.copy()
    quantile_groups_tmp_possible_lofs=quantile_groups_possible_lofs.copy()
    quantile_groups_tmp_indels=quantile_groups_indels.copy()
    quantile_groups_tmp_snv=quantile_groups_snv.copy()
    quantile_groups_tmp_simulated_genes=quantile_groups_simulated_genes.copy()
    for i in range(1,n_comp+1):
        quantile_groups_tmp["pc"+str(i)]=list(pc_decomposition[i])
        quantile_groups_tmp_clinvar["pc"+str(i)]=list(pc_decomposition_clinvar[i])
        quantile_groups_tmp_clinvar_no_indel["pc"+str(i)]=list(pc_decomposition_clinvar_no_indel[i])
        quantile_groups_tmp_possible_lofs["pc"+str(i)]=list(pc_decomposition_possible_lofs[i])
        quantile_groups_tmp_indels["pc"+str(i)]=list(pc_decomposition_indels[i])
        quantile_groups_tmp_snv["pc"+str(i)]=list(pc_decomposition_snv[i])
        quantile_groups_tmp_simulated_genes["pc"+str(i)]=list(pc_decomposition_simulated_genes[i])
    cols_to_keep=[col for col in quantile_groups_tmp.columns if 'pc' in str(col)]
    cols_to_keep_clinvar=[col for col in quantile_groups_tmp_clinvar.columns if 'pc' in str(col)]
    cols_to_keep_clinvar_no_indel=[col for col in quantile_groups_tmp_clinvar_no_indel.columns if 'pc' in str(col)]
    cols_to_keep_possible_lofs=[col for col in quantile_groups_tmp_possible_lofs.columns if 'pc' in str(col)]
    cols_to_keep_indels=[col for col in quantile_groups_tmp_indels.columns if 'pc' in str(col)]
    cols_to_keep_snv=[col for col in quantile_groups_tmp_snv.columns if 'pc' in str(col)]
    cols_to_keep_simulated_genes=[col for col in quantile_groups_tmp_simulated_genes.columns if 'pc' in str(col)]
    return [quantile_groups_tmp[cols_to_keep],quantile_groups_tmp_clinvar[cols_to_keep_clinvar],quantile_groups_tmp_clinvar_no_indel[cols_to_keep_clinvar_no_indel],quantile_groups_tmp_possible_lofs[cols_to_keep_possible_lofs],quantile_groups_tmp_indels[cols_to_keep_indels],quantile_groups_tmp_snv[cols_to_keep_snv],quantile_groups_tmp_simulated_genes[cols_to_keep_simulated_genes]]

def cluster_gaussian_mixture(quantile_groups_table,quantile_groups_table_clinvar,quantile_groups_table_clinvar_no_indel,quantile_groups_table_possible_lofs,quantile_groups_table_indels,quantile_groups_table_snv,quantile_groups_table_simulated,variants,output_prefix,n_clust_min,n_clust_max,annotation_groups,single=False,gene_dict=genes,cols_to_output=["transcript","annotation_grouped"],run_sims=0):
    bic=[]
    if single:
        variants["label"]=variants.transcript
    else:
        variants["label"]=variants.transcript+"_"+variants.annotation_grouped   # build label
    for n_clust in range(n_clust_min,n_clust_max+1):
        #for cv_type in ["spherical","tied","diag","full"]:
        for cv_type in ["full"]:
            clusters=GaussianMixture(n_components=n_clust,covariance_type=cv_type).fit(quantile_groups_table)
            clusters.labels_=clusters.predict(quantile_groups_table)
            clusters_proba=clusters.predict_proba(quantile_groups_table)
            # Predict clusters for ClinVar variants
            clinvar_labels=clusters.predict(quantile_groups_table_clinvar)
            quantile_groups_table_clinvar_tmp=quantile_groups_table_clinvar.reset_index()[["transcript"]]
            quantile_groups_table_clinvar_tmp["cluster"]=clinvar_labels
            # add in probabilities of each cluster
            clinvar_proba=clusters.predict_proba(quantile_groups_table_clinvar)
            for i in range(0,n_clust):
                quantile_groups_table_clinvar_tmp["prob"+str(i)]=clinvar_proba[:,i]
            quantile_groups_table_clinvar_tmp.to_csv(output_prefix+"_clusters_clinvar_"+str(n_clust)+".txt")
            # Clinvar without indels
            clinvar_labels_no_indel=clusters.predict(quantile_groups_table_clinvar_no_indel)
            quantile_groups_table_clinvar_tmp_no_indel=quantile_groups_table_clinvar_no_indel.reset_index()[["transcript"]]
            quantile_groups_table_clinvar_tmp_no_indel["cluster"]=clinvar_labels_no_indel
            clinvar_no_indel_proba=clusters.predict_proba(quantile_groups_table_clinvar_no_indel)
            for i in range(0,n_clust):
                quantile_groups_table_clinvar_tmp_no_indel["prob"+str(i)]=clinvar_no_indel_proba[:,i]
            quantile_groups_table_clinvar_tmp_no_indel.to_csv(output_prefix+"_clusters_clinvar_no_indels_"+str(n_clust)+".txt")
            # locations of possible lof variants
            possible_lof_labels=clusters.predict(quantile_groups_table_possible_lofs)
            quantile_groups_table_possible_lofs_tmp=quantile_groups_table_possible_lofs.reset_index()[["transcript"]]
            quantile_groups_table_possible_lofs_tmp["cluster"]=possible_lof_labels
            possible_lof_proba=clusters.predict_proba(quantile_groups_table_possible_lofs)
            for i in range(0,n_clust):
                quantile_groups_table_possible_lofs_tmp["prob"+str(i)]=possible_lof_proba[:,i]
            quantile_groups_table_possible_lofs_tmp.to_csv(output_prefix+"_clusters_possible_lofs_"+str(n_clust)+".txt")
            # UKBB LoF indels
            indels_labels=clusters.predict(quantile_groups_table_indels)
            quantile_groups_table_indels_tmp=quantile_groups_table_indels.reset_index()[["transcript"]]
            quantile_groups_table_indels_tmp["cluster"]=indels_labels
            indels_proba=clusters.predict_proba(quantile_groups_table_indels)
            for i in range(0,n_clust):
                quantile_groups_table_indels_tmp["prob"+str(i)]=indels_proba[:,i]
            quantile_groups_table_indels_tmp.to_csv(output_prefix+"_clusters_indels_"+str(n_clust)+".txt")
            # ukbb lof indels and lof snvs combined
            snv_labels=clusters.predict(quantile_groups_table_snv)
            quantile_groups_table_snv_tmp=quantile_groups_table_snv.reset_index()[["transcript"]]
            quantile_groups_table_snv_tmp["cluster"]=snv_labels
            snv_proba=clusters.predict_proba(quantile_groups_table_snv)
            for i in range(0,n_clust):
                quantile_groups_table_snv_tmp["prob"+str(i)]=snv_proba[:,i]
            quantile_groups_table_snv_tmp.to_csv(output_prefix+"_clusters_lof_snv_"+str(n_clust)+".txt")
            # simulated genes
            if run_sims:
                simulated_labels=clusters.predict(quantile_groups_table_simulated)
                quantile_groups_table_simulated_tmp=quantile_groups_table_simulated.reset_index()[["transcript"]]
                quantile_groups_table_simulated_tmp["cluster"]=simulated_labels
                simulated_proba=clusters.predict_proba(quantile_groups_table_simulated)
                for i in range(0,n_clust):
                    quantile_groups_table_simulated_tmp["prob"+str(i)]=simulated_proba[:,i]
                quantile_groups_table_simulated_tmp.to_csv(output_prefix+"_clusters_simulated_genes_"+str(n_clust)+".txt")
            #clusters.labels_=clusters.predict_proba(quantile_groups_table)
            bic.append(clusters.bic(quantile_groups_table))
            plot_cluster(clusters,variants,n_clust,quantile_groups_table,output_prefix+"_gaussian_mixture_"+cv_type,single,gene_dict,cols_to_output,annotation_groups,clusters_proba,False)
            create_3d_scatter(output_prefix+"_gaussian_mixture_"+cv_type+"_n_"+str(n_clust),quantile_groups_table,output_prefix+"_gaussian_mixture_"+cv_type+"_n_"+str(n_clust)+"_cluster_names.txt",single)
    return bic

def plot_cluster(clusters,variants,n_clust,quantile_groups_table,output_prefix,single,gene_dict,cols_to_output,annotation_groups,clusters_proba,kmeans=False):
    # initialise array to hold cluster labels
    colours=["#1F77B4","#FF7F0E","#2CA02C","#D62728"]
    if len(annotation_groups)<4:
        colours=colours[0:len(annotation_groups)]
    cluster_genes=pd.DataFrame()
    quantile_groups_table.columns=[col for col in quantile_groups_table.columns]
    quantile_groups_table=quantile_groups_table.reset_index()
    quantile_groups_table["label"]=clusters.labels_
    for i in range(0,n_clust):
        quantile_groups_table["prob"+str(i)]=clusters_proba[:,i]
    quantile_groups_table.to_csv(output_prefix+"_n_"+str(n_clust)+"_cluster_pcs_labels.txt")
    # plot clusters
    for i in range(n_clust):
        cluster_0=quantile_groups_table[cols_to_output][clusters.labels_==i]   # get IDs of clusters
        if single:
            labels=list(cluster_0.transcript)
        else:
            cluster_0["label"]=cluster_0.transcript+"_"+cluster_0.annotation_grouped
            labels=list(cluster_0.label) # concatenate labels
        variants_cluster=variants[variants.label.isin(labels)]# filter out from variants
        if len(variants_cluster)>0:
            annotations_cluster=list(set(variants_cluster.annotation_grouped))
            annots=np.array(annotation_groups)[np.isin(annotation_groups,annotations_cluster)].tolist()
            # check that there are more than 2 variants in each annotation group, else drop them for plotting
            annots_tmp=[]
            for anno in annots:
                if len(variants_cluster[variants_cluster.annotation_grouped.isin([anno])])>2:
                    annots_tmp.append(anno)
            annots=annots_tmp
            if len(annots)>0:
                cols=np.array(colours)[np.isin(annotation_groups,annotations_cluster)].tolist()
                for_fig=[variants_cluster.position[variants_cluster.annotation_grouped.isin([x])] for x in annots]
                fig_size=len(cluster_0)
                fig=ff.create_distplot(for_fig,annots,show_hist=False,show_rug=False,colors=cols)#,template="simple_white")    # plot
                fig.write_image(output_prefix+"_n_"+str(n_clust)+"_cluster_"+str(i)+"_"+str(fig_size)+".png")
        # save list of genes in each cluster to output
        if kmeans:
            # calculate distance to cluster centre
            # first get cluster centre coords for current cluster
            centres=clusters.cluster_centers_[i]
            quantile_groups_table["dist"]=0
            for j in range(len(centres)):
                quantile_groups_table.dist=quantile_groups_table.dist+(quantile_groups_table.iloc[:,j+len(cols_to_output)]-centres[j])**2
            quantile_groups_table.dist=quantile_groups_table.dist.apply(lambda x: math.sqrt(x))
            cluster_0=quantile_groups_table[cols_to_output+["dist"]][clusters.labels_==i].rename(columns={"annotation_grouped":"annotation"})
        else:
            cluster_0=quantile_groups_table[cols_to_output][clusters.labels_==i].rename(columns={"annotation_grouped":"annotation"})
        cluster_0["gene"]=cluster_0["transcript"].map(gene_dict)
        cluster_0["cluster"]=i
        cluster_genes=cluster_genes.append(cluster_0)
    # write out table of genes in each cluster
    cluster_genes.to_csv(output_prefix+"_n_"+str(n_clust)+"_cluster_names.txt",sep="\t")

def run_clustering_all(variants,output_prefix,annotation_groups,n_break,clinvar,possible_lofs,variants_indels,variants_snv,simulated_genes,run_sims):
    # count number of variants within each transcript for each variant type
    quantile_groups=variants.groupby(by=["transcript","annotation_grouped","quantile"]).size().reset_index().rename(columns={0: "counts"})
    # count total number of variants per transcipt and annotation
    quantile_groups_transcript_annotation=count_variants(variants,["transcript","annotation_grouped"],quantile_groups,["transcript","annotation_grouped"],["quantile"])
    # count clinvar variants in each quantile
    quantile_groups_clinvar=clinvar.groupby(by=["transcript","annotation_grouped","quantile"]).size().reset_index().rename(columns={0: "counts"})
    quantile_groups_clinvar_no_indel=clinvar.loc[clinvar.indel!=1].groupby(by=["transcript","annotation_grouped","quantile"]).size().reset_index().rename(columns={0: "counts"})
    # count total number of variants per transcript in clinvar
    quantile_groups_transcript_annotation_clinvar=count_variants(clinvar,["transcript","annotation_grouped"],quantile_groups_clinvar,["transcript","annotation_grouped"],["quantile"])
    quantile_groups_transcript_annotation_clinvar_no_indel=count_variants(clinvar.loc[clinvar.indel!=1],["transcript","annotation_grouped"],quantile_groups_clinvar_no_indel,["transcript","annotation_grouped"],["quantile"])
    # count possible lof variants in each quantile
    quantile_groups_possible_lof=possible_lofs.groupby(by=["transcript","annotation_grouped","quantile"]).size().reset_index().rename(columns={0: "counts"})
    quantile_groups_transcript_annotation_possible_lofs=count_variants(possible_lofs,["transcript","annotation_grouped"],quantile_groups_possible_lof,["transcript","annotation_grouped"],["quantile"])
    # count indels
    quantile_groups_indels=variants_indels.groupby(by=["transcript","annotation_grouped","quantile"]).size().reset_index().rename(columns={0: "counts"})
    quantile_groups_transcript_annotation_indels=count_variants(variants_indels,["transcript","annotation_grouped"],quantile_groups_indels,["transcript","annotation_grouped"],["quantile"])
    # count LoF SNVs
    quantile_groups_snvs=variants_snv.groupby(by=["transcript","annotation_grouped","quantile"]).size().reset_index().rename(columns={0: "counts"})
    quantile_groups_transcript_annotation_snvs=count_variants(variants_snv,["transcript","annotation_grouped"],quantile_groups_snvs,["transcript","annotation_grouped"],["quantile"])
    # simulated genes
    if run_sims:
        simulated_genes=simulated_genes.rename(columns={"gene":"transcript","quantile1":"1","quantile2":"2","quantile3":"3","quantile4":"4","quantile5":"5"})
        simulated_genes=simulated_genes.set_index("transcript")
        simulated_genes["n_var"]=simulated_genes.sum(axis=1)    # count number of variants
        for i in range(0,n_break):
            simulated_genes.iloc[:,i]=simulated_genes.iloc[:,i]/simulated_genes.n_var # generate proportions
        simulated_genes=simulated_genes.drop(["n_var"],axis=1)
    else:
        simulated_genes=quantile_groups_transcript_annotation   # if we don't want to run simulations, copy main table to pass to PCs so it gets some data
    ################## PCA
    quantile_groups_pcs_arr=pca_decomposition(quantile_groups_transcript_annotation,quantile_groups_transcript_annotation_clinvar,quantile_groups_transcript_annotation_clinvar_no_indel,quantile_groups_transcript_annotation_possible_lofs,quantile_groups_transcript_annotation_indels,quantile_groups_transcript_annotation_snvs,simulated_genes,n_comp=n_break)
    quantile_groups_pcs=pd.DataFrame(quantile_groups_pcs_arr[0])
    quantile_groups_pcs_clinvar=pd.DataFrame(quantile_groups_pcs_arr[1])
    quantile_groups_pcs_clinvar_no_indel=pd.DataFrame(quantile_groups_pcs_arr[2])
    quantile_groups_pcs_possible_lofs=pd.DataFrame(quantile_groups_pcs_arr[3])
    quantile_groups_pcs_indels=pd.DataFrame(quantile_groups_pcs_arr[4])
    quantile_groups_pcs_snvs=pd.DataFrame(quantile_groups_pcs_arr[5])
    quantile_groups_pcs_simulated=pd.DataFrame(quantile_groups_pcs_arr[6])
    bic_clust=[str(n_break)+" raw"]
    bic_clust=bic_clust+cluster_gaussian_mixture(quantile_groups_pcs,quantile_groups_pcs_clinvar,quantile_groups_pcs_clinvar_no_indel,quantile_groups_pcs_possible_lofs,quantile_groups_pcs_indels,quantile_groups_pcs_snvs,quantile_groups_pcs_simulated,variants,output_prefix+"_pca",n_clust_min,n_clust_max,annotation_groups,run_sims)
    return([bic_clust])

def create_3d_scatter(output_prefix,quantile_groups_transcript_annotation,clusters_file,single): # create 3d scatter_plot showing 6 clusters coloured
    clusters=pd.read_csv(clusters_file,sep="\t")
    clusters=clusters.drop(clusters.columns[0],axis=1)
    quantile_groups_transcript_annotation.columns=[col for col in quantile_groups_transcript_annotation.columns]
    quantile_groups_transcript_annotation=quantile_groups_transcript_annotation.reset_index()
    if single:
        cols_to_add=[quantile_groups_transcript_annotation.columns[0]]
        for i in range(1,len(quantile_groups_transcript_annotation.columns)):
            cols_to_add.append("pc"+str(i))
        quantile_groups_transcript_annotation.columns=[col for col in cols_to_add]
    # then merge in clusters
    for_scatter=pd.merge(quantile_groups_transcript_annotation,clusters)
    for_scatter["cluster"]=for_scatter["cluster"].astype(str)
    i=1
    j=2
    k=3
    fig=px.scatter_3d(for_scatter,x="pc"+str(i),y="pc"+str(j),z="pc"+str(k),color="cluster")
    fig.write_html(output_prefix+"_3D_scatter_"+str(i)+"_"+str(j)+"_"+str(k)+".html")

np.random.seed(42)
# run for MANE select transcripts
if run_mane:
    bic_df=pd.DataFrame()
    for n_break in range(n_break_min,n_break_max):
        tmp=read_variants(variants_file,n_break)	# Read variant data
        variants=tmp[0]
        variants_indels=tmp[1]
        variants=filter_mane(variants,mane_transcripts,["Ensembl_nuc","MANE_status"])   # filter non-MANE transcripts
        variants_indels=filter_mane(variants_indels,mane_transcripts,["Ensembl_nuc","MANE_status"])   # filter non-MANE transcripts
        if indel_thresh:
            variants_indels=filter_indel_length(variants_indels)
        if group_lof:
            variants.loc[variants["annotation_grouped"]=="lof_snv",["annotation_grouped"]]="lof"
            variants_snv=filter_lof_indel_length(variants)
        clinvar=read_clinvar_variants(clinvar_file,n_break) # read in the clinvar variants positions
        possible_lofs=read_possible_lof_variants(possible_lof_file,n_break) # read in the clinvar variants positions
        annotation_groups=list(set(variants.annotation_grouped))
        if run_sims:
            simulated_genes=read_simulated_genes(simulated_genes_file_prefix+str(n_break),n_break)
        else:
            simulated_genes=variants    # if we don't want to run simulations, just copy variants so we have some data to pass to clustering
        bic_arr=run_clustering_all(variants,output_prefix+"_MANE_transcript_n_break_"+str(n_break),annotation_groups,n_break,clinvar,possible_lofs,variants_indels,variants_snv,simulated_genes,run_sims)
        bic_df=bic_df.append(pd.DataFrame(bic_arr[0]).T)
        #bic_df=bic_df.append(pd.DataFrame(bic_arr[1]).T)
        #bic_df=bic_df.append(pd.DataFrame(bic_arr[2]).T)
    bic_df.to_csv(output_prefix+"_MANE_transcript_bic_summary.txt",header=None,index=None,sep="\t")

# single pseudotranscript
if run_pseudo:
    bic_df=pd.DataFrame()
    for n_break in range(n_break_min,n_break_max):
        variants_pseudotranscript=read_variants(variants_file_pseudotranscript,n_break)	# Read variant data
        variants_pseudotranscript=variants_pseudotranscript.rename(columns={"gene": "transcript"})
        if group_lof:
            variants_pseudotranscript.loc[variants_pseudotranscript["annotation_grouped"]=="lof_snv",["annotation_grouped"]]="lof"
        clinvar=read_clinvar_variants(clinvar_file_pseudotranscript,n_break) # read in the clinvar variants positions
        clinvar["transcript"]=clinvar["gene"]
        annotation_groups=list(set(variants_pseudotranscript.annotation_grouped))
        bic_arr=run_clustering_all(variants_pseudotranscript,output_prefix+"_pseudotranscript_n_break_"+str(n_break),annotation_groups,n_break,clinvar)
        #bic_df=bic_df.append(pd.DataFrame(bic_arr[0]).T)
        #bic_df=bic_df.append(pd.DataFrame(bic_arr[1]).T)
        #bic_df=bic_df.append(pd.DataFrame(bic_arr[2]).T)
    #bic_df.to_csv(output_prefix+"_pseudotranscript_bic_summary.txt",header=None,index=None,sep="\t")
