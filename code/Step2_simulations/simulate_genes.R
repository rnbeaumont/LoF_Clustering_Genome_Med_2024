#!/usr/bin/env Rscript
# R script for performing simulations to determine the number of genes we'd expect in each cluster
library(tidyverse)
library(magrittr)
library(rjson)

# set number of genes to generate and number of quantiles
n_gene <- 50000
n_quantile <- 5

# read in file locations
args<-fromJSON(file="file_locations")
output_prefix<-args$VARIANTS_FILE_PREFIX

# set random seed
set.seed(42)

# generate a set of randomly generated genes which I can then include in the cluster step
for(n_quantile in seq(4,10)){
    for(n_var in c(seq(5,20),30,50,100)){
        # sample from a uniform distribution between 0 and 1 n_var times
        genes_df <- seq(1,n_gene)   # gene names
        for(i in seq(1,n_var)){
            genes_df %<>% cbind(runif(n_gene))
        }
        genes_df %<>% as_tibble(.name_repair="universal")   # convert to a df
        names(genes_df) <- paste0("v",seq(0,n_var))
        # calculate number of variants in each quantile for each gene
        genes_df %<>% pivot_longer(cols=paste0("v",seq(1,n_var))) %>% 
            mutate(quantile=cut(value,breaks=seq(0,1,length.out=n_quantile+1),labels=seq(1,n_quantile),include.lowest=T)) %>% 
            group_by(v0) %>% 
            count(quantile) %>% 
            pivot_wider(names_from=quantile,values_from=n,names_prefix="quantile") %>% 
            rename(gene=v0)
        genes_df[is.na(genes_df)]=0 # replace NA with 0
        # write out as table
        write.table(genes_df,file=paste0(output_prefix,"_simulated_genes_n_quantile_",n_quantile,"_n_var_",n_var),row.names=F,quote=F)
    }
}

# generate based on UKBB profile
variants_file <- "/slade/home/rnb203/Projects/Rare_Variant_Projects/Variant_Clustering/results/test_output_variant_relative_positions_all_variants"
mane_transcripts <- "/slade/home/rnb203/Datasets/MANE_transcripts/MANE.GRCh38.v0.95.summary.txt.gz"
# first get the profile of the 10000 genes we include in the analysis
super_groups<-c(
    "LoF_HC"="lof",
    "LoF_HC_50bp_fail"="lof",
    "LoF_HC_50bp_pass"="lof",
    "LoF_LC"="lof",
    "LoF_LC_50bp_fail"="lof",
    "LoF_LC_50bp_pass"="lof",
    "LoF_SNV"="lof_snv",
    "LoF_SNV_50bp_fail"="lof_snv",
    "LoF_SNV_50bp_pass"="lof_snv",
    "missense_cadd24less_variant"="missense",
    "missense_cadd25plus_variant"="missense",
    "missense_cadd30plus_variant"="missense",
    "inframe_indel_cadd24less_variant"="indel",
    "inframe_indel_cadd25plus_variant"="indel",
    "inframe_indel_cadd30plus_variant"="indel",
    "synonymous"="synonymous"
)
variants <- read_delim(variants_file,delim=" ") %>% 
    mutate(annotation_grouped=super_groups[annotation]) %>%
    filter(!(annotation_grouped %in% c("indel")))
ensembl<-read_tsv(mane_transcripts,col_select=all_of(c("Ensembl_nuc","MANE_status"))) %>%
    rename(transcript=c("Ensembl_nuc","MANE_status")[1],MANE_status=c("Ensembl_nuc","MANE_status")[2]) %>%
    mutate(transcript=str_split_fixed(transcript,"[.]",2)[,1]) %>% 
    filter(MANE_status %in% c("MANE Select"))
variants %<>% filter(transcript %in% ensembl$transcript)
variants %<>% mutate(annotation_grouped=ifelse(annotation_grouped=="lof_snv","lof",annotation_grouped))
var_counts<-variants %>% 
    count(gene,transcript,annotation_grouped) %>% 
    pivot_wider(id_cols=c(gene,transcript),names_from=annotation_grouped,values_from=n,names_prefix="count_")
lengths<-variants %>% select(gene,transcript,transcript_length) %>% unique()
data<-merge(var_counts,lengths) %>% as_tibble %>% 
    filter(transcript_length>1000) %>%
    filter(count_synonymous >= 5) %>% 
    filter(count_missense >= 5) %>% 
    filter(count_lof >= 5)

n_sim <- 10000
for(n_quantile in seq(4,10)){
    # Do I need to filter to only those where Synonymous and Missense are in c 0/1?
    clusters<-read_tsv(paste0("/slade/home/rnb203/Projects/Rare_Variant_Projects/Variant_Clustering/plots/test_bayesian/test_output_clinvar_possible_lof_MANE_transcript_n_break_",n_quantile,"_pca_gaussian_mixture_full_n_7_cluster_names.txt")) %>% 
        pivot_wider(id_cols=c(gene,transcript),names_from=annotation,values_from=cluster) %>% 
        merge(data) %>% 
        filter((synonymous==0 | synonymous==1) & (missense==0 | missense==1)) %>%
        as_tibble
    for(type in c("synonymous","missense","lof","uniform")){
        # sample from synonymous (and missense?) variants in these genes rather than uniform?
        eval(parse(text=paste0("positions <- variants %>% filter(gene %in% clusters$gene & annotation_grouped %in% c(\"",type,"\")) %>% ",
            "select(position)")))
        gene_num=1
        print(type)
        for(n_var in clusters$count_lof){
            # sample from a uniform distribution between 0 and 1 n_var times
            genes_df <- seq(1,n_sim)
            genes_df <- paste0("sim_",genes_df,"_gene_",gene_num,"_",type)   # gene names
            for(i in seq(1,n_var)){
                if(type=="uniform"){
                    genes_df %<>% cbind(runif(n_sim))
                }else{
                    genes_df %<>% cbind(sample(positions$position,size=n_sim,replace=T))
                }
            }
            genes_df %<>% as_tibble(.name_repair="universal")   # convert to a df
            names(genes_df) <- paste0("v",seq(0,n_var))
            # calculate number of variants in each quantile for each gene
            genes_df %<>% pivot_longer(cols=paste0("v",seq(1,n_var))) %>% 
                mutate(value=as.numeric(value)) %>% 
                mutate(quantile=cut(value,breaks=seq(0,1,length.out=n_quantile+1),labels=seq(1,n_quantile),include.lowest=T)) %>% 
                group_by(v0) %>% 
                count(quantile) %>% 
                pivot_wider(names_from=quantile,values_from=n,names_prefix="quantile") %>% 
                rename(gene=v0)
            genes_df[is.na(genes_df)]=0 # replace NA with 0
            # write out as table
            if(gene_num==1 && type=="synonymous"){
                write.table(genes_df,file=paste0(output_prefix,"_simulated_genes_n_quantile_",n_quantile,"_ukb_profile"),row.names=F,quote=F)
            }else{
                write.table(genes_df,file=paste0(output_prefix,"_simulated_genes_n_quantile_",n_quantile,"_ukb_profile"),row.names=F,quote=F,append=T,col.names=F)
            }
            gene_num=gene_num+1
        }
    }
}