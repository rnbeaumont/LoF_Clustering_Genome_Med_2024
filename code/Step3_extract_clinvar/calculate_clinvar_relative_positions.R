#!/usr/bin/env Rscript
# R script for reading exon and clinvar variant positions and generating regional plots
library(tidyverse)
library(magrittr)
library(vroom)
library(glue)
library(ggplot2)
theme_set(theme_bw())
library(Cairo)
library(rjson)
library(GenomicRanges)

args<-fromJSON(file="file_locations")
exon_file<-args$EXON_POS
clinvar_pathogenic<-args$CLINVAR_FORMATTED
freq_file<-args$FREQ_FILE
output_prefix<-args$VARIANTS_FILE_PREFIX

read_exons<-function(exon_file){	# function to read the positions of the exons
	chrs_to_get<-c(seq(1,23),"X")
	data<-vroom(exon_file,show_col_types=F) %>% # read in the file
		select(`Gene name`,`Transcript stable ID`,`Exon stable ID`,`Chromosome/scaffold name`,`Genomic coding start`,`Genomic coding end`,`Strand`) %>% # pull out fields we want
		dplyr::rename(gene=`Gene name`,transcript=`Transcript stable ID`,exonID=`Exon stable ID`,chr=`Chromosome/scaffold name`,start=`Genomic coding start`,end=`Genomic coding end`,strand=Strand) %>% # rename columns
		filter(chr %in% chrs_to_get & !is.na(start) & !is.na(end) & !is.na(gene) & !is.na(transcript) & !is.na(strand)) %>% # filter only to complete rows on autosomes (no random patches)
		arrange(gene,transcript,start)
	data
}

calculate_exon_boundaries<-function(data){	# function to calculate exon boundaries relative to coding sequence
	exon_counts<-data %>%	# count the number of exons per transcript
		count(transcript)
	data %<>% mutate(relative_start=-1,relative_end=-1,exon_length=end-start) %>%	# calculate exon length
		group_by(transcript) %>% mutate(relative_end=cumsum(exon_length)) %>%	# calculate relative end position
		mutate(relative_start=relative_end-exon_length+1)	# then back-calculate relative start position
	# then pull out the total length of the transcript and merge into exon_counts
	exon_counts<-data %>% group_by(transcript) %>% 
		slice_max(relative_end) %>% select(transcript,relative_end,strand) %>%
		dplyr::rename(transcript_length=relative_end) %>% base::merge(exon_counts,.)
	list(data, exon_counts)
}

calculate_exon_boundaries_single_pseudogene<-function(data){	# function to calculate exon boundaries relative to coding sequence
	data %<>% arrange(chr,gene,start,end) %>%	# sort exons to be in order of genomic position
		ungroup() %>%	# remove grouping
		select(-c(transcript,relative_start,relative_end,exonID)) %>%	# drop transcript as we don't need this
		distinct() %>%	# remove duplicate rows
		mutate(relative_start=-1,relative_end=-1,exon_length=end-start)	# calculat[e exon length
	# Now we need to merge exons which overlap to form single (larger) exons. This may impact n_exons making it inaccurate!!! Work out how to count best
	prev_len<-nrow(data)
	new_len<-nrow(data)-1
	i=0 # check how many times it takes to stop removing more - just for info. TODO remove this
	while(prev_len>new_len){ # iteratively remove exons until we don't remove any further exons
		prev_len<-nrow(data)	# save current length
		data %<>% data.table	# convert to data.table to compare with prev and next rows
		data[,prev_start:=c(NA,start[-.N]),by=gene]	# generate temporary column of previous exons's start position for gene (first one will be NA)
		data[,prev_end:=c(NA,end[-.N]),by=gene]	# generate temporary column of previous exons's end position for gene (first one will be NA)
		data[,prev_exon_length:=c(NA,exon_length[-.N]),by=gene] # save previous exon length to check if start falls outside previous exon
		data %<>% arrange(chr,gene,start,-end) %>% data.table # reorder to get next exon's start and end position
		data[,next_start:=c(NA,start[-.N]),by=gene]	# generate temporary column of next exons's start position for gene
		data[,next_end:=c(NA,end[-.N]),by=gene]	# generate temporary column of next exons's end position for gene
		data %<>% as_tibble %>% # convert back to tibble
			mutate(diff_start=start-prev_start,diff_end=end-prev_end,diff_start_next=next_start-start) %>% # calculate difference between previous start and end positions
			mutate(start=ifelse(is.na(prev_start),start,ifelse(diff_end<=0 | diff_start<prev_exon_length,prev_start,start)),
				diff_start=start-prev_start,
				diff_start_next=next_start-start) %>% # check for overlapping exons, then update start position if necessary
			mutate(end=ifelse(!is.na(next_end) & !is.na(prev_end),	# need to check whether they overlap as well as checking if starts are equal
					ifelse(diff_start==0 & diff_start_next==0,
						pmax(end,prev_end,next_end),
						ifelse(diff_start==0,
							pmax(end,prev_end),
							ifelse(diff_start_next==0,
								pmax(end,next_end),
								end))),
					ifelse(is.na(next_end),
						ifelse(!is.na(prev_end),
							ifelse(diff_start==0,
								pmax(end,prev_end),
								end),
							end),
						ifelse(diff_start_next==0,
							pmax(end,next_end),
							end)))) %>% # check for overlapping exons, then update end position if necessary
			mutate(exon_length=end-start) %>% # update the exon length
			select(-c(prev_start,prev_end,next_end,diff_start,diff_end,prev_exon_length,next_start,diff_start_next)) %>% # and drop the previous start and end columns
			distinct %>% arrange(chr,gene,start,end)	# filter out now duplicate rows. TODO hopefully this will fix the issue, but maybe should iterate until the number of rows is the same before and after
		new_len<-nrow(data) # save new length to check if we're still removing any
		i=i+1
	}
	print(paste("number of iterations required",i))

	exon_counts<-data %>% count(gene)	# count the number of exons per gene
	data %<>% group_by(gene) %>% mutate(relative_end=cumsum(exon_length)) %>%	# calculate relative end position
		mutate(relative_start=relative_end-exon_length+1)	# then back-calculate relative start position
	# then pull out the total length of the transcript and merge into exon_counts
	exon_counts<-data %>% group_by(gene) %>% 
		dplyr::slice(which.max(relative_end)) %>% select(gene,relative_end,strand) %>%
		dplyr::rename(gene_length=relative_end) %>% base::merge(exon_counts,.)
	list(data, exon_counts)
}

read_masks<-function(annotations,exon_positions,exon_counts){	# function to read the annotations of variants and append them to a data frame, calculating their relative positions
	# read in Clinvar pathogenic variants
	variants_tmp <- vroom(annotations,show_col_types=F,delim="\t") %>%	# read in the files
		select(variantID) %>% # drop row number column
		mutate(chr=str_split(variantID,":",simplify=T)[,1]) %>% # split out chr
		mutate(pos=str_split(variantID,":",simplify=T)[,2]) %>% # split out position
		mutate(ref=str_split(variantID,":",simplify=T)[,3]) %>% # split out reference allele
		mutate(alt=str_split(variantID,":",simplify=T)[,4]) # split out other allele
	# add in gene name to each variant
	variants <- data.frame()
	for(i in c(1:22,"X")){
		variants_chr <- variants_tmp %>% filter(chr==i)
		exon_positions_tmp <- exon_positions %>% filter(chr==i) %>% 
			arrange(start)
		variants_chr <- variants_chr %$% GRanges(seqnames=paste0("chr",chr),
				IRanges(start=as.numeric(pos),end=as.numeric(pos)),
				ref=ref,
				alt=alt,
				chr=chr,
				variantID=variantID)
		exon_positions_tmp_grange <- exon_positions_tmp %$% GRanges(seqnames=paste0("chr",chr),
				IRanges(start=start,end=end),
				gene=gene)
		variants_chr <- mergeByOverlaps(variants_chr,exon_positions_tmp_grange) %>%
			as_tibble() %>% 
			mutate(pos=str_split(variantID,":",simplify=T)[,2]) %>% 
			select(variantID,chr,pos,ref,alt,gene) %>% distinct
		variants_chr %<>% #transform(transcript=exon_positions_tmp$transcript[findInterval(pos,exon_positions_tmp$start)]) %>% 
			merge(exon_positions_tmp) %>% 
			mutate(pos=as.numeric(pos)) %>% 
			filter(start<pos & end>pos) %>% # keep exons the variant is within
			mutate(pos=pos-start+relative_start) %>%	# get the position within the gene
			mutate(position=pos) %>% select(-pos)
		variants <- rbind(variants,variants_chr)
	}
	variants <- exon_counts %>% select(transcript,n,transcript_length,strand) %>% 
		dplyr::rename(n_exons=n) %>% merge(variants) %>%
		mutate(position=position/transcript_length) %>% 
		mutate(position=ifelse(strand==-1,1-position,position)) %>%	# calculate relative position as proportion of coding sequence
		select(variantID,gene,transcript,position,transcript_length,n_exons,strand)	# get columns we want and combine with variant df
	# return the final data frame
	variants
}

# first calculate boundaries based on transcripts
exon_positions<-read_exons(exon_file)
tmp<-calculate_exon_boundaries(exon_positions)
exon_positions<-tmp[[1]]
exon_counts<-tmp[[2]]
variants<-read_masks(clinvar_pathogenic,exon_positions,exon_counts)
write.table(variants,file=glue(output_prefix,"_clinvar_variants"),quote=F,row.names=F)

rm(variants)
gc()

# then create a single pseudo-transcript for each gene
tmp<-calculate_exon_boundaries_single_pseudogene(exon_positions)
exon_positions<-tmp[[1]] %>% mutate(transcript=gene)
exon_counts<-tmp[[2]] %>% mutate(transcript=gene) %>% dplyr::rename(transcript_length=gene_length)
variants<-read_masks(clinvar_pathogenic,exon_positions,exon_counts) %>% 
	dplyr::rename(gene_length=transcript_length) %>% 
	select(-transcript)
write.table(variants,file=glue(output_prefix,"_single_pseudo_transcript_clinvar_variants"),quote=F,row.names=F)