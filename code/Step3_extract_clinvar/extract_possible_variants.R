# R script to extract possible variant positions from Adam's database

# database location
# /slade/home/rnb203/Projects/Rare_Variant_Projects/Variant_Clustering/data/adams_possible_variants/ensembl_exons_codon_3.050820.tbl.gz

# want fields 21,22,23 (these correspond to the actual codon rather than the possible proximal codons)
# need to loop through each and create a set of conditionals (switch-case?) to determine the impact
# possible codons
all_codons <- c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT",
                "CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT",
                "GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
                "TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT")
# LoF consequences we're interested in
# Start Loss: change from ATG (TODO are there any others?)
start_codon <- c("ATG")
start_codon_loss <- all_codons[!(all_codons %in% start_codon)]
# Stop Gain: Change to TAG, TAA or TGA
stop_codon <- c("TAG","TAA","TGA")
stop_codon_loss <- all_codons[!(all_codons %in% stop_codon)]
# Splice Donor:     TODO add these in
splice_donor <- c()
splice_donor_loss <- all_codons[!(all_codons %in% splice_donor)]
# Splice Acceptor: 
splice_acceptor <- c()
splice_acceptor_loss <- all_codons[!(all_codons %in% splice_acceptor)]

library(tidyverse)
library(magrittr)
library(rjson)
library(GenomicRanges)
library(vroom)

args<-fromJSON(file="file_locations")
exon_file<-args$EXON_POS
output_prefix<-args$VARIANTS_FILE_PREFIX

data<-read_delim("/slade/home/rnb203/Projects/Rare_Variant_Projects/Variant_Clustering/data/adams_possible_variants/ensembl_exons_codon_3.050820.tbl.gz",col_names=F,delim="\t")
# determine if any of the variants are LoF
# function to test if a line has any mutations of a given consequence
compare_codons <- function(data,ref_type,consequence_type){
    data %<>% mutate(lof_var1=ifelse(X12 %in% ref_type & X21 %in% consequence_type,1,lof_var1)) %>%
        mutate(lof_var1=ifelse(X12 %in% ref_type & X22 %in% consequence_type,1,lof_var1)) %>%
        mutate(lof_var1=ifelse(X12 %in% ref_type & X23 %in% consequence_type,1,lof_var1)) %>%
        mutate(lof_var2=ifelse(X12 %in% ref_type & X24 %in% consequence_type,1,lof_var2)) %>%
        mutate(lof_var2=ifelse(X12 %in% ref_type & X25 %in% consequence_type,1,lof_var2)) %>%
        mutate(lof_var2=ifelse(X12 %in% ref_type & X26 %in% consequence_type,1,lof_var2)) %>%
        mutate(lof_var3=ifelse(X12 %in% ref_type & X27 %in% consequence_type,1,lof_var3)) %>%
        mutate(lof_var3=ifelse(X12 %in% ref_type & X28 %in% consequence_type,1,lof_var3)) %>%
        mutate(lof_var3=ifelse(X12 %in% ref_type & X29 %in% consequence_type,1,lof_var3))
    data
}
# first split out codons from fields 21-29
data %<>% mutate(X21=str_split(X21,":",simplify=T)[,1]) %>% 
    mutate(X22=str_split(X22,":",simplify=T)[,1]) %>% 
    mutate(X23=str_split(X23,":",simplify=T)[,1]) %>% 
    mutate(X24=str_split(X24,":",simplify=T)[,1]) %>% 
    mutate(X25=str_split(X25,":",simplify=T)[,1]) %>% 
    mutate(X26=str_split(X26,":",simplify=T)[,1]) %>% 
    mutate(X27=str_split(X27,":",simplify=T)[,1]) %>% 
    mutate(X28=str_split(X28,":",simplify=T)[,1]) %>% 
    mutate(X29=str_split(X29,":",simplify=T)[,1]) %>% 
    mutate(X2=str_split(X2,"[.]",simplify=T)[,1]) %>% # split out the transcript number
    mutate(lof_var1=0,lof_var2=0,lof_var3=0) %>% # start with assuming they're not, then change for each type of LoF variant
    compare_codons(start_codon,start_codon_loss) %>% # start loss
    compare_codons(stop_codon_loss,stop_codon) %>% # stop gain
    compare_codons(splice_donor,splice_donor_loss) %>% # splice donor loss  TODO do there need to compare the next/previous base pair?
    compare_codons(splice_donor_loss,splice_donor) %>% # splice donor gain
    compare_codons(splice_acceptor,splice_acceptor_loss) %>% # splice acceptor loss
    compare_codons(splice_acceptor_loss,splice_acceptor) # splice acceptor gain

lof_positions <- data %>% filter(lof_var1==1) %>% 
    select(X2,X7,X8)
lof_positions <- data %>% filter(lof_var2==1) %>% 
    select(X2,X7,X9) %>%
    rbind(lof_positions,.)
lof_positions <- data %>% filter(lof_var3==1) %>% 
    select(X2,X7,X10) %>%
    rbind(lof_positions,.)
names(lof_positions) <- c("transcript","chr","position")

# calculate relative positions TODO
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

read_masks<-function(lof_positions,exon_positions,exon_counts){	# function to read the annotations of variants and append them to a data frame, calculating their relative positions
    # merge lof_variants with exon positions
    variants <- merge(lof_positions,exon_positions) %>% 
        filter(position >= start & position <= end) %>% 
        mutate(position=position-start+relative_start) %>% 
        as_tibble()
	# add in gene name to each variant
	variants <- exon_counts %>% select(transcript,n,transcript_length,strand) %>% 
		dplyr::rename(n_exons=n) %>% merge(variants) %>%
		mutate(position=position/transcript_length) %>% 
		mutate(position=ifelse(strand==-1,1-position,position)) %>%	# calculate relative position as proportion of coding sequence
		select(gene,transcript,position,transcript_length,n_exons,strand) %>% 	# get columns we want and combine with variant df
        as_tibble()
	# return the final data frame
	variants
}

exon_positions<-read_exons(exon_file)
tmp<-calculate_exon_boundaries(exon_positions)
exon_positions<-tmp[[1]]
exon_counts<-tmp[[2]]
read_masks(lof_positions,exon_positions,exon_counts)
write.table(variants,file=paste0(output_prefix,"_possible_lof_positions.txt"),quote=F,row.names=F)