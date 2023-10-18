suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

####MG>=1####
output_dir <- "/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/qc_OG_taxa_count"
if (!dir.exists(output_dir)) {dir.create(output_dir)}
setwd(output_dir)
count_sum <- read.delim("~/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGanyr1_noSalBongori_taxa_stratified.tsv",header = T,sep = "\t",stringsAsFactors = F)
dna_taxa <- fread("~/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/ruvg_edgeR_cfu_covariates/dna_ruvg_edgeR_da_results.tsv",header = T,sep = "\t",stringsAsFactors = F,select = c("defined_taxa"))
rna_taxa <- fread("~/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/ruvg_edgeR_cfu_covariates/rna_ruvg_edgeR_da_results.tsv",header = T,sep = "\t",stringsAsFactors = F,select = c("defined_taxa"))

taxas <- unique(c(dna_taxa$defined_taxa,rna_taxa$defined_taxa))
count_sum_taxa <- dplyr::filter(count_sum, gtdb_taxa %in% taxas)
count_sum_taxa$defined_taxa <- str_split_fixed(count_sum_taxa$gtdb_taxa,"\\|",n=2)[,1]

####Reading stratified taxonomy counts and adding taxonomy information####
taxa_meta_data <- read.delim("~/Documents/microbiome/markerGeneDatabase/gtdb_r89_bac_arch_mags_metaData.tsv",header = T,sep = "\t",stringsAsFactors = F)
taxa_mag_meta_data <- read.delim("~/Documents/microbiome/gtdbtk_analysis/gtdbtk_v1.4.0_results/gtdbtk.bac120.summary.tsv",header = T,sep = "\t",stringsAsFactors = F)
taxa_mag_meta_data <- separate(taxa_mag_meta_data,classification,into = c("kingdom","phylum","class","order","family","genus","species"),sep=";",remove = F,convert = T)
df_mag_taxa <- taxa_mag_meta_data[,c(1:2)]
df_mag_taxa$defined_taxa <- sapply(1:nrow(taxa_mag_meta_data),FUN = function(x){
  if(taxa_mag_meta_data$species[x]!="s__"){
    return(taxa_mag_meta_data$species[x])
  } else if(taxa_mag_meta_data$genus[x]!="g__"){
    return(taxa_mag_meta_data$genus[x])
  } else if(taxa_mag_meta_data$family[x]!="f__"){
    return(taxa_mag_meta_data$family[x])
  } else if(taxa_mag_meta_data$order[x]!="o__"){
    return(taxa_mag_meta_data$order[x])
  } else if(taxa_mag_meta_data$class[x]!="c__"){
    return(taxa_mag_meta_data$class[x])
  } else if(taxa_mag_meta_data$phylum[x]!="p__"){
    return(taxa_mag_meta_data$phylum[x])
  } else {
    return(taxa_mag_meta_data$kingdom[x])
  }
})
df_mag_taxa$mags_taxa <- paste0(df_mag_taxa$user_genome,"|",df_mag_taxa$defined_taxa)
df_taxa <- taxa_meta_data[,c(1,4,14)]
df_taxa$ncbi_taxonomy[which(is.na(df_taxa$ncbi_taxonomy))] <- "none"
df_taxa$defined_taxa <- sapply(1:nrow(taxa_meta_data),FUN = function(x){
  if(taxa_meta_data$species[x]!="s__"){
    return(taxa_meta_data$species[x])
  } else if(taxa_meta_data$genus[x]!="g__"){
    return(taxa_meta_data$genus[x])
  } else if(taxa_meta_data$family[x]!="f__"){
    return(taxa_meta_data$family[x])
  } else if(taxa_meta_data$order[x]!="o__"){
    return(taxa_meta_data$order[x])
  } else if(taxa_meta_data$class[x]!="c__"){
    return(taxa_meta_data$class[x])
  } else if(taxa_meta_data$phylum[x]!="p__"){
    return(taxa_meta_data$phylum[x])
  } else {
    return(taxa_meta_data$kingdom[x])
  }
})

rowNums <- grep("iMGMC|ihMAG",df_taxa$accession)
df_taxa$defined_taxa[rowNums] <- paste0(df_taxa$accession[rowNums],"|",df_taxa$defined_taxa[rowNums])
count_grouped_defined_species <- left_join(count_sum_taxa,df_taxa[,c(2,4)],by=c("gtdb_taxa"="defined_taxa"))

###Filtering Taxa based on read counts...removing rows with 0 counts in all columns####
keep <- rowSums(count_grouped_defined_species[,8:31] != 0) >= 1
count_removing_zeros <- count_grouped_defined_species[keep,]

###marking mags
count_removing_zeros <- count_removing_zeros %>%
  mutate(mags = case_when(startsWith(defined_taxa, "iMGMC") ~ "iMGMC",
                          startsWith(defined_taxa, "ihMAG") ~ "ihMAGs"))
mags_with_full_species <- taxa_mag_meta_data$species[which(taxa_mag_meta_data$species != "s__")]
count_removing_zeros$mags[which(count_removing_zeros$defined_taxa %in% gsub(" ","_",mags_with_full_species))] <- "Pangenomes containing MAGs"
count_removing_zeros$mags[which(is.na(count_removing_zeros$mags))] <- "Pangenomes"
count_removing_zeros <- separate(count_removing_zeros,gtdb_taxonomy,
                                 into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep=";",
                                 remove=T,convert=T)

write.table(count_removing_zeros,file = paste0(output_dir,"/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGanyr1_noSalBongori_selected_taxa_stratified.tsv"),sep = "\t",row.names = F)

column_names <- colnames(count_removing_zeros)[-c(1:7,32:40)]
stra_og <- count_removing_zeros %>% group_by(associated_OG) %>% 
  dplyr::summarise(genes_fasta_headers = paste0(genes_fasta_headers,collapse = ";"),
                   max_annot_lvl_taxa = paste0(unique(max_annot_lvl_taxa),collapse = ";"),
                   COG_category = paste0(unique(COG_category),collapse = ";"),
                   Description = paste0(unique(Description),collapse = ";"),
                   defined_taxa = paste0(unique(gtdb_taxa),collapse = ";"),
                   Kingdom = paste0(unique(Kingdom),collapse = ";"),
                   Phylum = paste0(unique(Phylum),collapse = ";"),
                   Class = paste0(unique(Class),collapse = ";"),
                   Order = paste0(unique(Order),collapse = ";"),
                   Family = paste0(unique(Family),collapse = ";"),
                   Genus = paste0(unique(Genus),collapse = ";"),
                   Species = paste0(unique(Species),collapse = ";"),
                   mags = paste0(unique(mags),collapse = ";"),
                   across(all_of(column_names), ~ sum(.x, na.rm = TRUE)))
write.table(stra_og,file = paste0(output_dir,"/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGanyr1_noSalBongori_selected_OG_stratified.tsv"),sep = "\t",row.names = F)

####MG>=40%####
output_dir <- "/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/qc_OG_taxa_count"
if (!dir.exists(output_dir)) {dir.create(output_dir)}
setwd(output_dir)
count_sum <- fread("~/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGperc40_taxa_stratified.tsv",header = T,sep = "\t",stringsAsFactors = F)
dna_taxa <- fread("~/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/ruvg_edgeR_cfu_covariates/dna_ruvg_edgeR_da_results.tsv",header = T,sep = "\t",stringsAsFactors = F,select = c("defined_taxa"))
rna_taxa <- fread("~/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/ruvg_edgeR_cfu_covariates/rna_ruvg_edgeR_da_results.tsv",header = T,sep = "\t",stringsAsFactors = F,select = c("defined_taxa"))

taxas <- unique(c(dna_taxa$defined_taxa,rna_taxa$defined_taxa))
count_sum_taxa <- dplyr::filter(count_sum, gtdb_taxa %in% taxas)
count_sum_taxa$defined_taxa <- str_split_fixed(count_sum_taxa$gtdb_taxa,"\\|",n=2)[,1]

####Reading stratified taxonomy counts and adding taxonomy information####
taxa_meta_data <- read.delim("~/Documents/microbiome/markerGeneDatabase/gtdb_r89_bac_arch_mags_metaData.tsv",header = T,sep = "\t",stringsAsFactors = F)
taxa_mag_meta_data <- read.delim("~/Documents/microbiome/gtdbtk_analysis/gtdbtk_v1.4.0_results/gtdbtk.bac120.summary.tsv",header = T,sep = "\t",stringsAsFactors = F)
taxa_mag_meta_data <- separate(taxa_mag_meta_data,classification,into = c("kingdom","phylum","class","order","family","genus","species"),sep=";",remove = F,convert = T)
df_mag_taxa <- taxa_mag_meta_data[,c(1:2)]
df_mag_taxa$defined_taxa <- sapply(1:nrow(taxa_mag_meta_data),FUN = function(x){
  if(taxa_mag_meta_data$species[x]!="s__"){
    return(taxa_mag_meta_data$species[x])
  } else if(taxa_mag_meta_data$genus[x]!="g__"){
    return(taxa_mag_meta_data$genus[x])
  } else if(taxa_mag_meta_data$family[x]!="f__"){
    return(taxa_mag_meta_data$family[x])
  } else if(taxa_mag_meta_data$order[x]!="o__"){
    return(taxa_mag_meta_data$order[x])
  } else if(taxa_mag_meta_data$class[x]!="c__"){
    return(taxa_mag_meta_data$class[x])
  } else if(taxa_mag_meta_data$phylum[x]!="p__"){
    return(taxa_mag_meta_data$phylum[x])
  } else {
    return(taxa_mag_meta_data$kingdom[x])
  }
})
df_mag_taxa$mags_taxa <- paste0(df_mag_taxa$user_genome,"|",df_mag_taxa$defined_taxa)
df_taxa <- taxa_meta_data[,c(1,4,14)]
df_taxa$ncbi_taxonomy[which(is.na(df_taxa$ncbi_taxonomy))] <- "none"
df_taxa$defined_taxa <- sapply(1:nrow(taxa_meta_data),FUN = function(x){
  if(taxa_meta_data$species[x]!="s__"){
    return(taxa_meta_data$species[x])
  } else if(taxa_meta_data$genus[x]!="g__"){
    return(taxa_meta_data$genus[x])
  } else if(taxa_meta_data$family[x]!="f__"){
    return(taxa_meta_data$family[x])
  } else if(taxa_meta_data$order[x]!="o__"){
    return(taxa_meta_data$order[x])
  } else if(taxa_meta_data$class[x]!="c__"){
    return(taxa_meta_data$class[x])
  } else if(taxa_meta_data$phylum[x]!="p__"){
    return(taxa_meta_data$phylum[x])
  } else {
    return(taxa_meta_data$kingdom[x])
  }
})

rowNums <- grep("iMGMC|ihMAG",df_taxa$accession)
df_taxa$defined_taxa[rowNums] <- paste0(df_taxa$accession[rowNums],"|",df_taxa$defined_taxa[rowNums])
count_grouped_defined_species <- left_join(count_sum_taxa,df_taxa[,c(2,4)],by=c("gtdb_taxa"="defined_taxa"))

###Filtering Taxa based on read counts...removing rows with 0 counts in all columns####
keep <- rowSums(count_grouped_defined_species[,8:31] != 0) >= 1
count_removing_zeros <- count_grouped_defined_species[keep,]

###marking mags
count_removing_zeros <- count_removing_zeros %>%
  mutate(mags = case_when(startsWith(defined_taxa, "iMGMC") ~ "iMGMC",
                          startsWith(defined_taxa, "ihMAG") ~ "ihMAGs"))
mags_with_full_species <- taxa_mag_meta_data$species[which(taxa_mag_meta_data$species != "s__")]
count_removing_zeros$mags[which(count_removing_zeros$defined_taxa %in% gsub(" ","_",mags_with_full_species))] <- "Pangenomes containing MAGs"
count_removing_zeros$mags[which(is.na(count_removing_zeros$mags))] <- "Pangenomes"
count_removing_zeros <- separate(count_removing_zeros,gtdb_taxonomy,
                                 into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep=";",
                                 remove=T,convert=T)

write.table(count_removing_zeros,file = paste0(output_dir,"/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGperc40_selected_taxa_stratified.tsv"),sep = "\t",row.names = F)

column_names <- colnames(count_removing_zeros)[-c(1:7,32:40)]
stra_og <- count_removing_zeros %>% group_by(associated_OG) %>% 
  dplyr::summarise(genes_fasta_headers = paste0(genes_fasta_headers,collapse = ";"),
                   max_annot_lvl_taxa = paste0(unique(max_annot_lvl_taxa),collapse = ";"),
                   COG_category = paste0(unique(COG_category),collapse = ";"),
                   Description = paste0(unique(Description),collapse = ";"),
                   defined_taxa = paste0(unique(gtdb_taxa),collapse = ";"),
                   Kingdom = paste0(unique(Kingdom),collapse = ";"),
                   Phylum = paste0(unique(Phylum),collapse = ";"),
                   Class = paste0(unique(Class),collapse = ";"),
                   Order = paste0(unique(Order),collapse = ";"),
                   Family = paste0(unique(Family),collapse = ";"),
                   Genus = paste0(unique(Genus),collapse = ";"),
                   Species = paste0(unique(Species),collapse = ";"),
                   mags = paste0(unique(mags),collapse = ";"),
                   across(all_of(column_names), ~ sum(.x, na.rm = TRUE)))
write.table(stra_og,file = paste0(output_dir,"/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGperc40_selected_OG_stratified.tsv"),sep = "\t",row.names = F)
