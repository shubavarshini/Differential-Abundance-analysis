suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(edgeR))
suppressMessages(library(RUVSeq))
suppressMessages(library(pheatmap))
suppressMessages(library(sva))

setwd("/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/qc_taxa_count")
# setwd("/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/qc_taxa_count")
####Reading stratified taxonomy counts and adding taxonomy information####
count_sum <- read.delim("~/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGanyr1_noSalBongori_species_summed.tsv",header = T,sep = "\t",stringsAsFactors = F)
# count_sum <- read.delim("~/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGperc40_species_summed.tsv",header = T,sep = "\t",stringsAsFactors = F)
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
count_sum <- left_join(count_sum,df_mag_taxa[,c(1,3,4)],by=c("species"="user_genome"))
count_sum <- count_sum %>% mutate(defined_taxa = coalesce(mags_taxa,species))
count_grouped_defined_species <- left_join(count_sum,df_taxa[,2:4],by=c("defined_taxa"))
count_grouped_defined_species <- count_grouped_defined_species %>% group_by(defined_taxa,gtdb_taxonomy,ncbi_taxonomy) %>% 
  summarise_if(is.numeric, sum, na.rm = TRUE)

###Filtering Taxa based on read counts...removing rows with 0 counts in all columns####
keep <- rowSums(count_grouped_defined_species[,4:27] != 0) >= 1
count_removing_zeros <- count_grouped_defined_species[keep,]

###marking mags
count_removing_zeros <- count_removing_zeros %>%
  mutate(mags = case_when(startsWith(defined_taxa, "iMGMC") ~ "iMGMC",
                          startsWith(defined_taxa, "ihMAG") ~ "ihMAGs"))
mags_with_full_species <- taxa_mag_meta_data$species[which(taxa_mag_meta_data$species != "s__")]
count_removing_zeros$mags[which(count_removing_zeros$defined_taxa %in% gsub(" ","_",mags_with_full_species))] <- "Pangenomes with species-MAGs"
count_removing_zeros$mags[which(is.na(count_removing_zeros$mags))] <- "FiRe Pangenomes"

#SummarizedExperiment easy to handle all the information at once
colnames(count_removing_zeros) <- gsub("MOCK_","Uninfected_",colnames(count_removing_zeros))
colnames(count_removing_zeros) <- gsub("WT_","Infected_",colnames(count_removing_zeros))
count_removing_zeros <- count_removing_zeros[,c(1,2,3,4,6,8,10,12,14,16,18,20,22,24,26,5,7,9,11,13,15,17,19,21,23,25,27,28)]
taxa_df <- separate(count_removing_zeros[,c(1,2,3,28)],gtdb_taxonomy,
                    into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep=";",
                    remove=T,convert=T)
taxa_df$Species <- sapply(1:nrow(taxa_df),function(x){
  if(taxa_df$Species[x] == "s__"){
    return(taxa_df$defined_taxa[x])
  } else {
    return(taxa_df$Species[x])
  }
})
taxa <- as.matrix(taxa_df,row.names=count_removing_zeros[,1])


data <- SummarizedExperiment(assays = list(counts=as.matrix(count_removing_zeros[,-c(1,2,3,28)],
                                                            row.names=count_removing_zeros[,c(1)])),
                             rowData = taxa,colData = DataFrame(Samples=colnames(count_removing_zeros)[-c(1,2,3,28)], 
                                                                Condition=gsub("_\\d+","",colnames(count_removing_zeros)[-c(1,2,3,28)]),
                                                                Type=str_split_fixed(colnames(count_removing_zeros)[-c(1,2,3,28)],"_",n=3)[,3],
                                                                Replicate=str_split_fixed(colnames(count_removing_zeros)[-c(1,2,3,28)],"_",n=3)[,2],
                                                                infectionRate=scales::rescale(c(rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000),rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000)))))

y <- edgeR::DGEList(counts=assay(data),group = data$Condition,samples =colData(data),genes = rowData(data))

####Metagenome####
####Filtering MG_cpm_ge64_3samples####
x_dna <- edgeR::DGEList(counts=assay(data[,1:12]),group = data$Condition[1:12],samples =colData(data[,1:12]),
                        genes = rowData(data))
CPM <- cpm(x_dna,lib.size=colSums(x_dna$counts))
if(getwd()=="/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/qc_taxa_count"){
  CPM.Cutoff <- 64
} else {
  CPM.Cutoff <- 8
}

tol <- 1e-14
MinSampleSize <- 3
keep.CPM <- rowSums(CPM >= CPM.Cutoff) >= (MinSampleSize - tol)

x_dna <- x_dna[keep.CPM,,keep.lib.sizes=FALSE]
####RUVg####
design_ruvg <- model.matrix(~0+group, data = data.frame(cbind(samples=x_dna$samples$Samples,group=x_dna$samples$group),
                                                        row.names = x_dna$samples$Samples,stringsAsFactors = F))
design_ruvg[,"group1"] <- scales::rescale(c(rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000)))
colnames(design_ruvg) <- c("wt","mock")
all_contrasts <- makeContrasts(dna=wt-mock,levels=colnames(design_ruvg))

normalize_method <- "TMM"
x_dna_ruvg <- calcNormFactors(x_dna, method=normalize_method)
x_dna_ruvg <- estimateGLMCommonDisp(x_dna_ruvg, design_ruvg)
x_dna_ruvg <- estimateGLMTagwiseDisp(x_dna_ruvg, design_ruvg)

fit <- glmFit(x_dna_ruvg, design_ruvg)
lrt <- glmLRT(fit, contrast=all_contrasts)

top <- topTags(lrt, n=nrow(x_dna))$table
empirical <- x_dna_ruvg$genes$defined_taxa[which(!(x_dna_ruvg$genes$defined_taxa %in% top$defined_taxa[1:250]))]
####empirical RUVg####
a <- floor(x_dna_ruvg$counts)
row.names(a) <- x_dna_ruvg$genes$defined_taxa
set <- newSeqExpressionSet(a,phenoData = data.frame(x_dna_ruvg$samples$group, row.names=colnames(x_dna_ruvg$counts)))
set2 <- RUVg(set, empirical, k=2)
i=2
####Plotting counts####
cairo_pdf(file="./MG_infectionCovariate/logCPM_distribution_pca_rle_afterRUVg.pdf",bg="transparent",
          width = 11,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 5))
a_logCPM <- data.frame(cbind("defined_taxa" =rownames(a),edgeR::cpm.default(a,log=T)), stringsAsFactors = FALSE)
m_a <- reshape2::melt(a_logCPM,value.name = "log_cpm",id.vars=c("defined_taxa"))
m_a$Condition <- gsub("_\\d+","",m_a$variable)
m_a$Tpe=str_split_fixed(m_a$variable,"_",n=3)[,3]
m_a$Replicate=str_split_fixed(m_a$variable,"_",n=3)[,2]
m_a$log_cpm <- as.numeric(m_a$log_cpm)

vioplot::vioplot(m_a$log_cpm ~ m_a$variable, col=c(rep("blue",6),rep("red",6)), cex=0.4,
                 main=paste0("Log CPM violin plot\nN= ",nrow(a)), xlab="Replicates",ylab="Log CPM",xaxt = "n")

normalisedCounts <- data.frame(cbind("defined_taxa" =rownames(normCounts(set2)),edgeR::cpm.default(normCounts(set2),log=T)), stringsAsFactors = FALSE)
m_normalisedCounts <- reshape2::melt(normalisedCounts,value.name = "normalised_count",id.vars=c("defined_taxa"))
m_normalisedCounts$Condition <- gsub("_\\d+","",m_normalisedCounts$variable)
m_normalisedCounts$Tpe=str_split_fixed(m_normalisedCounts$variable,"_",n=3)[,3]
m_normalisedCounts$Replicate=str_split_fixed(m_normalisedCounts$variable,"_",n=3)[,2]
m_normalisedCounts$normalised_count <- as.numeric(m_normalisedCounts$normalised_count)

vioplot::vioplot(m_normalisedCounts$normalised_count ~ m_normalisedCounts$variable, col=c(rep("blue",6),rep("red",6)), cex=0.4,
                 main=paste0("Normalised Count from RUVg\nLog CPM violin plot k=",i,"\nN= ",nrow(normalisedCounts)), xlab="Replicates",ylab="Log CPM",xaxt = "n")

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=c(rep("blue",6),rep("red",6)),main=paste0("RLE plot k=",i),xaxt = "n", xlab="Replicates",ylab="RLE")

plotPCA(set2, col=c(rep("blue",6),rep("red",6)), cex=0.4,main=paste0("PCA plot k=",i))

corrected <- removeBatchEffect(cpm(x_dna, log=TRUE), design=design_ruvg, covariates=as.matrix(pData(set2)[,c("W_1","W_2")]))
corrected <- data.frame(cbind("defined_taxa" =x_dna$genes$defined_taxa,corrected), stringsAsFactors = FALSE)
m_corrected <- reshape2::melt(corrected,value.name = "corrected_logCPM",id.vars=c("defined_taxa"))
m_corrected$Condition <- gsub("_\\d+","",m_corrected$variable)
m_corrected$Tpe=str_split_fixed(m_corrected$variable,"_",n=3)[,3]
m_corrected$Replicate=str_split_fixed(m_corrected$variable,"_",n=3)[,2]
m_corrected$corrected_logCPM <- as.numeric(m_corrected$corrected_logCPM)
vioplot::vioplot(m_corrected$corrected_logCPM ~ m_corrected$variable, col=c(rep("blue",6),rep("red",6)), cex=0.4,
                 main=paste0("Corrected counts after RUVg\nLog CPM violin plot k=",i,"\nN= ",nrow(corrected)), xlab="Replicates",ylab="Log CPM",xaxt = "n")
dev.off()
####Calculating normalising factors####
x_dna <- calcNormFactors(x_dna,method = normalize_method)
####limma no RUVg####
cairo_pdf(file="./MG_infectionCovariate/pvalue_limma_noRUVg.pdf",bg="transparent",
          width = 6,height = 3.5,family = "Times New Roman")
par(mfrow = c(1, 2))
logCPM_afterNormFac_dna <- edgeR::cpm(x_dna, log=TRUE)
row.names(logCPM_afterNormFac_dna) <- x_dna$genes$defined_taxa

fit <- lmFit(logCPM_afterNormFac_dna , design_ruvg)
fit <- contrasts.fit(fit, contrasts=all_contrasts)
fit <- eBayes(fit, trend=TRUE)
dna_top_table_noruvg <- topTable(fit, coef=1,number = nrow(x_dna$counts))
hist(dna_top_table_noruvg$adj.P.Val,breaks=20,main = paste("Histogram of adj.P-value\nLimma (No RUVg) N=" , nrow(x_dna$counts)),
     xlab = "adj.P-value")
hist(dna_top_table_noruvg$logFC,breaks=50,main = paste("Histogram of logFC\nLimma N=" , nrow(x_dna$counts)),
     xlab = "logFC")
dev.off()
####Heatmap
cairo_pdf(file="./MG_infectionCovariate/limma_noRUVg_sigHeatmap.pdf",bg="transparent",
          width = 5,height = 4,family = "Times New Roman")
dge_limma <- subset(dna_top_table_noruvg,dna_top_table_noruvg$adj.P.Val<=0.1 & 
                      (dna_top_table_noruvg$logFC>=2|dna_top_table_noruvg$logFC<=-2))
dge_limma <- left_join(rownames_to_column(dge_limma,var = "defined_taxa"),a_logCPM,by=c("defined_taxa"))
dge_limma <- dge_limma %>% mutate(across(colnames(dge_limma)[-1],as.numeric))
dge_limma <- dge_limma[order(dge_limma$logFC,decreasing = T), , drop = FALSE]
mat_dna <- as.matrix(dge_limma[,8:19])
rownames(mat_dna) <- dge_limma[,1]
print(pheatmap(mat_dna, cluster_rows = F, cluster_cols = T,fontsize_row = 4,
               main = "Bi-clustering using hclust of log CPMs\nSignificant after limma (No RUVg)"))
dev.off()
####logFC and logCPM bar plots
cairo_pdf(file="./MG_infectionCovariate/limma_noRUVg_sig_logFC_logCPM_barPlots.pdf",bg="transparent",
          width = 4,height = 4,family = "Times New Roman")
logFC_barPlot <- ggplot(dge_limma,aes(x=reorder(defined_taxa,-logFC),y=logFC)) + 
  geom_col() + xlab(element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

m_mat_dna <- reshape2::melt(mat_dna,value.name = "log_CPM",id.vars=c("defined_taxa"))
m_mat_dna$Condition <- str_split_fixed(m_mat_dna$Var2,"_",n=3)[,1]
m_mat_dna$Type=str_split_fixed(m_mat_dna$Var2,"_",n=3)[,3]
m_mat_dna$Replicate=str_split_fixed(m_mat_dna$Var2,"_",n=3)[,2]
m_mat_dna$log_CPM <- as.numeric(m_mat_dna$log_CPM)
logCPM_boxPlot <- ggplot(m_mat_dna,aes(x=Var1,y=log_CPM,fill=fct_relevel(Condition,"Uninfected","Infected"))) + 
  geom_boxplot(outlier.size = 0.5) + xlab(element_blank()) + ylab("log CPM") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = "italic"),
        legend.position = "bottom",legend.title = element_blank()) + 
  scale_fill_manual(values = c("blue","red"))

theme_set(theme_bw(base_size = 12,base_family = "Times New Roman") +
            theme(strip.background = element_rect(fill = "transparent",colour = "NA"),
                  panel.background = element_rect(fill='transparent',colour = "NA"), #transparent panel bg
                  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                  # panel.grid.major = element_blank(), #remove major gridlines
                  panel.grid.minor = element_blank(), #remove minor gridlines
                  legend.background = element_rect(fill='transparent',colour = "NA"), #transparent legend bg
                  legend.box.background = element_rect(fill='transparent',colour = "NA"),
                  legend.key = element_rect(fill = "transparent",colour = "NA")))
print(ggpubr::ggarrange(logFC_barPlot,logCPM_boxPlot,nrow = 2,ncol = 1,heights = c(1,5),align = "v"))
dev.off()

####limma-voom no RUVg####
cairo_pdf(file="./MG_infectionCovariate/meanVariance_limmaVoom_noRUVg.pdf",bg="transparent",
          width = 8,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 3))
logCPM_afterNormFac_dna <- edgeR::cpm(x_dna, log=TRUE)
row.names(logCPM_afterNormFac_dna) <- x_dna$genes$defined_taxa
v <- voom(x_dna,design_ruvg,plot = T)
fit <- lmFit(v , design_ruvg)
fit <- contrasts.fit(fit, contrasts=all_contrasts)
fit <- eBayes(fit, trend=TRUE)
dna_top_table_limmaVoom_noruvg <- topTable(fit, coef=1,number = nrow(x_dna$counts))
hist(dna_top_table_limmaVoom_noruvg$adj.P.Val,breaks=20,main = paste("Histogram of adj.P-value\nLimma (No RUVg) N=" , nrow(x_dna$counts)),
     xlab = "adj.P-value")
hist(dna_top_table_limmaVoom_noruvg$logFC,breaks=50,main = paste("Histogram of logFC\nLimma N=" , nrow(x_dna$counts)),
     xlab = "logFC")
dev.off()
####Heatmap
cairo_pdf(file="./MG_infectionCovariate/limmaVoom_noRUVg_sigHeatmap.pdf",bg="transparent",
          width = 5,height = 4,family = "Times New Roman")
dge_limma <- subset(dna_top_table_limmaVoom_noruvg,dna_top_table_limmaVoom_noruvg$adj.P.Val<=0.1 & 
                      (dna_top_table_limmaVoom_noruvg$logFC>=2|dna_top_table_limmaVoom_noruvg$logFC<=-2))
dge_limma <- left_join(dge_limma,a_logCPM,by=c("defined_taxa"))
dge_limma <- dge_limma %>% mutate(across(colnames(dge_limma)[-1],as.numeric))
dge_limma <- dge_limma[order(dge_limma$logFC,decreasing = T), , drop = FALSE]
mat_dna <- as.matrix(dge_limma[,17:28])
rownames(mat_dna) <- dge_limma[,1]
print(pheatmap(mat_dna, cluster_rows = F, cluster_cols = T,fontsize_row = 4,
               main = "Bi-clustering using hclust of log CPMs\nSignificant after limma (No RUVg)"))
dev.off()
####logFC and logCPM bar plots
cairo_pdf(file="./MG_infectionCovariate/limmaVoom_noRUVg_sig_logFC_logCPM_barPlots.pdf",bg="transparent",
          width = 4,height = 4,family = "Times New Roman")
logFC_barPlot <- ggplot(dge_limma,aes(x=reorder(defined_taxa,-logFC),y=logFC)) + 
  geom_col() + xlab(element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

m_mat_dna <- reshape2::melt(mat_dna,value.name = "log_CPM",id.vars=c("defined_taxa"))
m_mat_dna$Condition <- str_split_fixed(m_mat_dna$Var2,"_",n=3)[,1]
m_mat_dna$Type=str_split_fixed(m_mat_dna$Var2,"_",n=3)[,3]
m_mat_dna$Replicate=str_split_fixed(m_mat_dna$Var2,"_",n=3)[,2]
m_mat_dna$log_CPM <- as.numeric(m_mat_dna$log_CPM)
logCPM_boxPlot <- ggplot(m_mat_dna,aes(x=Var1,y=log_CPM,fill=fct_relevel(Condition,"Uninfected","Infected"))) + 
  geom_boxplot(outlier.size = 0.5) + xlab(element_blank()) + ylab("log CPM") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = "italic"),
        legend.position = "bottom",legend.title = element_blank()) + 
  scale_fill_manual(values = c("blue","red"))

theme_set(theme_bw(base_size = 12,base_family = "Times New Roman") +
            theme(strip.background = element_rect(fill = "transparent",colour = "NA"),
                  panel.background = element_rect(fill='transparent',colour = "NA"), #transparent panel bg
                  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                  # panel.grid.major = element_blank(), #remove major gridlines
                  panel.grid.minor = element_blank(), #remove minor gridlines
                  legend.background = element_rect(fill='transparent',colour = "NA"), #transparent legend bg
                  legend.box.background = element_rect(fill='transparent',colour = "NA"),
                  legend.key = element_rect(fill = "transparent",colour = "NA")))
print(ggpubr::ggarrange(logFC_barPlot,logCPM_boxPlot,nrow = 2,ncol = 1,heights = c(1,5),align = "v"))
dev.off()

####voom limma####
cairo_pdf(file="./MG_infectionCovariate/meanVariance_limmaVoom_afterRUVg.pdf",bg="transparent",
          width = 8,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 3))
design_voom <- cbind(design_ruvg,set2$W_1,set2$W_2)
colnames(design_voom) <- c("wt","mock","W_1","W_2")
v <- voom(x_dna,design_voom,plot = T)
fit <- lmFit(v , design_ruvg)
fit <- contrasts.fit(fit, contrasts=all_contrasts)
fit <- eBayes(fit, trend=TRUE)
dna_top_table_limmaVoom <- topTable(fit, coef=1,number = nrow(x_dna$counts))
hist(dna_top_table_limmaVoom$adj.P.Val,breaks=20,main = paste("Histogram of adj.P-value\nLimma-voom N=" , nrow(x_dna$counts)),
     xlab = "adj.P-value")
hist(dna_top_table_limmaVoom$logFC,breaks=50,main = paste("Histogram of logFC\nLimma-voom N=" , nrow(x_dna$counts)),
     xlab = "logFC")
dim(subset(dna_top_table_limmaVoom,dna_top_table_limmaVoom$adj.P.Val<=0.05 & 
             (dna_top_table_limmaVoom$logFC>=2|dna_top_table_limmaVoom$logFC<=-2)))
dev.off()
####limma####
cairo_pdf(file="./MG_infectionCovariate/pvalue_limma_afterRUVg.pdf",bg="transparent",
          width = 6,height = 3.5,family = "Times New Roman")
par(mfrow = c(1, 2))
log_counts_corrected <- sapply(corrected[,-1], as.numeric)
rownames(log_counts_corrected) <- corrected[,1]
fit <- lmFit(log_counts_corrected , design_ruvg)
fit <- contrasts.fit(fit, contrasts=all_contrasts)
fit <- eBayes(fit, trend=TRUE)
dna_top_table_limma <- topTable(fit, coef=1,number = nrow(x_dna$counts))
hist(dna_top_table_limma$adj.P.Val,breaks=20,main = paste("Histogram of adj.P-value\nLimma N=" , nrow(x_dna$counts)),
     xlab = "adj.P-value")
hist(dna_top_table_limma$logFC,breaks=50,main = paste("Histogram of logFC\nLimma N=" , nrow(x_dna$counts)),
     xlab = "logFC")
dev.off()
####Heatmap
cairo_pdf(file="./MG_infectionCovariate/limma_afterRUVg_sigHeatmap.pdf",bg="transparent",
          width = 5,height = 4,family = "Times New Roman")
dge_limma <- subset(dna_top_table_limma,dna_top_table_limma$adj.P.Val<=0.05 & 
                      (dna_top_table_limma$logFC>=2|dna_top_table_limma$logFC<=-2))
dge_limma <- left_join(rownames_to_column(dge_limma,var = "defined_taxa"),a_logCPM,by=c("defined_taxa"))
dge_limma <- dge_limma %>% mutate(across(colnames(dge_limma)[-1],as.numeric))
dge_limma <- dge_limma[order(dge_limma$logFC,decreasing = T), , drop = FALSE]
mat_dna <- as.matrix(dge_limma[,8:19])
rownames(mat_dna) <- dge_limma[,1]
print(pheatmap(mat_dna, cluster_rows = F, cluster_cols = T,fontsize_row = 4,
               main = "Bi-clustering using hclust of log CPMs\nSignificant after limma"))
dev.off()
####logFC and logCPM bar plots
cairo_pdf(file="./MG_infectionCovariate/limma_afterRUVg_sig_logFC_logCPM_barPlots.pdf",bg="transparent",
          width = 4,height = 4,family = "Times New Roman")
logFC_barPlot <- ggplot(dge_limma,aes(x=reorder(defined_taxa,-logFC),y=logFC)) + 
  geom_col() + xlab(element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

m_mat_dna <- reshape2::melt(mat_dna,value.name = "log_CPM",id.vars=c("defined_taxa"))
m_mat_dna$Condition <- str_split_fixed(m_mat_dna$Var2,"_",n=3)[,1]
m_mat_dna$Type=str_split_fixed(m_mat_dna$Var2,"_",n=3)[,3]
m_mat_dna$Replicate=str_split_fixed(m_mat_dna$Var2,"_",n=3)[,2]
m_mat_dna$log_CPM <- as.numeric(m_mat_dna$log_CPM)
logCPM_boxPlot <- ggplot(m_mat_dna,aes(x=Var1,y=log_CPM,fill=fct_relevel(Condition,"Uninfected","Infected"))) + 
  geom_boxplot(outlier.size = 0.5) + xlab(element_blank()) + ylab("log CPM") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = "italic"),
        legend.position = "bottom",legend.title = element_blank()) + 
  scale_fill_manual(values = c("blue","red"))

theme_set(theme_bw(base_size = 12,base_family = "Times New Roman") +
            theme(strip.background = element_rect(fill = "transparent",colour = "NA"),
                  panel.background = element_rect(fill='transparent',colour = "NA"), #transparent panel bg
                  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                  # panel.grid.major = element_blank(), #remove major gridlines
                  panel.grid.minor = element_blank(), #remove minor gridlines
                  legend.background = element_rect(fill='transparent',colour = "NA"), #transparent legend bg
                  legend.box.background = element_rect(fill='transparent',colour = "NA"),
                  legend.key = element_rect(fill = "transparent",colour = "NA")))
print(ggpubr::ggarrange(logFC_barPlot,logCPM_boxPlot,nrow = 2,ncol = 1,heights = c(1,5),align = "v"))
dev.off()

####edgeR####
cairo_pdf(file="./MG_infectionCovariate/pvalue_fdr_edgeR_afterRUVg.pdf",bg="transparent",
          width = 10,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 4))
design_edgeR <- model.matrix(~x_dna_ruvg.samples.group + W_1 + W_2, data=pData(set2))
design_edgeR[,"x_dna_ruvg.samples.groupUninfected_DNA"] <- c(rep(0,6),rep(1,6))
y <- DGEList(counts=counts(set2), group=x_dna$samples$group)
y <- calcNormFactors(y, method="TMM")
y <- estimateDisp(y, design_edgeR)

plotBCV(y, xlab="Average log CPM", ylab="Biological coefficient of variation",
        pch=16, cex=0.2, col.common="red", col.trend="blue", col.tagwise="black")
fit <- glmFit(y, design_edgeR)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
hist(top$PValue,breaks=20,main = paste("Histogram of P-value\nedgeR N=" , nrow(x_dna$counts)),
     xlab = "P-value")
hist(top$FDR,breaks=20,main = paste("Histogram of FDR\nedgeR N=" , nrow(x_dna$counts)),
     xlab = "FDR")
plot(top$logCPM, top$logFC,main = paste("MA plot\nedgeR N=" , nrow(x_dna$counts)),
     xlab="Average log CPM", ylab="logFC",pch=20,cex=0.5)
points(subset(top$logCPM,top$PValue<=0.05 & top$FDR <= 0.05 & (top$logFC>=2|top$logFC<=-2)),
       subset(top$logFC,top$PValue<=0.05 & top$FDR <= 0.05 & (top$logFC>=2|top$logFC<=-2)),pch=20,cex=0.5,col="red")
dev.off()

####Heatmap
cairo_pdf(file="./MG_infectionCovariate/edgeR_afterRUVg_sigHeatmap.pdf",bg="transparent",
          width = 5,height = 5,family = "Times New Roman")
dge_edgeR <- subset(top,top$PValue<=0.05 & top$FDR <= 0.05 & (top$logFC>=2|top$logFC<=-2))
dge_edgeR <- left_join(rownames_to_column(dge_edgeR,var = "defined_taxa"),a_logCPM,by=c("defined_taxa"))
dge_edgeR <- dge_edgeR %>% mutate(across(colnames(dge_edgeR)[-1],as.numeric))
dge_edgeR <- dge_edgeR[order(dge_edgeR$logFC,decreasing = T), , drop = FALSE]
mat_dna <- as.matrix(dge_edgeR[,7:18])
rownames(mat_dna) <- dge_edgeR[,1]
print(pheatmap(mat_dna, cluster_rows = F, cluster_cols = T,fontsize_row = 4,
               main = "Bi-clustering using hclust of log CPMs\nSignificant after edgeR"))
dev.off()
####logFC and logCPM bar plots
cairo_pdf(file="./MG_infectionCovariate/edgeR_afterRUVg_sig_logFC_logCPM_barPlots.pdf",bg="transparent",
          width = 4,height = 4,family = "Times New Roman")
logFC_barPlot <- ggplot(dge_edgeR,aes(x=reorder(defined_taxa,-logFC),y=logFC)) + 
  geom_col() + xlab(element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

m_mat_dna <- reshape2::melt(mat_dna,value.name = "log_CPM",id.vars=c("defined_taxa"))
m_mat_dna$Condition <- str_split_fixed(m_mat_dna$Var2,"_",n=3)[,1]
m_mat_dna$Type=str_split_fixed(m_mat_dna$Var2,"_",n=3)[,3]
m_mat_dna$Replicate=str_split_fixed(m_mat_dna$Var2,"_",n=3)[,2]
m_mat_dna$log_CPM <- as.numeric(m_mat_dna$log_CPM)
logCPM_boxPlot <- ggplot(m_mat_dna,aes(x=Var1,y=log_CPM,fill=fct_relevel(Condition,"Uninfected","Infected"))) + 
  geom_boxplot(outlier.size = 0.5) + xlab(element_blank()) + ylab("log CPM") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = "italic"),
        legend.position = "bottom",legend.title = element_blank()) + 
  scale_fill_manual(values = c("blue","red"))

theme_set(theme_bw(base_size = 12,base_family = "Times New Roman") +
            theme(strip.background = element_rect(fill = "transparent",colour = "NA"),
                  panel.background = element_rect(fill='transparent',colour = "NA"), #transparent panel bg
                  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                  # panel.grid.major = element_blank(), #remove major gridlines
                  panel.grid.minor = element_blank(), #remove minor gridlines
                  legend.background = element_rect(fill='transparent',colour = "NA"), #transparent legend bg
                  legend.box.background = element_rect(fill='transparent',colour = "NA"),
                  legend.key = element_rect(fill = "transparent",colour = "NA")))
print(ggpubr::ggarrange(logFC_barPlot,logCPM_boxPlot,nrow = 2,ncol = 1,heights = c(1,5),align = "v"))
dev.off()

####Metatranscriptome####
####Filtering MT_cpm_ge64_3samples####
x_rna <- edgeR::DGEList(counts=assay(data[,13:24]),group = data$Condition[13:24],samples =colData(data[,13:24]),
                        genes = rowData(data))
CPM <- cpm(x_rna,lib.size=colSums(x_rna$counts))
if(getwd()=="/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/qc_taxa_count"){
  CPM.Cutoff <- 64
} else {
  CPM.Cutoff <- 8
}

tol <- 1e-14
MinSampleSize <- 3
keep.CPM <- rowSums(CPM >= CPM.Cutoff) >= (MinSampleSize - tol)

x_rna <- x_rna[keep.CPM,,keep.lib.sizes=FALSE]
####RUVg####
design_ruvg <- model.matrix(~0+group, data = data.frame(cbind(samples=x_rna$samples$Samples,group=x_rna$samples$group),
                                                        row.names = x_rna$samples$Samples,stringsAsFactors = F))
design_ruvg[,"group1"] <- scales::rescale(c(rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000)))
colnames(design_ruvg) <- c("wt","mock")
all_contrasts <- makeContrasts(rna=wt-mock,levels=colnames(design_ruvg))

normalize_method <- "TMM"
x_rna_ruvg <- calcNormFactors(x_rna, method=normalize_method)
x_rna_ruvg <- estimateGLMCommonDisp(x_rna_ruvg, design_ruvg)
x_rna_ruvg <- estimateGLMTagwiseDisp(x_rna_ruvg, design_ruvg)

fit <- glmFit(x_rna_ruvg, design_ruvg)
lrt <- glmLRT(fit, contrast=all_contrasts)

top <- topTags(lrt, n=nrow(x_rna))$table
empirical <- x_rna_ruvg$genes$defined_taxa[which(!(x_rna_ruvg$genes$defined_taxa %in% top$defined_taxa[1:250]))]
####empirical RUVg####
a <- floor(x_rna_ruvg$counts)
row.names(a) <- x_rna_ruvg$genes$defined_taxa
set <- newSeqExpressionSet(a,phenoData = data.frame(x_rna_ruvg$samples$group, row.names=colnames(x_rna_ruvg$counts)))
set2 <- RUVg(set, empirical, k=2)
i=2
####Plotting counts####
cairo_pdf(file="./MT_infectionCovariate/logCPM_distribution_pca_rle_afterRUVg.pdf",bg="transparent",
          width = 11,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 5))
a_logCPM <- data.frame(cbind("defined_taxa" =rownames(a),edgeR::cpm.default(a,log=T)), stringsAsFactors = FALSE)
m_a <- reshape2::melt(a_logCPM,value.name = "log_cpm",id.vars=c("defined_taxa"))
m_a$Condition <- gsub("_\\d+","",m_a$variable)
m_a$Tpe=str_split_fixed(m_a$variable,"_",n=3)[,3]
m_a$Replicate=str_split_fixed(m_a$variable,"_",n=3)[,2]
m_a$log_cpm <- as.numeric(m_a$log_cpm)

vioplot::vioplot(m_a$log_cpm ~ m_a$variable, col=c(rep("blue",6),rep("red",6)), cex=0.4,
                 main=paste0("Log CPM violin plot\nN= ",nrow(a)), xlab="Replicates",ylab="Log CPM",xaxt = "n")

normalisedCounts <- data.frame(cbind("defined_taxa" =rownames(normCounts(set2)),edgeR::cpm.default(normCounts(set2),log=T)), stringsAsFactors = FALSE)
m_normalisedCounts <- reshape2::melt(normalisedCounts,value.name = "normalised_count",id.vars=c("defined_taxa"))
m_normalisedCounts$Condition <- gsub("_\\d+","",m_normalisedCounts$variable)
m_normalisedCounts$Tpe=str_split_fixed(m_normalisedCounts$variable,"_",n=3)[,3]
m_normalisedCounts$Replicate=str_split_fixed(m_normalisedCounts$variable,"_",n=3)[,2]
m_normalisedCounts$normalised_count <- as.numeric(m_normalisedCounts$normalised_count)

vioplot::vioplot(m_normalisedCounts$normalised_count ~ m_normalisedCounts$variable, col=c(rep("blue",6),rep("red",6)), cex=0.4,
                 main=paste0("Normalised Count from RUVg\nLog CPM violin plot k=",i,"\nN= ",nrow(normalisedCounts)), xlab="Replicates",ylab="Log CPM",xaxt = "n")

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=c(rep("blue",6),rep("red",6)),main=paste0("RLE plot k=",i),xaxt = "n", xlab="Replicates",ylab="RLE")

plotPCA(set2, col=c(rep("blue",6),rep("red",6)), cex=0.4,main=paste0("PCA plot k=",i),xaxt = "n")

corrected <- removeBatchEffect(cpm(x_rna, log=TRUE), design=design_ruvg, covariates=as.matrix(pData(set2)[,c("W_1","W_2")]))
corrected <- data.frame(cbind("defined_taxa" =x_rna$genes$defined_taxa,corrected), stringsAsFactors = FALSE)
m_corrected <- reshape2::melt(corrected,value.name = "corrected_logCPM",id.vars=c("defined_taxa"))
m_corrected$Condition <- gsub("_\\d+","",m_corrected$variable)
m_corrected$Tpe=str_split_fixed(m_corrected$variable,"_",n=3)[,3]
m_corrected$Replicate=str_split_fixed(m_corrected$variable,"_",n=3)[,2]
m_corrected$corrected_logCPM <- as.numeric(m_corrected$corrected_logCPM)
vioplot::vioplot(m_corrected$corrected_logCPM ~ m_corrected$variable, col=c(rep("blue",6),rep("red",6)), cex=0.4,
                 main=paste0("Corrected counts after RUVg\nLog CPM violin plot k=",i,"\nN= ",nrow(corrected)), xlab="Replicates",ylab="Log CPM",xaxt = "n")
dev.off()
####Calculating normalising factors####
x_rna <- calcNormFactors(x_rna,method = normalize_method)
####limma no RUVg####
cairo_pdf(file="./MT_infectionCovariate/pvalue_limma_noRUVg.pdf",bg="transparent",
          width = 6,height = 3.5,family = "Times New Roman")
par(mfrow = c(1, 2))
logCPM_afterNormFac_rna <- edgeR::cpm(x_rna, log=TRUE)
row.names(logCPM_afterNormFac_rna) <- x_rna$genes$defined_taxa

fit <- lmFit(logCPM_afterNormFac_rna , design_ruvg)
fit <- contrasts.fit(fit, contrasts=all_contrasts)
fit <- eBayes(fit, trend=TRUE)
rna_top_table_noruvg <- topTable(fit, coef=1,number = nrow(x_rna$counts))
hist(rna_top_table_noruvg$adj.P.Val,breaks=20,main = paste("Histogram of adj.P-value\nLimma (No RUVg) N=" , nrow(x_rna$counts)),
     xlab = "adj.P-value")
hist(rna_top_table_noruvg$logFC,breaks=50,main = paste("Histogram of logFC\nLimma N=" , nrow(x_rna$counts)),
     xlab = "logFC")
dev.off()

####Heatmap
cairo_pdf(file="./MT_infectionCovariate/limma_noRUVg_sigHeatmap.pdf",bg="transparent",
          width = 6,height = 5,family = "Times New Roman")
dge_limma <- subset(rna_top_table_noruvg,rna_top_table_noruvg$adj.P.Val<=0.05 & 
                      (rna_top_table_noruvg$logFC>=2|rna_top_table_noruvg$logFC<=-2))
dge_limma <- left_join(rownames_to_column(dge_limma,var = "defined_taxa"),a_logCPM,by=c("defined_taxa"))
dge_limma <- dge_limma %>% mutate(across(colnames(dge_limma)[-1],as.numeric))
dge_limma <- dge_limma[order(dge_limma$logFC,decreasing = T), , drop = FALSE]
mat_rna <- as.matrix(dge_limma[,8:19])
rownames(mat_rna) <- dge_limma[,1]
print(pheatmap(mat_rna, cluster_rows = F, cluster_cols = T,fontsize_row = 4,
               main = "Bi-clustering using hclust of log CPMs\nSignificant after limma (No RUVg)"))
dev.off()
####logFC and logCPM bar plots
cairo_pdf(file="./MT_infectionCovariate/limma_noRUVg_sig_logFC_logCPM_barPlots.pdf",bg="transparent",
          width = 5,height = 5,family = "Times New Roman")
logFC_barPlot <- ggplot(dge_limma,aes(x=reorder(defined_taxa,-logFC),y=logFC)) + 
  geom_col() + xlab(element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

m_mat_rna <- reshape2::melt(mat_rna,value.name = "log_CPM",id.vars=c("defined_taxa"))
m_mat_rna$Condition <- str_split_fixed(m_mat_rna$Var2,"_",n=3)[,1]
m_mat_rna$Type=str_split_fixed(m_mat_rna$Var2,"_",n=3)[,3]
m_mat_rna$Replicate=str_split_fixed(m_mat_rna$Var2,"_",n=3)[,2]
m_mat_rna$log_CPM <- as.numeric(m_mat_rna$log_CPM)
logCPM_boxPlot <- ggplot(m_mat_rna,aes(x=Var1,y=log_CPM,fill=fct_relevel(Condition,"Uninfected","Infected"))) + 
  geom_boxplot(outlier.size = 0.5) + xlab(element_blank()) + ylab("log CPM") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = "italic"),
        legend.position = "bottom",legend.title = element_blank()) + 
  scale_fill_manual(values = c("blue","red"))

theme_set(theme_bw(base_size = 12,base_family = "Times New Roman") +
            theme(strip.background = element_rect(fill = "transparent",colour = "NA"),
                  panel.background = element_rect(fill='transparent',colour = "NA"), #transparent panel bg
                  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                  # panel.grid.major = element_blank(), #remove major gridlines
                  panel.grid.minor = element_blank(), #remove minor gridlines
                  legend.background = element_rect(fill='transparent',colour = "NA"), #transparent legend bg
                  legend.box.background = element_rect(fill='transparent',colour = "NA"),
                  legend.key = element_rect(fill = "transparent",colour = "NA")))
print(ggpubr::ggarrange(logFC_barPlot,logCPM_boxPlot,nrow = 2,ncol = 1,heights = c(1,5),align = "v"))
dev.off()

####limma-voom no RUVg####
cairo_pdf(file="./MT_infectionCovariate/meanVariance_limmaVoom_noRUVg.pdf",bg="transparent",
          width = 8,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 3))
logCPM_afterNormFac_rna <- edgeR::cpm(x_rna, log=TRUE)
row.names(logCPM_afterNormFac_rna) <- x_rna$genes$defined_taxa
v <- voom(x_rna,design_ruvg,plot = T)
fit <- lmFit(v , design_ruvg)
fit <- contrasts.fit(fit, contrasts=all_contrasts)
fit <- eBayes(fit, trend=TRUE)
rna_top_table_limmaVoom_noruvg <- topTable(fit, coef=1,number = nrow(x_rna$counts))
hist(rna_top_table_limmaVoom_noruvg$adj.P.Val,breaks=50,main = paste("Histogram of adj.P-value\nLimma (No RUVg) N=" , nrow(x_rna$counts)),
     xlab = "adj.P-value")
hist(rna_top_table_limmaVoom_noruvg$logFC,breaks=50,main = paste("Histogram of logFC\nLimma N=" , nrow(x_rna$counts)),
     xlab = "logFC")
dev.off()
####Heatmap
cairo_pdf(file="./MT_infectionCovariate/limmaVoom_noRUVg_sigHeatmap.pdf",bg="transparent",
          width = 6,height = 5,family = "Times New Roman")
dge_limma <- subset(rna_top_table_limmaVoom_noruvg,rna_top_table_limmaVoom_noruvg$adj.P.Val<=0.05 & 
                      (rna_top_table_limmaVoom_noruvg$logFC>=2|rna_top_table_limmaVoom_noruvg$logFC<=-2))
dge_limma <- left_join(dge_limma,a_logCPM,by=c("defined_taxa"))
dge_limma <- dge_limma %>% mutate(across(colnames(dge_limma)[-1],as.numeric))
dge_limma <- dge_limma[order(dge_limma$logFC,decreasing = T), , drop = FALSE]
mat_rna <- as.matrix(dge_limma[,17:28])
rownames(mat_rna) <- dge_limma[,1]
print(pheatmap(mat_rna, cluster_rows = F, cluster_cols = T,fontsize_row = 4,
               main = "Bi-clustering using hclust of log CPMs\nSignificant after limma (No RUVg)"))
dev.off()
####logFC and logCPM bar plots
cairo_pdf(file="./MT_infectionCovariate/limmaVoom_noRUVg_sig_logFC_logCPM_barPlot.pdf",bg="transparent",
          width = 6,height = 5,family = "Times New Roman")
logFC_barPlot <- ggplot(dge_limma,aes(x=reorder(defined_taxa,-logFC),y=logFC)) + 
  geom_col() + xlab(element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

m_mat_rna <- reshape2::melt(mat_rna,value.name = "log_CPM",id.vars=c("defined_taxa"))
m_mat_rna$Condition <- str_split_fixed(m_mat_rna$Var2,"_",n=3)[,1]
m_mat_rna$Type=str_split_fixed(m_mat_rna$Var2,"_",n=3)[,3]
m_mat_rna$Replicate=str_split_fixed(m_mat_rna$Var2,"_",n=3)[,2]
m_mat_rna$log_CPM <- as.numeric(m_mat_rna$log_CPM)
logCPM_boxPlot <- ggplot(m_mat_rna,aes(x=Var1,y=log_CPM,fill=fct_relevel(Condition,"Uninfected","Infected"))) + 
  geom_boxplot(outlier.size = 0.5) + xlab(element_blank()) + ylab("log CPM") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = "italic"),
        legend.position = "bottom",legend.title = element_blank()) + 
  scale_fill_manual(values = c("blue","red"))

theme_set(theme_bw(base_size = 12,base_family = "Times New Roman") +
            theme(strip.background = element_rect(fill = "transparent",colour = "NA"),
                  panel.background = element_rect(fill='transparent',colour = "NA"), #transparent panel bg
                  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                  # panel.grid.major = element_blank(), #remove major gridlines
                  panel.grid.minor = element_blank(), #remove minor gridlines
                  legend.background = element_rect(fill='transparent',colour = "NA"), #transparent legend bg
                  legend.box.background = element_rect(fill='transparent',colour = "NA"),
                  legend.key = element_rect(fill = "transparent",colour = "NA")))
print(ggpubr::ggarrange(logFC_barPlot,logCPM_boxPlot,nrow = 2,ncol = 1,heights = c(1,5),align = "v"))
dev.off()
####voom limma####
cairo_pdf(file="./MT_infectionCovariate/meanVariance_limmaVoom_afterRUVg.pdf",bg="transparent",
          width = 8,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 3))
design_voom <- cbind(design_ruvg,set2$W_1,set2$W_2)
colnames(design_voom) <- c("wt","mock","W_1","W_2")
v <- voom(x_rna,design_voom,plot = T)
fit <- lmFit(v , design_ruvg)
fit <- contrasts.fit(fit, contrasts=all_contrasts)
fit <- eBayes(fit, trend=TRUE)
rna_top_table_limmaVoom <- topTable(fit, coef=1,number = nrow(x_rna$counts))
hist(rna_top_table_limmaVoom$adj.P.Val,breaks=20,main = paste("Histogram of adj.P-value\nLimma-voom N=" , nrow(x_rna$counts)),
     xlab = "adj.P-value")
hist(rna_top_table_limmaVoom$logFC,breaks=50,main = paste("Histogram of logFC\nLimma-voom N=" , nrow(x_rna$counts)),
     xlab = "logFC")
dim(subset(rna_top_table_limmaVoom,rna_top_table_limmaVoom$adj.P.Val<=0.05 & 
             (rna_top_table_limmaVoom$logFC>=2|rna_top_table_limmaVoom$logFC<=-2)))
dev.off()
####Heatmap
cairo_pdf(file="./MT_infectionCovariate/limmaVoom_afterRUVg_sigHeatmap.pdf",bg="transparent",
          width = 5,height = 4,family = "Times New Roman")
dge_limma <- subset(rna_top_table_limmaVoom,rna_top_table_limmaVoom$adj.P.Val<=0.05 & 
                      (rna_top_table_limmaVoom$logFC>=2|rna_top_table_limmaVoom$logFC<=-2))
dge_limma <- left_join(dge_limma,a_logCPM,by=c("defined_taxa"))
dge_limma <- dge_limma %>% mutate(across(colnames(dge_limma)[-1],as.numeric))
dge_limma <- dge_limma[order(dge_limma$logFC,decreasing = T), , drop = FALSE]
mat_rna <- as.matrix(dge_limma[,17:28])
rownames(mat_rna) <- dge_limma[,1]
print(pheatmap(mat_rna, cluster_rows = F, cluster_cols = T,fontsize_row = 4,
               main = "Bi-clustering using hclust of log CPMs\nSignificant after limma"))
dev.off()
####logFC and logCPM bar plots
cairo_pdf(file="./MT_infectionCovariate/limmaVoom_afterRUVg_sig_logFC_logCPM_barPlots.pdf",bg="transparent",
          width = 4,height = 4,family = "Times New Roman")
logFC_barPlot <- ggplot(dge_limma,aes(x=reorder(defined_taxa,-logFC),y=logFC)) + 
  geom_col() + xlab(element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

m_mat_rna <- reshape2::melt(mat_rna,value.name = "log_CPM",id.vars=c("defined_taxa"))
m_mat_rna$Condition <- str_split_fixed(m_mat_rna$Var2,"_",n=3)[,1]
m_mat_rna$Type=str_split_fixed(m_mat_rna$Var2,"_",n=3)[,3]
m_mat_rna$Replicate=str_split_fixed(m_mat_rna$Var2,"_",n=3)[,2]
m_mat_rna$log_CPM <- as.numeric(m_mat_rna$log_CPM)
logCPM_boxPlot <- ggplot(m_mat_rna,aes(x=Var1,y=log_CPM,fill=fct_relevel(Condition,"Uninfected","Infected"))) + 
  geom_boxplot(outlier.size = 0.5) + xlab(element_blank()) + ylab("log CPM") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = "italic"),
        legend.position = "bottom",legend.title = element_blank()) + 
  scale_fill_manual(values = c("blue","red"))

theme_set(theme_bw(base_size = 12,base_family = "Times New Roman") +
            theme(strip.background = element_rect(fill = "transparent",colour = "NA"),
                  panel.background = element_rect(fill='transparent',colour = "NA"), #transparent panel bg
                  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                  # panel.grid.major = element_blank(), #remove major gridlines
                  panel.grid.minor = element_blank(), #remove minor gridlines
                  legend.background = element_rect(fill='transparent',colour = "NA"), #transparent legend bg
                  legend.box.background = element_rect(fill='transparent',colour = "NA"),
                  legend.key = element_rect(fill = "transparent",colour = "NA")))
print(ggpubr::ggarrange(logFC_barPlot,logCPM_boxPlot,nrow = 2,ncol = 1,heights = c(1,5),align = "v"))
dev.off()

####limma####
cairo_pdf(file="./MT_infectionCovariate/pvalue_limma_afterRUVg.pdf",bg="transparent",
          width = 6,height = 3.5,family = "Times New Roman")
par(mfrow = c(1, 2))
log_counts_corrected <- sapply(corrected[,-1], as.numeric)
rownames(log_counts_corrected) <- corrected[,1]
fit <- lmFit(log_counts_corrected , design_ruvg)
fit <- contrasts.fit(fit, contrasts=all_contrasts)
fit <- eBayes(fit, trend=TRUE)
rna_top_table_limma <- topTable(fit, coef=1,number = nrow(x_rna$counts))
hist(rna_top_table_limma$adj.P.Val,breaks=20,main = paste("Histogram of adj.P-value\nLimma N=" , nrow(x_rna$counts)),
     xlab = "adj.P-value")
hist(rna_top_table_limma$logFC,breaks=50,main = paste("Histogram of logFC\nLimma N=" , nrow(x_rna$counts)),
     xlab = "logFC")
dev.off()
####Heatmap
cairo_pdf(file="./MT_infectionCovariate/limma_afterRUVg_sigHeatmap.pdf",bg="transparent",
          width = 5,height = 4,family = "Times New Roman")
dge_limma <- subset(rna_top_table_limma,rna_top_table_limma$adj.P.Val<=0.05 & 
                      (rna_top_table_limma$logFC>=2|rna_top_table_limma$logFC<=-2))
dge_limma <- left_join(rownames_to_column(dge_limma,var = "defined_taxa"),a_logCPM,by=c("defined_taxa"))
dge_limma <- dge_limma %>% mutate(across(colnames(dge_limma)[-1],as.numeric))
dge_limma <- dge_limma[order(dge_limma$logFC,decreasing = T), , drop = FALSE]
mat_rna <- as.matrix(dge_limma[,8:19])
rownames(mat_rna) <- dge_limma[,1]
print(pheatmap(mat_rna, cluster_rows = F, cluster_cols = T,fontsize_row = 4,
               main = "Bi-clustering using hclust of log CPMs\nSignificant after limma"))
dev.off()
####logFC and logCPM bar plots
cairo_pdf(file="./MT_infectionCovariate/limma_afterRUVg_sig_logFC_logCPM_barPlots.pdf",bg="transparent",
          width = 4,height = 4,family = "Times New Roman")
logFC_barPlot <- ggplot(dge_limma,aes(x=reorder(defined_taxa,-logFC),y=logFC)) + 
  geom_col() + xlab(element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

m_mat_rna <- reshape2::melt(mat_rna,value.name = "log_CPM",id.vars=c("defined_taxa"))
m_mat_rna$Condition <- str_split_fixed(m_mat_rna$Var2,"_",n=3)[,1]
m_mat_rna$Type=str_split_fixed(m_mat_rna$Var2,"_",n=3)[,3]
m_mat_rna$Replicate=str_split_fixed(m_mat_rna$Var2,"_",n=3)[,2]
m_mat_rna$log_CPM <- as.numeric(m_mat_rna$log_CPM)
logCPM_boxPlot <- ggplot(m_mat_rna,aes(x=Var1,y=log_CPM,fill=fct_relevel(Condition,"Uninfected","Infected"))) + 
  geom_boxplot(outlier.size = 0.5) + xlab(element_blank()) + ylab("log CPM") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = "italic"),
        legend.position = "bottom",legend.title = element_blank()) + 
  scale_fill_manual(values = c("blue","red"))

theme_set(theme_bw(base_size = 12,base_family = "Times New Roman") +
            theme(strip.background = element_rect(fill = "transparent",colour = "NA"),
                  panel.background = element_rect(fill='transparent',colour = "NA"), #transparent panel bg
                  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                  # panel.grid.major = element_blank(), #remove major gridlines
                  panel.grid.minor = element_blank(), #remove minor gridlines
                  legend.background = element_rect(fill='transparent',colour = "NA"), #transparent legend bg
                  legend.box.background = element_rect(fill='transparent',colour = "NA"),
                  legend.key = element_rect(fill = "transparent",colour = "NA")))
print(ggpubr::ggarrange(logFC_barPlot,logCPM_boxPlot,nrow = 2,ncol = 1,heights = c(1,5),align = "v"))
dev.off()

####edgeR####
cairo_pdf(file="./MT_infectionCovariate/pvalue_fdr_edgeR_afterRUVg.pdf",bg="transparent",
          width = 10,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 4))
design_edgeR <- model.matrix(~x_rna_ruvg.samples.group + W_1 + W_2, data=pData(set2))
design_edgeR[,"x_rna_ruvg.samples.groupUninfected_RNA"] <- c(rep(0,6),rep(1,6))
y <- DGEList(counts=counts(set2), group=x_rna$samples$group)
y <- calcNormFactors(y, method="TMM")
y <- estimateDisp(y, design_edgeR)

plotBCV(y, xlab="Average log CPM", ylab="Biological coefficient of variation",
        pch=16, cex=0.2, col.common="red", col.trend="blue", col.tagwise="black")
fit <- glmFit(y, design_edgeR)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
hist(top$PValue,breaks=20,main = paste("Histogram of P-value\nedgeR N=" , nrow(x_rna$counts)),
     xlab = "P-value")
hist(top$FDR,breaks=20,main = paste("Histogram of FDR\nedgeR N=" , nrow(x_rna$counts)),
     xlab = "FDR")
plot(top$logCPM, top$logFC,main = paste("MA plot\nedgeR N=" , nrow(x_dna$counts)),
     xlab="Average log CPM", ylab="logFC",pch=20,cex=0.5)
points(subset(top$logCPM,top$PValue<=0.05 & top$FDR <= 0.05 & (top$logFC>=2|top$logFC<=-2)),
       subset(top$logFC,top$PValue<=0.05 & top$FDR <= 0.05 & (top$logFC>=2|top$logFC<=-2)),pch=20,cex=0.5,col="red")
dev.off()

####Heatmap
cairo_pdf(file="./MT_infectionCovariate/edgeR_afterRUVg_sigHeatmap.pdf",bg="transparent",
          width = 5,height = 6,family = "Times New Roman")
dge_edgeR <- dge_edgeR <- subset(top,top$PValue<=0.05 & top$FDR <= 0.05 & (top$logFC>=2|top$logFC<=-2))
dge_edgeR <- left_join(rownames_to_column(dge_edgeR,var = "defined_taxa"),a_logCPM,by=c("defined_taxa"))
dge_edgeR <- dge_edgeR %>% mutate(across(colnames(dge_edgeR)[-1],as.numeric))
dge_edgeR <- dge_edgeR[order(dge_edgeR$logFC,decreasing = T), , drop = FALSE]
mat_rna <- as.matrix(dge_edgeR[,7:18])
rownames(mat_rna) <- dge_edgeR[,1]
print(pheatmap(mat_rna, cluster_rows = F, cluster_cols = T,fontsize_row = 4,
               main = "Bi-clustering using hclust of log CPMs\nSignificant after edgeR"))
dev.off()
####logFC and logCPM bar plots
cairo_pdf(file="./MT_infectionCovariate/edgeR_afterRUVg_sig_logFC_logCPM_barPlots.pdf",bg="transparent",
          width = 6,height = 5,family = "Times New Roman")
logFC_barPlot <- ggplot(dge_edgeR,aes(x=reorder(defined_taxa,-logFC),y=logFC)) + 
  geom_col() + xlab(element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

m_mat_rna <- reshape2::melt(mat_rna,value.name = "log_CPM",id.vars=c("defined_taxa"))
m_mat_rna$Condition <- str_split_fixed(m_mat_rna$Var2,"_",n=3)[,1]
m_mat_rna$Type=str_split_fixed(m_mat_rna$Var2,"_",n=3)[,3]
m_mat_rna$Replicate=str_split_fixed(m_mat_rna$Var2,"_",n=3)[,2]
m_mat_rna$log_CPM <- as.numeric(m_mat_rna$log_CPM)
logCPM_boxPlot <- ggplot(m_mat_rna,aes(x=Var1,y=log_CPM,fill=fct_relevel(Condition,"Uninfected","Infected"))) + 
  geom_boxplot(outlier.size = 0.5) + xlab(element_blank()) + ylab("log CPM") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = "italic"),
        legend.position = "bottom",legend.title = element_blank()) + 
  scale_fill_manual(values = c("blue","red"))

theme_set(theme_bw(base_size = 12,base_family = "Times New Roman") +
            theme(strip.background = element_rect(fill = "transparent",colour = "NA"),
                  panel.background = element_rect(fill='transparent',colour = "NA"), #transparent panel bg
                  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                  # panel.grid.major = element_blank(), #remove major gridlines
                  panel.grid.minor = element_blank(), #remove minor gridlines
                  legend.background = element_rect(fill='transparent',colour = "NA"), #transparent legend bg
                  legend.box.background = element_rect(fill='transparent',colour = "NA"),
                  legend.key = element_rect(fill = "transparent",colour = "NA")))
print(ggpubr::ggarrange(logFC_barPlot,logCPM_boxPlot,nrow = 2,ncol = 1,heights = c(1,5),align = "v"))
dev.off()

