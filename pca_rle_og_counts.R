suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(edgeR))
suppressMessages(library(RUVSeq))
suppressMessages(library(pheatmap))

####MG>=1####
output_dir <- "/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/qc_OG_taxa_count"
if (!dir.exists(output_dir)) {dir.create(output_dir)}
setwd(output_dir)
count_sum <- read.delim("~/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/qc_OG_taxa_count/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGanyr1_noSalBongori_selected_OG_stratified.tsv",header = T,sep = "\t",stringsAsFactors = F)
# count_sum <- read.delim("~/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/qc_OG_taxa_count/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGanyr1_noSalBongori_selected_taxa_stratified.tsv",header = T,sep = "\t",stringsAsFactors = F)

####MG>=40%####
# output_dir <- "/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/qc_OG_taxa_count"
# if (!dir.exists(output_dir)) {dir.create(output_dir)}
# setwd(output_dir)
# count_sum <- read.delim("/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/qc_OG_taxa_count/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGperc40_selected_OG_stratified.tsv",header = T,sep = "\t",stringsAsFactors = F)
# count_sum <- read.delim("/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/qc_OG_taxa_count/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGperc40_selected_taxa_stratified.tsv",header = T,sep = "\t",stringsAsFactors = F)

####Function to plot log CPM dist, PCA and RLE########
logcpm_dist_pca_rle <- function(degObject,outfile_prefix,folder){
  ####Log CPM distribution####
  log_counts <- edgeR::cpm.default(degObject,log=T)
  row.names(log_counts) <- degObject$genes$associated_OG
  log_counts_df <- rownames_to_column(data.frame(log_counts,stringsAsFactors = F),var = "associated_OG")
  m_log_counts_df <- reshape2::melt(log_counts_df,value.name = "log_CPM",id.vars=c("associated_OG"))
  m_log_counts_df$Condition <- gsub("_\\d+","",m_log_counts_df$variable)
  m_log_counts_df$Type=str_split_fixed(m_log_counts_df$variable,"_",n=3)[,3]
  m_log_counts_df$Replicate=str_split_fixed(m_log_counts_df$variable,"_",n=3)[,2]
  m_log_counts_df$log_CPM <- as.numeric(m_log_counts_df$log_CPM)
  log2_count_violin <- ggplot(m_log_counts_df,aes(x=variable,y=log_CPM,group=variable,
                                                  color=fct_relevel(Condition,"Uninfected_DNA","Infected_DNA","Uninfected_RNA","Infected_RNA"))) +
    geom_violin() + 
    ggtitle("Log CPM") + xlab("Replicates") + ylab("log CPM") + 
    theme(legend.title = element_blank(), legend.position = "bottom",axis.text = element_text(color = "black")) +
    guides(color=guide_legend(nrow=2)) + 
    scale_x_discrete(labels=setNames(degObject$samples$Replicate, degObject$samples$Samples))
  
  log2_count_violin_low <- subset(m_log_counts_df,m_log_counts_df$log_CPM<=5) %>% ggplot(aes(x=variable,y=log_CPM,group=variable,
                                                  color=fct_relevel(Condition,"Uninfected_DNA","Infected_DNA","Uninfected_RNA","Infected_RNA"))) +
    geom_violin() + 
    ggtitle("Log CPM <= 5") + xlab("Replicates") + ylab("log CPM") + 
    theme(legend.title = element_blank(), legend.position = "bottom",axis.text = element_text(color = "black")) +
    guides(color=guide_legend(nrow=2)) + 
    scale_x_discrete(labels=setNames(degObject$samples$Replicate, degObject$samples$Samples)) +
    scale_y_continuous(limits = c(-3,5))
  
  log2_count_violin_high <- subset(m_log_counts_df,m_log_counts_df$log_CPM>5) %>% ggplot(aes(x=variable,y=log_CPM,group=variable,
                                                  color=fct_relevel(Condition,"Uninfected_DNA","Infected_DNA","Uninfected_RNA","Infected_RNA"))) +
    geom_violin() + 
    ggtitle("Log CPM > 5") + xlab("Replicates") + ylab("log CPM") + 
    theme(legend.title = element_blank(), legend.position = "bottom",axis.text = element_text(color = "black")) +
    guides(color=guide_legend(nrow=2)) + 
    scale_x_discrete(labels=setNames(degObject$samples$Replicate, degObject$samples$Samples))+
    scale_y_continuous(limits = c(5,18))
  
  ####Perform PCA for both MG and MT####
  pca_data <- prcomp(t(log_counts))
  # Calculate the percentage of rows contributing to each PC
  variance_explained <- round(100 * pca_data$sdev^2 / sum(pca_data$sdev^2),2)
  # Create a PCA plot with labeled and colored points
  pca_data_df <- data.frame(rownames_to_column(data.frame(pca_data$x,stringsAsFactors = f),var = "Sample"),stringsAsFactors = F)
  pca_data_df$Condition <- gsub("_\\d+","",pca_data_df$Sample)
  pca_data_df$Type=str_split_fixed(pca_data_df$Sample,"_",n=3)[,3]
  pca_data_df$Replicate=str_split_fixed(pca_data_df$Sample,"_",n=3)[,2]
  pca <- ggplot(pca_data_df,aes(x=PC1,y=PC2,color=fct_relevel(Condition,"Uninfected_DNA","Infected_DNA","Uninfected_RNA","Infected_RNA"))) + 
    geom_point(size=1) + geom_text(aes(label = Replicate),size = 3.5,vjust = 0, nudge_y = 0.5) +
    ggtitle("PCA Plot") + xlab(paste0("PC1", ": ", variance_explained[1], "%")) +
    ylab(paste0("PC2", ": ", variance_explained[2], "%")) + 
    theme(legend.title = element_blank(), legend.position = "bottom",axis.text = element_text(color = "black")) +
    guides(color=guide_legend(nrow=2))
  
  ####Calculate RLE values####
  # Convert the log_counts columns to numeric
  log_count_columns <- colnames(log_counts_df)[2:length(log_counts_df)]
  log_counts_df[log_count_columns] <- lapply(log_counts_df[log_count_columns], as.numeric)
  
  # Perform median subtraction by group (defined_taxa)
  # Calculate the row-wise median for all columns except associated_OG
  medians <- log_counts_df %>%
    select(-associated_OG) %>%
    mutate(across(everything(), median)) 
  
  # Subtract the row-wise median from all columns except associated_OG
  rle_df <- log_counts_df %>%
    select(associated_OG) %>%
    bind_cols(log_counts_df %>%
                select(-associated_OG) %>%
                mutate(across(everything(), ~ . - medians[[cur_column()]])))
  
  m_rle_df <- reshape2::melt(rle_df,value.name = "rle",id.vars=c("associated_OG"))
  m_rle_df$Condition <- gsub("_\\d+","",m_rle_df$variable)
  m_rle_df$Type=str_split_fixed(m_rle_df$variable,"_",n=3)[,3]
  m_rle_df$Replicate=str_split_fixed(m_rle_df$variable,"_",n=3)[,2]
  
  rle_boxplot <- ggplot(m_rle_df,aes(x=variable,y=rle,color=fct_relevel(Condition,"Uninfected_DNA","Infected_DNA","Uninfected_RNA","Infected_RNA"))) +
    geom_boxplot(outlier.shape = NA) +
    ggtitle("RLE Plot") + xlab("Replicates") + ylab("Relative log CPM expression") + 
    theme(legend.title = element_blank(), legend.position = "bottom",axis.text = element_text(color = "black")) +
    guides(color=guide_legend(nrow=2)) + 
    scale_x_discrete(labels=setNames(degObject$samples$Replicate, degObject$samples$Samples)) + scale_y_continuous(limits = c(-2.5,2.5))
  
  ####Plotting logCPM dist, PCA, RLE to files####
  cairo_pdf(file=paste0("./",folder,"/","logCPM_distribution_pca_rle_",outfile_prefix,".pdf"),bg="transparent",
            width = 12,height = 10,family = "Times New Roman")
  theme_set(theme_bw(base_size = 12,base_family = "Times New Roman") +
              theme(strip.background = element_rect(fill = "transparent",colour = "NA"),
                    panel.background = element_rect(fill='transparent',colour = "NA"), #transparent panel bg
                    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                    # panel.grid.major = element_blank(), #remove major gridlines
                    panel.grid.minor = element_blank(), #remove minor gridlines
                    legend.background = element_rect(fill='transparent',colour = "NA"), #transparent legend bg
                    legend.box.background = element_rect(fill='transparent',colour = "NA"),
                    legend.key = element_rect(fill = "transparent",colour = "NA")))
  print(ggpubr::ggarrange(log2_count_violin,pca,rle_boxplot,
                          log2_count_violin_low,log2_count_violin_high, ncol = 3,nrow = 2,
                          common.legend = T,legend = "bottom"))
  dev.off()
  # cairo_pdf(file=paste0("./",folder,"/","heatmap_biclustering_hclust_",outfile_prefix,".pdf"),bg="transparent",
  #           width = 6,height = 12,family = "Times New Roman")
  # print(pheatmap(log_counts, cluster_rows = T, cluster_cols = T,fontsize_row = 1,main = "Bi-clustering using hclust of log CPMs"))
  # dev.off()
}

ruvg_function <- function(dgeObject,normalize_method,outfile_prefix,folder){
  design <- model.matrix(~0+group, data = data.frame(cbind(samples=dgeObject$samples$Samples,group=dgeObject$samples$group),
                                                     row.names = dgeObject$samples$Samples,stringsAsFactors = F))
  design[,"group1"] <- scales::rescale(c(rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000)))
  colnames(design) <- c("wt","mock")
  all_contrasts <- makeContrasts(dna=wt-mock,levels=colnames(design))
  
  dgeObject <- calcNormFactors(dgeObject, method=normalize_method)
  dgeObject <- estimateGLMCommonDisp(dgeObject, design)
  dgeObject <- estimateGLMTagwiseDisp(dgeObject, design)
  
  fit <- glmFit(dgeObject, design)
  lrt <- glmLRT(fit, contrast=all_contrasts)
  
  top <- topTags(lrt, n=nrow(dgeObject),sort.by = "PValue")$table
  # empirical <- dgeObject$genes$associated_OG[which(!(dgeObject$genes$associated_OG %in% top$associated_OG[1:250]))]
  empirical <- top$associated_OG[which(top$PValue>0.8 & top$FDR>0.1)]
  ## ----emp_ruvg-----------------------------------------------------------------
  a <- floor(dgeObject$counts)
  row.names(a) <- dgeObject$genes$associated_OG
  set <- newSeqExpressionSet(a,phenoData = data.frame(dgeObject$samples$group, row.names=colnames(dgeObject$counts)))
  cairo_pdf(file=paste0("./",folder,"/","ruvg_",normalize_method,"_",outfile_prefix,".pdf"),bg="transparent",
            width = 10,height = 3.5,family = "Times New Roman",onefile = T)
  # par(mfrow = c(5, 2))
  for (i in 1:5) {
    set2 <- RUVg(set, empirical, k=i)
    normalisedCounts <- data.frame(cbind("associated_OG" =rownames(normCounts(set2)),edgeR::cpm.default(normCounts(set2),log=T)), stringsAsFactors = FALSE)
    m_normalisedCounts <- reshape2::melt(normalisedCounts,value.name = "normalised_count",id.vars=c("associated_OG"))
    m_normalisedCounts$Condition <- gsub("_\\d+","",m_normalisedCounts$variable)
    m_normalisedCounts$Tpe=str_split_fixed(m_normalisedCounts$variable,"_",n=3)[,3]
    m_normalisedCounts$Replicate=str_split_fixed(m_normalisedCounts$variable,"_",n=3)[,2]
    m_normalisedCounts$normalised_count <- as.numeric(m_normalisedCounts$normalised_count)
    
    par(mfrow = c(1, 3))
    vioplot::vioplot(m_normalisedCounts$normalised_count ~ m_normalisedCounts$variable, col=c(rep("blue",6),rep("red",6)), cex=0.4,
                     main=paste0("Normalised Count from RUVg\nLog CPM violin plot k=",i,"\nN= ",nrow(normalisedCounts)), xlab="Replicates",ylab="Log CPM",xaxt = "n")
    plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=c(rep("blue",6),rep("red",6)),main=paste0("RLE plot k=",i),xaxt = "n", xlab="Replicates",ylab="RLE")
    plotPCA(set2, col=c(rep("blue",6),rep("red",6)), cex=0.4,main=paste0("PCA plot k=",i),xaxt = "n")
    par(mfrow = c(1, 1))
  }
  dev.off()
}

###Filtering Taxa based on read counts...removing rows with 0 counts in all columns####
keep <- rowSums(count_sum[,7:30] != 0) >= 1
count_removing_zeros <- count_sum[keep,]
count_removing_zeros <- count_removing_zeros[which(count_removing_zeros$associated_OG!="UNKNOWN"),]

#SummarizedExperiment easy to handle all the information at once
colnames(count_removing_zeros) <- gsub("MOCK_","Uninfected_",colnames(count_removing_zeros))
colnames(count_removing_zeros) <- gsub("WT_","Infected_",colnames(count_removing_zeros))

taxa <- as.matrix(count_removing_zeros[,c(1:6)],row.names=count_removing_zeros[,1])


data <- SummarizedExperiment(assays = list(counts=as.matrix(count_removing_zeros[,-c(1:6)],
                                                            row.names=count_removing_zeros[,c(1)])),
                             rowData = taxa,colData = DataFrame(Samples=colnames(count_removing_zeros)[-c(1:6)], 
                                                                Condition=gsub("_\\d+","",colnames(count_removing_zeros)[-c(1:6)]),
                                                                Type=str_split_fixed(colnames(count_removing_zeros)[-c(1:6)],"_",n=3)[,3],
                                                                Replicate=str_split_fixed(colnames(count_removing_zeros)[-c(1:6)],"_",n=3)[,2],
                                                                infectionRate=scales::rescale(c(rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000),rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000)))))

y <- edgeR::DGEList(counts=assay(data),group = data$Condition,samples =colData(data),genes = rowData(data))
logcpm_dist_pca_rle(y,"countsOverZero",NULL)
# ruvg_function(y,"TMM","countsOverZero",NULL)
# ruvg_function(y,"upperquartile","countsOverZero",NULL)

####Metagenome####
if (!dir.exists(paste0(output_dir,"/MG_infectionCovariate"))) {dir.create(paste0(output_dir,"/MG_infectionCovariate"))}
####MG_simple_filterbyExpr####
keep.exprs <- edgeR::filterByExpr(y, group=y$samples$group,min.count=20,min.total.count=25)

###MG_simple_filterbyExpr TMM
x_dna <- edgeR::DGEList(counts=assay(data[keep.exprs,1:12]),group = data$Condition[1:12],samples =colData(data[,1:12]),
                        genes = rowData(data[keep.exprs,]))
ruvg_function(x_dna,"TMM","MG_simple_filterbyExpr","MG_infectionCovariate")

x_dna$samples$infection <- round(scales::rescale(c(rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000))),4)
x_dna <- calcNormFactors(x_dna)
logcpm_dist_pca_rle(x_dna,"MG_simple_filterbyExpr","MG_infectionCovariate")

####MG_simple_filterbyExpr_upperquartile
x_dna <- edgeR::DGEList(counts=assay(data[keep.exprs,1:12]),group = data$Condition[1:12],samples =colData(data[,1:12]),
                        genes = rowData(data[keep.exprs,]))
ruvg_function(x_dna,"upperquartile","MG_simple_filterbyExpr","MG_infectionCovariate")

x_dna$samples$infection <- round(scales::rescale(c(rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000))),4)
x_dna <- calcNormFactors(x_dna,method="upperquartile")
logcpm_dist_pca_rle(x_dna,"MG_simple_filterbyExpr_upperquartile","MG_infectionCovariate")


####MG_cpm_ge32_3samples####
x_dna <- edgeR::DGEList(counts=assay(data[,1:12]),group = data$Condition[1:12],samples =colData(data[,1:12]),
                        genes = rowData(data))
CPM <- cpm(x_dna,lib.size=colSums(x_dna$counts))
if(getwd()=="/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/qc_OG_taxa_count"){
  CPM.Cutoff <- 32
} else {
  CPM.Cutoff <- 32
}

tol <- 1e-14
MinSampleSize <- 3
keep.CPM <- rowSums(CPM >= CPM.Cutoff) >= (MinSampleSize - tol)

x_dna <- x_dna[keep.CPM,,keep.lib.sizes=FALSE]
logcpm_dist_pca_rle(x_dna,"MG_cpm_ge32_3samples","MG_infectionCovariate")
ruvg_function(x_dna,"TMM","MG_cpm_ge32_3samples","MG_infectionCovariate")
ruvg_function(x_dna,"upperquartile","MG_cpm_ge32_3samples","MG_infectionCovariate")

####Metatranscriptome####
if (!dir.exists(paste0(output_dir,"/MT_infectionCovariate"))) {dir.create(paste0(output_dir,"/MT_infectionCovariate"))}
####MT_simple_filterbyExpr####
keep.exprs <- edgeR::filterByExpr(y, group=y$samples$group,min.count=20,min.total.count=25)

###MT_simple_filterbyExpr TMM
x_rna <- edgeR::DGEList(counts=assay(data[keep.exprs,13:24]),group = data$Condition[13:24],samples =colData(data[,13:24]),
                        genes = rowData(data[keep.exprs,]))
ruvg_function(x_rna,"TMM","MT_simple_filterbyExpr","MT_infectionCovariate")

x_rna$samples$infection <- round(scales::rescale(c(rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000))),4)
x_rna <- calcNormFactors(x_rna)
logcpm_dist_pca_rle(x_rna,"MT_simple_filterbyExpr","MT_infectionCovariate")

####MT_simple_filterbyExpr_upperquartile
x_rna <- edgeR::DGEList(counts=assay(data[keep.exprs,13:24]),group = data$Condition[13:24],samples =colData(data[,13:24]),
                        genes = rowData(data[keep.exprs,]))
ruvg_function(x_rna,"upperquartile","MT_simple_filterbyExpr","MT_infectionCovariate")

x_rna$samples$infection <- round(scales::rescale(c(rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000))),4)
x_rna <- calcNormFactors(x_rna,method="upperquartile")
logcpm_dist_pca_rle(x_rna,"MT_simple_filterbyExpr_upperquartile","MT_infectionCovariate")

####MT_cpm_ge32_3samples####
x_rna <- edgeR::DGEList(counts=assay(data[,13:24]),group = data$Condition[13:24],samples =colData(data[,13:24]),
                        genes = rowData(data))
CPM <- cpm(x_rna,lib.size=colSums(x_rna$counts))
if(getwd()=="/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/qc_OG_taxa_count/"){
  CPM.Cutoff <- 32
} else {
  CPM.Cutoff <- 32
}

tol <- 1e-14
MinSampleSize <- 3
keep.CPM <- rowSums(CPM >= CPM.Cutoff) >= (MinSampleSize - tol)

x_rna <- x_rna[keep.CPM,,keep.lib.sizes=FALSE]
logcpm_dist_pca_rle(x_rna,"MT_cpm_ge32_3samples","MT_infectionCovariate")
ruvg_function(x_rna,"TMM","MT_cpm_ge32_3samples","MT_infectionCovariate")
ruvg_function(x_rna,"upperquartile","MT_cpm_ge32_3samples","MT_infectionCovariate")
