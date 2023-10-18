suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(edgeR))
suppressMessages(library(RUVSeq))
suppressMessages(library(pheatmap))

####MG>=1####
output_dir <- "/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/ruvg_edgeR_cfu_covariates"
if (!dir.exists(output_dir)) {dir.create(output_dir)}
setwd(output_dir)
count_sum <- read.delim("~/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGanyr1_noSalBongori_species_summed.tsv",header = T,sep = "\t",stringsAsFactors = F)

####MG>=40%####
# output_dir <- "/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/ruvg_edgeR_cfu_covariates"
# if (!dir.exists(output_dir)) {dir.create(output_dir)}
# setwd(output_dir)
# count_sum <- read.delim("~/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGperc40_species_summed.tsv",header = T,sep = "\t",stringsAsFactors = F)

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
count_removing_zeros$mags[which(count_removing_zeros$defined_taxa %in% gsub(" ","_",mags_with_full_species))] <- "Pangenomes containing MAGs"
count_removing_zeros$mags[which(is.na(count_removing_zeros$mags))] <- "Pangenomes"

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
if(length(grep("/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/",getwd()))==1){
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
cairo_pdf(file="dna_logCPM_distribution_pca_rle_afterRUVg.pdf",bg="transparent",
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
####edgeR####
cairo_pdf(file="dna_pvalue_fdr_edgeR_afterRUVg.pdf",bg="transparent",
          width = 10,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 4),bg="transparent",family="Times New Roman")
design_edgeR <- model.matrix(~x_dna_ruvg.samples.group + W_1 + W_2, data=pData(set2))
design_edgeR[,"x_dna_ruvg.samples.groupUninfected_DNA"] <- c(rep(0,6),rep(1,6))
y_dna <- DGEList(counts=counts(set2), group=x_dna$samples$group)
y_dna <- calcNormFactors(y_dna, method="TMM")
y_dna <- estimateDisp(y_dna, design_edgeR)

plotBCV(y_dna, xlab="Average log CPM", ylab="Biological coefficient of variation",
        pch=16, cex=0.2, col.common="red", col.trend="blue", col.tagwise="black",
        main = paste0("BCV plot\nN=" , nrow(x_dna$counts)))
fit_dna <- glmFit(y_dna, design_edgeR)
lrt_dna <- glmLRT(fit_dna, coef=2)
top_dna <- topTags(lrt_dna, n=nrow(set))$table
hist(top_dna$PValue,breaks=20,main = paste0("Histogram of P-value\nN=" , nrow(x_dna$counts)),
     xlab = "P-value")
hist(top_dna$FDR,breaks=20,main = paste0("Histogram of FDR\nN=" , nrow(x_dna$counts)),
     xlab = "FDR")
plot(top_dna$logCPM, top$logFC,main = paste0("MA plot\nN=" , nrow(x_dna$counts)),
     xlab="Average log CPM", ylab="logFC",pch=20,cex=0.5)
points(subset(top_dna$logCPM,top_dna$PValue<=0.05 & top_dna$FDR <= 0.05 & (top_dna$logFC>=2|top_dna$logFC<=-2)),
       subset(top_dna$logFC,top_dna$PValue<=0.05 & top_dna$FDR <= 0.05 & (top_dna$logFC>=2|top_dna$logFC<=-2)),
       pch=20,cex=0.5,col="red")
dev.off()

top_dna <- left_join(rownames_to_column(top_dna,var = "defined_taxa"),x_dna$genes,by="defined_taxa")
top_dna <- left_join(top_dna,a_logCPM,by="defined_taxa")
top_dna <- top_dna %>% mutate_at(.vars = colnames(a_logCPM)[-1],as.numeric)
write.table(top_dna,"dna_ruvg_edgeR_da_results.tsv",sep = "\t",row.names = F)
####Heatmap Figure 7####
suppressMessages(library(ComplexHeatmap))

##Significant rows dna
dna_limma_daa_sig <- subset(top_dna,top_dna$PValue<=0.05 & top_dna$FDR <= 0.05 & (top_dna$logFC>=2|top_dna$logFC<=-2))

##Renaming MAGs without defined genus or species
# dna_limma_daa_sig$Species[which(dna_limma_daa_sig$Species =="s__")] <- paste0("Unclassified|",sapply(1:length(which(dna_limma_daa_sig$Species =="s__")),FUN=function(y){
#   str_split(dna_limma_daa_sig$defined_taxa[which(dna_limma_daa_sig$Species =="s__")][y],"\\|")[[1]][2]}))
dna_limma_daa_sig$Genus[which(dna_limma_daa_sig$Genus =="g__")] <- paste0("Unclassified|",sapply(1:length(which(dna_limma_daa_sig$Genus =="g__")),FUN=function(y){
  str_split(dna_limma_daa_sig$defined_taxa[which(dna_limma_daa_sig$Genus =="g__")][y],"\\|")[[1]][2]}))
dna_limma_daa_sig <- dna_limma_daa_sig[order(dna_limma_daa_sig$logFC,decreasing = T),]

####Significant species/taxa heatmaps dna and rna####
mat_dna = as.matrix(dna_limma_daa_sig[,16:27])
mat_dna = apply(mat_dna, 2 ,as.numeric)
Sample = gsub("_\\d+_DNA|RNA", "", colnames(mat_dna))
# palette <- randomcoloR::distinctColorPalette(length(unique(c(dna_limma_daa_sig$Class,rna_limma_daa_sig$Class))))
# names(palette) <- as.character(unique(c(dna_limma_daa_sig$Class,rna_limma_daa_sig$Class)))
palette <- c("c__Gammaproteobacteria"= "#7D7CD4",
             "c__Clostridia" = "#D5CDD1",
             "c__Bacteroidia" = "#D49FD5",
             "c__Bacilli" = "#DCE15A",
             "c__Actinobacteria" = "#AE4BE3",
             "c__Verrucomicrobiae" = "#7EBDDB",
             "c__Coriobacteriia" = "#D49FD5",
             "c__Alphaproteobacteria"="#D1D098",
             "c__Bradimonadia"="#7BDBC1",
             "c__Campylobacteria"="#D98471",
             "c__Deinococci"="#7CE574")

ht_dna = Heatmap(dna_limma_daa_sig$Class, name = "Class",width = unit(5, "mm"),
                 rect_gp = gpar(col = "white", lwd = 1),col = palette) +
  Heatmap(mat_dna, name = "log CPM", cluster_rows = T, 
          col=circlize::colorRamp2(c(range(mat_dna)[1], 0, range(mat_dna)[2]), c("#023FA5", "white", "#8E063B")),
          top_annotation = HeatmapAnnotation(Sample = fct_relevel(Sample,"Uninfected","Infected"), 
                                             CFUs = round(scales::rescale(c(rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000))),4),
                                             annotation_name_side = "left",
                                             col = list(Sample = c("Infected"="lightpink2","Uninfected"="olivedrab"),
                                                        CFUs = circlize::colorRamp2(c(0, 0.0002, 1), c("lightgrey","mistyrose","hotpink")))),
          cluster_columns = F, show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE,
          rect_gp = gpar(col = "white", lwd = 1)) +
  rowAnnotation(logFC = anno_barplot(dna_limma_daa_sig$logFC,baseline = 0, width = unit(1, "cm"))) +
  rowAnnotation(Pvalue = anno_barplot(dna_limma_daa_sig$PValue,baseline = 0, width = unit(1, "cm"))) +
  rowAnnotation(FDR = anno_barplot(dna_limma_daa_sig$FDR,baseline = 0, width = unit(1, "cm")))
cairo_pdf("dna_heatmap_logCPMs_logFC_class.pdf",width=6,height=2.5,family = "Times New Roman")
draw(ht_dna,merge_legend = TRUE, heatmap_legend_side = "right",annotation_legend_side = "right")
dev.off()

####Significant Salmonella barplot####
df_sal <- dna_limma_daa_sig[grep("Salmonella",dna_limma_daa_sig$Genus),]
m_df_sal <- reshape2::melt(df_sal,id.vars=c(colnames(df_sal)[1:15]))
m_df_sal$value <- as.numeric(m_df_sal$value)
m_df_sal$condition <- gsub("_\\d+_DNA","",m_df_sal$variable)
m_df_sal$replicate <-  gsub("_DNA","",m_df_sal$variable)
m_df_sal$replicate <-  gsub("Uninfected_|Infected_","",m_df_sal$replicate)
m_df_sal$replicate <- paste0("R",m_df_sal$replicate)
pc_dna <- ggplot(m_df_sal,aes(x=fct_relevel(condition,"Uninfected","Infected"),y=value)) + 
  geom_boxplot(outlier.size = 1) + geom_point(aes(color=replicate)) + 
  theme_bw(base_size = 12,base_family = "Times New Roman") + 
  ylab("log2 CPM") + xlab(element_blank()) +
  # facet_wrap(~defined_taxa,scales = "free",nrow = 1) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(6,"Dark2")) + 
  annotate("text", x=1, y=4, label= paste0("P-value = ",round(unique(m_df_sal$PValue),4),"\nFDR = ",round(unique(m_df_sal$FDR),4))) + 
  theme(axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        legend.position = "right",panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent',colour = NA), #transparent legend bg
        legend.box.background = element_rect(fill='transparent',colour = NA) #transparent legend panel
  ) + 
  ggtitle("Metagenome") + scale_y_continuous(trans = scales::log2_trans(),
                                             breaks = scales::trans_breaks("log2", function(x) 2^x),
                                             labels = scales::trans_format("log2", scales::math_format(2^.x)))

cairo_pdf(filename = "dna_s__Salmonella_enterica_da_boxplot_log2CPM.pdf",width = 4.5,height = 2.5,
          family = "Times New Roman",onefile = T)
print(pc_dna)
dev.off()
