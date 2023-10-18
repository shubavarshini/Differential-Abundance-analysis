suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(edgeR))
suppressMessages(library(RUVSeq))
suppressMessages(library(pheatmap))

####MG>=1####
output_dir <- "/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/qc_OG_taxa_count/ruvg_edgeR_cfu_covariates/MT_infectionCovariate"
if (!dir.exists(output_dir)) {dir.create(output_dir,recursive = T)}
setwd(output_dir)
count_sum <- read.delim("/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_MGanyr1_noSalBongori_salmonQuant/qc_OG_taxa_count/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGanyr1_noSalBongori_selected_OG_stratified.tsv",header = T,sep = "\t",stringsAsFactors = F)
normalize_method <- "TMM"
i=3
####MG>=40%####
# output_dir <- "/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/qc_OG_taxa_count/ruvg_edgeR_cfu_covariates/MT_infectionCovariate"
# if (!dir.exists(output_dir)) {dir.create(output_dir,recursive = T)}
# setwd(output_dir)
# count_sum <- read.delim("/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/salmon_GTDB_iMGMC_ihMAG_markerGenes/mockWT_GTDB_iMGMC_ihMAG_MGperc40_salmonQuant/qc_OG_taxa_count/salmonCount_mockWT_GTDB_iMGMC_ihMAG_MGperc40_selected_OG_stratified.tsv",header = T,sep = "\t",stringsAsFactors = F)
# normalize_method <- "upperquartile"
# i=2

###Filtering Taxa based on read counts...removing rows with 0 counts in all columns####
keep <- rowSums(count_sum[,15:38] != 0) >= 1
count_removing_zeros <- count_sum[keep,]
count_removing_zeros <- count_removing_zeros[which(count_removing_zeros$associated_OG!="UNKNOWN"),]

#SummarizedExperiment easy to handle all the information at once
colnames(count_removing_zeros) <- gsub("MOCK_","Uninfected_",colnames(count_removing_zeros))
colnames(count_removing_zeros) <- gsub("WT_","Infected_",colnames(count_removing_zeros))

taxa <- as.matrix(count_removing_zeros[,c(1:14)],row.names=count_removing_zeros[,1])

data <- SummarizedExperiment(assays = list(counts=as.matrix(count_removing_zeros[,-c(1:14)],
                                                            row.names=count_removing_zeros[,c(1)])),
                             rowData = taxa,colData = DataFrame(Samples=colnames(count_removing_zeros)[-c(1:14)], 
                                                                Condition=gsub("_\\d+","",colnames(count_removing_zeros)[-c(1:14)]),
                                                                Type=str_split_fixed(colnames(count_removing_zeros)[-c(1:14)],"_",n=3)[,3],
                                                                Replicate=str_split_fixed(colnames(count_removing_zeros)[-c(1:14)],"_",n=3)[,2],
                                                                infectionRate=scales::rescale(c(rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000),rep(0,6),c(219342105,193750000,37179,1781250,62500,15250000)))))

y <- edgeR::DGEList(counts=assay(data),group = data$Condition,samples =colData(data),genes = rowData(data))

####Metatranscriptome####
####Filtering MG_cpm_ge32_3samples####
x_rna <- edgeR::DGEList(counts=assay(data[,13:24]),group = data$Condition[13:24],samples =colData(data[,13:24]),
                        genes = rowData(data))
CPM <- cpm(x_rna,lib.size=colSums(x_rna$counts))
CPM.Cutoff <- 32

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

x_rna_ruvg <- calcNormFactors(x_rna, method=normalize_method)
x_rna_ruvg <- estimateGLMCommonDisp(x_rna_ruvg, design_ruvg)
x_rna_ruvg <- estimateGLMTagwiseDisp(x_rna_ruvg, design_ruvg)

fit <- glmFit(x_rna_ruvg, design_ruvg)
lrt <- glmLRT(fit, contrast=all_contrasts)

top <- topTags(lrt, n=nrow(x_rna))$table
empirical <- top$associated_OG[which(top$PValue>0.8 & top$FDR>0.1)]
####empirical RUVg####
a <- floor(x_rna_ruvg$counts)
row.names(a) <- x_rna_ruvg$genes$associated_OG
set <- newSeqExpressionSet(a,phenoData = data.frame(x_rna_ruvg$samples$group, row.names=colnames(x_rna_ruvg$counts)))
set2 <- RUVg(set, empirical, k=i)
####Plotting counts####
cairo_pdf(file="rna_logCPM_distribution_pca_rle_afterRUVg.pdf",bg="transparent",
          width = 11,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 5))
a_logCPM <- data.frame(cbind("associated_OG" =rownames(a),edgeR::cpm.default(a,log=T)), stringsAsFactors = FALSE)
m_a <- reshape2::melt(a_logCPM,value.name = "log_cpm",id.vars=c("associated_OG"))
m_a$Condition <- gsub("_\\d+","",m_a$variable)
m_a$Tpe=str_split_fixed(m_a$variable,"_",n=3)[,3]
m_a$Replicate=str_split_fixed(m_a$variable,"_",n=3)[,2]
m_a$log_cpm <- as.numeric(m_a$log_cpm)

vioplot::vioplot(m_a$log_cpm ~ m_a$variable, col=c(rep("blue",6),rep("red",6)), cex=0.4,
                 main=paste0("Log CPM violin plot\nN= ",nrow(a)), xlab="Replicates",ylab="Log CPM",xaxt = "n")

normalisedCounts <- data.frame(cbind("associated_OG" =rownames(normCounts(set2)),edgeR::cpm.default(normCounts(set2),log=T)), stringsAsFactors = FALSE)
m_normalisedCounts <- reshape2::melt(normalisedCounts,value.name = "normalised_count",id.vars=c("associated_OG"))
m_normalisedCounts$Condition <- gsub("_\\d+","",m_normalisedCounts$variable)
m_normalisedCounts$Tpe=str_split_fixed(m_normalisedCounts$variable,"_",n=3)[,3]
m_normalisedCounts$Replicate=str_split_fixed(m_normalisedCounts$variable,"_",n=3)[,2]
m_normalisedCounts$normalised_count <- as.numeric(m_normalisedCounts$normalised_count)

vioplot::vioplot(m_normalisedCounts$normalised_count ~ m_normalisedCounts$variable, col=c(rep("blue",6),rep("red",6)), cex=0.4,
                 main=paste0("Normalised Count from RUVg\nLog CPM violin plot k=",i,"\nN= ",nrow(normalisedCounts)), xlab="Replicates",ylab="Log CPM",xaxt = "n")

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=c(rep("blue",6),rep("red",6)),main=paste0("RLE plot k=",i),xaxt = "n", xlab="Replicates",ylab="RLE")

plotPCA(set2, col=c(rep("blue",6),rep("red",6)), cex=0.4,main=paste0("PCA plot k=",i))

corrected <- removeBatchEffect(cpm(x_rna, log=TRUE), design=design_ruvg, covariates=as.matrix(pData(set2)[,c("W_1","W_2")]))
corrected <- data.frame(cbind("associated_OG" =x_rna$genes$associated_OG,corrected), stringsAsFactors = FALSE)
m_corrected <- reshape2::melt(corrected,value.name = "corrected_logCPM",id.vars=c("associated_OG"))
m_corrected$Condition <- gsub("_\\d+","",m_corrected$variable)
m_corrected$Tpe=str_split_fixed(m_corrected$variable,"_",n=3)[,3]
m_corrected$Replicate=str_split_fixed(m_corrected$variable,"_",n=3)[,2]
m_corrected$corrected_logCPM <- as.numeric(m_corrected$corrected_logCPM)
vioplot::vioplot(m_corrected$corrected_logCPM ~ m_corrected$variable, col=c(rep("blue",6),rep("red",6)), cex=0.4,
                 main=paste0("Corrected counts after RUVg\nLog CPM violin plot k=",i,"\nN= ",nrow(corrected)), xlab="Replicates",ylab="Log CPM",xaxt = "n")
dev.off()
####Calculating normalising factors####
x_rna <- calcNormFactors(x_rna,method = normalize_method)
####edgeR####
cairo_pdf(file="rna_pvalue_fdr_edgeR_afterRUVg.pdf",bg="transparent",
          width = 10,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 4),bg="transparent",family="Times New Roman")
design_edgeR <- model.matrix(~x_rna_ruvg.samples.group + W_1 + W_2, data=pData(set2))
design_edgeR[,"x_rna_ruvg.samples.groupUninfected_RNA"] <- c(rep(0,6),rep(1,6))
y_rna <- DGEList(counts=counts(set2), group=x_rna$samples$group,genes = x_rna$genes)
y_rna <- calcNormFactors(y_rna, method=normalize_method)
y_rna <- estimateDisp(y_rna, design_edgeR)

plotBCV(y_rna, xlab="Average log CPM", ylab="Biological coefficient of variation",
        pch=16, cex=0.2, col.common="red", col.trend="blue", col.tagwise="black",
        main = paste0("BCV plot\nN=" , nrow(x_rna$counts)))
fit_rna <- glmFit(y_rna, design_edgeR)
lrt_rna <- glmLRT(fit_rna, coef=2)
top_rna <- topTags(lrt_rna, n=nrow(set))$table
hist(top_rna$PValue,breaks=20,main = paste0("Histogram of P-value\nN=" , nrow(y_rna$counts)),
     xlab = "P-value")
hist(top_rna$FDR,breaks=20,main = paste0("Histogram of FDR\nN=" , nrow(y_rna$counts)),
     xlab = "FDR")
plot(top_rna$logCPM, top$logFC,main = paste0("MA plot\nN=" , nrow(y_rna$counts)),
     xlab="Average log CPM", ylab="logFC",pch=20,cex=0.5)
points(subset(top_rna$logCPM,top_rna$PValue<=0.05 & top_rna$FDR <= 0.05 & (top_rna$logFC>=2|top_rna$logFC<=-2)),
       subset(top_rna$logFC,top_rna$PValue<=0.05 & top_rna$FDR <= 0.05 & (top_rna$logFC>=2|top_rna$logFC<=-2)),
       pch=20,cex=0.5,col="red")
dev.off()

top_rna <- left_join(top_rna,a_logCPM,by="associated_OG")
top_rna <- top_rna %>% mutate_at(.vars = colnames(a_logCPM)[-1],as.numeric)
write.table(top_rna,"rna_ruvg_edgeR_da_results.tsv",sep = "\t",row.names = F)

####Fry GO####
h_path <- read.table("/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/mapping_eggnog_OGs_GOs_forFry.tsv", sep="\t", stringsAsFactors=F, header=T, quote="\"")
colnames(h_path) <- c("path_id","gene_name","Species","Source","Pathway_description")

get_h_path_indices_forFry <- function(x, gene_names){
  path_list <- list()
  for (i in as.character(unique(x$path_id))){
    if(length(which(gene_names %in% x[x$path_id == i, 2])) < 100 && length(which(gene_names %in% x[x$path_id == i, 2])) > 15){
      path_list[[i]] = gene_names[gene_names %in% x[x$path_id == i, 2]]
    }
  }
  return(path_list)
}
path_indices_forFry <- get_h_path_indices_forFry(h_path, y_rna$genes$associated_OG)

fry_rna <- fry(y_rna, path_indices_forFry, design=design_edgeR,
               geneid=y_rna$genes$associated_OG,contrast=2)

reformat_h_path <- function(path_list, h_path){
  path_names <- c()
  source <- c()
  geneNames <- c()
  for (i in rownames(path_list)){
    index <- which(h_path$path_id == i)[1]
    path_names <- c(path_names, h_path[index, 5])
    source <- c(source, h_path[index, 4])
    index <- which(h_path$path_id == i)
    geneNames <- c(geneNames, paste(h_path[index,2],sep="",collapse=";"))
  }
  new_list <- cbind(path_list, "path_name"=path_names, "source"=source, "gene_name"=geneNames)
  return(new_list)
}
f_rna <- reformat_h_path(fry_rna, h_path)

significant_gene_list <- function(reformatted_fry_results,lmfit_de){
  nSigGene <- c()
  sig_geneExpr <- c()
  for(i in 1:nrow(reformatted_fry_results)){
    gene_names <- strsplit(as.character(reformatted_fry_results$gene_name[i]),";")[[1]]
    gnames <- subset(lmfit_de$associated_OG,lmfit_de$associated_OG %in% gene_names & (lmfit_de$logFC >= 2 | lmfit_de$logFC <= -2) & lmfit_de$FDR <= 0.05 & lmfit_de$PValue <= 0.05)
    nSigGene <- c(nSigGene,length(gnames))
    sig_geneExpr <- c(sig_geneExpr, paste(gnames,sep="",collapse=";"))
  }
  new_list <- cbind(reformatted_fry_results,"NSignificantGenes"=nSigGene,"NameSignificantGenes"=sig_geneExpr)
  return(new_list)
}

sig_f_rna <- significant_gene_list(f_rna,top_rna)
write.table(rownames_to_column(sig_f_rna,var = "GO"),file="rna_fry_results_withGO.tsv",row.names = F,quote = T,sep="\t")

cairo_pdf(file="rna_pvalue_fdr_fryGO.pdf",bg="transparent",
          width = 10,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 4),bg="transparent",family="Times New Roman")
hist(sig_f_rna$PValue,breaks=20,main = paste0("Histogram of P-value\nN=" , nrow(y_rna$counts)),
     xlab = "P-value")
hist(sig_f_rna$FDR,breaks=20,main = paste0("Histogram of FDR\nN=" , nrow(y_rna$counts)),
     xlab = "FDR")
hist(sig_f_rna$PValue.Mixed,breaks=20,main = paste0("Histogram of non-directional P-value\nN=" , nrow(y_rna$counts)),
     xlab = "P-value")
hist(sig_f_rna$FDR.Mixed,breaks=20,main = paste0("Histogram of non-directional FDR\nN=" , nrow(y_rna$counts)),
     xlab = "FDR")
dev.off()

sig_f_rna <- sig_f_rna %>% mutate(NGenes = ifelse(Direction == "Down", -NGenes, NGenes))
# Create the bar plot
pathway_bar_plot <- ggplot(subset(sig_f_rna,!is.na(sig_f_rna$path_name)), aes(x = NGenes, y = reorder(path_name, NGenes))) +
  geom_bar(stat = "identity", aes(color = PValue, fill=FDR)) + 
  labs(x = "Number of OGs", y = "Pathway", title = "GO terms (after FRY)") +
  scale_colour_gradient2(trans = "sqrt", name = "PValue",low = "brown", mid = "white", high = "grey", midpoint = .25) +
  scale_fill_gradient2(trans = "sqrt", name = "FDR",low = "brown", mid = "white", high = "grey", midpoint = .25) +
  theme_minimal(base_size = 12,base_family = "Times New Roman") + 
  theme(axis.text.y = element_text(size=3,colour = "black"),
        axis.text.x = element_text(colour = "black"),
        legend.position = "right",panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent',colour = NA), #transparent legend bg
        legend.box.background = element_rect(fill='transparent',colour = NA) #transparent legend panel
  )
cairo_pdf("rna_all_GO_pathway_bar_plot.pdf",width=6,height=10,family = "Times New Roman",bg = "transparent")
print(pathway_bar_plot)
dev.off()

sig_pathways <- subset(sig_f_rna,sig_f_rna$PValue <= 0.05 & sig_f_rna$FDR <= 0.05 & 
                         sig_f_rna$NSignificantGenes > 0 & !is.na(sig_f_rna$path_name))
if(nrow(sig_pathways)==0){
  sig_pathways <- subset(sig_f_rna,sig_f_rna$PValue.Mixed <= 0.05 & sig_f_rna$FDR.Mixed <= 0.05 & 
                           sig_f_rna$NSignificantGenes > 0 & !is.na(sig_f_rna$path_name))
}
# Create the bar plot
sig_pathway_bar_plot <- ggplot(sig_pathways, aes(x = NGenes, y = reorder(path_name, NGenes))) +
  geom_bar(stat = "identity", aes(fill = Direction)) +
  labs(x = "Number of OGs", y = "Pathway", title = "Significant GO terms (after FRY)") +
  scale_fill_manual(values = c("Up" = "blue", "Down" = "red")) +
  theme_minimal(base_size = 12,base_family = "Times New Roman") + 
  theme(axis.text.y = element_text(size = 5,colour = "black"),
        axis.text.x = element_text(colour = "black"),
        legend.position = "right",panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent',colour = NA), #transparent legend bg
        legend.box.background = element_rect(fill='transparent',colour = NA) #transparent legend panel
  )
cairo_pdf("rna_sig_GO_pathway_bar_plot.pdf",width=6,height=6,family = "Times New Roman",bg = "transparent")
print(sig_pathway_bar_plot)
dev.off()

####FRY KEGG####
h_path <- read.table("/home/shuba/Documents/microbiome/wildtype_MG_MT_nextflow/newFiles/mapping_eggnog_OGs_KEGG_pathway_forFry_noDescriptions_removedRepeats.tsv", sep="\t", stringsAsFactors=F, header=T, quote="\"")
colnames(h_path) <- c("path_id","gene_name","Species","Source","Pathway_description")

path_indices_forFry <- get_h_path_indices_forFry(h_path, y_rna$genes$associated_OG)

fry_rna <- fry(y_rna, path_indices_forFry, design=design_edgeR,
               geneid=y_rna$genes$associated_OG,contrast=2)

f_rna <- reformat_h_path(fry_rna, h_path)

if(length(which(f_rna$source == "KEGG_Pathway"))!=0){
  row_nos <- which(f_rna$source == "KEGG_Pathway")
  f_rna$path_name[row_nos] <- sapply(row_nos,function(x){
    tryCatch(KEGGREST::keggGet(rownames(f_rna)[x])[[1]]$NAME,error=function(e){NA})
  })
}
#ko00471= "D-Glutamine and D-glutamate metabolism"

sig_f_rna <- significant_gene_list(f_rna,top_rna)
sig_f_rna$path_name[which(rownames(sig_f_rna)=="ko00471")] <- "D-Glutamine and D-glutamate metabolism"
write.table(rownames_to_column(sig_f_rna,var = "KEGG"),file="rna_fry_results_withKEGG.tsv",row.names = F,quote = T,sep="\t")

cairo_pdf(file="rna_pvalue_fdr_fryKEGG.pdf",bg="transparent",
          width = 10,height = 2.5,family = "Times New Roman")
par(mfrow = c(1, 4),bg="transparent",family="Times New Roman")
hist(sig_f_rna$PValue,breaks=20,main = paste0("Histogram of P-value\nN=" , nrow(y_rna$counts)),
     xlab = "P-value")
hist(sig_f_rna$FDR,breaks=20,main = paste0("Histogram of FDR\nN=" , nrow(y_rna$counts)),
     xlab = "FDR")
hist(sig_f_rna$PValue.Mixed,breaks=20,main = paste0("Histogram of non-directional P-value\nN=" , nrow(y_rna$counts)),
     xlab = "P-value")
hist(sig_f_rna$FDR.Mixed,breaks=20,main = paste0("Histogram of non-directional FDR\nN=" , nrow(y_rna$counts)),
     xlab = "FDR")
dev.off()

sig_f_rna <- sig_f_rna %>%
  mutate(NGenes = ifelse(Direction == "Down", -NGenes, NGenes))
# Create the bar plot
pathway_bar_plot <- ggplot(subset(sig_f_rna,!is.na(sig_f_rna$path_name)), aes(x = NGenes, y = reorder(path_name, NGenes))) +
  geom_bar(stat = "identity", aes(color = PValue, fill=FDR)) + 
  labs(x = "Number of OGs", y = "Pathway", title = "KEGG pathways (after FRY)") +
  scale_colour_gradient2(trans = "sqrt", name = "PValue",low = "brown", mid = "white", high = "grey", midpoint = .25) +
  scale_fill_gradient2(trans = "sqrt", name = "FDR",low = "brown", mid = "white", high = "grey", midpoint = .25) +
  theme_minimal(base_size = 12,base_family = "Times New Roman") + 
  theme(axis.text.y = element_text(size=5,colour = "black"),
        axis.text.x = element_text(colour = "black"),
        legend.position = "right",panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent',colour = NA), #transparent legend bg
        legend.box.background = element_rect(fill='transparent',colour = NA) #transparent legend panel
  )
cairo_pdf("rna_all_KEGG_pathway_bar_plot.pdf",width=6,height=10,family = "Times New Roman",bg = "transparent")
print(pathway_bar_plot)
dev.off()

sig_pathways <- subset(sig_f_rna,sig_f_rna$PValue <= 0.05 & sig_f_rna$FDR <= 0.05 &
                         sig_f_rna$NSignificantGenes > 0 & !is.na(sig_f_rna$path_name))
# Create the bar plot
sig_pathway_bar_plot <- ggplot(sig_pathways, aes(x = NGenes, y = reorder(path_name, NGenes))) +
  geom_bar(stat = "identity", aes(fill = Direction)) +
  labs(x = "Number of OGs", y = "Pathway", title = "Significant KEGG pathways (after FRY)") +
  scale_fill_manual(values = c("Up" = "blue", "Down" = "red")) +
  theme_minimal(base_size = 12,base_family = "Times New Roman") + 
  theme(axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        legend.position = "right",panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent',colour = NA), #transparent legend bg
        legend.box.background = element_rect(fill='transparent',colour = NA) #transparent legend panel
  )
cairo_pdf("rna_sig_KEGG_pathway_bar_plot.pdf",width=6,height=6,family = "Times New Roman",bg = "transparent")
print(sig_pathway_bar_plot)
dev.off()


####Heatmap Figure 7####
suppressMessages(library(ComplexHeatmap))

##Significant rows rna
rna_limma_daa_sig <- subset(top_rna,top_rna$PValue<=0.05 & top_rna$FDR <= 0.05 & (top_rna$logFC>=2|top_rna$logFC<=-2))
# Split the "Class" column into separate rows
df <- top_rna %>% separate_rows(Class, sep = ";")
# palette <- randomcoloR::distinctColorPalette(length(unique(c(df$Class,df$Class))))
# names(palette) <- as.character(unique(c(df$Class,df$Class)))
# palette <- c("c__Gammaproteobacteria"= "#7D7CD4",
#              "c__Clostridia" = "#D5CDD1",
#              "c__Bacteroidia" = "#D49FD5",
#              "c__Bacilli" = "#DCE15A",
#              "c__Bacilli_A" = "#DCE15A",
#              "c__Actinobacteria" = "#AE4BE3",
#              "c__Verrucomicrobiae" = "#7EBDDB",
#              "c__Coriobacteriia" = "tomato1",
#              "c__Alphaproteobacteria"="#D1D098",
#              "c__Bradimonadia"="#7BDBC1",
#              "c__Campylobacteria"="#D98471",
#              "c__Deinococci"="#7CE574",
#              "c__Polyangia"="#D61600",
#              "c__Gemmatimonadetes"="#FFF400",
#              "c__Myxococcia" = "royalblue1")
palette <- c("c__Clostridia" = "#DF7C5F","c__Bacteroidia" = "#DE77BD", "c__Verrucomicrobiae" = "#89E0A6", 
             "c__Bacilli" = "#DCD668", "c__Actinobacteria" = "#83D3DA","c__Coriobacteriia" = "#8AE85F", 
             "c__Gammaproteobacteria" = "#7987D9", "c__Alphaproteobacteria" = "#D7D0B1", "c__Bacilli_A" = "#B047DA", 
             "c__Myxococcia" = "#C5B1CF","c__Polyangia" = "#D61600","c__Methanococci" = "#D98471",
             "c__Thermoleophilia" = "#FFF400")

####Volcano plot####
rna_class <- ggplot(df, aes(x = logCPM, y = logFC,color=Class, fill=Class)) +
  ggrastr::geom_point_rast(size=1,alpha=0.2,raster.dpi = 150) + 
  theme_bw(base_size = 12,base_family = "Times New Roman") + theme(legend.position = "right") +
  scale_color_manual(values = palette) + theme(axis.text.y = element_text(colour = "black"),
                                               axis.text.x = element_text(colour = "black"),
                                               legend.position = "right",panel.background = element_rect(fill='transparent'), #transparent panel bg
                                               plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                                               # panel.grid.major = element_blank(), #remove major gridlines
                                               panel.grid.minor = element_blank(), #remove minor gridlines
                                               legend.background = element_rect(fill='transparent',colour = NA), #transparent legend bg
                                               legend.box.background = element_rect(fill='transparent',colour = NA) #transparent legend panel
  )

# Create the scatter plot with color scale based on pvalue and fdr
rna_pval_fdr <- ggplot(df, aes(x = logCPM, y = logFC, color = PValue, fill = FDR)) +
  ggrastr::geom_point_rast(size=1,alpha=0.2,raster.dpi = 150) +  
  scale_colour_gradient2(trans = "sqrt", name = "PValue",low = "brown", mid = "white", high = "grey", midpoint = .25) +
  scale_fill_gradient2(trans = "sqrt", name = "FDR",low = "brown", mid = "white", high = "grey", midpoint = .25) +
  theme_bw(base_size = 12,base_family = "Times New Roman") + theme(legend.position = "right") +
  theme(axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        legend.position = "right",panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent',colour = NA), #transparent legend bg
        legend.box.background = element_rect(fill='transparent',colour = NA) #transparent legend panel
  )

cairo_pdf(filename = "rna_og_ruvg_edgeR_MAplot.pdf",width = 5,height = 3.75,
          family = "Times New Roman",onefile = T)
print(rna_class)
print(rna_pval_fdr)
dev.off()
