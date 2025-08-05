# install.packages("BiocManager", dependencies=TRUE)
# BiocManager::install("DESeq2")

#step1: uploading  libraries
library("DESeq2")
library("tidyverse")
library("kableExtra")

#step2: Preparing data
data=read.table('/home/vaishnavi-goyal/rna_analysis/R-yeast/counts.txt', sep='\t', header=TRUE)
# getwd()  for curent directory
# setwd("/home/vaishnavi-goyal/rna_analysis/R-yeast") for setting directory

data=as.data.frame(data)
row.names(data)=data$Geneid
data=data[,-c(1:6)]
head(data)
colnames(data)
columns= c('batch1', 'batch2', 'batch3', 'chem1', 'chem2', 'chem3')
colnames(data)=columns
colnames(data)
data=data[,sort(colnames(data))]
head(data)
colData=data.frame(row.names=colnames(data),group=c(rep("batch",3),rep("chemostat",3)))
colData

#step3: creating DESeq2 Object
dds = DESeqDataSetFromMatrix(countData = data,
                             colData = colData,
                             design = ~group)
dds = DESeq(dds)
nrow(dds)
dds = dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)
dds

# BiocManager::install("biomaRt")
library(biomart)

# Connect to Ensembl for S. cerevisiae
ensembl <- useMart("ensembl", dataset = "scerevisiae_gene_ensembl")


# Get gene annotation information
yeast_genes <- getBM(
  attributes = c(
    "ensembl_gene_id", 
    "external_gene_name", 
    "chromosome_name", 
    "start_position", 
    "end_position", 
    "strand", 
    "gene_biotype"
  ),
  mart = ensembl
)

# Save to CSV file
write.csv(yeast_genes, "saccharomyces_cerevisiae_gene_annotation.csv", row.names = FALSE)

#step4: add annotation file
annotation <- read.csv("/home/vaishnavi-goyal/rna_analysis/R-yeast/saccharomyces_cerevisiae_gene_annotation.csv", header = T, stringsAsFactors = F)

#install lib for visualisation
library("RcolorBrewer")
library("grDevices")
# library("pheatmap")
# library("ggrepel")
# library("edgeR")

# Step5: Visualisation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
plotDists=function(vsd.obj){
  sampleDists <-dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix)<- paste(vsd.obj$group)
  colour <- colorRampPalette(rev(brewer.pal(9,"Blues")))(225)
  pheatmap::pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colour)
}
plotDists(vsd)


#Variable Genes Heatmap
variable_gene_heatmap <- function (vsd.obj, num_genes = 500, annotation, title = "") {
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette( RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  stabilized_counts <- assay(vsd.obj)
  row_variances <- rowVars(stabilized_counts)
  top_variable_genes <- stabilized_counts[order(row_variances, decreasing=T)[1:num_genes],]
  top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm=T)
  gene_names <- annotation$Gene.name[match(rownames(top_variable_genes), annotation$Gene.stable.ID)]
  rownames(top_variable_genes) <- gene_names
  coldata <- as.data.frame(vsd.obj@colData)
  coldata$sizeFactor <- NULL
  pheatmap::pheatmap(top_variable_genes, color = mr, annotation_col = coldata, fontsize_col = 8, fontsize_row = 250/num_genes, border_color = NA, main = title)
}

variable_gene_heatmap(vsd, num_genes = 40, annotation = annotation)


# pca graph
plot_PCA = function (vsd.obj) {
  pcaData <- plotPCA(vsd.obj,  intgroup = c("group"), returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=group)) +
    geom_point(size=3) +
    labs(x = paste0("PC1: ",percentVar[1],"% variance"),
         y = paste0("PC2: ",percentVar[2],"% variance"),
         title = "PCA Plot colored by group") +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
}
plot_PCA(vsd)


# Step6: Differential gene analysis 
dds <- DESeq(dds)
res <- results(dds)
kable(head(res)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")

res <- results(dds, alpha = 0.5, lfcThreshold=0.01)
summary(res)

# sort genes by fold change
res <- res[order(abs( res$log2FoldChange), decreasing=TRUE),]
kable(head(res)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")

# MA plot
DESeq2::plotMA(res,  ylim = c(-5, 5))

# volcano plot
res1 = as.data.frame(res)
res1 = mutate(res1, sig=ifelse(res1$padj<0.05, "FDR<0.05", "Not Sig"))
res1[which(abs(res1$log2FoldChange)<0.5),'sig'] <- "Not Sig"

p = ggplot(res1, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black"))
p

# Step7: Differtial Gene Annotation and Summary ----
generate_DE_results <- function (dds, comparisons, padjcutoff = 0.001, log2cutoff = 0.5, cpmcutoff = 2) {
  raw_counts <- counts(dds, normalized = F)
  cpms <- enframe(rowMeans(edgeR::cpm(raw_counts)))
  colnames(cpms) <- c("ensembl_id", "avg_cpm")
  res <- results(dds)
  res <- as_tibble(res, rownames = "ensembl_id")
  my_annotation <- read.csv("/home/vaishnavi-goyal/rna_analysis/R-yeast/saccharomyces_cerevisiae_gene_annotation.csv", header = T, stringsAsFactors = F)
  res <- left_join(res, my_annotation, by = c("ensembl_id" = "ensembl_gene_id"))
  res <- left_join(res, cpms, by = c("ensembl_id" = "ensembl_id"))
  normalized_counts <- round(counts(dds, normalized = TRUE),3)
  pattern <- str_c(comparisons[1], "|", comparisons[2])
  combined_data <- as_tibble(cbind(res, normalized_counts[,grep(pattern, colnames(normalized_counts))] ))
  combined_data <- combined_data[order(combined_data$log2FoldChange, decreasing = T),]
  res_prot <- res[which(res$Gene.type == "protein_coding"),]
  res_prot_ranked <- res_prot[order(res_prot$log2FoldChange, decreasing = T),c("Gene.name"="external_gene_name", "log2FoldChange")]
  res_prot_ranked <- na.omit(res_prot_ranked)
  res_prot_ranked$Gene.name <- str_to_upper(res_prot_ranked$Gene.name)
  res <- res[order(res$log2FoldChange, decreasing=TRUE ),]
  de_genes_padj <- res[which(res$padj < padjcutoff),]
  de_genes_log2f <- res[which(abs(res$log2FoldChange) > log2cutoff & res$padj < padjcutoff),]
  de_genes_cpm <- res[which(res$avg_cpm > cpmcutoff & res$padj < padjcutoff),]
  write.csv (combined_data, file = paste0(comparisons[1], "_vs_", comparisons[2], "_allgenes.csv"), row.names =F)
  writeLines( paste0("For the comparison: ", comparisons[1], "_vs_", comparisons[2], ", out of ", nrow(combined_data), " genes, there were: \n", 
                     nrow(de_genes_padj), " genes below padj ", padjcutoff, "\n",
                     nrow(de_genes_log2f), " genes below padj ", padjcutoff, " and above a log2FoldChange of ", log2cutoff, "\n",
                     nrow(de_genes_cpm), " genes below padj ", padjcutoff, " and above an avg cpm of ", cpmcutoff, "\n",
                     "Gene lists ordered by log2fchange with the cutoffs above have been generated.") )
  gene_count <- tibble (cutoff_parameter = c("padj", "log2fc", "avg_cpm" ), 
                        cutoff_value = c(padjcutoff, log2cutoff, cpmcutoff), 
                        signif_genes = c(nrow(de_genes_padj), nrow(de_genes_log2f), nrow(de_genes_cpm)))
  invisible(gene_count)
}
dds_result_summary= generate_DE_results (dds, c("batch","chemostat"))
# Step8: Further analysis ----
res = read.csv("/home/vaishnavi-goyal/rna_analysis/R-yeast/batch_vs_chemostat_allgenes.csv", header = T)
head(res)

# get the significant genes
resSig = as.data.frame(subset(res,padj<0.5) )
resSig = resSig[order(resSig$log2FoldChange,decreasing=TRUE),]
head(resSig)

# save file with significant genes
write.csv(resSig,"yeast_SigGenes.csv")

# volcano plot 
plot_volcano <- function (res, padj_cutoff, nlabel = 10, label.by = "padj"){
  res <- mutate(res, significance=ifelse(res$padj<padj_cutoff, paste0("padj < ", padj_cutoff), paste0("padj > ", padj_cutoff)))
  res = res[!is.na(res$significance),]
  significant_genes <- res %>% filter(significance == paste0("padj < ", padj_cutoff))
  if (label.by == "padj") {
    top_genes <- significant_genes %>% arrange(padj) %>% head(nlabel)
    bottom_genes <- significant_genes %>% filter (log2FoldChange < 0) %>% arrange(padj) %>% head (nlabel)
  } else if (label.by == "log2FoldChange") {
    top_genes <- head(arrange(significant_genes, desc(log2FoldChange)),nlabel)
    bottom_genes <- head(arrange(significant_genes, log2FoldChange),nlabel)
  } else
    stop ("Invalid label.by argument. Choose either padj or log2FoldChange.")
  
  ggplot(res, aes(log2FoldChange, -log(padj))) +
    geom_point(aes(col=significance)) + 
    scale_color_manual(values=c("red", "black")) + 
    ggrepel::geom_text_repel(data=top_genes, aes(label=head(external_gene_name,nlabel)), size = 3)+
    ggrepel::geom_text_repel(data=bottom_genes, aes(label=head(external_gene_name,nlabel)), color = "#619CFF", size = 3)+
    labs ( x = "Log2FoldChange", y = "-(Log normalized p-value)")+
    geom_vline(xintercept = 0, linetype = "dotted")+
    theme_minimal()
}

plot_volcano(res, 0.0005, nlabel = 15, label.by = "padj")
