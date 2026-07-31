# last updated: 12/25/2025
# load libraries
require(biomaRt)
library(devtools)
library(edgeR)
library(limma)
library(Glimma)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(msigdbr)
library(WriteXLS)
library(scales)
library(wesanderson)
library(gdata)
library(ggpubr)
library(factoextra)
library(tidyr)
library(readxl)
library(pheatmap)
library(openxlsx)
library(dplyr)
library(magrittr)
library(DT)
library(org.Ce.eg.db)
library(gplots)
library(rstudioapi)
library(stringr)
library(data.table)
library(RRHO)
library(gridExtra)
library(dendextend)
library(GO.db)
library(tidyverse)
library(ggrepel)
library(kableExtra)
library(knitr)

# create raw counts matrix from STAR output
# use count_result_star folder as directory for download_raw_reads

Download_raw_reads <- function(featurecounts_dir){
  setwd(featurecounts_dir)
  
  temp_data <- read.csv(dir()[!grepl(".summary$",dir())][1],header=T,sep='\t',skip = 1)
  counts_star <- data.frame(ID=temp_data$Geneid)
  rownames(counts_star) <- counts_star$ID
  for (i in dir()[!grepl(".summary$",dir())]){
    temp_name <- strsplit(i,".count")[[1]][1]
    temp_name <- gsub("-","_",temp_name)
    temp_data <- read.csv(i,header=T,sep='\t',skip = 1)
    temp_counts <- temp_data[,7]
    counts_star[temp_name] <- temp_counts
  }
  rm(temp_data)
  counts_star$ID <- rownames(counts_star)
  counts_star <- counts_star[,-1]
  return(counts_star)
}

count_matrix <- Download_raw_reads('/home/wmitchell/greer_RNAseq_celegans/count_result_star')

write.csv(as.data.frame(counts_matrix), 'raw_genecounts.csv')

########### excluding transcripts with low expression ##########

keep <- edgeR::filterByExpr(count_matrix, design = NULL)
counts.keep <- count_matrix[keep,]

write.csv(as.data.frame(counts.keep), 'filtered_genecounts.csv')

# normalize input counts with edgeR

dgeObj <- DGEList(counts.keep)
barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2)

png('boxplot_rawcounts.png', width = 6*400, height= 8*400, pointsize = 8, res = 400)
logcounts <- cpm(dgeObj,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million")
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalized)")
dev.off()

dgeObj <- calcNormFactors(dgeObj, method = 'RLE')
logcounts <- cpm(dgeObj,log=T)

boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (normalized)")

write.csv(as.data.frame(logcounts), 'normalized_genecounts_hsp1_inputs.csv')

############ import counts and perform PCA #####################################################
# set wd
setwd("/home/wmitchell/greer_new/data")
# import filtered normalized genecounts (filtered IP counts divided by input counts)

IP <- read.csv("filtered_IP.csv", header = TRUE, row.names = 1)
input <- read.csv("filtered_input.csv", header = TRUE, row.names = 1)
te <- log2((IP+1)/(input+1))
te <- scale(te)
write.csv(as.data.frame(te), 'te_hsp1.csv')

# read metadata
meta <- read.csv("meta_hsp1_IP.csv")

# pca
pca <- prcomp(t(te))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PCAsummary <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = meta$gene) 

g <- ggplot(data = PCAsummary, aes_string(x = "PC1", y = "PC2", color = "group")) + 
  geom_point(size = 3) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linetype = "solid", fill = NA, size = 0.75), plot.margin = margin(1,1,1,1,"mm")) +
  theme(axis.ticks.length = unit(1, "mm"), axis.ticks = element_line(size = 0.5), axis.text = element_text(size = rel(0.75))) +
  theme(axis.title = element_text(size = rel(0.75)))+
  theme(legend.text = element_text(size = rel(0.65)), legend.title = element_blank(), legend.key.size = unit(15, "pt"), legend.margin = margin(0,0,0,0,"pt"), legend.box.margin=margin(0,0,0,-10, "pt"))

ggsave(filename="PCA_IP_hsp1.pdf", plot=g, device="pdf", units="in", width=4, height=4)
g

############## differential expression analysis #####################################
design <- model.matrix(~factor(meta$gene))
colnames(te[,colnames(te)%in% meta$sample])==meta$sample

#### prepare DEGs
fit <- lmFit(te, design, method = 'robust') #input method = robust
ebayes <- eBayes(fit)
dt_fdr <- decideTests(fit)
summary(dt_fdr)

tab <- topTable(ebayes, coef=2, adjust="fdr", n=nrow(counts))
tab$gene = rownames(tab)

####### check the sign of diff expression 
gene1 = tab$gene[tab$logFC>0][1]
counts[grep(gene1, row.names(counts)),]

#### if the sign is reversed, run this
tab$logFC = tab$logFC * (-1)

write.csv(as.data.frame(tab), paste('DEGs_hsp1_IP','.csv', sep = ''), row.names = T)

################# heatmap of differentially expressed genes (FDR < 0.05) ############################
colors <- colorRampPalette(c("blue","black","yellow"))(256)
plotData <- te

# keep only genes that are significantly DE by FDR < 0.05
DE_genes <- tab[tab$adj.P.Val <= 0.05,] %>% row.names()

plotData <- plotData[DE_genes,] # %>% scale() only use for inputs
distCol <- t(plotData) %>% dist()
hclustCol <- hclust(distCol, method = 'complete') %>% rotate(., order = colnames(plotData)) 


p <- pheatmap(plotData, clustering_distance_rows = 'correlation', clustering_method = "average", scale = 'row', cluster_cols = hclustCol,
              color = colors,  border_color = NA, show_rownames = FALSE, fontsize = 8, silent = TRUE
              
)

png(filename="IP_heatmap_hsp1.png", res = 300, height = 4, width = 4, units = 'in')
p
invisible(dev.off())

################# GO enrichment analysis #########################################
tab = read.csv("DEGs_hsp1_IP.csv", header = T, row.names = 1)
tab$rank = ifelse(tab$logFC>0, 1, -1)*(-log10(tab$P.Value))
tab = tab[order(tab$rank, decreasing = T), ]
geneList = tab$rank
names(geneList)= tab$gene
head(geneList)

gsea_cc <- gseGO(geneList = geneList, OrgDb = org.Ce.eg.db, keyType = "WORMBASE", ont = "CC", exponent = 0, eps = 1e-100, minGSSize = 20, maxGSSize = 500, pvalueCutoff = 0.05, verbose  = FALSE)
gsea_bp <- gseGO(geneList = geneList, OrgDb = org.Ce.eg.db, keyType = "WORMBASE", ont = "BP", exponent = 0, eps = 1e-100, minGSSize = 20, maxGSSize = 500, pvalueCutoff = 0.05, verbose  = FALSE)
gsea_mf <- gseGO(geneList = geneList, OrgDb = org.Ce.eg.db, keyType = "WORMBASE", ont = "MF", exponent = 0, eps = 1e-100, minGSSize = 20, maxGSSize = 500, pvalueCutoff = 0.05, verbose  = FALSE)
write.csv(gsea_mf@result,  paste('GSEA_GOmf_hsp1_IP_redone','.csv', sep = ''))

png('GO_hsp1_IP.png')
dotplot(gse_GO, showCategory = 15, split = '.sign', font.size = 12)+ 
  facet_grid(.~.sign)+ 
  theme(strip.text.x = element_text(size = 16), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))+
  theme_bw()
dev.off()

####### convert WORMBASE IDs to SYMBOL #########################
eg = bitr(unique(tab$gene), fromType="WORMBASE", toType="SYMBOL", OrgDb="org.Ce.eg.db")
names(eg)[1] <- "gene"
tab = merge(tab, eg, by="gene", all = T)
write.csv(tab, "DEGs_symbol_IP_hsp1.csv")

############### volcano plot ##################################
# function for outputting tables 
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}

# import data
data <- read.csv("DEGs_symbol_IP_hsp1.csv")

head(data) %>% 
  knitr_table()

# simple volcano plot
p1 <- ggplot(data, aes(logFC, -log10(adj.P.Val)))+
  geom_point(size = 4/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))
p1

# add color to plot
data <- data %>% 
  mutate(
    Change = case_when(logFC >= 1 & -log10(adj.P.Val)	 >= 2 ~ "upregulated",
                       logFC <= -1 & -log10(adj.P.Val)	 >= 2 ~ "downregulated",
                       TRUE ~ "Unchanged")
  )

head(data) %>% 
  knitr_table()

# add FDR labels 
data <- data %>% 
  mutate(
    Significance = case_when(
      abs(logFC) >= 1 & adj.P.Val <= 0.05 & adj.P.Val > 0.01 ~ "FDR 0.05", 
      abs(logFC) >= 1 & adj.P.Val <= 0.01 & adj.P.Val > 0.001 ~ "FDR 0.01",
      abs(logFC) >= 1 & adj.P.Val <= 0.001 ~ "FDR 0.001", 
      TRUE ~ "Unchanged")
  )

# plot
png('input_hsp1_volcano.png', width = 3.5*400, height= 2.6*400, pointsize = 8, res = 400)
p2 <- ggplot(data, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color = Significance), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  theme_bw()
dev.off()
p2

# label top 25 proteins upregulated and downregulated
top <- 25
top_genes <- bind_rows(
  data %>% 
    filter(Change == 'upregulated') %>% 
    arrange(adj.P.Val	, desc(abs(logFC))) %>% 
    head(top),
  data %>% 
    filter(Change == 'downregulated') %>% 
    arrange(adj.P.Val	, desc(abs(logFC))) %>% 
    head(top),
)
top_genes %>% 
  knitr_table()

# plot
png('IP_hsp1_volcano.png', width = 8*400, height= 12*400, pointsize = 8, res = 400)
p3 <-  p2 +
  geom_label_repel(data = top_genes, max.overlaps = Inf,
                   mapping = aes(logFC, -log10(adj.P.Val)	, label = SYMBOL),
                  size = 2)
p3
dev.off()

################################## end of script ################################################################################################

