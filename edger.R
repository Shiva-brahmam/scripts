#load required packages
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(colorRamps)
library(dplyr)
library(ggrepel)

#set working directory
setwd("D:/Documents/RNAseq_projects/")
list.files()

# Read count data
count_data <- read.csv('counts.csv', header = TRUE, row.names = "Geneid")
head(count_data)

# Read the coldata (metadata)
coldata <- read.csv('coldata.csv', header = TRUE, row.names = 1)
coldata

# Ensure that the rownames of coldata match the colnames of count_data
all(colnames(count_data) %in% rownames(coldata))  # Should return TRUE

# Ensure that the column names of count data match the row names of coldata
all(colnames(count_data) == rownames(coldata))  # Should return TRUE

coldata$Condition <- relevel(factor(coldata$Condition), ref = "LNNA")

# Create a DGEList object for edgeR normalization
dge <- DGEList(counts = count_data, group = coldata$Condition)
dge
dim(dge)

#remove the low count genes(keeps at least 10 reads count)
keep <- filterByExpr(y = dge)

dge <- dge[keep, , keep.lib.sizes=FALSE]

#TMM normalization
dge <- calcNormFactors(object = dge)
dge_norm <- dge
#plotMDS(dge)

#log transformed counts
log_counts <- cpm(dge, log = TRUE)
write.csv(log_counts, file = "log_transformed_counts.csv", row.names = TRUE)

# Save normalized counts
normalized_counts <- cpm(dge_norm, normalized.lib.sizes = TRUE)
write.csv(normalized_counts, file = "normalized_counts.csv", row.names = TRUE)

# Create a boxplot to visualize distribution of normalized counts across samples
boxplot(normalized_counts, las = 2, main = "Normalized Counts (CPM)",
        ylab = "CPM", xlab = "Samples", col = "lightblue", outline = FALSE)

# Plot density for normalized counts
plot(density(normalized_counts[,1]), col = "blue", lwd = 2, main = "Density Plot of Normalized Counts")
for (i in 2:ncol(normalized_counts)) {
  lines(density(normalized_counts[,i]), col = i, lwd = 2)
}
legend("topright", legend = colnames(normalized_counts), col = 1:ncol(normalized_counts), lwd = 2)

#estimate both common and tagwise dispersion in one run
dge <- estimateDisp(y = dge)

#testing dge(differential gene expression)
et <- exactTest(object = dge)
et

#extract the table with adjusted p values(FDR)
top_degs = topTags(object = et, n = "Inf")
top_degs
#top_degs
dim(top_degs)

#summary (logfoldchange min 1 and padj <0.05)
summary(decideTests(object = et))

plotMD(object = et)#,main = "",xlab = "",ylab="")
#write into csv
write.csv(as.data.frame(top_degs), file="edger_res.csv")

#reading the edge results
setwd("D:/Documents/RNAseq_projects/") 

edge_res <- read.csv("edger_res.csv", header = T, row.names = 1)
dim(edge_res)

#filtering results 
edge <- edge_res %>% filter(edge_res$PValue <= 0.05)
edge_fdr <- edge_res %>% filter(edge_res$FDR <= 0.05)

#segregating the genes based upon the log fold change and Pvalue
edge_up <- edge %>% filter(edge$logFC >= 1)
edge_up2 <- edge %>% filter(edge$logFC >= 2)

#write the output into the new file 
write.csv(edge_up,"PValue_up.csv")
write.csv(edge_up2,"PValue_up2.csv")

edge_down <- edge %>% filter(edge$logFC <= -1)
edge_down2 <- edge %>% filter(edge$logFC <= -2)

#write the output into new file
write.csv(edge_down,"PValue_down.csv")
write.csv(edge_down2,"PValue_down2.csv")

#segregating the genes based upon the log fold change and FDR
edge_up <- edge_fdr %>% filter(edge_fdr$logFC >= 1)
edge_up2 <- edge_fdr %>% filter(edge_fdr$logFC >= 2)

#write the output into the new file 
write.csv(edge_up,"FDR_up.csv")
write.csv(edge_up2,"FDR_up2.csv")

edge_down <- edge_fdr %>% filter(edge_fdr$logFC <= -1)
edge_down2 <- edge_fdr %>% filter(edge_fdr$logFC <= -2)

#write the output into new file
write.csv(edge_down,"FDR_down.csv")
write.csv(edge_down2,"FDR_down2.csv")

##################################################
#extracting log transformed values for the DEGs
setwd("D:/Documents/")

log_counts <- read.csv("log_transformed_counts.csv",header = T, row.names = 1)

# for up1 (logFC >= 1)
up1 <- read.csv("FDR_up.csv", header = T, row.names = 1)
up1 <- row.names(up1)
up1 <- log_counts[up1, ]
write.csv(up1, "up1_log2_counts.csv")

#for down1
down1 <- read.csv("FDR_down.csv", header = T, row.names = 1)
down1 <- row.names(down1)
down1 <- log_counts[down1, ]
write.csv(down1, "down1_log2_counts.csv")

# for up1 (logFC >= 2)
up2 <- read.csv("FDR_up2.csv", header = T, row.names = 1)
up2 <- row.names(up2)
up2 <- log_counts[up2, ]
write.csv(up2, "up2_log2_counts.csv")

#for down1
down2 <- read.csv("FDR_down2.csv", header = T, row.names = 1)
down2 <- row.names(down2)
down2 <- log_counts[down2, ]
write.csv(down2, "down2_log2_counts.csv")

################ volcano plot #####################

#labeling genes for giving color
df <- read.csv('D:/Documents/RNAseq_projects/edger_res.csv',
               header =T, row.names =1)

df$diffexpressed <- 'NO'
df$diffexpressed[df$logFC > 2 & df$FDR <= 0.05] <- 'UP'
df$diffexpressed[df$logFC < -2 & df$FDR <= 0.05] <- 'DOWN'
dim(df)

# select top genes 
topGenes <- rownames(df)[df$FDR <= sort(df$FDR)[10] &!is.na(df$FDR)]
topGenes

#top20genes <- head(df[order(df$PValue), 'gene_symbol'], 20)
#top10genes <- head(rownames(df)[order(df$PValue)], 10)
#top10genes

#add new column to dataframe with top genes
df$delabel <- ifelse(rownames(df) %in% topGenes, rownames(df), NA)
#df$delabel <- ifelse(df$gene_symbol %in% topGenes, df$gene_symbol, NA)

## volcano plot
ggplot(data = df, aes(x = logFC, y = -log10(FDR),col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-2,2), col = 'gray', linetype = "dashed")+
  geom_hline(yintercept = c(-log10(0.05)), col = 'gray', linetype = "dashed")+
  geom_point()+
  #scale_color_manual(values = c('blue','black','red'))
  scale_color_manual(values = c("#00AFBB","gray","#bb0c00"),
                     labels = c("Down","Not sig","Up")
  )+
  #geom_text_repel(max.overlaps = Inf)+ ## remove from comment if you top genes to be labelled
  theme_test()+
  #coord_cartesian(xlim = c(-10, 10))+
  #scale_x_continuous(breaks = seq(-10,10,2))+
  labs( 
    color = 'Genes',
    x = expression('log'[2]*'FC'), y = expression('-log'[10]*'FDR'))+
  theme(axis.title = element_text(size = 12),  # Increase axis label size
        axis.text = element_text(size = 12))+   # Increase axis text size
  ggtitle("SN vs LN")+
  theme(
    legend.position = "bottom",
    #legend.position.inside = c(0.85, 0.85),  # Position the legend in the top-right corner
    legend.justification = c(1, 1),    # Justify the legend to the top-right corner
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5) # Center the plot title
  )+ 
  #theme(legend.position = "") # to remove the legend
  theme(plot.title = element_text(hjust = 0.5))#to centre the title
  #geom_text_repel(max.overlaps = Inf)


#########################################################
###################  HEATMAP OF GENES ########################

###  Heatmap of top 50 and bottom 50 degs stacked together ###
#setwd("D:/Documents/RNAseq_projects/RNAseq_shreya/DGEA/SNNA_vs_LNNA")  #LN_vs_LNNA SN_vs_SNNA

#read log counts 
setwd("D:/Documents/RNAseq_projects/SN_vs_LN")
log_counts <- read.csv('log_transformed_counts.csv',header =T, row.names =1)
head(log_counts)

#read up genes
up_genes <- read.csv('FDR_up2.csv',header = T, row.names = 1) 
#up_genes <- read.csv('Pvalue_up.csv',header = T, row.names = 1)

#up_hits <- up_genes[order(up_genes$logFC,decreasing =T), ][1:25, ]
up_hits <- up_genes[order(up_genes$FDR,decreasing =F), ][1:20, ]
up_hits

up_hits <- row.names(up_hits)
up_hits

#read downgenes
down_genes <- read.csv('FDR_down.csv',header = T, row.names = 1)
#down_genes <- read.csv('Pvalue_down.csv',header = T, row.names = 1)
#bottom_hits <- down_genes[order(down_genes$logFC,decreasing = F), ][1:25, ]
bottom_hits <- down_genes[order(down_genes$FDR,decreasing = F), ][1:20, ]
bottom_hits

bottom_hits <- row.names(bottom_hits)
bottom_hits

#keep the colour palette
rainbow_palette <- colorRamps::blue2red(100)

# Calculate average expression for the first two samples
SN <- rowMeans(log_counts[, 1:3])

# Calculate average expression for the next two samples
LN <- rowMeans(log_counts[, 4:6])

# Combine average expression values into a matrix
avg_counts <- cbind(SN, LN)
head(avg_counts)

deg <- avg_counts[bottom_hits, ]

#deg <- read.csv('D:/Documents/dissertation/deg_data/pathway genes/4.csv',header =T, row.names =1)
#setwd("D:/Documents/dissertation/deg_data")
#deg <- read.csv("Phenypropanoid.csv", header = T, row.names = 1)
#deg <- t(deg)
pheatmap(deg, 
         cluster_rows = F, 
         show_rownames = T,
         cluster_cols = F,
         legend = T,
         display_numbers = T,
         fontsize_number = 10,
         #border_color = "black",
         number_color = "black",
         #legend_breaks = c(""),
         #legend_labels = 
         #legend = T  ,   # TRUE or FALSE def true
         color = rainbow_palette,
         #color = colorRampPalette(c("#FF0000", "#FF4000", "#FF8000", "#FFBF00", "#FFFF00"))(100),
         #color = colorRampPalette(c("#FF0000", "white","blue"))(100),
         #color = cm.colors(100),
         angle_col = 0,
         #main = "Down SER vs GAR"
         )

deg <- avg_counts[up_hits, ]
pheatmap(deg, 
         cluster_rows = T, 
         show_rownames = T,
         cluster_cols = F,
         legend = T,
         color = rainbow_palette,
         angle_col = 0
         #main = "Up SER vs GAR"
         )

print("Completed")
############################################################
