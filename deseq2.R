#load required libraries
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)

#setup working directory
setwd("D:/Documents/project")
list.files()

#read counts data
counts_data <- read.csv('counts.csv',header = TRUE,row.names = 1)
head(counts_data)

#read in sample info
colData <- read.csv('coldata.csv', header = TRUE,row.names = 1)
colData

#make sure the row names in colData matches to column names in counts_data
all(colnames((counts_data) %in% rownames(colData)))

#are they in same order
all(colnames(counts_data) == rownames(colData))

#set factor levels
colData$condition <- factor(colData$condition)

#construct a DESeqDataSet object ------
dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = colData, 
                              design = ~condition )
dds

#pre-filtering: removing rows with low gene counts
#keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#perform the statistical tests to identify differentially expressed genes
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds)
res

#explore the results
summary(res)

#variance stablizing transformation
vsd <- vst(dds,blind = FALSE )

#use the transformed values to generate a pca plot
plotPCA(vsd,intgroup = "condition")

#MAplot
plotMA(res)

#volcano plot
res_dataframe <- as.data.frame(res)
EnhancedVolcano(res_dataframe,
                lab = rownames(res_dataframe),
                x = 'log2FoldChange',
                y = 'pvalue')

#heatmap
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds))
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=df)
