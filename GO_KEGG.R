#load the required libraries
library(clusterProfiler)
library(org.At.tair.db)   # Annotation database for Arabidopsis thaliana
library(AnnotationHub)
library(biomaRt)
library(KEGGREST)
library(ggplot2)
library(dplyr)
library(stringr)

setwd("D:/Documents/RNAseq_projects/RNAseq_shreya/DGEA/SN_vs_LN")
genes <- read.csv('FDR_up2.csv',header =T)
head(genes)

genes <- genes$X

#Gene list
gene_list <- genes
#gene_list <- c(" AT5G47340","AT5G49180","AT3G09540","AT5G10160","AT3G49360","AT2G41850")

# Use biomaRt to convert TAIR IDs to Entrez IDs
mart <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "https://plants.ensembl.org")

# Connect to Ensembl and the human dataset (Homo sapiens)
#mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# List all attributes for Arabidopsis thaliana in Ensembl Plants
attributes <- listAttributes(mart)
head(attributes, 20)  # Check the first 20 attributes

genes_entrez <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
                      filters = "ensembl_gene_id", 
                      values = gene_list, 
                      mart = mart)

head(genes_entrez)

#performe GO

go_enrichment <- enrichGO(
  gene          = genes_entrez$entrezgene_id,  # Use Entrez IDs
  OrgDb         = org.At.tair.db,              # Arabidopsis thaliana database
  keyType       = "ENTREZID",                  # Key type
  ont           = "ALL",                       # You can specify "BP", "MF", or "CC"
  #ont           = "BP",
  pAdjustMethod = "BH",                        # Adjust p-values for multiple testing
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# View the GO results
head(go_enrichment)

# Convert the enrichment result to a data frame for visualization
go_enrichment_df <- as.data.frame(go_enrichment)
write.csv(go_enrichment_df, "UP_GO.csv")

# Plot the top 10 GO terms, coloring by p.adjust value
ggplot(go_enrichment_df[1:10, ], aes(x = reorder(Description, -Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip coordinates for a horizontal bar plot
  scale_fill_gradient(low = "red", high = "blue") +  # Color gradient based on p.adjust
  xlab("GO Terms") +
  ylab("Gene Count") +
  #ggtitle("Top") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 10)
  )

# Get the data from the dotplot
dot_data <- dotplot(go_enrichment, showCategory = 10)$data

# Create a custom dotplot with gene counts displayed on the plot
ggplot(dot_data, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
  geom_point() +
  #geom_text(aes(label = Count), vjust = -1.5) +  # Add gene count labels
  scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
  scale_size(range = c(2, 4)) +  # Adjust the size range of the dots
  xlab("Gene Ratio") +
  ylab("GO Term") +
  ggtitle("GO Enrichment (Biological Process)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  )


# for shinygo like plot
# Calculate Fold Enrichment and -log10(FDR)
data <- go_enrichment_df %>%
  mutate(
    GeneRatio_num = as.numeric(str_split(GeneRatio, "/", simplify = TRUE)[,1]),
    GeneRatio_den = as.numeric(str_split(GeneRatio, "/", simplify = TRUE)[,2]),
    BgRatio_num = as.numeric(str_split(BgRatio, "/", simplify = TRUE)[,1]),
    BgRatio_den = as.numeric(str_split(BgRatio, "/", simplify = TRUE)[,2]),
    # Fold Enrichment based on the corrected understanding
    Fold_Enrichment = (GeneRatio_num / GeneRatio_den) / (BgRatio_num / BgRatio_den),
    #Fold_Enrichment = (GeneRatio_num / GeneRatio_den) / (BgRatio_num / BgRatio_den),
    Log10FDR = -log10(p.adjust)
  )

# Select the top rows based on Log10FDR and then Fold Enrichment
top_data <- data %>%
  arrange(desc(Log10FDR), desc(Fold_Enrichment)) %>%  # Sort by Log10FDR (descending) and then Fold Enrichment
  head(50)  # specify the number of top rows

# Define min and max for Log10FDR to set dynamic limits for the color scale
log10FDR_min <- min(data$Log10FDR, na.rm = TRUE)  # use only when points dont show colors
log10FDR_max <- max(data$Log10FDR, na.rm = TRUE)

# Plot
ggplot(top_data, aes(x = Fold_Enrichment, y = reorder(Description, Fold_Enrichment))) +
  geom_segment(aes(x = 0, xend = Fold_Enrichment, yend = reorder(Description, Fold_Enrichment), color = Log10FDR),
               size = 0.8) + # linewidth = 0.8  for 3.4.0 ggplot # Color of line follows Log10FDR
  geom_point(aes(size = Count, color = Log10FDR)) +
  scale_color_gradientn(
    #colors = cm.colors(100),
    colors = c("blue", "purple","red"),
                        #limits = c(log10FDR_min, log10FDR_max),  # Dynamic color scale limits
                        oob = scales::squish,  # Ensures values outside range are squished to nearest color
                        guide = guide_colorbar(barwidth = 1.5, barheight = 10)) +
  scale_size_continuous(range = c(1, 4)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 9, colour = "black"),
    axis.text.y = element_text(size = 9, colour = "black", hjust = 1),
    axis.title.x = element_text(size = 11),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2))

#for viusalization of GO results
# Bar plot of GO terms
barplot(go_enrichment, showCategory = 20)

# Dot plot
dotplot(go_enrichment, showCategory = 20)

# GO graph plot
plotGOgraph(go_enrichment)

#perform KEGG
kegg_enrichment <- enrichKEGG(
  gene         = gene_list,
  organism     = 'ath',                       # Arabidopsis thaliana KEGG organism code
  pvalueCutoff = 0.05
)

kegg_data <- as.data.frame(kegg_enrichment@result)
head(kegg_data)
write.csv(kegg_data, "up kegg.csv")

# View the KEGG results
head(kegg_enrichment)

#visualize Kegg
# Bar plot of KEGG pathways
barplot(kegg_enrichment, showCategory = 20)

# Dot plot
dotplot(kegg_enrichment, showCategory = 20)


#custom plot 
# Calculate -log10 of the p-values for better visualization
kegg_data$p_log10 = -log10(kegg_data$pvalue)
# Remove " - Arabidopsis thaliana (thale cress)" from Description column
kegg_data$Description <- str_replace(kegg_data$Description, " - Arabidopsis thaliana \\(thale cress\\)", "")

# Convert GeneRatio to numeric (assuming format is "3/5", "1/5", etc.)
kegg_data$GeneRatio_numeric <- sapply(kegg_data$GeneRatio, function(x) {
  ratio <- as.numeric(str_split(x, "/")[[1]])
  return(ratio[1] / ratio[2])  # Return the numeric value of the GeneRatio
})

# Create the dot plot
ggplot(kegg_data[1:50,], aes(x = GeneRatio_numeric, y = reorder(Description, GeneRatio_numeric), 
                      size = Count, color = p.adjust)) +
  geom_point() +
  scale_size_continuous(range = c(1.5, 4)) +  # Size of dots based on Gene Count
  scale_color_gradient(low = "blue", high = "red") +  # Color based on p.adjust values
  #scale_color_gradient(low = "lightgreen", high = "darkred") +
  #scale_color_gradient(low = "blue", high = "yellow")+
  #scale_color_gradient(low = "purple", high = "orange")+
  #scale_color_gradient(low = "lightblue", high = "darkblue")+
  #scale_color_gradient(low = "purple", high = "yellow")+
  #scale_color_gradientn(colors = colorRamps::green2red(100))+ 
  labs(x = "Gene Ratio", y = "KEGG Pathway", 
      # title = "KEGG Enrichment An",
       #subtitle = "Gene Ratio, Size = Gene Count, Color = Adjusted p-value",
       color = "P.Adjust", size = "Gene Count") +
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(size = 9, colour = "black", hjust = 1)
    )

# bar plot
# Create the bar plot with colorRamps::blue2red gradient
ggplot(kegg_data[1:50,], aes(x = GeneRatio_numeric, y = reorder(Description, GeneRatio_numeric), 
                      fill = p.adjust)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), 
            hjust = -0.1,  # Position text slightly outside the bar
            size = 3.5) +    # Adjust text size as needed
  #scale_fill_gradientn(colors = colorRamps::blue2red(100)) +  
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Gene Ratio", y = "KEGG Pathway", 
       #title = "KEGG Enrichment Analysis Bar Plot",
       #subtitle = "Gene Ratio, Color = Adjusted p-value",
       fill = "P.Adjust") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(size = 9, colour = "black")#, hjust = 1)
    )

# Pathway visualization (using KEGG)
# This requires pathview package for detailed visualization
#BiocManager::install("pathview")
library(pathview)

# Visualize a specific KEGG pathway (e.g., pathway ID "ath00010" for Glycolysis / Gluconeogenesis)
pathview(gene.data = genes_entrez$entrezgene_id, pathway.id = "ath01200", species = "ath")

