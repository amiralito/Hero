# RNA-Seq pipeline - Poisson Distance between samples

Load required libraries:
```R
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library(RColorBrewer)
```
Input data from DeSeq2:
```R
poisd <- PoissonDistance(t(counts(dds)))
```

Generate Matrix:
```R
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$name)
colnames(samplePoisDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(
  samplePoisDistMatrix, 
  clustering_distance_rows = poisd$dd, 
  clustering_distance_cols = poisd$dd, 
  col = colors, 
  fontsize = 15,          # Adjust font size
  legend = TRUE,          # Ensure the legend is shown
  legend.title = "Title", # Set a title for the legend if desired
  cellwidth = 50,         # Increase width of each cell
  cellheight = 50,        # Increase height of each cell
  width = 50,             # Increase overall plot width
  height = 50              # Increase overall plot height
)
```

#RNA-Seq pipeline - PCA of samples

Transform data from DeSeq2:
```R
vsd <- vst(dds, blind = FALSE)
```

Generate PCA plot:
```R
PCA <- plotPCA(vsd, intgroup = c("tissue", "name"))
PCA + 
  coord_fixed(ratio = 20, xlim = NULL, ylim = NULL) + 
  theme_minimal() +
  theme(
    text = element_text(size = 15),  # Increase base font size
    plot.title = element_text(size = 20, face = "bold"),  # Title font size
    axis.title = element_text(size = 16),  # Axis titles font size
    axis.text = element_text(size = 14)   # Axis text font size
  ) +
  geom_point(size = 6)  # Increase point size
```
  
#Plot normalized counts in same plot from list

Load required libraries:
```R
library(tidyr)
```

Input gene IDs for normalized counts data:
```R
gene_data <- read.csv("path/to/your/gene_data.csv")
gene_ids <- unique(gene_data$ID)
combined_data <- data.frame()
```

Generate figure for normalized counts data in one plot:
```R
plot_data <- tryCatch({
    plotCounts(dds, gene = gene_id, intgroup = "tissue", returnData = TRUE)
  }, error = function(e) {
    message(paste("Error for gene:", gene_id, " - ", e$message))
    return(NULL)  # Return NULL if there’s an error
  })
  if (!is.null(plot_data)) {
    plot_data$gene_id <- gene_id
     combined_data <- rbind(combined_data, plot_data)
  }
}
if (nrow(combined_data) == 0) {
  stop("No valid data was found for the specified gene IDs.")
}
combined_data$gene_id <- factor(combined_data$gene_id, levels = gene_ids)
ggplot(combined_data, aes(x = gene_id, y = count, color = tissue)) +
  geom_jitter(width = 0.2, size = 4, shape = 1, stroke = 3) +  
  theme_bw() +
  scale_y_log10(limits = c(1, 10000000), oob = scales::squish) +
  labs(x = "Gene ID", y = "Log10 Count", title = "Gene Expression") +
  theme(
    legend.title = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
```
