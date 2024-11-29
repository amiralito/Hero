# Heatmap of expression data

Load required libraries:
```R
library(ggplot2)
library(pheatmap)
```
Input data file:
```R
data <- read.csv("/Users/luedke/Desktop/Heinz_expression.csv", row.names = 1)
```

Convert to matrix:
```R
expression_matrix <- as.matrix(data)
```

Apply log2 transformation:
```R
log_expression_matrix <- log2(expression_matrix + 1)
```

Create the heatmap with log-transformed values displayed:
```R
pheatmap(
  log_expression_matrix,
  scale = "none",            # Do not standardize, since we're using log-transformed data
  clustering_distance_rows = "euclidean",  # Distance measure for rows
  clustering_distance_cols = "euclidean",  # Distance measure for columns
  clustering_method = "complete",          # Clustering method
  color = colorRampPalette(c("blue", "white", "red"))(50),  # Custom color palette
    fontsize_row = 8,         # Adjust row font size
  fontsize_col = 10,        # Adjust column font size
  main = "Gene Expression Heatmap (Log2 Scale)"
)
```