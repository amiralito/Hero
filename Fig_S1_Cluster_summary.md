# Bar chart for cluster analysis

```R
# Load required libraries
library(ggplot2)
library(dplyr)
library(readxl)  # For reading Excel files

# Replace 'your_file.xlsx' with the path to your Excel file
data <- read_excel("/path/to/folder/Cluster_frq.xlsx")

# Define colors
color_homo <- "#1f77b4"  # Blue for Homo
color_hetero <- "#ff7f0e"  # Orange for Hetero
color_hetero_helper <- "#b30000"  # Dark Red for Hetero with Helper

# Calculate counts for each cluster
count_data <- data %>%
  group_by(ClusterSize, ClusterType) %>%
  summarize(count = n())

# Create a stacked bar chart
plot <- ggplot(count_data, aes(x = ClusterSize, y = count, fill = ClusterType)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = c("Homo" = color_homo, "Hetero" = color_hetero, "HeteroHelper" = color_hetero_helper)) +
  labs(x = "Number of Genes in Cluster", y = "Frequency") +
  theme_minimal()

# Add numbers to the bars
plot + geom_text(aes(label = count), position = position_stack(vjust = 0.5))
```
# Upset plot for cluster analysis
```R
# Load required libraries
library(UpSetR)

# Input data
expressionInput <- c(
  Rx = 9,
  H = 1,
  SD = 5,
  CCO = 8,
  CCG10 = 3,
  CCR = 0,
  TIR = 4,
  'Rx&SD' = 1,
  'Rx&CCO&TIR' = 1,
  'SD&H' = 1,
  'Rx&H' = 2,
  'CCO&H' = 1,
  'TIR&CCO' = 1
)
 

# Create the UpSet-style plot
upset(
  fromExpression(expressionInput),
  sets = c("Rx", "CCO", "SD", "TIR", "CCG10", "H", "CCR"),
  point.size = 3.5,
  line.size = 1, 
  mainbar.y.label = "Intersections",
  sets.x.label = "Set Size", 
  text.scale = c(2, 2, 0.1, 2, 2, 2),
  keep.order = TRUE
  )
```
