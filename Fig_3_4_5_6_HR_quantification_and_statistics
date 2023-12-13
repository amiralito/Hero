#HR quantification

```R
# Load required libraries
library(tidyverse)
library(besthr)
library(ggplot2)
library(patchwork)

# Read in data

# Replicate = replicate number: 1, 2, 3, n
# Date = date experiment started
# DPI = days post-infiltration scored
# Leaf = leaf image filename
# Genotype = genotype name
# Genotype_code = genotype seed name
# NLR1 = First infiltration variable. Name can be changed to reflect experiment
# NLR2 = Second infiltration variable. Name can be changed to reflect experiment. E.g. Mutation
# Effector = Third infiltration variable. Name can be changed to reflect experiment. E.g. Mutation
# Score = numeric HR score: 1, 2, 3, ..., 7
# Treatment = the combination of Genotype, Genotype_code, NLR1, NLR2, and Effector relevant to your experiment

df <- readxl::read_excel("/path/to/folder/NRC_HR_quantification_example_Hero_A.xlsx") %>% # Add your .xlsx here
  filter(!is.na(Score)) %>% # Remove missing values
  mutate(Replicate = factor(Replicate, levels = c("1", "2", "3", "4", "5", "6")),
         Genotype = factor(Genotype, levels = c("wt")),
         Genotype_code = factor(Genotype_code, levels = c("wt")),
         NLR1 = factor(NLR1, levels = c("NRC6a", "EV")), # Add appropriate levels
         NLR2 = factor(NLR2, levels = c("Hero_A", "EV")), # Add appropriate levels
         Effector = factor(Effector, levels = c("WT", "DV")), # Add appropriate levels
         Score = factor(Score, levels = c("0", "1", "2", "3", "4", "5", "6", "7")),
         Treatment = paste0(NLR2, NLR1, Effector)) # Add appropriate treatment interaction for your experiment

# Always manually check your dataframe to make sure NAs were not introduced
head(df)

# Group data by Score, NLR1, NLR2, Effector, and count number of occurrences 
df_grouped <- df %>% 
  group_by(Score, Genotype, Genotype_code, NLR1, NLR2, Effector, Treatment) %>%
  summarise(count = n())

# Some statistics that are probably not very meaningful but might give an approximate idea of differences
df_stats <- df %>%
  mutate(Score = as.numeric(as.character(Score))) %>%
  group_by(Genotype, Genotype_code, NLR1, NLR2, Effector, Treatment) %>%
  summarise(N = length(Score), 
            Min = signif(min(Score), digits = 2),
            Max = signif(max(Score), digits = 2),
            Median = signif(median(Score), digits = 2),
            Std = signif(sd(Score), digits = 2),
            Mean = mean(Score)) %>%
  mutate(SE = signif(Std / sqrt(N), digits = 2),
         Lower_ci = signif(Mean - qt(1 - (0.05 / 2), N - 1) * SE, digits = 2),
         Upper_ci = signif(Mean + qt(1 - (0.05 / 2), N - 1) * SE, digits = 2)) %>%
  ungroup()

################################################################
################################################################
################################################################

# Assign colours, Change as appropriate to the experiment
dfcols <- c("WT" = "#8faad0",
            "DV" = "#c9a9ff")

# Plot setup with Effector on x-axis and Score on y-axis
df_plot <- ggplot(df_grouped, aes(x = Treatment, y = Score)) +
  geom_count(
    aes(size = count, fill = Effector), 
    shape = 21, 
    stroke = 1,
    show.legend = TRUE
  ) +
  geom_jitter(
    data = df, 
    aes(color = Replicate), 
    size = 1.5, 
    width = 0.3, 
    height = 0.2
  ) +
  scale_fill_manual("", values = dfcols, breaks = NULL) +
  scale_y_discrete(limits = c("0", "1", "2", "3", "4", "5", "6", "7")) +
  scale_x_discrete(limits = c("Hero_AEVWT", "Hero_AEVDV", "Hero_ANRC6aWT", "Hero_ANRC6aDV", "EVNRC6aWT", "EVNRC6aDV")) +
  labs(x = "", y = "Cell death index", size = "Count") +
  scale_color_manual(values = c ("#d0800b", "#00aa00","#F8F300", "#004488")) +
  scale_size_continuous(range = c(1, 12)) +
  #facet_wrap( ~ NLR2, nrow = 1, scales = "free_x") + # If you want to facet wrap along one of the variables. Otherwise comment this line out
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "lightgray"),
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text.x  = element_text(size = 16, colour = "black"),
    strip.background = element_blank(), 
    strip.placement = "outside",
    legend.key.size = unit(0.35,'cm')
  )


# Display the plot
df_plot

ggsave("/path/to/folder/NRC_HR_quantification_Hero_A.pdf", width = 5, height = 5) # Change name to appropriate filename
```

# HR statistical analysis

# Define the desired order of treatments
desired_order <- c(
  "Hero_AEVWT", "Hero_AEVDV", "Hero_ANRC6aWT", 
  "Hero_ANRC6aDV", "EVNRC6aWT", "EVNRC6aDV"
)

# Modify the Treatment variable levels
df$Treatment <- factor(df$Treatment, levels = desired_order)

hr_est_combined_Sample <- estimate(df, Score, Treatment, Replicate, Control = "Hero_ANRC6aWT", nits = 1000) # Change control to appropriate control

p <- plot(hr_est_combined_Sample)
p[[1]] <- p[[1]] + 
  scale_colour_manual(values = c("#8faad0", "#c9a9ff", "#8faad0", "#c9a9ff", "#8faad0", "#c9a9ff")) 
p[[2]] <- p[[2]] + scale_fill_manual(
  values = c("#5f5fffff", "#FF5F5FFF", "lightgray"),
  name = "bootstrap percentile", labels=c("<0.025", ">0.975", "0.025−0.975"),
  guide = guide_legend(reverse=TRUE)
)
plot(p)


ggsave("/path/to/folder/Hero_A/NRC_HR_stats_Hero_A.pdf", width = 10, height = 3.5) # Change name to appropriate filename
```
