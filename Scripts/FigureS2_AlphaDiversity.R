
# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(readxl)
library(stringr)
library(vegan)
library("emmeans")
library(patchwork)

#####
#
#1) Reading in data
#
#####

#****ABUNDANCE DATA

#read in abundance data
abund0 <- read.table("novogene_analysis/featureTable.sample.total.absolute.txt", header = TRUE, comment.char = "")

#remove names and phylogeny
abund1 <- abund0[, -c(1, NCOL(abund0))]

abund2 <- apply(abund1, MARGIN = 2, function(x) round(100000 * (x / sum(x))))

#final set
abund <- abund2

#****METADATA

#read in metadata
meta0 <- data.frame(read_excel("~/ColiShield/Pearo_16S/Pearo_16S/Meta/All stool samples adj.xlsx"))

#add names to metadata
rownames(meta0) <- meta0$sample

#final set
meta <- meta0

#****TAXONOMY

#last column in abundance data has taxonomy
#split into a 7-col matrix
taxa0 <- str_split_fixed(string = abund0[, NCOL(abund0)], pattern = ";", n = 7)

#remove all "P__"-like patterns
taxa1 <- gsub(pattern = ".__", replacement = "", x = taxa0)

#name correctly
colnames(taxa1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#final set
taxa <- taxa1

#####
#2) Build phyloseq object
#####

ps <- phyloseq(otu_table(abund, taxa_are_rows = TRUE), 
               sample_data(meta), 
               tax_table(taxa))

table(meta$study.day)

ps.sub <- subset_samples(ps, groupe %in% c(1, 2, 3) & study.day %in% c(1, 8, 15, 29, 54) & Lysis %in% c("FP"))

# Calculate Shannon diversity for each sample
shannon_diversity <- estimate_richness(ps.sub, measures = "Shannon")

# Combine Shannon index with sample metadata
sample_data_df <- data.frame(sample_data(ps.sub))
shannon_data <- cbind(sample_data_df, Shannon = shannon_diversity$Shannon)

# Ensure 'groupe' and 'study.day' are treated as factors
shannon_data$groupe <- as.factor(shannon_data$groupe)
shannon_data$study.day <- as.factor(shannon_data$study.day)

LM=lm(Shannon ~ groupe*study.day, data=shannon_data)

which(abs(rstandard(LM))>4)

emm=emmeans(LM, ~ groupe | study.day, adjust="sidak")

contrast(emm, method = "tukey")

plot(LM)


# Shared theme for ALL plots
theme_myplot <- theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(
      angle  = 45,
      hjust  = 1,
      colour = "black",
      face   = "bold"
    ),
    axis.text.y = element_text(
      colour = "black",
      face   = "bold"
    ),
    axis.line = element_line(
      linetype = "solid",
      linewidth = 1.0
    ),
    plot.title = element_text(
      hjust = 0.5,
      face  = "bold"
    ),
    legend.position = "right",
    panel.background = element_blank(),
    plot.background  = element_blank()
  )

# Day labels / levels
day_levels <- c("1", "8", "15", "29", "54")
day_labels <- setNames(
  c("Day 1", "Day 8", "Day 15", "Day 29", "Day 54"),
  day_levels
)

# BOX PLOT SHANNON
s_plot <- ggplot(shannon_data, aes(x = as.factor(study.day), y = Shannon, fill = groupe)) +
  geom_boxplot( position = position_dodge(width = 0.75)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 2, alpha = 0.4) +  
  scale_fill_manual(
    values = c("#2E78DB", "#FF9505", "#DB402E"),  # Custom colors
    labels = c("Control", "VHH1", "VHH1+2")          # Custom group names
  ) +  scale_x_discrete(labels = day_labels) +
  labs(x = " ", y = "Shannon Diversity Index", fill = "Group") +
  theme_myplot
print(s_plot)

###diversity

# Calculate various diversity indices for each sample
diversity_indices <- estimate_richness(ps.sub, measures = c("Observed"))

# Combine diversity indices with sample metadata
sample_data_df <- data.frame(sample_data(ps.sub))
diversity_data <- cbind(sample_data_df, diversity_indices)

# Ensure 'groupe' is treated as a factor
diversity_data$groupe <- as.factor(diversity_data$groupe)

# Print the first few rows of the combined data
head(diversity_data)

# Ensure 'groupe' and 'study.day' are treated as factors
diversity_data$groupe <- as.factor(diversity_data$groupe)
diversity_data$study.day <- as.factor(diversity_data$study.day)

LM=lm(Observed ~ groupe*study.day, data=diversity_data)

which(abs(rstandard(LM))>4)

emm=emmeans(LM, ~ groupe | study.day, adjust="sidak")

contrast(emm, method = "tukey")

plot(LM)

# BOX PLOT OBSERVED SPECIES
d_plot  <- ggplot(diversity_data, aes(x = as.factor(study.day), y = Observed, fill = groupe)) +
   geom_boxplot( position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 2, alpha = 0.4) +
  scale_fill_manual(
    values = c("#2E78DB", "#FF9505", "#DB402E"),  # Custom colors
    labels = c("Control", "VHH1", "VHH1+2")          # Custom group names
  )  + scale_x_discrete(labels = day_labels) +
  labs(x = "Study Day", y = "Observed Species", fill = "Group") +
  theme_myplot

print(d_plot)
# Combine the two plots vertically
combined_plot <- s_plot / d_plot

# Display the combined plot
print(combined_plot)

# Combine the two plots vertically  
combined_plot <- s_plot / d_plot
ggsave("plots/s_plot.pdf",  s_plot, width = 12, height = 10)
ggsave("plots/d_plot.pdf", d_plot, width = 12, height = 10)
ggsave("plots/combined_plot.pdf", combined_plot, width = 12, height = 14)





