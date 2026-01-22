# Load necessary libraries
library(phyloseq)
library(dplyr)
library(stringr)
library(readxl)
library(openxlsx)
library(writexl)
library(tidyr)

# Clear workspace
rm(list=ls()); gc()
writeTable=F
# 1) Reading in data

# Abundance data
abund0 <- read.table("novogene_analysis/featureTable.sample.total.absolute.txt", header=T, comment.char = "")
abund1 <- abund0[,-c(1,NCOL(abund0))]
abund2= apply(abund1, MARGIN = 2, function(x) round(100000*(x/sum(x))))
#final set
abund=abund2



# Metadata
meta0 <- data.frame(read_excel("~/ColiShield/Pearo_16S/Pearo_16S/Meta/All stool samples adj.xlsx"))
rownames(meta0) <- meta0$sample
meta <- meta0

# Taxonomy
taxa0 <- str_split_fixed(string = abund0[,NCOL(abund0)], pattern = ";", n = 7)
taxa1 <- gsub(pattern = ".__", replacement = "", x = taxa0)
colnames(taxa1) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa <- taxa1

# 2) Build phyloseq object
ps <- phyloseq(otu_table(abund, taxa_are_rows=T), sample_data(meta), tax_table(taxa))
ps.sub <- subset_samples(ps, groupe %in% c(1,2,3) & study.day %in% c(1,8,15,29,54) & Lysis %in% c("FP"))
ps.sub.genus <- tax_glom(ps.sub, taxrank = "Genus")

# Ensure 'groupe' is treated as a factor
ps.sub.genus@sam_data$groupe <- as.factor(ps.sub.genus@sam_data$groupe)

# 3) Data processing

# Melt the phyloseq object to a data frame
psdf_genus <- psmelt(ps.sub.genus)

# Calculate relative abundance
psdf_genus <- psdf_genus %>%
  group_by(Sample) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

# Filter the top 20 genera based on mean relative abundance
top_genera <- psdf_genus %>%
  group_by(Genus) %>%
  summarize(MeanRelativeAbundance = mean(RelativeAbundance)) %>%
  arrange(desc(MeanRelativeAbundance)) %>%
  slice(1:5) %>%
  pull(Genus)

# Filter the data for the top 20 genera
filtered_genus <- psdf_genus %>%
  filter(Genus %in% top_genera)

# 4) Calculate summary statistics

summary_stats <- filtered_genus %>%
  group_by(groupe, study.day, Genus) %>%
  summarize(
    AV = mean(RelativeAbundance),
    SD = sd(RelativeAbundance),
    SE = SD / sqrt(n()),
    MIN = min(RelativeAbundance),
    MAX = max(RelativeAbundance)
  ) %>%
  ungroup()

# Pivot the table for better readability
summary_table <- summary_stats %>%
  pivot_wider(names_from = Genus, values_from = c(AV, SD, SE, MIN, MAX))

# Print the summary table
print(summary_table)

# Write the summary table to an Excel file
write_xlsx(summary_table, "genus_summary_statistics.xlsx")

# 5) Create and save tables for specific genera of interest

# List of genera of interest
genera_of_interest <- c("Escherichia-Shigella", "Lactobacillus", "Prevotella", "Prevotella_9", "Ruminococcus", "Clostridium_sensu_stricto_1", "Blautia", "Alloprevotella")

# Function to create summary table for a given genus
create_genus_table <- function(genus_name, data) {
  genus_data <- data %>%
    filter(Genus == genus_name) %>%
    group_by(groupe, study.day) %>%
    summarize(
      MeanRelativeAbundance = mean(RelativeAbundance),
      SD = sd(RelativeAbundance),
      SE = SD / sqrt(n()),
      MIN = min(RelativeAbundance),
      MAX = max(RelativeAbundance)
    ) %>%
    ungroup()
  
  return(genus_data)
}

# Create and save tables for each genus
for (genus in genera_of_interest) {
  genus_table <- create_genus_table(genus, psdf_genus)
  write_xlsx(genus_table, paste0(genus, "_summary_statistics.xlsx"))
}

# Print the tables to the console (optional) and save to a single Excel file
all_genus_tables <- list()

for (genus in genera_of_interest) {
  genus_table <- create_genus_table(genus, psdf_genus)
  print(genus)
  print(genus_table)
  all_genus_tables[[genus]] <- genus_table
}

# Write all tables to a single Excel file with different sheets
write.xlsx(all_genus_tables, file = "all_genus_summary_statistics.xlsx")

# 6) Family-level analysis

# Subset the phyloseq object at the family level
ps.sub.fam <- tax_glom(ps.sub, taxrank = "Family")

# Melt the phyloseq object to a data frame
psdf_family <- psmelt(ps.sub.fam)

# Calculate relative abundance
psdf_family <- psdf_family %>%
  group_by(Sample) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

if(writeTable==T) {
  # Filter the top 20 families based on mean relative abundance
  top_families <- psdf_family %>%
    group_by(Family) %>%
    summarize(MeanRelativeAbundance = mean(RelativeAbundance)) %>%
    arrange(desc(MeanRelativeAbundance)) %>%
    slice(1:20) %>%
    pull(Family)
  
  # Filter the data for the top 20 families
  filtered_family <- psdf_family %>%
    filter(Family %in% top_families)
  
  # Calculate summary statistics
  summary_stats_family <- filtered_family %>%
    group_by(groupe, study.day, Family) %>%
    summarize(
      AV = mean(RelativeAbundance),
      SD = sd(RelativeAbundance),
      SE = SD / sqrt(n()),
      MIN = min(RelativeAbundance),
      MAX = max(RelativeAbundance)
    ) %>%
    ungroup()
  
  # Pivot the table for better readability
  summary_table_family <- summary_stats_family %>%
    pivot_wider(names_from = Family, values_from = c(AV, SD, SE, MIN, MAX))
  
  # Print the summary table
  print(summary_table_family)
  
  # Write the summary table to an Excel file
  write_xlsx(summary_table_family, "family_summary_statistics.xlsx")
  
  
  
  # Pivot the table for better readability
  summary_table_family <- summary_stats_family %>%
    pivot_wider(names_from = Family, values_from = c(AV, SD, SE, MIN, MAX))
  
  # Print the summary table
  print(summary_table_family)
  
  # Write each top family to a different sheet in an Excel file
  family_tables <- list()
  
  for (family in top_families) {
    family_table <- summary_stats_family %>%
      filter(Family == family) %>%
      pivot_wider(names_from = Family, values_from = c(AV, SD, SE, MIN, MAX))
    family_tables[[family]] <- family_table
  }
  
  write.xlsx(family_tables, file = "family_summary_statistics.xlsx")
  
} 

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Ensure 'groupe' is treated as a factor in psdf_family
psdf_family$groupe <- as.factor(psdf_family$groupe)
psdf_family$logAbundance = log10(psdf_family$Abundance+1)

# Filter for Enterobacteriaceae
plot_data <- psdf_family %>%
  filter(Family == "Enterobacteriaceae")

# Define custom labels
day_labels <- c("Day 1", "Day 8", "Day 15", "Day 29", "Day 54")
names(day_labels) <- c("1", "8", "15", "29", "54")

group_labels <- c("Control", "Ablacto+1", "Ablacto+2")
names(group_labels) <- c("1", "2", "3")

# Create plot for Enterobacteriaceae
e <- ggplot(plot_data, aes(x = as.factor(study.day), y = logAbundance, color = groupe, group = groupe)) +
  geom_point(position = position_dodge(width = 1), size = 2, alpha = 0.6) +  # Adjust position dodge width
  geom_smooth(method = "loess", se = TRUE, aes(group = groupe), linetype = "dashed") +  # Add trend lines with error bands
  scale_x_discrete(labels = day_labels) +
  scale_color_manual(values = c("1" = "#2E78DB", "2" = "#FF9505", "3" = "#DB402E"), labels = group_labels) +
  labs(title = "Dynamic Changes in Enterobacteriaceae",
       x = "Day",
       y = "Relative Abundance",
       color = "Group") +
  #ylim(0, 0.2) +  # Set y-axis limits
  theme(
    axis.text.x = element_text(size = 22, colour = "black", vjust = 0.25, angle = 45, face = 'plain'),
    axis.text.y = element_text(size = 22, colour = "black", vjust = 0.25, face = 'plain'),
    strip.text = element_text(size = 10, face = "bold"),
    axis.line = element_line(linetype = "solid", linewidth = 1.0),
    axis.title = element_text(size = 22),
    plot.title = element_text(size = 20),
    legend.position = "none",
    panel.background = element_blank(),  # Remove panel background
    plot.background = element_blank(),   # Remove plot background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )

# Print the plot
print(e)


#### calculate percentages of top 5 of each?##

# Top 5 Phyla
ps.sub.phyla <- tax_glom(ps.sub, taxrank = "Phylum")
psdf_phyla <- psmelt(ps.sub.phyla) %>%
  group_by(Sample) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

top_phyla <- psdf_phyla %>%
  group_by(Phylum) %>%
  summarize(MeanRelativeAbundance = mean(RelativeAbundance)) %>%
  arrange(desc(MeanRelativeAbundance)) %>%
  slice(1:5) %>%
  pull(Phylum)

phyla_percentage <- psdf_phyla %>%
  filter(Phylum %in% top_phyla) %>%
  group_by(Phylum) %>%
  summarize(Percentage = mean(RelativeAbundance) * 100) %>%
  arrange(desc(Percentage))  # Sort by Percentage in descending order


# Top 5 Families
ps.sub.family <- tax_glom(ps.sub, taxrank = "Family")
psdf_family <- psmelt(ps.sub.family) %>%
  group_by(Sample) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

top_families <- psdf_family %>%
  group_by(Family) %>%
  summarize(MeanRelativeAbundance = mean(RelativeAbundance)) %>%
  arrange(desc(MeanRelativeAbundance)) %>%
  slice(1:5) %>%
  pull(Family)

family_percentage <- psdf_family %>%
  filter(Family %in% top_families) %>%
  group_by(Family) %>%
  summarize(Percentage = mean(RelativeAbundance) * 100) %>%
  arrange(desc(Percentage))  # Sort by Percentage in descending order


# Top 5 Genera
ps.sub.genus <- tax_glom(ps.sub, taxrank = "Genus")
psdf_genus <- psmelt(ps.sub.genus) %>%
  group_by(Sample) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

top_genera <- psdf_genus %>%
  group_by(Genus) %>%
  summarize(MeanRelativeAbundance = mean(RelativeAbundance)) %>%
  arrange(desc(MeanRelativeAbundance)) %>%
  slice(1:5) %>%
  pull(Genus)

genus_percentage <- psdf_genus %>%
  filter(Genus %in% top_genera) %>%
  group_by(Genus) %>%
  summarize(Percentage = mean(RelativeAbundance) * 100)  %>%
  arrange(desc(Percentage))  # Sort by Percentage in descending order

# Print the results
print("Top 5 Phyla Percentage:")
print(phyla_percentage)

print("Top 5 Families Percentage:")
print(family_percentage)

print("Top 5 Genera Percentage:")
print(genus_percentage)