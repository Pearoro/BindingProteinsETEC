## Load libraries ----
library(phyloseq)
library(dplyr)
library(ggplot2)
library(microshades)
library(forcats)
library(stringr)
library(readxl)
library(openxlsx)
library(writexl)

## Custom palette ----
palette <- c(
  Other                 = "#404E4D",  # Very Dark Gray/Black
  Bariatricus           = "#5E0511",  # Deep Red
  Blautia               = "#8D0819",  # Soft Red-Orange
  Christensenella       = "#BB0A21",  # Beige
  Clostridium           = "#016FB9",  # Deep Blue
  Dysosmobacter         = "#8EA604",  # Medium Gray
  Lactobacillus         = "#FF9505",
  Mediterraneibacter    = "#7E3267",  # Rust
  Oscillibacter         = "#5F468A",  # Dark Brown
  Phascolarctobacterium = "#7878BA",  # Saddle Brown
  Terrisporobacter      = "#4059AD"   # Royal Blue
)

## Read the Excel files ---- differnt 16S run in nanopore
abundance_df1 <- read_excel(
  "~/ColiShield/Pearo_16S/Pearo_16S/Nano/Nano1/Emu-Species-plot.xlsx"
)
abundance_df2 <- read_excel(
  "~/ColiShield/Pearo_16S/Pearo_16S/Nano/Nano1/Nano run 2025/emu-combined-species.run22025_plot.xlsx"
)

## Merge the dataframes by the common columns ----
merge_columns <- c("species", "genus", "family", "order", "class", "phylum", "superkingdom")

abundance <- full_join(abundance_df1, abundance_df2, by = merge_columns)

writexl::write_xlsx(
  abundance,
  path = "~/ColiShield/Pearo_16S/Pearo_16S/abundance.xlsx"
)

## Separate table in two (OTU and taxa) ----
# Identify the columns for taxa levels
taxa_columns <- c(
  "species", "genus", "family",
  "order", "class", "phylum", "superkingdom"
)

# Identify the columns for samples (those starting with 'S')
sample_columns <- grep("^S", names(abundance), value = TRUE)

# OTU table
otu  <- abundance[, sample_columns]

# Taxa table
taxa <- abundance[, taxa_columns]

## Name taxonomy columns correctly ----
colnames(taxa) <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Superkingdom")

# Re-order as the other file
taxa1 <- taxa[, c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]

## Ensure that the OTU table contains only numeric values ----
otu <- otu %>%
  mutate(across(everything(), as.numeric))

## Read in metadata ----
metadata_path <- "~/ColiShield/Pearo_16S/Pearo_16S/Nano/Nano1/meta16full.xlsx"

# Check if the file exists
if (!file.exists(metadata_path)) {
  stop("Metadata file does not exist at the specified path: ", metadata_path)
}

# Read metadata
meta <- data.frame(read_excel(metadata_path))

# Print the first few rows of the metadata to check if it has been read correctly
print(head(meta))

# Add names to metadata
rownames(meta) <- meta$sample

# Final set
meta <- meta

## Name taxonomy columns correctly (again, if needed) ----
colnames(taxa) <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Superkingdom")

## Convert taxa_df to a matrix ----
taxa_matrix <- as.matrix(taxa1)

# Check the structure of the matrix
str(taxa_matrix)

## Create phyloseq object with the converted matrix ----
ps <- phyloseq(
  otu_table(otu, taxa_are_rows = TRUE),
  sample_data(meta),
  tax_table(taxa_matrix)
)

ps_melt <- psmelt(ps)

## Remove taxa with NA and zero counts ----
# Remove taxa with NA in any taxonomic rank
ps <- prune_taxa(!apply(is.na(tax_table(ps)), 1, any), ps)

# Remove taxa with zero counts across all samples
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

## Subset samples by specific study days ----
ps.sub <- subset_samples(ps, study.day %in% c(1, 8, 15, 29, 54))

## Load fantaxtic ----
library(fantaxtic)

## Get nested top taxa ----
top_nested <- nested_top_taxa(
  ps,
  top_tax_level   = "Genus",
  nested_tax_level = "Species",
  n_top_taxa      = 10,
  n_nested_taxa   = 5
)

## Plot nested bar plot ----
p <- plot_nested_bar(
  ps_obj       = top_nested$ps_obj,
  top_level    = "Genus",
  nested_level = "Species",
  palette      = palette
)

p

## ---- Optional: labels for plotting facets/scales later ----
day_labels   <- c("1" = "Day 1", "8" = "Day 8", "15" = "Day 15",
                  "29" = "Day 29", "54" = "Day 54")
group_labels <- c("1" = "Control", "2" = "VHH1", "3" = "VHH1+2")

p + facet_wrap(~ study.day + groupe, nrow = 1, scales = "free_x", labeller = labeller(study.day = day_labels, groupe = group_labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom") +
  labs(title = "Top 5 Taxa", x = "Taxa", y = "Relative Abundance")


## top 15

top_nested <- nested_top_taxa(ps,
                              top_tax_level = "Genus",
                              nested_tax_level = "Species",
                              n_top_taxa = 15, 
                              n_nested_taxa = 5)



p15<- plot_nested_bar(ps_obj = top_nested$ps_obj,
                      top_level = "Genus",
                      nested_level = "Species",
                      palette = palette15)
p15
p15 + facet_wrap(~  study.day + groupe  , nrow = 1, scales = "free_x", labeller = labeller(study.day = day_labels, groupe = group_labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom") +
  labs(title = "Top 5 Taxa", x = "Taxa", y = "Relative Abundance")



# Define your custom palette
palette15 <- c(
  Other = "#404E4D",  
  Anaerobutyricum = "#5E0511",
  Bariatricus =  "#8D0819", 
  Christensenella =  "#BB0A21",  
  Clostridium = "#016FB9",
  Blautia =  "#8ecae6", 
  Dysosmobacter =  "#58b4d1",  
  Eubacterium =  "#219ebc",
  Lacrimispora = "#023047",
  Lactobacillus = "#ffb703",  
  Mediterraneibacter =  "#fb8500",  
  Sporobacter = "#999300",
  Oscillibacter =  "#7E3267",  
  Phascolarctobacterium = "#5F468A",  
  Roseburia =  "#7878BA", 
  Terrisporobacter = "#4059AD"
)
