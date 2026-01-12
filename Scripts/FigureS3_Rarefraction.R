rm(list=ls()); gc()

library(phyloseq)
library(stringr)
library(readxl)
library(vegan)



#####
#
#1) Reading in data
#
#####

#****ABUNDANCE DATA

#read in abundance data
abund0=read.table("novogene_analysis/featureTable.sample.total.absolute.txt", header=T, comment.char = "")

#remove names and phylogeny

abund1=abund0[,-c(1,NCOL(abund0))]
#normalization

abund2= apply(abund1, MARGIN = 2, function(x) round(100000*(x/sum(x))))


#final set
abund=abund1

#****METADATA

#read in metadata
meta0=data.frame(read_excel("~/ColiShield/Pearo_16S/Pearo_16S/Meta/All stool samples adj.xlsx"))

#add names to metadata
rownames(meta0)=meta0$sample

#final set
meta=meta0

#****TAXONOMY

#last column in abundance data has taxonomy
#split into a 7-col matrix
taxa0=str_split_fixed(string = abund0[,NCOL(abund0)], pattern = ";", n = 7)

#remove all "P__"-like patterns
taxa1=gsub(pattern = ".__", replacement = "",x = taxa0)

#name correctly
colnames(taxa1) = c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")

#final set
taxa=taxa1

#####
#2) Build phyloseq object
#####

ps <- phyloseq(otu_table(abund, taxa_are_rows=T), 
               sample_data(meta), 
               tax_table(taxa))

table(meta$study.day)

ps.sub=subset_samples(ps, groupe %in% c(1,2,3) & study.day %in% c(1,8,15,29,54) & Lysis %in% c("FP"))


###RAREFRACTION
# Extract abundance matrix from phyloseq object
otu_matrix <- as(otu_table(ps.sub), "matrix")

# If OTUs are in columns, transpose it
if (taxa_are_rows(ps.sub)) {
  otu_matrix <- t(otu_matrix)
}

# Extract group information from sample metadata (replace 'Group' with your grouping variable name)
group_info <- sample_data(ps)$groupe

# Define your preferred colors for the three groups
color_palette <- c("#2E78DB", "#FF9505", "#DB402E")
group_colors <- color_palette[as.factor(group_info)]

# Plot rarefaction curve with specified group colors
r<- rarecurve(otu_matrix, step = 1000, 
          col = group_colors, cex = 0.6, xlab = "Sample size", ylab = "Species richness")
ggsave("plots/rare.pdf",  f, width = 8, height = 6)

