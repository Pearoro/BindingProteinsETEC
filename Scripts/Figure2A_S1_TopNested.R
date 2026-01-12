rm(list=ls()); gc()

library(phyloseq)
library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)
library(writexl)
library(fantaxtic)
library(ggtext)

#####
# 1) Reading in data
#####

# ****ABUNDANCE DATA
abund0 = read.table("novogene_analysis/featureTable.sample.total.absolute.txt",
                    header = TRUE, comment.char = "", check.names = FALSE)

# remove names and phylogeny
abund1 = abund0[, -c(1, NCOL(abund0))]

# scale to 100000 per sample
abund2 = apply(abund1, MARGIN = 2, function(x) round(100000 * (x / sum(x))))

# final set
abund = abund2


# ****METADATA
meta0 = data.frame(read_excel("~/ColiShield/Pearo_16S/Pearo_16S/Meta/All stool samples adj.xlsx"))
rownames(meta0) = meta0$sample
meta = meta0


# ****TAXONOMY
taxa0 = str_split_fixed(string = abund0[, NCOL(abund0)], pattern = ";", n = 7)
taxa1 = gsub(pattern = ".__", replacement = "", x = taxa0)
colnames(taxa1) = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
taxa = taxa1


#####
# 2) Build phyloseq object
#####

ps <- phyloseq(
  otu_table(abund, taxa_are_rows = TRUE),
  sample_data(meta),
  tax_table(taxa)
)

table(meta$study.day)

table(StudyDay = meta$study.day,
      Group    = meta$groupe)

ps.sub = subset_samples(ps,
                        groupe %in% c(1,2,3) &
                          study.day %in% c(1,8,15,29,54) &
                          Lysis %in% c("FP"))

ps.sub.genus = tax_glom(ps.sub, taxrank = "Genus")
Genus_ps_sub = psmelt(ps.sub.genus)


#####
# 3) Nested top taxa + plots
#####

# Nested top taxa (Phylum -> Family)
top_nested <- nested_top_taxa(ps.sub.genus,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Family",
                              n_top_taxa = 5,
                              n_nested_taxa = 3)

palette <- c(
  Bacteroidota     = "#8D0819",
  Euryarchaeota    = "#999300",
  Firmicutes       = "#016FB9",
  Proteobacteria   = "#7E3267",
  Spirochaetota    = "#7D3200",
  Actinobacteriota = "black"
)

p <- plot_nested_bar(ps_obj = top_nested$ps_obj,
                     top_level = "Phylum",
                     nested_level = "Family",
                     palette = palette)

day_labels <- c("Day 1", "Day 8", "Day 15", "Day 29", "Day 54")
names(day_labels) <- c("1", "8", "15", "29", "54")

group_labels <- c("1" = "Control", "2" = "VHH1", "3" = "VHH1+2")

f <- p +
  facet_wrap(~ study.day + groupe,
             nrow = 1,
             scales = "free_x",
             labeller = labeller(study.day = day_labels, groupe = group_labels)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_markdown(size = 22, colour = "black", vjust = 0.25, face = 'plain'),
        strip.text = element_text(size = 10, face = "bold"),
        axis.line = element_line(linetype = "solid", linewidth = 1.0),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22),
        legend.position = "bottom") +
  labs(title = "Top 5 Taxa", x = "Taxa", y = "Relative Abundance")
print(f)

ggsave("fig2barplotlegend.pdf", plot = f, width = 8.78, height = 5.18, units = "in")


# Nested top taxa (Family -> Genus)
top_nested_g <- nested_top_taxa(ps.sub.genus,
                                top_tax_level = "Family",
                                nested_tax_level = "Genus",
                                n_top_taxa = 5,
                                n_nested_taxa = 3)

palette <- c(
  Prevotellaceae   = "#A71D31",
  Ruminococcaceae  = "#3A405A",
  Clostridiaceae   = "#3F88C5",
  Lachnospiraceae  = "#5603AD",
  Lactobacillaceae = "#7CFEF0"
)

g <- plot_nested_bar(ps_obj = top_nested_g$ps_obj,
                     top_level = "Family",
                     nested_level = "Genus",
                     palette = palette)

g + facet_wrap(~ study.day + groupe,
               nrow = 1,
               scales = "free_x",
               labeller = labeller(study.day = day_labels, groupe = group_labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom") +
  labs(title = "Top 5 Taxa", x = "Taxa", y = "Relative Abundance")



