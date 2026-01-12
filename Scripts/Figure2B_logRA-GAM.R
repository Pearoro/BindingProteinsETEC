##############################
## Setup
##############################

# Load necessary libraries
library(phyloseq)
library(dplyr)
library(stringr)
library(readxl)
library(openxlsx)
library(writexl)
library(tidyr)
library(ggplot2)
library(mgcv)  # 

# Clear workspace
rm(list = ls())
gc()

##############################
## 1) Read data
##############################

# Abundance data
abund0 <- read.table(
  "novogene_analysis/featureTable.sample.total.absolute.txt",
  header       = TRUE,
  comment.char = ""
)

# Normalize/RA
abund1 <- abund0[, -c(1, NCOL(abund0))]
abund2 <- apply(abund1, MARGIN = 2, function(x) round(100000 * (x / sum(x))))
abund  <- abund2  # final OTU table

# Metadata
meta0 <- read_excel("~/ColiShield/Pearo_16S/Pearo_16S/Meta/All stool samples adj.xlsx") |>
  as.data.frame()
rownames(meta0) <- meta0$sample
meta <- meta0

# Taxonomy
taxa0 <- str_split_fixed(
  string  = abund0[, NCOL(abund0)],
  pattern = ";",
  n       = 7
)
taxa1 <- gsub(pattern = ".__", replacement = "", x = taxa0)
colnames(taxa1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa <- taxa1

##############################
## 2) Build phyloseq object
##############################

ps <- phyloseq(
  otu_table(abund, taxa_are_rows = TRUE),
  sample_data(meta),
  tax_table(taxa)
)

ps.sub <- subset_samples(
  ps,
  groupe    %in% c(1, 2, 3) &
    study.day %in% c(1, 8, 15, 29, 54) &
    Lysis     %in% c("FP")
)



##############################
## 3) Global labels, colors, theme
##############################

# Day labels / levels
day_levels <- c("1", "8", "15", "29", "54")
day_breaks <- as.numeric(day_levels)
day_labels <- setNames(
  c("1", "8", "15", "29", "54"),
  day_breaks
)

# Group labels & colors
group_labels <- setNames(
  c("Control", "VHH1", "VHH1+2"),
  c("1", "2", "3")
)

group_colors <- c(
  "1" = "#2E78DB",
  "2" = "#FF9505",
  "3" = "#DB402E"
)

# Y-limits (adjust if needed)
y_limits_family <- c(0, 6)
y_limits_genus  <- c(0, 5)
y_limits_phyla  <- c(0, 6)

# Common y-axis label
y_label_RA <- expression(log[10] * "(RA)")

# Shared theme for ALL plots
theme_myplot <- theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(
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
    legend.position = "none",
    panel.background = element_blank(),
    plot.background  = element_blank()
  )

##############################
## 4) Prepare phyloseq melts
##############################

# Subset at different taxonomic levels
ps.sub.fam <- tax_glom(ps.sub, taxrank = "Family")
ps.sub.ph  <- tax_glom(ps.sub, taxrank = "Phylum")
ps.sub.gen <- tax_glom(ps.sub, taxrank = "Genus")

# Melt to data frames
psdf_family <- psmelt(ps.sub.fam)
psdf_phyla  <- psmelt(ps.sub.ph)
psdf_genus  <- psmelt(ps.sub.gen)

# Common point position (jitter+dodge)
pos <- position_jitterdodge(
  dodge.width  = 4,
  jitter.width = 0.25,
  jitter.height = 0
)

# Ensure factors, numeric day, log-abundance are consistent

## Family
psdf_family <- psdf_family |>
  mutate(
    groupe        = factor(groupe, levels = c("1", "2", "3")),
    study.day     = factor(as.character(study.day), levels = day_levels),
    study.day.num = as.numeric(as.character(study.day)),  # numeric for spacing
    logAbundance  = log10(Abundance + 1)
  )

## Genus
psdf_genus <- psdf_genus |>
  mutate(
    groupe        = factor(groupe, levels = c("1", "2", "3")),
    study.day     = factor(as.character(study.day), levels = day_levels),
    study.day.num = as.numeric(as.character(study.day)),
    logAbundance  = log10(Abundance + 1)
  )

## Phyla
psdf_phyla <- psdf_phyla |>
  mutate(
    groupe        = factor(groupe, levels = c("1", "2", "3")),
    study.day     = factor(as.character(study.day), levels = day_levels),
    study.day.num = as.numeric(as.character(study.day)),
    logAbundance  = log10(Abundance + 1)
  )

##############################
## 5) FAMILY plots (GAM + continuous days)
##############################

df_family <- psdf_family

plot_family <- function(fam) {
  sub <- df_family |> filter(Family == fam)
  
  ggplot(sub, aes(
    x     = study.day.num,
    y     = logAbundance,
    color = groupe,
    group = groupe
  )) +
    geom_point(
      position = pos,
      size     = 3,
      alpha    = 0.6
    ) +
    geom_smooth(
      method  = "gam",
      formula = y ~ s(x, k = 5),
      se      = TRUE,
      linetype = "dashed",
      na.rm   = TRUE,
      fill    = "grey70",
      alpha   = 0.3
    ) +
    coord_cartesian(ylim = y_limits_family) +
    scale_x_continuous(
      breaks = day_breaks,
      labels = day_labels,
      expand = expansion(add = c(0.5, 0.5))
    ) +
    scale_color_manual(values = group_colors, labels = group_labels) +
    labs(
      title = paste("Dynamic Changes in", fam),
      x     = "",
      y     = y_label_RA,
      color = "Group"
    ) +
    theme_myplot
}

families <- c("Prevotellaceae", "Lactobacillaceae", "Enterobacteriaceae")
plots_family <- setNames(lapply(families, plot_family), families)

##############################
## 6) GENUS plots (GAM + continuous days)
##############################

df_genus <- psdf_genus

plot_genus <- function(gen) {
  sub <- df_genus |> filter(Genus == gen)
  
  ggplot(sub, aes(
    x     = study.day.num,
    y     = logAbundance,
    color = groupe,
    group = groupe
  )) +
    geom_point(
      position = pos,
      size     = 3,
      alpha    = 0.6
    ) +
    geom_smooth(
      method  = "gam",
      formula = y ~ s(x, k = 5),
      se      = TRUE,
      linetype = "dashed",
      na.rm   = TRUE,
      fill    = "grey70",
      alpha   = 0.3
    ) +
    coord_cartesian(ylim = y_limits_genus) +
    scale_x_continuous(
      breaks = day_breaks,
      labels = day_labels,
      expand = expansion(add = c(0.5, 0.5))
    ) +
    scale_color_manual(values = group_colors, labels = group_labels) +
    labs(
      title = paste("Dynamic Changes in", gen),
      x     = "",
      y     = y_label_RA,
      color = "Group"
    ) +
    theme_myplot
}

genera <- c(
  "Escherichia-Shigella",
  "Prevotella_9",
  "Lactobacillus",
  "Limosilactobacillus",
  "Clostridium_sensu_stricto_1",
  "Ligilactobacillus",
  "Blautia",
  "Ruminococcus",
  "Coprococcus",
  "Treponema",
  "Roseburia",
  "Campylobacter"
)

# Example prints
print(plot_genus("Escherichia-Shigella"))
print(plot_genus("Prevotella_9"))


plots_genus <- setNames(lapply(genera, plot_genus), genera)

##############################
## 7) PHYLA plots (GAM + continuous days)
##############################

# Firmicutes
plot_data_phyla_F <- psdf_phyla |> filter(Phylum == "Firmicutes")

f <- ggplot(plot_data_phyla_F, aes(
  x     = study.day.num,
  y     = logAbundance,
  color = groupe,
  group = groupe
)) +
  geom_point(
    position = pos,
    size     = 3,
    alpha    = 0.6
  ) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    linetype = "dashed",
    na.rm   = TRUE,
    fill    = "grey70",
    alpha   = 0.3
  ) +
  coord_cartesian(ylim = y_limits_phyla) +
  scale_x_continuous(
    breaks = day_breaks,
    labels = day_labels,
    expand = expansion(add = c(0.5, 0.5))
  ) +
  scale_color_manual(values = group_colors, labels = group_labels) +
  labs(
    title = "Dynamic Changes in Firmicutes",
    x     = "",
    y     = y_label_RA,
    color = "Group"
  ) +
  theme_myplot

# Bacteroidota
plot_data_phyla_B <- psdf_phyla |> filter(Phylum == "Bacteroidota")

b <- ggplot(plot_data_phyla_B, aes(
  x     = study.day.num,
  y     = logAbundance,
  color = groupe,
  group = groupe
)) +
  geom_point(
    position = pos,
    size     = 3,
    alpha    = 0.6
  ) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    linetype = "dashed",
    na.rm   = TRUE,
    fill    = "grey70",
    alpha   = 0.3
  ) +
  coord_cartesian(ylim = y_limits_phyla) +
  scale_x_continuous(
    breaks = day_breaks,
    labels = day_labels,
    expand = expansion(add = c(0.5, 0.5))
  ) +
  scale_color_manual(values = group_colors, labels = group_labels) +
  labs(
    title = "Dynamic Changes in Bacteroidota",
    x     = "",
    y     = y_label_RA,
    color = "Group"
  ) +
  theme_myplot

##############################
## 8) Save plots as separate PDFs
##############################

dir.create("plots", showWarnings = FALSE)

## Families – one PDF per family
for (nm in names(plots_family)) {
  ggsave(
    filename = file.path("plots", paste0("family_", nm, ".pdf")),
    plot     = plots_family[[nm]],
    width    = 8,
    height   = 8
  )
}

## Genera – one PDF per genus
for (nm in names(plots_genus)) {
  ggsave(
    filename = file.path("plots", paste0("genus_", nm, ".pdf")),
    plot     = plots_genus[[nm]],
    width    = 8,
    height   = 8
  )
}

## Phyla – one PDF per phylum plot
ggsave("plots/phyla_Firmicutes.pdf",  f, width = 8, height = 6)
ggsave("plots/phyla_Bacteroidota.pdf", b, width = 8, height = 6)

##############################
## 9) Optional sanity checks
##############################

unique(psdf_genus$Genus)
sum(is.na(psdf_genus$groupe))
