rm(list=ls()); gc()

library(phyloseq)
library(stringr)
library(readxl)
library(vegan)

# Define color gradient function
color.gradient <- function(x, colors=c("#2E78DB","#B9DB2E","#DB402E"), colsteps=100) {
  return(colorRampPalette(colors)(colsteps)[findInterval(x, seq(min(x), max(x), length.out=colsteps))])
}

#####
# 1) Reading in data
#####

# Abundance Data
abund0 <- read.table("novogene_analysis/featureTable.sample.total.absolute.txt", header=TRUE, comment.char="")
abund1 <- abund0[,-c(1, NCOL(abund0))]
abund2 <- apply(abund1, MARGIN=2, function(x) round(100000 * (x / sum(x))))
abund <- abund2

# Metadata
meta0 <- data.frame(read_excel("~/ColiShield/Pearo_16S/Pearo_16S/Meta/All stool samples adj.xlsx"))
rownames(meta0) <- meta0$sample
meta <- meta0

# Taxonomy
taxa0 <- str_split_fixed(string=abund0[,NCOL(abund0)], pattern=";", n=7)
taxa1 <- gsub(pattern=".__", replacement="", x=taxa0)
colnames(taxa1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa <- taxa1

#####
# 2) Build phyloseq object
#####

ps <- phyloseq(otu_table(abund, taxa_are_rows=TRUE), sample_data(meta), tax_table(taxa))
table(meta$study.day)

ps.sub <- subset_samples(ps, groupe %in% c(1,2,3) & study.day %in% c(1,8,15,29,54) & Lysis %in% c("FP"))
ps.sub.filt <- filter_taxa(ps.sub, function(x) mean(x) > 1, TRUE)

#####
# 3) Make ordination
#####

ps.sub.ord <- ordinate(ps.sub.filt, method="NMDS")

par(mar=c(2.5, 2.5, 2.5, 0.5), mgp=c(1.4, 0.4, 0))
plot(ps.sub.ord$points, pch=ps.sub.filt@sam_data$groupe + 21, bg=color.gradient(log(ps.sub.filt@sam_data$study.day)), cex=1.5, xlim=c(-1.5, 1.5))

# Create a gradient of colors matching the values used in the plot
study.days <- log(ps.sub.filt@sam_data$study.day)
unique_days <- sort(unique(study.days))
color_values <- color.gradient(unique_days)

# Add a legend for the color gradient
legend("topright", 
       legend=round(exp(unique_days), 1),  # Convert back from log scale to original study day values
       fill=color_values, 
       title="Study Day", 
       cex=1, 
       border=NA, 
       bty="n")

#####
# 4) Make PERMANOVA
#####

otu.filt <- wisconsin(sqrt(t(otu_table(ps.sub.filt))))
DIST <- vegdist(otu.filt)
DAY <- ps.sub.filt@sam_data$study.day
GROUPE <- ps.sub.filt@sam_data$groupe

# Include the by argument to test terms separately
ado <- adonis2(formula=DIST ~ DAY * GROUPE, permutations=9999, by="terms")
print(ado)

#####
# 5) Make pairwise
#####

DAYS <- unique(ps.sub.filt@sam_data$study.day)
for (i in DAYS) {
  ps.sub.filt.day <- subset_samples(ps.sub.filt, study.day == i)
  ps.sub.ord.day <- ordinate(ps.sub.filt.day, method="NMDS", verbose=0)
  
  # Define colors for each group
  group_colors <- c("#2E78DB", "#B9DB2E", "#DB402E")  # Adjust colors as needed
  
  plot(ps.sub.ord.day$points, 
       pch=ps.sub.filt.day@sam_data$groupe + 21, 
       col="black",  # Outline color
       bg=group_colors[ps.sub.filt.day@sam_data$groupe],  # Fill color
       cex=2,  # Increase the size of the symbols
       main=i)
  
  otu.filt.day <- wisconsin(sqrt(t(otu_table(ps.sub.filt.day))))
  DIST <- vegdist(otu.filt.day)
  DAY <- ps.sub.filt.day@sam_data$study.day
  GROUPE <- ps.sub.filt.day@sam_data$groupe
  
  print(i)
  ado <- adonis2(DIST ~ GROUPE, permutations=9999, by="terms")
  print(ado)
  
  # Extract p-value and R²
  p_value <- ado$aov.tab$`Pr(>F)`[1]
  r_squared <- ado$aov.tab$R2[1]
  
  # Add legend with p-value and R²
  legend("topright", 
         legend=c(paste("p-value:", format(p_value, digits=3)), 
                  paste("R²:", format(r_squared, digits=3))),
         bty="n")
}
