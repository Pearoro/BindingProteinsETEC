library(dplyr)
library(phyloseq)
library(writexl)

psdf <- psmelt(top_nested$ps_obj)

family_stats <- psdf %>%
  group_by(Sample, study.day, groupe) %>%
  mutate(RelAb = Abundance / sum(Abundance)) %>%   # per-sample relative abundance (includes Other)
  ungroup() %>%
  group_by(study.day, groupe, Family) %>%
  summarise(
    Percentage = mean(RelAb) * 100,
    SD = sd(RelAb) * 100,
    SE = (sd(RelAb) / sqrt(n())) * 100,
    MIN = min(RelAb) * 100,
    MAX = max(RelAb) * 100,
    n = n(),
    .groups = "drop"
  )

print(family_stats)
write_xlsx(family_stats, "family_stats_meanRA_.xlsx")

psdf <- psmelt(top_nested_g$ps_obj)

genus_stats <- psdf %>%
  group_by(Sample, study.day, groupe) %>%
  mutate(RelAb = Abundance / sum(Abundance)) %>%   # per-sample relative abundance (includes Other)
  ungroup() %>%
  group_by(study.day, groupe, Genus) %>%
  summarise(
    Percentage = mean(RelAb) * 100,
    SD = sd(RelAb) * 100,
    SE = (sd(RelAb) / sqrt(n())) * 100,
    MIN = min(RelAb) * 100,
    MAX = max(RelAb) * 100,
    n = n(),
    .groups = "drop"
  )

print(genus_stats)
write_xlsx(genus_stats, "genus_stats_meanRA_.xlsx")

