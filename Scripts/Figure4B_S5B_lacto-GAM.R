library(readxl)
library(dplyr)
library(emmeans)
library(ggplot2)
library(writexl)
library(tidyr)
library(mgcv)

A <- read_excel("C:\\Users\\pearoro\\OneDrive - Danmarks Tekniske Universitet\\Documents\\ColiShield\\Pearo_16S\\Pearo_16S\\Lactod.xlsx")

A$Sample <- gsub("-1[0]+","", A$Sample)

B <- read.csv("Coli_16s.csv", header = TRUE, sep = ";")

intersect(A$Sample, B$Sample)
setdiff(A$Sample, B$Sample)
setdiff(B$Sample, A$Sample)

C <- merge(A, B, by.x ="Sample", by.y = "Sample", all = FALSE)

C %>% count(Sample) %>% filter(n > 4) 
C %>% count(Sample) %>% filter(n > 5) 
C %>% count(Sample) %>% filter(n < 4) 

C$Conc     <- as.numeric(C$Conc)
C$theo_16S <- as.numeric(C$Conc)

# copies / mg 
C$Conc_mg <- C$Conc / 2

### PLOTS ###

plot_data <- C

# global day settings
day_levels <- c("1", "8", "15", "29", "54")
day_breaks <- as.numeric(day_levels)
day_labels <- setNames(
  c("1", "8", "15", "29", "54"),
  day_breaks
)

group_labels <- c("1" = "Control", "2" = "VHH1", "3" = "VHH1+2")
group_colors <- c("1" = "#2E78DB", "2" = "#FF9505", "3" = "#DB402E")

# Ensure columns are in the correct format AND LOG for visualization
plot_data <- plot_data %>%
  mutate(
    study.day.y   = factor(as.character(study.day.y), levels = day_levels),
    study.day.num = as.numeric(as.character(study.day.y)),
    logConc_mg    = log10(Conc_mg + 1),
    logtheo_mg    = log10(theo_mg + 1),
    groupe.y      = factor(groupe.y, levels = c("1","2","3"))
  )

pos          <- position_jitterdodge(dodge.width = 5,
                                     jitter.width = 0.25,
                                     jitter.height = 0)
y_limits     <- c(0, 8)
y_label_expr <- expression(log[10]*"(cp/mg)")

################
# 16S (E. coli)
################

df16S <- subset(plot_data, Target == "16S")

e <- ggplot(df16S, aes(x = study.day.num, y = logConc_mg,
                       color = groupe.y, group = groupe.y)) +
  geom_point(position = pos, size = 3, alpha = 0.6) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    linetype = "dashed",
    na.rm   = TRUE,
    fill    = "grey70",
    alpha   = 0.3
  ) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels,
                     expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(breaks = 0:7, expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = y_limits) +
  scale_color_manual(values = group_colors, breaks = names(group_colors)) +
  labs(title = "16S", x = "", y = y_label_expr, color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,
                               colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background  = element_blank()
  )

print(e)
ggsave("plots/EdPCR.pdf",  e , width = 8, height = 8)

# with legend
el <- ggplot(df16S, aes(x = study.day.num, y = logConc_mg,
                        color = groupe.y, group = groupe.y)) +
  geom_point(position = pos, size = 3, alpha = 0.6) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    linetype = "dashed",
    na.rm   = TRUE,
    fill    = "grey70",
    alpha   = 0.3
  ) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels,
                     expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(breaks = 0:7, expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = y_limits) +
  scale_color_manual(values = group_colors,
                     breaks = names(group_labels),
                     labels = group_labels) +
  labs(title = expression(italic('E. coli')~"(16S rRNA gene)"),
       x     = "",
       y     = y_label_expr,
       color = "Group")+
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,
                               colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    plot.background  = element_blank()
  )
print(el)
ggsave("plots/EdPCR_legend.pdf",  el , width = 8, height = 8)

############
# ybbw
############

dfyb <- subset(plot_data, Target == "ybbw")

y <- ggplot(dfyb, aes(x = study.day.num, y = logConc_mg,
                      color = groupe.y, group = groupe.y)) +
  geom_point(position = pos, size = 3, alpha = 0.6) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    linetype = "dashed",
    na.rm   = TRUE,
    fill    = "grey70",
    alpha   = 0.3
  ) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels,
                     expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(breaks = 0:7, expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = y_limits) +
  scale_color_manual(values = group_colors, breaks = names(group_colors)) +
  labs(title = "ybbw", x = "", y = y_label_expr, color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,
                               colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background  = element_blank()
  )
print(y)
ggsave("plots/EybbwdPCR.pdf",  y, width = 8, height = 8)

# ybbw with legend
yl <- ggplot(dfyb, aes(x = study.day.num, y = logConc_mg,
                       color = groupe.y, group = groupe.y)) +
  geom_point(position = pos, size = 3, alpha = 0.6) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    linetype = "dashed",
    na.rm   = TRUE,
    fill    = "grey70",
    alpha   = 0.3
  ) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels,
                     expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(breaks = 0:7, expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = y_limits) +
  scale_color_manual(values = group_colors,
                     breaks = names(group_labels),
                     labels = group_labels) +
  labs(
    title = expression(italic('E. coli (ybbw)')),
    x     = "",
    y     = y_label_expr,
    color = "Group"
  )+
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,
                               colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    plot.background  = element_blank()
  )
print(yl)
ggsave("plots/EybbwdPCR_legend.pdf", yl, width = 8, height = 8)

############
# 16L
############

df16L <- subset(plot_data, Target == "16L")

LS <- ggplot(df16L, aes(x = study.day.num, y = logConc_mg,
                        color = groupe.y, group = groupe.y)) +
  geom_point(position = pos, size = 3, alpha = 0.6) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    linetype = "dashed",
    na.rm   = TRUE,
    fill    = "grey70",
    alpha   = 0.3
  ) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels,
                     expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(breaks = 0:8, expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = y_limits) +
  scale_color_manual(values = group_colors, breaks = names(group_colors)) +
  labs(title = "16S", x = "", y = y_label_expr, color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,
                               colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background  = element_blank()
  )
print(LS)
ggsave("plots/EdPCR_Lreuteri.pdf",  LS , width = 8, height = 8)

LSl <- ggplot(df16L, aes(x = study.day.num, y = logConc_mg,
                         color = groupe.y, group = groupe.y)) +
  geom_point(position = pos, size = 3, alpha = 0.6) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    linetype = "dashed",
    na.rm   = TRUE,
    fill    = "grey70",
    alpha   = 0.3
  ) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels,
                     expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(breaks = 0:8, expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = y_limits) +
  scale_color_manual(values = group_colors,
                     breaks = names(group_labels),
                     labels = group_labels) +
  labs(title = expression(italic('L. reu')~"(16S rRNA gene)"),
       x     = "",
       y     = y_label_expr,
       color = "Group")+
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,
                               colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    plot.background  = element_blank()
  )
print(LSl)
ggsave("plots/EdPCR_legend_Lreuteri.pdf",  LSl , width = 8, height = 8)


# factors for modelling
df16LM <- df16L %>%
  mutate(
    group = factor(groupe.x),
    study.day = factor(study.day.x)
  )

LM_REU <- lm(logConc_mg ~ group * study.day, data = df16LM)
par(mfrow = c(2, 2))
plot(LM_REU)

emm_reu <- emmeans(LM_REU, ~ group | study.day, adjust = "sidak")
contrast(emm_reu, method = "tukey")
plot(emm_reu, comparisons = TRUE)



############
# purB
############

dfpur <- subset(plot_data, Target == "purB")

amy <- ggplot(dfpur, aes(x = study.day.num, y = logConc_mg,
                         color = groupe.y, group = groupe.y)) +
  geom_point(position = pos, size = 3, alpha = 0.6) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    linetype = "dashed",
    na.rm   = TRUE,
    fill    = "grey70",
    alpha   = 0.3
  ) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels,
                     expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(breaks = 0:8, expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = y_limits) +
  scale_color_manual(values = group_colors, breaks = names(group_colors)) +
  labs(title = "purB", x = "", y = y_label_expr, color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,
                               colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background  = element_blank()
  )
print(amy)
ggsave("plots/EdPCR_purB.pdf",  amy , width = 8, height = 8)

amyL <- ggplot(dfpur, aes(x = study.day.num, y = logConc_mg,
                          color = groupe.y, group = groupe.y)) +
  geom_point(position = pos, size = 3, alpha = 0.6) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    linetype = "dashed",
    na.rm   = TRUE,
    fill    = "grey70",
    alpha   = 0.3
  ) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels,
                     expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(breaks = 0:8, expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = y_limits) +
  scale_color_manual(
    values = group_colors,
    breaks = names(group_labels),
    labels = group_labels
  ) +
  labs(title = expression(italic('L. amy')~"(16S rRNA gene)"),
       x     = "",
       y     = y_label_expr,
       color = "Group")+
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,
                               colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    plot.background  = element_blank()
  )
print(amyL)
ggsave("plots/EdPCR_legend_purB.pdf",  amyL , width = 8, height = 8)



# factors for modelling
dfPURLM <- dfpur %>%
  mutate(
    group = factor(groupe.x),
    study.day = factor(study.day.x)
  )

LM_amy <- lm(logConc_mg ~ group * study.day, data = dfPURLM)
par(mfrow = c(2, 2))
plot(LM_amy)

emm_amy <- emmeans(LM_amy, ~ group | study.day, adjust = "sidak")
contrast(emm_amy, method = "tukey")
plot(emm_amy, comparisons = TRUE)

### trend figure 

plot_data_trend_lacto <- plot_data %>%
  filter(Target %in% c("16S", "ybbw", "purB", "16L")) %>%
  mutate(
    y_mg = if_else(
      Target %in% c("16S", "16L"),
      logtheo_mg,    # use this for 16S & 16L
      logConc_mg     # use this for ybbw & purB
    )
  )

# day and group settings 
day_levels <- c("1", "8", "15", "29", "54")
day_breaks <- as.numeric(day_levels)
day_labels <- setNames(c("1", "8", "15", "29", "54"),
                       day_breaks)
group_labels <- c("1" = "Control", "2" = "VHH1", "3" = "VHH1+2")

# colors and labels
target_colors <- c(
  "16S" = "#8E8070",
  "16L" = "#111D4A",
  "ybbw" = "#9278C3",
  "purB" = "#00AFB5"
)
target_labels <- c(
  "16S"  = expression(italic('E. coli')~"(16S rRNA gene)"),
  "16L"  = expression(italic('L. reuteri')~"(16S rRNA gene)"),
  "ybbw" = expression(italic('E. coli')~"(ybbw)"),
  "purB" = expression(italic('E. coli')~"(purB)")
)

plot_data_trend_lacto <- plot_data_trend_lacto %>%
  mutate(
    study.day.y   = factor(as.character(study.day.y), levels = day_levels),
    study.day.num = as.numeric(as.character(study.day.y)),   # numeric for spacing
    groupe.y      = factor(groupe.y, levels = c("1","2","3")),
    Target        = factor(Target, levels = c("16S", "ybbw", "purB", "16L"))
  )

# y-limits
y_limits <- c(0, 8)

lacto <- ggplot(plot_data_trend_lacto,
                aes(x = study.day.num,
                    y = y_mg,
                    color = Target,
                    fill  = Target,
                    group = Target)) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),   # GAM instead of loess
    se      = TRUE,
    size    = 1.5,
    alpha   = 0.25,
    na.rm   = TRUE
  ) +
  facet_grid(. ~ groupe.y, labeller = labeller(groupe.y = group_labels)) +
  scale_color_manual(values = target_colors, labels = target_labels, name = "Series") +
  scale_fill_manual(values = target_colors, labels = target_labels, guide = "none") +
  scale_x_continuous(
    breaks = day_breaks,
    labels = day_labels,
    expand = expansion(add = c(0.5, 0.5))
  ) +
  scale_y_continuous(limits = y_limits,
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Trends",
    x = "Day",
    y = y_label_expr   # mixed source, but still log10(cp/mg)
  ) +
  theme_classic(base_size = 30) +
  theme(
    axis.text.x = element_text(hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(face = "bold", size = 18),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none"
  )

print(lacto)
ggsave("plots/lacto_legend.pdf", lacto, width = 10, height = 8)




library(dplyr)

summary_tbl <- plot_data %>%
  group_by(Target, groupe.y, study.day.y) %>%
  summarise(
    n = sum(!is.na(logConc_mg)),
    mean_log = mean(logConc_mg, na.rm = TRUE),
    sd_log   = sd(logConc_mg, na.rm = TRUE),
    .groups = "drop"
  )

summary_tbl

library(writexl)

write_xlsx(summary_tbl, "summary_tbl.xlsx")

library(dplyr)

summary_tbltheo <- plot_data %>%
  group_by(Target, groupe.y, study.day.y) %>%
  summarise(
    n = sum(!is.na(logtheo_mg)),
    mean_log = mean(logtheo_mg, na.rm = TRUE),
    sd_log   = sd(logtheo_mg, na.rm = TRUE),
    .groups = "drop"
  )

summary_tbltheo

write_xlsx(summary_tbltheo, "summary_tblteho.xlsx")



