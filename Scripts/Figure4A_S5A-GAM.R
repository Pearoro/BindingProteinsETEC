library(readxl)
library(dplyr)
library(emmeans)
library(ggplot2)
library(writexl)
library(tidyr)
library(mgcv)

A <- read_excel("C:\\Users\\pearoro\\OneDrive - Danmarks Tekniske Universitet\\Documents\\ColiShield\\Pearo_16S\\Pearo_16S\\E COLId.xlsx")

A$Sample <- gsub("-1[0]+","", A$Sample)

B <- read.csv("Coli_16s.csv", header = TRUE, sep = ";")

intersect(A$Sample, B$Sample)
setdiff(A$Sample, B$Sample)
setdiff(B$Sample, A$Sample)

C <- merge(A, B, by.x ="Sample", by.y = "Sample", all = FALSE)

C %>% count(Sample) %>% filter(n > 5)

C$Conc <- as.numeric(C$Conc)
C$Conc_mg <- C$Conc / 2

### PLOTS ###

plot_data <- C

plot_data$logConc_mg <- log10(plot_data$Conc_mg + 1)
plot_data$groupe.y   <- as.factor(plot_data$groupe.y)

# day / group labels & colors
day_levels   <- c("1", "8", "15", "29", "54")
day_breaks   <- as.numeric(day_levels)
day_labels   <- setNames(c("1", "8", "15", "29", "54"), day_breaks)
group_colors <- c("1" = "#2E78DB", "2" = "#FF9505", "3" = "#DB402E")
group_labels <- c("1" = "Control", "2" = "VHH1", "3" = "VHH1+2")

y_limits     <- c(0, 8)
y_label_expr <- expression(log[10]*"(cp/mg)")

# make factor + numeric day for plotting
plot_data <- plot_data %>%
  mutate(
    study.day.y  = factor(study.day.y, levels = day_levels),
    study.day.num = as.numeric(as.character(study.day.y)),
    groupe.y     = factor(groupe.y, levels = c("1","2","3"))
  )

# point position
pos <- position_jitterdodge(
  dodge.width  = 5,
  jitter.width = 0.25,
  jitter.height = 0
)

#### 16S ####
df16S <- subset(plot_data, Target == "16S")

e <- ggplot(df16S,
            aes(x = study.day.num, y = logConc_mg,
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
  coord_cartesian(ylim = y_limits) +        # y_limits <- c(0, 8)
  scale_color_manual(values = group_colors,
                     breaks = names(group_colors),
                     drop   = FALSE) +
  labs(title = "16S", x = "", y = y_label_expr, color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(hjust = 1, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background  = element_blank()
  )



print(e)
ggsave("plots/EdPCR.pdf",  e , width = 8, height = 8)

# 16S with legend
el <- ggplot(df16S,
             aes(x = study.day.num, y = logConc_mg,
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
  scale_color_manual(
    values = group_colors,
    breaks = names(group_labels),
    labels = group_labels
  ) +
  labs(title = expression(italic('E. coli')~"(16S rRNA gene)"),
       x     = "",
       y     = y_label_expr,
       color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text( hjust = 1, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    plot.background  = element_blank()
  )


print(el)
ggsave("plots/EdPCR_legend.pdf",  el , width = 8, height = 8)

#### ybbw ####
dfyb <- subset(plot_data, Target == "ybbw")

y <- ggplot(dfyb,
            aes(x = study.day.num, y = logConc_mg,
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
  scale_color_manual(values = group_colors, breaks = names(group_colors), drop = FALSE) +
  labs(title = "ybbw", x = "", y = y_label_expr, color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(hjust = 1, colour = "black", face = "bold"),
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
yl <- ggplot(dfyb,
             aes(x = study.day.num, y = logConc_mg,
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
  scale_color_manual(
    values = group_colors,
    breaks = names(group_labels),
    labels = group_labels
  ) +
  labs(
    title = expression(italic('E. coli (ybbw)')),
    x     = "",
    y     = y_label_expr,
    color = "Group"
  ) +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(hjust = 1, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    plot.background  = element_blank()
  )

print(yl)
ggsave("plots/EybbwdPCR_legend.pdf",  yl, width = 8, height = 8)

#### LT ####
dfLT <- subset(plot_data, Target == "LT")

lt <- ggplot(dfLT,
             aes(x = study.day.num, y = logConc_mg,
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
  scale_color_manual(values = group_colors, breaks = names(group_colors), drop = FALSE) +
  labs(title = "LT", x = "", y = y_label_expr, color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(hjust = 1, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background  = element_blank()
  )
print(lt)
ggsave("plots/LTdPCR.pdf", lt , width = 8, height = 8)

# LT with legend
ltl <- ggplot(dfLT,
              aes(x = study.day.num, y = logConc_mg,
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
  scale_color_manual(
    values = group_colors,
    breaks = names(group_labels),
    labels = group_labels
  ) +
  labs(title = expression("ETEC ("*italic(lt)*")"),
       x     = "",
       y     = y_label_expr,
       color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    plot.background  = element_blank()
  )
print(ltl)
ggsave("plots/LT_legenddPCR.pdf", ltl , width = 8, height = 8)

# factors for modelling
dfLT <- dfLT %>%
  mutate(
    group = factor(groupe.x),
    study.day = factor(study.day.x)
  )

LM_LT <- lm(logConc_mg ~ group * study.day, data = dfLT)
par(mfrow = c(2, 2))
plot(LM_LT)

emm_LT <- emmeans(LM_LT, ~ group | study.day, adjust = "sidak")
contrast(emm_LT, method = "tukey")
plot(emm_LT, comparisons = TRUE)

#### F4 ####
dfF4 <- subset(plot_data, Target == "F4")

F4 <- ggplot(dfF4,
             aes(x = study.day.num, y = logConc_mg,
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
  scale_color_manual(values = group_colors, breaks = names(group_colors), drop = FALSE) +
  labs(title = "F4", x = "", y = y_label_expr, color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(hjust = 1, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background  = element_blank()
  )
print(F4)
ggsave("plots/F4dPCR.pdf", F4 , width = 8, height = 8)

F4l <- ggplot(dfF4,
              aes(x = study.day.num, y = logConc_mg,
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
  scale_color_manual(
    values = group_colors,
    breaks = names(group_labels),
    labels = group_labels
  ) +
  labs( title = expression("ETEC F4 ("*italic(faeG)*")"),
        x     = "",
        y     = y_label_expr,
        color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    plot.background  = element_blank()
  )
print(F4l)
ggsave("plots/F4_legenddPCR.pdf", F4l , width = 8, height = 8)

dfF4 <- dfF4 %>%
  mutate(
    group = factor(groupe.x),
    study.day = factor(study.day.x)
  )

LM_F4 <- lm(logConc_mg ~ group * study.day, data = dfF4)
par(mfrow = c(2, 2))
plot(LM_F4)

emm_F4 <- emmeans(LM_F4, ~ group | study.day, adjust = "sidak")
contrast(emm_F4, method = "tukey")
plot(emm_F4, comparisons = TRUE)

#### F18 ####
dfF18 <- subset(plot_data, Target == "F18")

f18 <- ggplot(dfF18,
              aes(x = study.day.num, y = logConc_mg,
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
  scale_color_manual(values = group_colors, breaks = names(group_colors), drop = FALSE) +
  labs(title = "F18", x = "", y = y_label_expr, color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background  = element_blank()
  )
print(f18)
ggsave("plots/F18.pdf", f18 , width = 8, height = 8)

ggplot(dfF18,
       aes(x = study.day.num, y = logConc_mg,
           color = groupe.y, group = groupe.y)) +
  geom_point(position = pos, size = 3, alpha = 0.6) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 4),
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
  scale_color_manual(
    values = group_colors,
    breaks = names(group_labels),
    labels = group_labels
  ) +
  labs(title = expression("ETEC F18 ("*italic(fedA)*")"),
       x     = "",
       y     = y_label_expr,
       color = "Group") +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black", face = "bold"),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    plot.background  = element_blank()
  )

dfF18 <- dfF18 %>%
  mutate(
    group = factor(groupe.x),
    study.day = factor(study.day.x)
  )

LM_F18 <- lm(logConc_mg ~ group * study.day, data = dfF18)
par(mfrow = c(2, 2))
plot(LM_F18)

emm_F18 <- emmeans(LM_F18, ~ group | study.day, adjust = "sidak")
contrast(emm_F18, method = "tukey")
plot(emm_F18, comparisons = TRUE)

#### PLOTTING TRENDS ####

plot_data_trend_all <- plot_data %>%
  filter(Target %in% c("16S", "ybbw","LT", "F4", "F18"))

target_colors <- c(
  "16S" = "#8E8070",
  "ybbw" = "#9278C3",
  "LT"   = "#3C965A",
  "F4"   = "#BA5624",
  "F18"  = "#F6CA83"
)
target_labels <- c(
  "16S"  = expression(italic('E. coli')~"(16S rRNA gene)"),
  "ybbw" = expression(italic('E. coli')~"(ybbw)"),
  "LT"   = expression("ETEC ("*italic(lt)*")"),
  "F4"   = expression("ETEC F4 ("*italic(faeG)*")"),
  "F18"  = expression("ETEC F18 ("*italic(fedA)*")")
)

plot_data_trend_all$Target <- factor(plot_data_trend_all$Target,
                                     levels = c("16S", "ybbw","LT", "F4", "F18"))

y_limits <- c(0, 8)

t1 <- ggplot(plot_data_trend_all,
             aes(x = study.day.num, y = logConc_mg,
                 color = Target, fill = Target, group = Target)) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    size    = 1.5,
    alpha   = 0.25,
    na.rm   = TRUE
  ) +
  facet_grid(. ~ groupe.y, labeller = labeller(groupe.y = group_labels)) +
  scale_color_manual(values = target_colors, labels = target_labels, name = "Series") +
  scale_fill_manual(values = target_colors, labels = target_labels, guide = "none") +
  scale_x_continuous(breaks = day_breaks, labels = day_labels,
                     expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Trends",
    x = "Day",
    y = y_label_expr    
  ) +
  theme_classic(base_size = 30) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold" ),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(face = "bold", size = 18),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "right"
  )
print(t1)
ggsave("plots/trend1.pdf", t1 , width = 10, height = 8)

# second trend with legend on left
t2 <- ggplot(plot_data_trend_all,
             aes(x = study.day.num, y = logConc_mg,
                 color = Target, fill = Target, group = Target)) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 5),
    se      = TRUE,
    size    = 1.5,
    alpha   = 0.25,
    na.rm   = TRUE
  ) +
  facet_grid(. ~ groupe.y, labeller = labeller(groupe.y = group_labels)) +
  scale_color_manual(values = target_colors, labels = target_labels, name = "Series") +
  scale_fill_manual(values = target_colors, labels = target_labels, guide = "none") +
  scale_x_continuous(breaks = day_breaks, labels = day_labels,
                     expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Trends",
    x = "Day",
    y = y_label_expr    
  ) +
  theme_classic(base_size = 30) +
  theme(
    axis.text.x = element_text( hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold" ),
    axis.line   = element_line(linetype = "solid", linewidth = 1.0),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(face = "bold", size = 18),
    strip.background = element_rect(fill = "gray80", color = "black"),
    legend.position = "none"
  )

print(t2)
ggsave("plots/trend2.pdf", t2 , width = 10, height = 8)
