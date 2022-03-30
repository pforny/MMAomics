#Overview all enzyme Activity Assays
### RELATIVE TO WT1

library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)
library(ggsci)
library(patchwork)


system("mkdir Figs/enzymatics")

fig_path <- c("Figs/enzymatics/")
fig_path_pdf <- c("Figs/v8/pdf/SuppFig_InteractionEnzymeAct/")



#############################
#OGDH plot

#load data
df1805 <- read.csv("Data/enzymatics/OGDH/180521correctOGDHNADHproduced_Florian.csv", header = T, sep=",")
df3105 <- read_excel("Data/enzymatics/OGDH/310521correctOGDHNADHproduced_Florian.xlsx")
df1506 <- read_excel("Data/enzymatics/OGDH/150621correctOGDHNADHproduced_Florian.xlsx")
df1805$date <- rep("18.05", 15)
df3105$date <- rep("31.05", 17)
df1506$date <- rep("15.06", 17)
df1805$group <- c("Control", "OGDH KO", "DLST KO", "DLST KO", "MMUT KO", "MMUT KO", "GLUD1 KO", "GLUD1 KO", "GLUD1 KO", "GLUD2 KO", "GLUD2 KO", "GLUD2 KO", "GLUD1/2 KO", "GLUD1/2 KO", "GLUD1/2 KO")
df3105$group <- c("Control","Control","Control", "OGDH KO", "DLST KO", "DLST KO", "MMUT KO", "MMUT KO", "GLUD1 KO", "GLUD1 KO", "GLUD1 KO", "GLUD2 KO", "GLUD2 KO", "GLUD2 KO", "GLUD1/2 KO", "GLUD1/2 KO", "GLUD1/2 KO")
df1506$group <- c("Control","Control","Control", "OGDH KO", "DLST KO", "DLST KO", "MMUT KO", "MMUT KO", "GLUD1 KO", "GLUD1 KO", "GLUD1 KO", "GLUD2 KO", "GLUD2 KO", "GLUD2 KO", "GLUD1/2 KO", "GLUD1/2 KO", "GLUD1/2 KO")

master.df <- bind_rows(df1805, df3105,  df1506)

master.df<- master.df %>% 
  mutate(cell.line = fct_relevel(cell.line, 
                                 "WT", "MTHFR A368L", "MTHFR A368G","OGDH.KO", "DLST.KO.A03", "DLST.KO.E10", "MMUT.KO.2T2", 
                                 "MMUT.KO.8T2", "GLUD1.KO.DO8", "GLUD1.KO.CO6", 
                                 "GLUD1.KO.BO3", "GLUD2.KO.HO7", "GLUD2.KO.F6",
                                 "GLUD2.KO.FO8", "GLUD1_2.KO.C12", "GLUD1_2.KO.G11", "GLUD1_2.KO.CO3"))
master.df <- master.df %>% 
  mutate(group = fct_relevel(group, 
                             "Control", "MMUT KO", "OGDH KO", "DLST KO",, "GLUD1 KO", "GLUD2 KO", "GLUD1/2 KO"))

# compare <- list(c("WT", "MMUT.KO.2T2", "MMUT.KO.8T2"), c("WT","GLUD1_2.KO.C12", "GLUD1_2.KO.G11", "GLUD1_2.KO.CO3"), c("WT", "GLUD2.KO.HO7", "GLUD2.KO.F6", "GLUD2.KO.FO8"), 
#                 c("WT",  "GLUD1.KO.DO8", "GLUD1.KO.CO6", "GLUD1.KO.BO3"))

# master.df1 <- copy(master.df)

# ggbarplot(master.df1, x = "cell.line", y = "norm.activity.to.WT", add = c("mean_sd","point"),  fill = "cell.line", 
#           palette =c('#C2C0C0','#C2C0C0','#C2C0C0','#7EB5D6','#BDD7E7', '#BDD7E7', '#E04A4E', '#E04A4E','#FFFFD4','#FFFFD4','#FFFFD4','#FEE391','#FEE391','#FEE391','#FACA69','#FACA69','#FACA69'),)+
#   geom_hline(yintercept = 1)+
#   theme(axis.text.x = element_text(face="bold", angle=90))+
#   stat_compare_means(method = "t.test", comparisons = compare, label = 'p.signif')

### RELATIVE TO MEAN OF ALL 3 WT'S
master.df2 <- copy(master.df)
levels(master.df2$cell.line) <- c("Control #1", "Control #2", "Control #3","OGDH KO", "DLST KO #1", "DLST KO #2", "MMUT KO #1", "MMUT KO #2", "GLUD1 KO #1", "GLUD1 KO #2", "GLUD1 KO #3", "GLUD2 KO #1", "GLUD2 KO #2", "GLUD2 KO #3", "GLUD1/2 KO #1", "GLUD1/2 KO #2", "GLUD1/2 KO #3")


plot_df <- master.df2 %>% 
  left_join(master.df2 %>%
              filter(group == "Control") %>% 
              group_by(date) %>% 
              summarize(mean.group.activity = mean(norm.activity, na.rm = TRUE))) %>% 
  mutate(normalized_to_mean_wt = norm.activity/mean.group.activity)

plot_df <- plot_df %>% 
  mutate(group = fct_relevel(group, 
                             "Control", "MMUT KO", "OGDH KO", "DLST KO", "GLUD1 KO", "GLUD2 KO", "GLUD1/2 KO"))

plot_df <- data.table(plot_df)


# t_tests <- plot_df %>% 
#   group_by(group) %>% 
#   summarise(P = t.test(normalized_to_mean_wt, mu = 1)$p.value,
#             Sig = ifelse(P < 0.05, "*", ""),
#             MaxWidth = max(normalized_to_mean_wt))

# rel.comb.WT<-ggbarplot(plot_df, "group", "normalized_to_mean_wt",
#                        fill = "cell.line",  palette = c('#C2C0C0','#C2C0C0','#C2C0C0','#7EB5D6','#BDD7E7', '#BDD7E7', '#E04A4E', '#E04A4E','#FFFFD4','#FFFFD4','#FFFFD4','#FEE391','#FEE391','#FEE391','#FACA69','#FACA69','#FACA69'),
#                        # position = position_dodge2(width = 0.3, preserve = "single"),
#                        position = position_dodge(0.9),
#                        label = FALSE,
#                        add=c("mean_sd", "point"), title ="GLUD activity assay")+
#   geom_text(aes(label = Sig, y = MaxWidth + 0.2), size = 10,
#             data = t_tests)+
#   geom_hline(yintercept =1)
# # ggsave("120721_OGDH_barplot_rel_combWT",rel.comb.WT)

plot_df_means <- plot_df %>% 
  group_by(cell.line) %>%
  summarize(mean = mean(normalized_to_mean_wt), sd = sd(normalized_to_mean_wt))
plot_df_means <- data.table(plot_df_means)
plot_df_means[, group := plot_df[match(plot_df_means$cell.line, plot_df$cell.line), ]$group]

ogdh_barplot <- 
ggplot(plot_df, aes(x = cell.line, y = normalized_to_mean_wt)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.6) +
  geom_errorbar(data = plot_df_means, aes(x = cell.line, ymin = mean-sd, ymax = mean+sd), inherit.aes = FALSE, size = 0.3, width = 0.3) +
  geom_col(data = plot_df_means, aes(x = cell.line, y = mean, fill = group), inherit.aes = FALSE) +
  geom_jitter(fill = "white", shape = 21, width = 0.1) +
  facet_grid(~group, scales = "free_x", space = "free")+
  ylab("KGDH enzymatic activity\n[normalized to control]") +
  xlab("Individual cell clones") +
  scale_fill_aaas() +
  labs(fill = "Clone type") +
  theme_pubr() +
  rotate_x_text(angle = 45) +
  theme(axis.text.x = element_blank(), legend.position = "right", strip.text=element_blank())

ggsave(paste(fig_path,"OGDHActivity_barplot.png", sep = ""), ogdh_barplot, device = png(), width = 7, height = 4)
dev.off()
ggsave(paste0(fig_path_pdf,"OGDHActivity_barplot.pdf"), ogdh_barplot, device = "pdf", width = 7, height = 2.5)
dev.off()


# same plot excluding GLUD KOs
ogdh_barplot_red <- 
ggplot(plot_df[group != "GLUD1 KO" & group != "GLUD2 KO" & group != "GLUD1/2 KO", ], aes(x = cell.line, y = normalized_to_mean_wt)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.6) +
  geom_errorbar(data = plot_df_means[group != "GLUD1 KO" & group != "GLUD2 KO" & group != "GLUD1/2 KO", ], aes(x = cell.line, ymin = mean-sd, ymax = mean+sd), inherit.aes = FALSE, size = 0.3, width = 0.3) +
  geom_col(data = plot_df_means[group != "GLUD1 KO" & group != "GLUD2 KO" & group != "GLUD1/2 KO", ], aes(x = cell.line, y = mean, fill = group), inherit.aes = FALSE) +
  geom_jitter(fill = "white", shape = 21, width = 0.1) +
  facet_grid(~group, scales = "free_x", space = "free")+
  ylab("KGDH enzymatic activity\n[normalized to control]") +
  xlab("Individual cell clones") +
  scale_fill_aaas() +
  labs(fill = "Clone type") +
  theme_pubr() +
  rotate_x_text(angle = 45) +
  theme(axis.text.x = element_blank(), legend.position = "right", strip.text=element_blank())

ggsave(paste0(fig_path_pdf,"OGDHActivity_barplot_red.pdf"), ogdh_barplot_red, device = "pdf", width = 5, height = 3)
dev.off()



#############################
#GLUD1 plot

#load data
df1805 <- read.csv("Data/enzymatics/GLUD1/180521correctGLUDNADHproduced_Florian.csv", header = T, sep=",")
df3105 <- read_excel("Data/enzymatics/GLUD1/310521correctGLUDNADHproduced_Florian.xlsx")
df1506 <- read_excel("Data/enzymatics/GLUD1/150621correctGLUDNADHproduced_Florian.xlsx")
df1805$date <- rep("18.05", 15)
df3105$date <- rep("31.05", 17)
df1506$date <- rep("15.06", 17)
df1805$group <- c("Control", "OGDH KO", "DLST KO", "DLST KO", "MMUT KO", "MMUT KO", "GLUD1 KO", "GLUD1 KO", "GLUD1 KO", "GLUD2 KO", "GLUD2 KO", "GLUD2 KO", "GLUD1/2 KO", "GLUD1/2 KO", "GLUD1/2 KO")
df3105$group <- c("Control","Control","Control", "OGDH KO", "DLST KO", "DLST KO", "MMUT KO", "MMUT KO", "GLUD1 KO", "GLUD1 KO", "GLUD1 KO", "GLUD2 KO", "GLUD2 KO", "GLUD2 KO", "GLUD1/2 KO", "GLUD1/2 KO", "GLUD1/2 KO")
df1506$group <- c("Control","Control","Control", "OGDH KO", "DLST KO", "DLST KO", "MMUT KO", "MMUT KO", "GLUD1 KO", "GLUD1 KO", "GLUD1 KO", "GLUD2 KO", "GLUD2 KO", "GLUD2 KO", "GLUD1/2 KO", "GLUD1/2 KO", "GLUD1/2 KO")

master.df <- bind_rows(df1805, df3105,  df1506)

master.df<- master.df %>% 
  mutate(cell.line = fct_relevel(cell.line, 
                                 "WT", "MTHFR A368L", "MTHFR A368G","OGDH.KO", "DLST.KO.A03", "DLST.KO.E10", "MMUT.KO.2T2", 
                                 "MMUT.KO.8T2", "GLUD1.KO.DO8", "GLUD1.KO.CO6", 
                                 "GLUD1.KO.BO3", "GLUD2.KO.HO7", "GLUD2.KO.F6",
                                 "GLUD2.KO.FO8", "GLUD1_2.KO.C12", "GLUD1_2.KO.G11", "GLUD1_2.KO.CO3"))

master.df <- master.df %>% 
  mutate(group = fct_relevel(group, 
                             "Control", "MMUT KO", "OGDH KO", "DLST KO", "GLUD1 KO", "GLUD2 KO", "GLUD1/2 KO"))

# compare <- list(c("WT", "MMUT.KO.2T2", "MMUT.KO.8T2"), c("WT","GLUD1_2.KO.C12", "GLUD1_2.KO.G11", "GLUD1_2.KO.CO3"), c("WT", "GLUD2.KO.HO7", "GLUD2.KO.F6", "GLUD2.KO.FO8"), 
#                 c("WT",  "GLUD1.KO.DO8", "GLUD1.KO.CO6", "GLUD1.KO.BO3"))

# master.df1 <- copy(master.df)

# ggbarplot(master.df1, x = "cell.line", y = "norm.activity.to.WT", add = c("mean_sd","point"),  fill = "cell.line", 
#           palette =c('#C2C0C0','#C2C0C0','#C2C0C0','#7EB5D6','#BDD7E7', '#BDD7E7', '#E04A4E', '#E04A4E','#FFFFD4','#FFFFD4','#FFFFD4','#FEE391','#FEE391','#FEE391','#FACA69','#FACA69','#FACA69'),)+
#   geom_hline(yintercept = 1)+
#   theme(axis.text.x = element_text(face="bold", angle=90))+
#   stat_compare_means(method = "t.test", comparisons = compare, label = 'p.signif')

### RELATIVE TO MEAN OF ALL 3 WT'S
master.df2 <- copy(master.df)
levels(master.df2$cell.line) <- c("Control #1", "Control #2", "Control #3","OGDH KO", "DLST KO #1", "DLST KO #2", "MMUT KO #1", "MMUT KO #2", "GLUD1 KO #1", "GLUD1 KO #2", "GLUD1 KO #3", "GLUD2 KO #1", "GLUD2 KO #2", "GLUD2 KO #3", "GLUD1/2 KO #1", "GLUD1/2 KO #2", "GLUD1/2 KO #3")


plot_df <- master.df2 %>% 
  left_join(master.df2 %>%
              filter(group == "Control") %>% 
              group_by(date) %>% 
              summarize(mean.group.activity = mean(norm.activity, na.rm = TRUE))) %>% 
  mutate(normalized_to_mean_wt = norm.activity/mean.group.activity)

plot_df <- plot_df %>% 
  mutate(group = fct_relevel(group, 
                             "Control", "MMUT KO","OGDH KO", "DLST KO", "GLUD1 KO", "GLUD2 KO", "GLUD1/2 KO"))

plot_df <- data.table(plot_df)


# t_tests <- plot_df %>% 
#   group_by(group) %>% 
#   summarise(P = t.test(normalized_to_mean_wt, mu = 1)$p.value,
#             Sig = ifelse(P < 0.05, "*", ""),
#             MaxWidth = max(normalized_to_mean_wt))

# rel.comb.WT<-ggbarplot(plot_df, "group", "normalized_to_mean_wt",
#           fill = "cell.line",  palette = c('#C2C0C0','#C2C0C0','#C2C0C0','#7EB5D6','#BDD7E7', '#BDD7E7', '#E04A4E', '#E04A4E','#FFFFD4','#FFFFD4','#FFFFD4','#FEE391','#FEE391','#FEE391','#FACA69','#FACA69','#FACA69'),
#           # position = position_dodge2(width = 0.3, preserve = "single"),
#           position = position_dodge(0.9),
#           label = FALSE,
#           add=c("mean_sd", "point"), title ="OGDH activity assay")+
#   geom_text(aes(label = Sig, y = MaxWidth + 0.2), size = 10,
#             data = t_tests)+
#   geom_hline(yintercept =1)
# # ggsave("120721_GDH_barplot_rel_combWT",rel.comb.WT)


plot_df_means <- plot_df %>% 
  group_by(cell.line) %>%
  summarize(mean = mean(normalized_to_mean_wt), sd = sd(normalized_to_mean_wt))
plot_df_means <- data.table(plot_df_means)
plot_df_means[, group := plot_df[match(plot_df_means$cell.line, plot_df$cell.line), ]$group]

glud_barplot <- 
ggplot(plot_df, aes(x = cell.line, y = normalized_to_mean_wt)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.6) +
  geom_errorbar(data = plot_df_means, aes(x = cell.line, ymin = mean-sd, ymax = mean+sd), inherit.aes = FALSE, size = 0.3, width = 0.3) +
  geom_col(data = plot_df_means, aes(x = cell.line, y = mean, fill = group), inherit.aes = FALSE) +
  geom_jitter(fill = "white", shape = 21, width = 0.1) +
  facet_grid(~group, scales = "free_x", space = "free")+
  ylab("GLUD1 enzymatic activity\n[normalized to control]") +
  xlab("Individual cell clones") +
  scale_fill_aaas() +
  labs(fill = "Clone type") +
  theme_pubr() +
  rotate_x_text(angle = 45) +
  theme(axis.text.x = element_blank(), legend.position = "right", strip.text=element_blank())

ggsave(paste(fig_path,"GLUDActivity_barplot.png", sep = ""), glud_barplot, device = png(), width = 7, height = 4)
dev.off()
ggsave(paste0(fig_path_pdf,"GLUDActivity_barplot.pdf"), glud_barplot, device = "pdf", width = 7, height = 2.5)
dev.off()






combo_plot <- 
ogdh_barplot + glud_barplot + plot_layout(nrow = 2, guides = "collect")

ggsave(paste0(fig_path_pdf,"ComboActivity_barplot.pdf"), combo_plot, device = "pdf", width = 7, height = 5)
dev.off()
