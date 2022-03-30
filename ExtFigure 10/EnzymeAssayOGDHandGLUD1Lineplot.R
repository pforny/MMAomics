#line plot NADH produced over time by OGDH
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(data.table)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
require(ggrepel)
require(ggsci)
require(patchwork)

system("mkdir Figs/enzymatics")

fig_path <- c("Figs/enzymatics/")
fig_path_pdf <- c("Figs/v8/pdf/SuppFig_InteractionEnzymeAct/")



##########################
# OGDH


#---WT
df1805 <- read.csv("Data/enzymatics/OGDH/1805_scatterdataWT.csv", header = T, sep=",")
df3105 <- read.csv("Data/enzymatics/OGDH/3105_scatterdataWT.csv", header = T, sep=",")
df1506 <- read.csv("Data/enzymatics/OGDH/1506_scatterdataWT.csv", header = T, sep=",")

df3105_mean<-df3105  %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1805_mean <- df1805 %>% 
  group_by(time) %>% 
  summarise(mean_NADH =mean(normNADH))

df1506_mean <- df1506 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.wt <- bind_rows(df1805_mean, df3105_mean, df1506_mean)
master.wt_mean <- master.wt %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "Control")

#----OGDH
df1805OGDH <- read.csv("Data/enzymatics/OGDH/1805_scatterdataOGDH.csv", header = T, sep=",")
df3105OGDH <- read.csv("Data/enzymatics/OGDH/3105_scatterdataOGDH.csv", header = T, sep=",")
df1506OGDH <- read.csv("Data/enzymatics/OGDH/1506_scatterdataOGDH.csv", header = T, sep=",")

df1805OGDH_mean <- df1805OGDH %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df3105OGDH_mean <- df3105OGDH %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1506OGDH_mean <- df1506OGDH %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.OGDH <- bind_rows(df1805OGDH_mean, df3105OGDH_mean, df1506OGDH_mean)
master.OGDH_mean <- master.OGDH %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "OGDH KO")

#----DLST
df1805DLST <- read.csv("Data/enzymatics/OGDH/1805_scatterdataDLST.csv", header = T, sep=",")
df3105DLST <- read.csv("Data/enzymatics/OGDH/3105_scatterdataDLST.csv", header = T, sep=",")
df1506DLST <- read.csv("Data/enzymatics/OGDH/1506_scatterdataDLST.csv", header = T, sep=",")

df1805DLST_mean <- df1805DLST %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df3105DLST_mean <- df3105DLST %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1506DLST_mean <- df1506DLST %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.DLST <- bind_rows(df1805DLST_mean, df3105DLST_mean, df1506DLST_mean)
master.DLST_mean <- master.DLST %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "DLST KO")

#----MMUT
df1805MMUT <- read.csv("Data/enzymatics/OGDH/1805_scatterdataMUT.csv", header = T, sep=",")
df3105MMUT <- read.csv("Data/enzymatics/OGDH/3105_scatterdataMUT.csv", header = T, sep=",")
df1506MMUT <- read.csv("Data/enzymatics/OGDH/1506_scatterdataMUT.csv", header = T, sep=",")

df1805MMUT_mean <- df1805MMUT %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df3105MMUT_mean <- df3105MMUT %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1506MMUT_mean <- df1506MMUT %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.MMUT <- bind_rows(df1805MMUT_mean, df3105MMUT_mean, df1506MMUT_mean)
master.MMUT_mean <- master.MMUT %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "MMUT KO")

#----GLUD1
df1805GLUD1 <- read.csv("Data/enzymatics/OGDH/1805_scatterdataGLUD1.csv", header = T, sep=",")
df3105GLUD1 <- read.csv("Data/enzymatics/OGDH/3105_scatterdataGLUD1.csv", header = T, sep=",")
df1506GLUD1 <- read.csv("Data/enzymatics/OGDH/1506_scatterdataGLUD1.csv", header = T, sep=",")

df1805GLUD1_mean <- df1805GLUD1 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df3105GLUD1_mean <- df3105GLUD1 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1506GLUD1_mean <- df1506GLUD1 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.GLUD1 <- bind_rows(df1805GLUD1_mean, df3105GLUD1_mean, df1506GLUD1_mean)
master.GLUD1_mean <- master.GLUD1 %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "GLUD1 KO")

#----GLUD2
df1805GLUD2 <- read.csv("Data/enzymatics/OGDH/1805_scatterdataGLUD2.csv", header = T, sep=",")
df3105GLUD2 <- read.csv("Data/enzymatics/OGDH/3105_scatterdataGLUD2.csv", header = T, sep=",")
df1506GLUD2 <- read.csv("Data/enzymatics/OGDH/1506_scatterdataGLUD2.csv", header = T, sep=",")

df1805GLUD2_mean <- df1805GLUD2 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df3105GLUD2_mean <- df3105GLUD2 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1506GLUD2_mean <- df1506GLUD2 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.GLUD2 <- bind_rows(df1805GLUD2_mean, df3105GLUD2_mean, df1506GLUD2_mean)
master.GLUD2_mean <- master.GLUD2 %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "GLUD2 KO")

#----GLUD1_2
df1805GLUD1_2 <- read.csv("Data/enzymatics/OGDH/1805_scatterdataGLUD1_2.csv", header = T, sep=",")
df3105GLUD1_2 <- read.csv("Data/enzymatics/OGDH/3105_scatterdataGLUD1_2.csv", header = T, sep=",")
df1506GLUD1_2 <- read.csv("Data/enzymatics/OGDH/1506_scatterdataGLUD1_2.csv", header = T, sep=",")

df1805GLUD1_2_mean <- df1805GLUD1_2 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df3105GLUD1_2_mean <- df3105GLUD1_2 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1506GLUD1_2_mean <- df1506GLUD1_2 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.GLUD1_2 <- bind_rows(df1805GLUD1_2_mean, df3105GLUD1_2_mean, df1506GLUD1_2_mean)
master.GLUD1_2_mean <- master.GLUD1_2 %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "GLUD1/2 KO")
#------
master.df <- bind_rows(master.wt_mean, master.OGDH_mean, master.DLST_mean, master.MMUT_mean, master.GLUD1_mean, master.GLUD2_mean, master.GLUD1_2_mean)

cumulative.plot <- master.df %>% 
  group_by(group) %>%
  mutate(cumulative.NADH =(cumsum(meanNADH)))

cumulative.plot <- cumulative.plot %>% 
  mutate(time2 = seq(0, 660, by = 60))

#--------
cumulative.plot$group <- fct_relevel(cumulative.plot$group, c("Control", "MMUT KO",  "OGDH KO", "DLST KO", "GLUD1 KO", "GLUD2 KO", "GLUD1/2 KO"))

cumulative.plot <- cumulative.plot %>% mutate(label = ifelse(time2 == max(time2), as.character(group), NA_character_))

cumulative.plot <- data.table(cumulative.plot)



ogdh_lineplot <- 
ggplot(cumulative.plot, aes(x = time2, y = cumulative.NADH, col = group)) +
  geom_label_repel(aes(label = label, fill = group), na.rm = TRUE, hjust = 0, nudge_x = 800, direction = "y", size = 3, color = "white", segment.color = "black", segment.size = 0.25) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = cumulative.NADH-sd, ymax = cumulative.NADH+sd), width = 0, size = 0.5) +
  geom_smooth(size = 0.5) +
  xlab("Time [s]")+
  ylab("Mean NADH produced\nby KGDH enzyme [nmol/mL]")+
  scale_color_aaas() +
  scale_fill_aaas() +
  scale_x_continuous(limits = c(0, 940), breaks = c(0,200, 400, 600)) +
  labs(color = "Clone type") +
  theme_pubr() +
  theme(legend.position = "none")

ggsave(paste(fig_path,"OGDHActivity_lineplot.png", sep = ""), ogdh_lineplot, device = png(), width = 4, height = 3)
dev.off()
ggsave(paste0(fig_path_pdf,"OGDHActivity_lineplot.pdf"), ogdh_lineplot, device = "pdf", width = 5, height = 3)
dev.off()



# same plot but excluding the GLUD KOs

ogdh_lineplot_red <- 
ggplot(cumulative.plot[group != "GLUD1 KO" & group != "GLUD2 KO" & group != "GLUD1/2 KO", ], aes(x = time2, y = cumulative.NADH, col = group)) +
  geom_label_repel(aes(label = label, fill = group), na.rm = TRUE, hjust = 0, nudge_x = 800, direction = "y", size = 3, color = "white", segment.color = "black", segment.size = 0.25) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = cumulative.NADH-sd, ymax = cumulative.NADH+sd), width = 0, size = 0.5) +
  geom_smooth(size = 0.5) +
  xlab("Time [s]")+
  ylab("Mean NADH produced\nby KGDH enzyme [nmol/mL]")+
  scale_color_aaas() +
  scale_fill_aaas() +
  scale_x_continuous(limits = c(0, 940), breaks = c(0,200, 400, 600)) +
  labs(color = "Clone type") +
  theme_pubr() +
  theme(legend.position = "none")

ggsave(paste0(fig_path_pdf,"OGDHActivity_lineplot_red.pdf"), ogdh_lineplot_red, device = "pdf", width = 4, height = 3)
dev.off()
  

##########################
# GLUD1



#---WT
df1805 <- read.csv("Data/enzymatics/GLUD1/1805_scatterdataWT_GLUD.csv", header = T, sep=",")
df3105 <- read.csv("Data/enzymatics/GLUD1/3105_scatterdataWT_GLUD.csv", header = T, sep=",")
df1506 <- read.csv("Data/enzymatics/GLUD1/1506_scatterdataWT_GLUD.csv", header = T, sep=",")

df3105_mean<-df3105  %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1805_mean <- df1805 %>% 
  group_by(time) %>% 
  summarise(mean_NADH =mean(normNADH))

df1506_mean <- df1506 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.wt <- bind_rows(df1805_mean, df3105_mean, df1506_mean)
master.wt_mean <- master.wt %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "Control")

#----OGDH
df1805OGDH <- read.csv("Data/enzymatics/GLUD1/1805_scatterdataOGDH_GLUD.csv", header = T, sep=",")
df3105OGDH <- read.csv("Data/enzymatics/GLUD1/3105_scatterdataOGDH_GLUD.csv", header = T, sep=",")
df1506OGDH <- read.csv("Data/enzymatics/GLUD1/1506_scatterdataOGDH_GLUD.csv", header = T, sep=",")

df1805OGDH_mean <- df1805OGDH %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df3105OGDH_mean <- df3105OGDH %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1506OGDH_mean <- df1506OGDH %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.OGDH <- bind_rows(df1805OGDH_mean, df3105OGDH_mean, df1506OGDH_mean)
master.OGDH_mean <- master.OGDH %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "OGDH KO")

#----DLST
df1805DLST <- read.csv("Data/enzymatics/GLUD1/1805_scatterdataDLST_GLUD.csv", header = T, sep=",")
df3105DLST <- read.csv("Data/enzymatics/GLUD1/3105_scatterdataDLST_GLUD.csv", header = T, sep=",")
df1506DLST <- read.csv("Data/enzymatics/GLUD1/1506_scatterdataDLST_GLUD.csv", header = T, sep=",")

df1805DLST_mean <- df1805DLST %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df3105DLST_mean <- df3105DLST %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1506DLST_mean <- df1506DLST %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.DLST <- bind_rows(df1805DLST_mean, df3105DLST_mean, df1506DLST_mean)
master.DLST_mean <- master.DLST %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "DLST KO")

#----MMUT
df1805MMUT <- read.csv("Data/enzymatics/GLUD1/1805_scatterdataMUT_GLUD.csv", header = T, sep=",")
df3105MMUT <- read.csv("Data/enzymatics/GLUD1/3105_scatterdataMUT_GLUD.csv", header = T, sep=",")
df1506MMUT <- read.csv("Data/enzymatics/GLUD1/1506_scatterdataMUT_GLUD.csv", header = T, sep=",")

df1805MMUT_mean <- df1805MMUT %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df3105MMUT_mean <- df3105MMUT %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1506MMUT_mean <- df1506MMUT %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.MMUT <- bind_rows(df1805MMUT_mean, df3105MMUT_mean, df1506MMUT_mean)
master.MMUT_mean <- master.MMUT %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "MMUT KO")

#----GLUD1
df1805GLUD1 <- read.csv("Data/enzymatics/GLUD1/1805_scatterdataGLUD1_GLUD.csv", header = T, sep=",")
df3105GLUD1 <- read.csv("Data/enzymatics/GLUD1/3105_scatterdataGLUD1_GLUD.csv", header = T, sep=",")
df1506GLUD1 <- read.csv("Data/enzymatics/GLUD1/1506_scatterdataGLUD1_GLUD.csv", header = T, sep=",")

df1805GLUD1_mean <- df1805GLUD1 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df3105GLUD1_mean <- df3105GLUD1 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1506GLUD1_mean <- df1506GLUD1 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.GLUD1 <- bind_rows(df1805GLUD1_mean, df3105GLUD1_mean, df1506GLUD1_mean)
master.GLUD1_mean <- master.GLUD1 %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "GLUD1 KO")

#----GLUD2
df1805GLUD2 <- read.csv("Data/enzymatics/GLUD1/1805_scatterdataGLUD2_GLUD.csv", header = T, sep=",")
df3105GLUD2 <- read.csv("Data/enzymatics/GLUD1/3105_scatterdataGLUD2_GLUD.csv", header = T, sep=",")
df1506GLUD2 <- read.csv("Data/enzymatics/GLUD1/1506_scatterdataGLUD2_GLUD.csv", header = T, sep=",")

df1805GLUD2_mean <- df1805GLUD2 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df3105GLUD2_mean <- df3105GLUD2 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1506GLUD2_mean <- df1506GLUD2 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.GLUD2 <- bind_rows(df1805GLUD2_mean, df3105GLUD2_mean, df1506GLUD2_mean)
master.GLUD2_mean <- master.GLUD2 %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "GLUD2 KO")

#----GLUD1_2
df1805GLUD1_2 <- read.csv("Data/enzymatics/GLUD1/1805_scatterdataGLUD1_2_GLUD.csv", header = T, sep=",")
df3105GLUD1_2 <- read.csv("Data/enzymatics/GLUD1/3105_scatterdataGLUD1_2_GLUD.csv", header = T, sep=",")
df1506GLUD1_2 <- read.csv("Data/enzymatics/GLUD1/1506_scatterdataGLUD1_2_GLUD.csv", header = T, sep=",")

df1805GLUD1_2_mean <- df1805GLUD1_2 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df3105GLUD1_2_mean <- df3105GLUD1_2 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

df1506GLUD1_2_mean <- df1506GLUD1_2 %>% 
  group_by(time) %>% 
  summarise(mean_NADH = mean(normNADH))

master.GLUD1_2 <- bind_rows(df1805GLUD1_2_mean, df3105GLUD1_2_mean, df1506GLUD1_2_mean)
master.GLUD1_2_mean <- master.GLUD1_2 %>% 
  group_by(time) %>% 
  summarise(meanNADH = mean(mean_NADH),
            SEM = sd(mean_NADH)/sqrt(3),
            sd = sd(mean_NADH)) %>% 
  mutate(group = "GLUD1/2 KO")
#------

#combine to master.df
master.df <- bind_rows(master.wt_mean, master.OGDH_mean, master.DLST_mean, master.MMUT_mean, master.GLUD1_mean, master.GLUD2_mean, master.GLUD1_2_mean)

cumulative.plot <- master.df %>% 
  group_by(group) %>%
  mutate(cumulative.NADH =(cumsum(meanNADH)))

cumulative.plot <- cumulative.plot %>% 
  mutate(time2 = seq(0, max(master.df$time)-60, by = 60))

#--------
cumulative.plot$group <- fct_relevel(cumulative.plot$group, c("Control", "MMUT KO", "OGDH KO", "DLST KO","GLUD1 KO", "GLUD2 KO", "GLUD1/2 KO"))
  
cumulative.plot <- cumulative.plot %>% mutate(label = ifelse(time2 == max(time2), as.character(group), NA_character_))


glud_lineplot <- 
ggplot(cumulative.plot, aes(x = time2, y = cumulative.NADH, color = group)) +
  geom_label_repel(aes(label = label, fill = group), na.rm = TRUE, hjust = 0, nudge_x = 4000, direction = "y", size = 3, color = "white", segment.color = "black", segment.size = 0.25) +
  geom_point(size = 1, alpha = 1, shape = 16) +
  xlab("Time [sec]")+
  ylab("Mean NADH produced\nby GLUD [nmol/mL]")+
  scale_color_aaas() +
  scale_fill_aaas() +
  scale_x_continuous(limits = c(0, 4000), breaks = c(0,1000, 2000, 3000)) +
  labs(color = "Clone type") +
  theme_pubr() +
  theme(legend.position = "none")

ggsave(paste(fig_path,"GLUDActivity_lineplot.png", sep = ""), glud_lineplot, device = png(), width = 4, height = 3)
dev.off()
ggsave(paste0(fig_path_pdf,"GLUDActivity_lineplot.pdf"), glud_lineplot, device = "pdf", width = 5, height = 3)
dev.off()







##### combined plot

combo_plot <- 
ogdh_lineplot + glud_lineplot + plot_layout(nrow = 2, guides = "collect")

ggsave(paste0(fig_path_pdf,"ComboActivity_lineplot.pdf"), combo_plot, device = "pdf", width = 5, height = 5)
dev.off()


