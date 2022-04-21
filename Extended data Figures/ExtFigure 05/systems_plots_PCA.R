# libraries
require(data.table)
require(ggplot2)
require(tidyverse)
require(ggpubr)
library(factoextra)
library(ggsci)
library(patchwork)


### load some helper functions
source("csrproject/Code/normalizationFunctions.R")


# create figures path
system("mkdir Figs/systems_plots")

fig_path <- c("Figs/systems_plots/")
fig_path_pdf <- c("Figs/v8/pdf/")



######################################################################
######################################################################
## PCA plot of protein dataset
######################################################################
######################################################################

load("interimData/prot_mat.RData")

prot_exp <- data.table(ENSG = row.names(prot_mat), prot_mat)
prot_exp[,PG.Qvalue:=NULL] 
prot_exp2 <- melt(prot_exp, id = c("ENSG"))

# length(unique(prot_exp2[variable == "MMA001",]$ENSG))
# mma1_genes <- prot_exp2[variable == "MMA001",]$ENSG
# mma1_genes[duplicated(mma1_genes)]
# prot_exp2[ENSG == mma1_genes[duplicated(mma1_genes)],]

# there is a duplicate gene, to be removed:
prot_exp3 <- prot_exp2 %>% distinct(ENSG, variable, .keep_all = TRUE)

prot_exp4 <- dcast(prot_exp3, variable ~ ENSG, value.var = c("value"))

row.names(prot_exp4) <- paste(prot_exp4$variable, row.names(prot_exp4), sep = "_")
prot_exp4$variable <- NULL

# normalize
#prot_exp5 <- scale(prot_exp4)
prot_exp5 <- scalerRound(prot_exp4, 7)
prot_exp5 <- data.table(prot_exp5)

hist(prot_exp4[[1244]], breaks = 30)
hist(prot_exp5[[1244]], breaks = 30)

# PCA plot
prot_exp_pca <- prcomp(prot_exp5)
plot(prot_exp_pca$x[,1], prot_exp_pca$x[,2])

df_out_prot <- as.data.frame(prot_exp_pca$x)
df_out_prot$no <- as.numeric(row.names(prot_exp5))

df_out_prot$type <- ifelse(df_out_prot$no < 151, "MMUT def.", ifelse(df_out_prot$no > 150 & df_out_prot$no < 211, "unknown", "unaffected"))
df_out_prot$type2 <- ifelse(df_out_prot$no < 151, "MMUT def.", "control")


#calculate percentages
perc_prot_pca <- vector("character", nrow(prot_exp_pca$x))

for(i in 1:nrow(prot_exp_pca$x)) {
	perc_prot_pca[i] <- round(((prot_exp_pca$sdev[i])^2) / sum(prot_exp_pca$sdev^2) * 100, 2)
	perc_prot_pca[i] <- paste(colnames(df_out_prot)[i], " (", paste(as.character(perc_prot_pca[i]), "%", ")", sep = ""), sep = "")
}

# scree plot with package factoextra
scree_plot_prot <- 
fviz_eig(prot_exp_pca) +
	ggtitle("Scree plot proteomics PCA") +
	ylab("% of explained variances") +
	theme_pubr() +
	theme(plot.title = element_text(size = 12))


g_prot_pca <-
ggplot(df_out_prot, aes(x = PC1, y = PC2, color = type)) +
	geom_point() +
	ggtitle("Proteomics PCA") +
	xlab(perc_prot_pca[1]) +
	ylab(perc_prot_pca[2]) +
	theme_bw()


g_prot_pca2 <-
ggplot(df_out_prot, aes(x = PC1, y = PC2, color = type2)) +
	geom_point() +
	ggtitle("Proteomics PCA") +
	xlab(perc_prot_pca[1]) +
	ylab(perc_prot_pca[2]) +
	scale_color_aaas() +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), legend.title = element_blank(), legend.position = "bottom")



df_out_prot <- data.table(df_out_prot)
pcs_and_type <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "type")

df_out_prot_melt <- melt(df_out_prot[, ..pcs_and_type])
setnames(df_out_prot_melt, c("type", "PC", "value"))

g_prot_pca_dens <- 
ggplot(df_out_prot_melt, aes(x = value, color = type, fill = type)) +
		geom_density(alpha = 0.5) +
		facet_wrap(.~PC) +
		ggtitle("Single PCs") +
		theme_bw()



pca_prot2 <- 
g_prot_pca2 / scree_plot_prot

ggsave("Figs/systems_plots/pca_plots_prot_scalerRound_simple.png", pca_prot2, device = png(), width = 4, height = 5)
dev.off()
ggsave(paste(fig_path_pdf, "SuppFig4/", "pca_plots_prot_scalerRound_simple.pdf", sep = ""), pca_prot2, device = "pdf", width = 3, height = 5, bg = "white")
dev.off()


######################################################################
######################################################################
## PCA plot of RNA dataset
######################################################################
######################################################################

load("interimData/rnafiltered_mat.RData")


rna_exp <- data.table(ENSG = row.names(rnafiltered_mat), rnafiltered_mat)
rna_exp2 <- melt(rna_exp, id = c("ENSG"))

length(unique(rna_exp2[variable == "MMA001",]$ENSG)) == length(rna_exp2[variable == "MMA001",]$ENSG)
# mma1_genes <- prot_exp2[variable == "MMA001",]$ENSG
# mma1_genes[duplicated(mma1_genes)]
# prot_exp2[ENSG == mma1_genes[duplicated(mma1_genes)],]

rna_exp3 <- dcast(rna_exp2, variable ~ ENSG, value.var = c("value"))

row.names(rna_exp3) <- paste(rna_exp3$variable, row.names(rna_exp3), sep = "_")
rna_exp3$variable <- NULL


# normalize
rna_exp4 <- scalerRound(rna_exp3, 7)
rna_exp4 <- scale(rna_exp3)
rna_exp4 <- data.table(rna_exp4)

hist(rna_exp3[[1244]], breaks = 30)
hist(rna_exp4[[1244]], breaks = 30)


# PCA plot
rna_exp_pca <- prcomp(rna_exp4)
plot(rna_exp_pca$x[,1], rna_exp_pca$x[,2])

df_out_rna <- as.data.frame(rna_exp_pca$x)
df_out_rna$group <- sapply(strsplit(as.character(row.names(rna_exp4)), "_"), "[[", 1 )
df_out_rna$no <- sub("MMA", "", df_out_rna$group)
df_out_rna$no <- as.numeric(df_out_rna$no)
df_out_rna$type <- ifelse(df_out_rna$no < 151, "MMUT def.", ifelse(df_out_rna$no > 150 & df_out_rna$no < 211, "unknown", "unaffected"))
df_out_rna$type2 <- ifelse(df_out_rna$no < 151, "MMUT def.", "control")

#calcluate percentages
perc_rna_pca <- vector("character", nrow(rna_exp_pca$x))

for(i in 1:nrow(rna_exp_pca$x)) {
	perc_rna_pca[i] <- round(((rna_exp_pca$sdev[i])^2) / sum(rna_exp_pca$sdev^2) * 100, 2)
	perc_rna_pca[i] <- paste(colnames(df_out_prot)[i], " (", paste(as.character(perc_rna_pca[i]), "%", ")", sep = ""), sep = "")
}

# scree plot with package factoextra
scree_plot_rna <- 
fviz_eig(rna_exp_pca) +
	ggtitle("Scree plot transcript. PCA") +
	ylab("% of explained variances") +
	theme_pubr() +
	theme(plot.title = element_text(size = 12))

g_rna_pca <-
ggplot(df_out_rna, aes(x = PC1, y = PC2, color = type)) +
	geom_point() +
	ggtitle("Transcriptomics PCA") +
	xlab(perc_rna_pca[1]) +
	ylab(perc_rna_pca[2]) +
	theme_bw()

g_rna_pca2 <-
ggplot(df_out_rna, aes(x = PC1, y = PC2, color = type2)) +
	geom_point() +
	ggtitle("Transcriptomics PCA") +
	xlab(perc_rna_pca[1]) +
	ylab(perc_rna_pca[2]) +
	scale_color_aaas() +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), legend.title = element_blank(), legend.position = "bottom")


df_out_rna <- data.table(df_out_rna)
pcs_and_type <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "type")

df_out_rna_melt <- melt(df_out_rna[, ..pcs_and_type])
setnames(df_out_rna_melt, c("type", "PC", "value"))

g_rna_pca_dens <- 
ggplot(df_out_rna_melt, aes(x = value, color = type, fill = type)) +
		geom_density(alpha = 0.5) +
		facet_wrap(.~PC) +
		ggtitle("Single PCs") +
		theme_bw()



pca_rna2 <- 
g_rna_pca2 / scree_plot_rna

ggsave("Figs/systems_plots/pca_plots_rna_scalerRound_simple.png", pca_rna2, device = png(), width = 4, height = 5)
dev.off()
ggsave(paste(fig_path_pdf, "SuppFig4/", "pca_plots_rna_scalerRound_simple.pdf", sep = ""), pca_rna2, device = "pdf", width = 3, height = 5, bg = "white")
dev.off()




pca_mix <- 
g_prot_pca2 / scree_plot_prot | g_rna_pca2 / scree_plot_rna

ggsave(paste(fig_path_pdf, "SuppFig4/", "pca_RNAandProtein.pdf", sep = ""), pca_mix, device = "pdf", width = 6, height = 5, bg = "white")
dev.off()

























