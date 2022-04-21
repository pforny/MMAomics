# libraries
require(data.table)
require(ggplot2)
require(tidyverse)
require(ggpubr)
require(ggcorrplot)
require(corrr)
library(ggsci)
library(ggrepel)
library(cowplot)
library(patchwork)
library(matrixTests)
library(viridis)



mypal <- pal_aaas("default", alpha = 1)(9)


### load some helper functions
# source("csrproject/Code/normalizationFunctions.R")


# create figures path
system("mkdir Figs/systems_plots")

fig_path <- c("Figs/systems_plots/")
fig_path_pdf <- c("Figs/v8/pdf/")



######################################################################
######################################################################
## overall correlation of transcript and protein data
######################################################################
######################################################################

##################################
# load protein matrix
load("interimData/prot_mat.RData")

prot_exp <- data.table(ENSG = row.names(prot_mat), prot_mat)
prot_exp <- prot_exp %>% distinct(ENSG, .keep_all = TRUE)
prot_exp_samples <- prot_exp[, -"ENSG"]
prot_exp_prots <- t(prot_exp_samples)

prot_exp_samples_norm <- scale(as.matrix(prot_exp_samples))
prot_exp_samples_norm <- data.table(prot_exp_samples_norm)

# hist(as.numeric(prot_exp_samples_norm[4,]))

# prot_exp2 <- melt(prot_exp, id = c("ENSG"))
# prot_exp2 <- prot_exp2 %>% distinct(ENSG, variable, .keep_all = TRUE)



##################################
# load transcript matrix
load("interimData/rnafiltered_mat.RData")

rna_exp <- data.table(ENSG = row.names(rnafiltered_mat), rnafiltered_mat)
rna_exp_samples <- rna_exp[, -"ENSG"]
rna_exp_rnas <- t(rna_exp_samples)







##################################
# simple analysis of expression, fold change and pvalue, generate ranking plots

# RNA

# calculate stats of expression levels
rna_exp_fold <- copy(rna_exp)
ctrl_names <- names(rna_exp_fold[, 145:222])
MMUT_names <- names(rna_exp_fold[, 2:144])

# calculate means per group
rna_exp_fold[, mean_ctrl := rowMeans(rna_exp_fold[,ctrl_names, with = FALSE], na.rm = TRUE)]
rna_exp_fold[, mean_MMUT := rowMeans(rna_exp_fold[,MMUT_names, with = FALSE], na.rm = TRUE)]

# exclude NAs
rna_exp_fold <- rna_exp_fold[!is.na(rna_exp_fold$mean_ctrl), ]
rna_exp_fold <- rna_exp_fold[!is.na(rna_exp_fold$mean_MMUT), ]

# calculate fold change
rna_exp_fold[, fold := mean_MMUT/mean_ctrl]

# calculate pvals
stats_tbl <- data.table(row_t_welch(rna_exp_fold[, MMUT_names, with = FALSE], rna_exp_fold[, ctrl_names, with = FALSE]))

# summary table
stats_tbl1 <- data.table(gene = rna_exp_fold$ENSG, fold = rna_exp_fold$fold, pval = stats_tbl$pvalue)
stats_tbl_pval <- stats_tbl1[order(pval), ]
stats_tbl_fold <- stats_tbl1[order(-fold), ]

pval_exp_rank_rna <- 
ggplot(stats_tbl_pval, aes(y=-log10(pval), x=reorder(gene, pval))) +
	geom_col(fill = "grey") +
	geom_point(data=stats_tbl_pval[gene == "ENSG00000146085"], size = 2) +
  geom_text_repel(data=stats_tbl_pval[gene == "ENSG00000146085"], label = "MMUT", nudge_x = 4000, nudge_y = 0, show.legend = FALSE) +
	geom_hline(yintercept=-log10(0.01), linetype = 2, alpha = 0.6) +
	theme_pubr() +
	labs(y="-log10(p value)", x="Transcript rank") +
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

fold_exp_rank_rna <- 
ggplot(stats_tbl_pval, aes(y=-log10(fold), x=reorder(gene, fold))) +
	geom_col(fill = "grey") +
	geom_point(data=stats_tbl_pval[gene == "ENSG00000146085"], size = 2) +
  geom_text_repel(data=stats_tbl_pval[gene == "ENSG00000146085"], label = "MMUT", nudge_x = 4000, nudge_y = 0.3, show.legend = FALSE) +
	theme_pubr() +
	labs(y="-log10(fold change)", x="Transcript rank") +
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())








# protein

prot_exp_fold <- copy(prot_exp)
ctrl_names_prot <- names(prot_exp_fold[, 145:222])
MMUT_names_prot <- names(prot_exp_fold[, 2:144])

# calculate means per group
prot_exp_fold[, mean_ctrl := rowMeans(prot_exp_fold[,ctrl_names_prot, with = FALSE], na.rm = TRUE)]
prot_exp_fold[, mean_MMUT := rowMeans(prot_exp_fold[,MMUT_names_prot, with = FALSE], na.rm = TRUE)]

# exclude NAs
prot_exp_fold <- prot_exp_fold[!is.na(prot_exp_fold$mean_ctrl), ]
prot_exp_fold <- prot_exp_fold[!is.na(prot_exp_fold$mean_MMUT), ]

# calculate fold change
prot_exp_fold[, fold := mean_MMUT/mean_ctrl]

# calculate pvals
stats_tbl_prot <- data.table(row_t_welch(prot_exp_fold[, MMUT_names_prot, with = FALSE], prot_exp_fold[, ctrl_names_prot, with = FALSE]))

# summary table
stats_tbl_prot1 <- data.table(gene = prot_exp_fold$ENSG, fold = prot_exp_fold$fold, pval = stats_tbl_prot$pvalue)
stats_tbl_pval_prot <- stats_tbl_prot1[order(pval), ]
stats_tbl_fold_prot <- stats_tbl_prot1[order(-fold), ]

pval_exp_rank_prot <- 
ggplot(stats_tbl_prot1, aes(y=-log10(pval), x=reorder(gene, pval))) +
	geom_col(fill = "grey") +
	geom_point(data=stats_tbl_prot1[gene == "ENSG00000146085"], size = 2) +
  geom_text_repel(data=stats_tbl_prot1[gene == "ENSG00000146085"], label = "MMUT", nudge_x = 1200, nudge_y = 0, show.legend = FALSE) +
	geom_hline(yintercept=-log10(0.01), linetype = 2, alpha = 0.6) +
	theme_pubr() +
	labs(y="-log10(p value)", x="Protein rank") +
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

fold_exp_rank_prot <- 
ggplot(stats_tbl_prot1, aes(y=-log10(fold), x=reorder(gene, fold))) +
	geom_col(fill = "grey") +
	geom_point(data=stats_tbl_prot1[gene == "ENSG00000146085"], size = 2) +
  geom_text_repel(data=stats_tbl_prot1[gene == "ENSG00000146085"], label = "MMUT", nudge_x = 1200, nudge_y = 0.3, show.legend = FALSE) +
	theme_pubr() +
	labs(y="-log10(fold change)", x="Protein rank") +
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())



rna_ranking <- 
pval_exp_rank_rna + fold_exp_rank_rna + plot_annotation(title = "Transcriptomics: gene ranking", theme = theme(plot.title = element_text(size = 12)))

prot_ranking <- 
pval_exp_rank_prot + fold_exp_rank_prot + plot_annotation(title = "Proteomics: gene ranking", theme = theme(plot.title = element_text(size = 12)))

gene_ranking <- 
ggarrange(rna_ranking, prot_ranking, nrow = 2)


ggsave(paste(fig_path,"GeneExpressionRanks.png", sep = ""), gene_ranking, device = png(), width = 5, height = 4)
dev.off()

pval_ranking <- 
pval_exp_rank_rna + pval_exp_rank_prot + plot_annotation(title = "Gene expression ranking", theme = theme(plot.title = element_text(size = 12)))

fold_ranking <- 
fold_exp_rank_rna / fold_exp_rank_prot + plot_annotation(title = "Gene expression ranking", theme = theme(plot.title = element_text(size = 12)))

ggsave(paste(fig_path_pdf,"Fig1/","GeneExpressionRanks_pvals.pdf", sep = ""), pval_ranking, device = "pdf", width = 5, height = 3)
dev.off()
ggsave(paste(fig_path_pdf,"SuppFig2/","GeneExpressionRanks_fold.pdf", sep = ""), fold_ranking, device = "pdf", width = 3, height = 5.5)
dev.off()


### CORRELOGRAM for proteomics

# correlation of samples among each other

p_correlo_dat <- cor(prot_exp_samples_norm, use = "complete.obs", method = "pearson")
rownames(p_correlo_dat) <- as.numeric(sub("MMA", "", rownames(p_correlo_dat)))
colnames(p_correlo_dat) <- as.numeric(sub("MMA", "", colnames(p_correlo_dat)))
# p_p_mat <- cor_pmat(prot_exp_samples)

p_correlo <- 
ggcorrplot(p_correlo_dat,
	outline.col = NA,
	hc.order = FALSE,
	type = "full",
	title = "prot_corr_samples_pearson")
p_correlo <- p_correlo +
	rotate_x_text(angle = 90) + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

ggsave(paste(fig_path,"prot_corr_samples.png", sep = ""), p_correlo, device = png(), width = 160, height = 160, units = c("mm"), dpi = 600)
dev.off()


# correlation of proteins among each other

p_correlo_dat1 <- cor(prot_exp_prots, use = "complete.obs", method = "pearson")
# p_p_mat1 <- cor_pmat(prot_exp_prots)

p_correlo1 <- 
ggcorrplot(p_correlo_dat1,
	outline.col = NA,
	hc.order = TRUE,
	type = "full",
	title = "prot_corr_proteins_ALL_pearson")  +
	rotate_x_text(angle = 90) + 
	theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

ggsave(paste(fig_path,"prot_corr_proteins_ALL.png", sep = ""), p_correlo1, device = png(), width = 160, height = 160, units = c("mm"), dpi = 300)
dev.off()


p_correlo_dat1_MMUT <- cor(prot_exp_prots[1:150, ], use = "complete.obs", method = "spearman")
p_correlo1_MMUT <- 
ggcorrplot(p_correlo_dat1_MMUT,
	outline.col = NA,
	hc.order = TRUE,
	type = "full",
	title = "Spearman correlation of proteins (MMUT def.)") +
	rotate_x_text(angle = 90) + 
	theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), plot.title = element_text(size = 10))

p_correlo_dat1_REST <- cor(prot_exp_prots[151:230, ], use = "complete.obs", method = "spearman")
p_correlo1_REST <- 
ggcorrplot(p_correlo_dat1_REST,
	outline.col = NA,
	hc.order = TRUE,
	type = "full",
	title = "Spearman correlation of proteins (control)") +
	rotate_x_text(angle = 90) + 
	theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), plot.title = element_text(size = 10))


p_correlo_MMUTvsREST <- ggarrange(p_correlo1_MMUT, p_correlo1_REST)

ggsave(paste(fig_path,"prot_corr_proteins_MMUTvsREST.png", sep = ""), p_correlo_MMUTvsREST, device = png(), width = 200, height = 100, units = c("mm"), dpi = 300, bg = "white")
dev.off()



### CORRELOGRAM for transcriptomics

# correlation of samples among each other

r_correlo_dat <- cor(rna_exp_samples, use = "complete.obs", method = "pearson")
nrow(r_correlo_dat)
# rownames(r_correlo_dat) <- sub("MMA", "", rownames(r_correlo_dat))
# colnames(r_correlo_dat) <- sub("MMA", "", colnames(r_correlo_dat))
r_p_mat <- cor_pmat(rna_exp_samples)

r_correlo <- 
ggcorrplot(r_correlo_dat,
	outline.col = NA,
	hc.order = FALSE,
	type = "full",
	title = "rna_corr_samples_pearson")
r_correlo <- r_correlo +
	rotate_x_text(angle = 90)

ggsave(paste(fig_path,"rna_corr_samples.png", sep = ""), r_correlo, device = png(), width = 160, height = 160, units = c("mm"), dpi = 600)
dev.off()


# correlation of transcripts among each other

r_correlo_dat1 <- cor(rna_exp_rnas, use = "complete.obs", method = "pearson")
# r_p_mat1 <- cor_pmat(rna_exp_rnas)

r_correlo1 <- 
ggcorrplot(r_correlo_dat1,
	outline.col = NA,
	hc.order = TRUE,
	type = "full",
	title = "rna_corr_rnas_ALL_pearson")

ggsave(paste(fig_path,"rna_corr_rnas_ALL.png", sep = ""), r_correlo1, device = png(), width = 160, height = 160, units = c("mm"), dpi = 300)
dev.off() # takes a lot of memory to compute/generate image file


r_correlo_dat1_MMUT <- cor(rna_exp_rnas[1:143, ], use = "complete.obs", method = "spearman")
r_correlo1_MMUT <- 
ggcorrplot(r_correlo_dat1_MMUT,
	outline.col = NA,
	hc.order = TRUE,
	type = "full",
	title = "Spearman correlation of transcripts (MMUT def.)") +
	rotate_x_text(angle = 90) + 
	theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), plot.title = element_text(size = 10))

ggsave(paste(fig_path,"rna_corr_rnas_MMUT.png", sep = ""), r_correlo1_MMUT, device = png(), width = 100, height = 100, units = c("mm"), dpi = 300, bg = "white")
dev.off()

r_correlo_dat1_REST <- cor(rna_exp_rnas[144:221, ], use = "complete.obs", method = "spearman")
r_correlo1_REST <- 
ggcorrplot(r_correlo_dat1_REST,
	outline.col = NA,
	hc.order = TRUE,
	type = "full",
	title = "Spearman correlation of transcripts (control)") +
	rotate_x_text(angle = 90) + 
	theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), plot.title = element_text(size = 10))

ggsave(paste(fig_path,"rna_corr_rnas_REST.png", sep = ""), r_correlo1_REST, device = png(), width = 100, height = 100, units = c("mm"), dpi = 300, bg = "white")
dev.off()

r_correlo_MMUTvsREST <- ggarrange(r_correlo1_MMUT, r_correlo1_REST)

ggsave(paste(fig_path,"rna_corr_rnas_MMUTvsREST.png", sep = ""), r_correlo_MMUTvsREST, device = png(), width = 200, height = 100, units = c("mm"), dpi = 300, bg = "white")
dev.off()




######################################################################
######################################################################
## spearman and pearson correlation histograms
######################################################################
######################################################################

# create protein and transcript merged data table, important: only include shared genes and shared samples, as protein dataset has less genes and transcript dataset lacks some samples

rna_exp_common_genes <- rna_exp[which(rna_exp[, ENSG] %in% prot_exp[, ENSG]), ]

keep_samples_in_prot <- (which(colnames(prot_exp)[-1] %in% colnames(rna_exp)[-1])+1)
prot_exp_common_genes <- prot_exp[, ..keep_samples_in_prot]
prot_exp_common_genes[, ENSG := prot_exp[, ENSG]]

ncol(rna_exp_common_genes) 
ncol(rna_exp_common_genes) == ncol(prot_exp_common_genes)
nrow(prot_exp_common_genes)
nrow(rna_exp_common_genes)

prot_rna_merge <- merge(prot_exp_common_genes, rna_exp_common_genes, by = "ENSG") # the merge function excludes automatically all the genes (ENSG) which are not shared

ncol(prot_rna_merge)
nrow(prot_rna_merge)
colnames(prot_rna_merge)

# columns 2:222 are protein data
# columns 223:443 are transcript data

#######################################
# correlation of every protein with every transcript of the same gene (only the genes which are shared among the datasets)
# create vectors with correlation data

p_cor_vector <- vector()
s_cor_vector <- vector()

for(i in 1:dim(prot_rna_merge)[1]) {

cat(i, "\n")
p_cor_vector <- c(p_cor_vector,  cor(unlist(prot_rna_merge[i, 2:222]), unlist(prot_rna_merge[i, 223:443]), use = "pairwise.complete.obs", method = "pearson"))
s_cor_vector <- c(s_cor_vector,  cor(unlist(prot_rna_merge[i, 2:222]), unlist(prot_rna_merge[i, 223:443]), use = "pairwise.complete.obs", method = "spearman"))

}

cor_pearson <- data.table(p_cor_vector)
cor_spearman <- data.table(s_cor_vector)

cor_spearman1 <- cor_spearman
cor_spearman1[, ENSG := prot_rna_merge[, ENSG]]
setorder(cor_spearman1, s_cor_vector)

# MMUT Ensembl: ENSG00000146085
# OGDH Ensembl: ENSG00000105953
# GLUD1 Ensembl: ENSG00000148672
genes_of_interest <- c("ENSG00000146085", "ENSG00000105953", "ENSG00000148672")
cor_spearman1[is.element(cor_spearman1$ENSG, genes_of_interest), ]
# cor_spearman1[ENSG == "ENSG00000146648"]


# plot histograms

# ggplot(cor_pearson, aes(x = p_cor_vector)) +
# 	geom_histogram(aes(y = ..density..), color = "black", fill = "white") +
# 	geom_density(alpha = 0.2, fill = "darkblue", color = NA) +
# 	theme_bw()

# ggplot(cor_spearman, aes(x = s_cor_vector)) +
# 	geom_histogram(aes(y = ..density..), color = "black", fill = "white") +
# 	geom_density(alpha = 0.2, fill = "darkblue", color = NA) +
# 	theme_bw()

mean(s_cor_vector)

prot_rna_corr_s_plot <- 
ggplot(cor_spearman, aes(x = s_cor_vector)) +
	geom_histogram(bins = 16, alpha = 1) +
	geom_vline(xintercept = mean(s_cor_vector), color = "black", linetype = "dashed", alpha = 0.6) +
	ggtitle("Transcript-protein expression:", subtitle = "Gene-wise correlation") +
	ylab("Gene count") +
	xlab("Spearman correlation") +
	scale_x_continuous(breaks = c(-.65,-.32,0,.32,.65)) +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 12))

ggsave(paste(fig_path,"ProtRNACorrGeneSpearman.png", sep = ""), prot_rna_corr_s_plot, device = png(), width = 4, height = 3)
dev.off()


prot_rna_corr_p_plot <- 
ggplot(cor_pearson, aes(x = p_cor_vector)) +
	geom_histogram(bins = 16, alpha = 1) +
	geom_vline(xintercept = mean(p_cor_vector), color = "black", linetype = "dashed", alpha = 0.6) +
	ggtitle("Transcript-protein expression:", subtitle = "Gene-wise correlation") +
	ylab("Gene count") +
	xlab("Pearson correlation") +
	scale_x_continuous(breaks = c(-.65,-.32,0,.32,.65)) +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 12))

ggsave(paste(fig_path,"ProtRNACorrGenePearson.png", sep = ""), prot_rna_corr_p_plot, device = png(), width = 4, height = 3)
dev.off()

# prot_rna_merge
# columns 2:222 are protein data
# 2:151 MMUT
# 152:222 REST

# columns 223:443 are transcript data
# 223:372 MMUT
# 373:443

#######################################
# correlation of every protein with every transcript of the same gene split in the MMUT and REST group

p_cor_vector_MMUT <- vector()
p_cor_vector_REST <- vector()
s_cor_vector_MMUT <- vector()
s_cor_vector_REST <- vector()

for(i in 1:dim(prot_rna_merge)[1]) {

cat(i, "\n")
p_cor_vector_MMUT <- c(p_cor_vector_MMUT,  cor(unlist(prot_rna_merge[i, 2:151]), unlist(prot_rna_merge[i, 223:372]), use = "pairwise.complete.obs", method = "pearson"))
p_cor_vector_REST <- c(p_cor_vector_REST,  cor(unlist(prot_rna_merge[i, 152:222]), unlist(prot_rna_merge[i, 373:443]), use = "pairwise.complete.obs", method = "pearson"))
s_cor_vector_MMUT <- c(s_cor_vector_MMUT,  cor(unlist(prot_rna_merge[i, 2:151]), unlist(prot_rna_merge[i, 223:372]), use = "pairwise.complete.obs", method = "spearman"))
s_cor_vector_REST <- c(s_cor_vector_REST,  cor(unlist(prot_rna_merge[i, 152:222]), unlist(prot_rna_merge[i, 373:443]), use = "pairwise.complete.obs", method = "spearman"))

}

cor_pearson_MMUT <- data.table(p_cor_vector_MMUT)
cor_pearson_REST <- data.table(p_cor_vector_REST)
cor_spearman_MMUT <- data.table(s_cor_vector_MMUT)
cor_spearman_REST <- data.table(s_cor_vector_REST)


cor_spearman_MMUT1 <- copy(cor_spearman_MMUT)
cor_spearman_MMUT1[, ENSG := prot_rna_merge[, ENSG]]
cor_spearman_REST1 <- copy(cor_spearman_REST)
cor_spearman_REST1[, ENSG := prot_rna_merge[, ENSG]]

cor_spearman_comp <- merge(cor_spearman_MMUT1, cor_spearman_REST1)
cor_spearman_comp[, mean_cor := rowMeans(cor_spearman_comp[,-1])]


ggplot(cor_spearman_comp, aes(x = s_cor_vector_REST, y = s_cor_vector_MMUT)) +
	geom_point(alpha = 0.4) +
	theme_test()


genetbl <- fread("Data/gene_sets/ENSG_HGNC_list_p.csv", header = TRUE)


TCA_PROTEINS_TBL_MAN = "Data/gene_sets/tca_geneset_manual.csv"
tca_tbl = fread(TCA_PROTEINS_TBL_MAN, header = FALSE)[,V1]
# tca_tbl_mod <- unique(rna_prot_merge_fnl[which(rna_prot_merge_fnl$ENSG %in% tca_tbl), ]$ENSG)
mito_tbl = fread("Data/gene_sets/MitoCarta3.0/MitoCarta_highconfidence.csv")[, c("EnsemblGeneID_mapping_version_20200130", "Symbol")]

genes_of_interest <- c("MMUT", "OGDH", "GLUD1", "SUCLA2", "CS", "GLS")
genesOI_tbl <- genetbl[which(genetbl$hgnc_symbol %in% genes_of_interest), ]
genes_of_interest2 <- c("OGDH", "SUCLA2")
genesOI_tbl2 <- genetbl[which(genetbl$hgnc_symbol %in% genes_of_interest2), ]

# add label for genes of interest
for(i in 1:dim(genesOI_tbl)[1]) {
  cor_spearman_comp[which(cor_spearman_comp$ENSG %in% genesOI_tbl[i, 1]), genesOI := as.character(genesOI_tbl[i, 2])]
}
cor_spearman_comp[!is.na(cor_spearman_comp$genesOI), ]

for(i in 1:dim(genetbl)[1]) {
  cor_spearman_comp[which(cor_spearman_comp$ENSG %in% genetbl[i, 1]), gene_symbol := as.character(genetbl[i, 2])]
}

for(i in 1:length(tca_tbl)) {
  cor_spearman_comp[which(cor_spearman_comp$ENSG %in% tca_tbl[i]), tca_genes := as.character(tca_tbl[i])]
}
for(i in 1:dim(mito_tbl)[1]) {
  cor_spearman_comp[which(cor_spearman_comp$ENSG %in% mito_tbl[i, 1]), mito_genes := as.character(mito_tbl[i, 2])]
}
for(i in 1:dim(genesOI_tbl2)[1]) {
  cor_spearman_comp[which(cor_spearman_comp$ENSG %in% genesOI_tbl2[i, 1]), genesOI2 := as.character(genesOI_tbl2[i, 2])]
}
# calculate difference from diagonal (euclidean distance )
cor_spearman_comp[, dist := abs(s_cor_vector_MMUT - s_cor_vector_REST) / sqrt(2)]
# Pythagorean Theorem
c_side <- sqrt(0.32^2+0.32^2)
dist_to_diagonal <- sqrt(0.32^2-(c_side/2)^2)
#check if mito genes are enriched beyond the euclidean distance cut-off
cor_spearman_comp[, dist_cutoff := 0]
cor_spearman_comp[dist > dist_to_diagonal, dist_cutoff := 1]
cor_spearman_comp[, mito_gene_disc := 0]
cor_spearman_comp[!is.na(mito_genes), mito_gene_disc := 1]
table(cor_spearman_comp[, c("dist_cutoff", "mito_gene_disc")])

# genes beyond cut-off
cut_off_tbl <- cor_spearman_comp[dist_cutoff == 1, ]
cut_off_tbl[order(-dist), ]$gene_symbol



# add grid lines to draw on plot
grid_lines <- c(-.65,-.32,.32,.65)

# cor_spearman_comp_col <- melt.data.table(cor_spearman_comp, id.vars = c("ENSG", "genesOI"))
# ggplot(cor_spearman_comp_col, aes(x=value, color = variable)) +
	# geom_density()

# check enrichment of tca genes
ggplot(cor_spearman_comp, aes(x = s_cor_vector_REST, y = s_cor_vector_MMUT, label = tca_genes)) +
	geom_point() +
	geom_point(data = cor_spearman_comp[!is.na(tca_genes), ], shape = 21, color = "white", fill = "blue", size = 3)

# check enrichment of mito genes
ggplot(cor_spearman_comp, aes(x = s_cor_vector_REST, y = s_cor_vector_MMUT, label = mito_genes)) +
	geom_point() +
	geom_point(data = cor_spearman_comp[!is.na(mito_genes), ], shape = 21, color = "white", fill = "blue", size = 3)



group_corr_plt_euclidean <- 
ggplot(cor_spearman_comp, aes(x = s_cor_vector_REST, y = s_cor_vector_MMUT, label = genesOI2, color = dist)) +
	geom_point() +
	geom_vline(xintercept = grid_lines, linetype = 3, alpha = 0.3) +
	geom_hline(yintercept = grid_lines, linetype = 3, alpha = 0.3) +
	geom_abline(slope = 1, intercept = 0, linetype = 2) +
	geom_abline(slope = 1, intercept = .32, linetype = 2) +
	geom_abline(slope = 1, intercept = -.32, linetype = 2) +
	geom_label_repel(size = 3, color = "black") +
	geom_point(data = cor_spearman_comp[!is.na(genesOI2), ], inherit.aes = FALSE, aes(x = s_cor_vector_REST, y = s_cor_vector_MMUT), shape = 21, color = "white", fill = "blue", size = 3) +
	scale_x_continuous(breaks = c(-.65,-.32,0,.32,.65)) +
	scale_y_continuous(breaks = c(-.65,-.32,0,.32,.65)) +
	scale_color_viridis(direction = -1, option = "rocket") +
	labs(y= "MMUT def. correlations", x= "Control correlations", color = "Euclidean\ndistance") +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), legend.position = "bottom", legend.title = element_text(size = 10))



ggsave(paste(fig_path,"ProtRnaCorr_MMUTvsREST_euclidean.png", sep = ""), group_corr_plt_euclidean, device = png(), width = 6, height = 4)
dev.off()
ggsave(paste(fig_path_pdf,"Fig4/","ProtRnaCorr_MMUTvsREST_euclidean.pdf", sep = ""), group_corr_plt_euclidean, device = "pdf", width = 3, height = 3.5)
dev.off()





group_corr_plt <- 
ggplot(cor_spearman_comp, aes(x = s_cor_vector_REST, y = s_cor_vector_MMUT, label = genesOI)) +
	geom_vline(xintercept = grid_lines, linetype = 3, alpha = 0.3) +
	geom_hline(yintercept = grid_lines, linetype = 3, alpha = 0.3) +
	geom_point() +
	geom_point(data=cor_spearman_comp[between(cor_spearman_comp$s_cor_vector_REST, -.32, .32) & between(cor_spearman_comp$s_cor_vector_MMUT, -.32, .32), ], alpha = 1, shape = 16, color = "grey") +
	geom_point(data=cor_spearman_comp[s_cor_vector_REST > .65 | s_cor_vector_MMUT > .65, ], alpha = 1, shape = 16, aes(color = "a")) +
	geom_point(data=cor_spearman_comp[between(cor_spearman_comp$s_cor_vector_REST, .32, .65) & s_cor_vector_MMUT < .65, ], alpha = 1, shape = 16, aes(color = "b")) +
	geom_point(data=cor_spearman_comp[between(cor_spearman_comp$s_cor_vector_MMUT, .32, .65) & s_cor_vector_REST < .32, ], alpha = 1, shape = 16, aes(color = "b")) +
	geom_point(data=cor_spearman_comp[s_cor_vector_REST < -.32 & s_cor_vector_MMUT < -.32, ], alpha = 1, shape = 16, aes(color = "d")) +
	geom_point(data=cor_spearman_comp[between(cor_spearman_comp$s_cor_vector_REST, -.32, 0) & s_cor_vector_MMUT < -.32, ], alpha = 1, shape = 16, aes(color = "c")) +
	geom_point(data=cor_spearman_comp[between(cor_spearman_comp$s_cor_vector_MMUT, -.32, .32) & s_cor_vector_REST < -.32, ], alpha = 1, shape = 16, aes(color = "c")) +
	geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(cor.coef.name = "rho", method = "spearman", r.accuracy = 0.01, size = 3) +
	geom_label_repel(size = 3) +
	geom_point(data = cor_spearman_comp[!is.na(genesOI), ], inherit.aes = FALSE, aes(x = s_cor_vector_REST, y = s_cor_vector_MMUT), color = "white", fill ="black", size = 3, shape = 21) +
	labs(y= "MMUT def. correlations", x= "Control correlations", color = "RNA-protein\ncorrelation") +
	scale_x_continuous(breaks = c(-.65,-.32,0,.32,.65)) +
	scale_y_continuous(breaks = c(-.65,-.32,0,.32,.65)) +
	scale_color_brewer(palette = "PRGn", labels = c("strong pos.", "positive", "negative", "strong neg."), guide = guide_legend(nrow = 2)) +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), legend.position = "bottom", legend.title = element_text(size = 10))


# what's the proportion of significant pairs (and how many)?
all_prots <- dim(cor_spearman_comp)[1]
non_sig_prots <- dim(cor_spearman_comp[between(cor_spearman_comp$s_cor_vector_REST, -.32, .32) & between(cor_spearman_comp$s_cor_vector_MMUT, -.32, .32), ])[1]
sig_prots <- all_prots-non_sig_prots
sig_prots
sig_prots/all_prots

combo_corr_plots <- group_corr_plt + group_corr_plt_euclidean + plot_annotation(title = "Transcript-protein correlation", theme = theme(plot.title = element_text(size = 12)))


ggsave(paste(fig_path,"ProtRnaCorr_MMUTvsREST.png", sep = ""), group_corr_plt, device = png(), width = 6, height = 4)
dev.off()
ggsave(paste(fig_path_pdf,"Fig4/","ProtRnaCorr_MMUTvsREST.pdf", sep = ""), group_corr_plt, device = "pdf", width = 4.5, height = 4.5)
dev.off()



# save combo plot
ggsave(paste(fig_path,"ProtRnaCorr_MMUTvsREST_combo.png", sep = ""), combo_corr_plots, device = png(), width = 6, height = 4)
dev.off()
ggsave(paste(fig_path_pdf,"Fig4/","ProtRnaCorr_MMUTvsREST_combo.pdf", sep = ""), combo_corr_plots, device = "pdf", width = 7.4, height = 4.3)
dev.off()





#######################################
# plot the top and bottom hits (genes) of the spearman correlation

# create protein and transcript merge table

prot_cols <- which(colnames(prot_exp) %in% colnames(rna_exp))
rna_rows <- which(unique(rna_exp$ENSG) %in% unique(prot_exp$ENSG))
prot_rows <- which(unique(prot_exp$ENSG) %in% unique(rna_exp$ENSG))

rna_prot_merge <- merge(prot_exp[prot_rows, ..prot_cols], rna_exp[rna_rows, ], by = "ENSG", all = TRUE)
rna_prot_merge_melt <- melt(rna_prot_merge, id.vars = c("ENSG"))
rna_prot_merge_melt[, molecule := sub("MMA\\d{3}\\.x", "prot_lvl", rna_prot_merge_melt$variable)]
rna_prot_merge_melt[, molecule := sub("MMA\\d{3}\\.y", "rna_lvl", rna_prot_merge_melt$molecule)]
rna_prot_merge_melt[, mma_id := sub("\\..+", "", rna_prot_merge_melt$variable)]
rna_prot_merge_fnl <- dcast(rna_prot_merge_melt, ENSG + mma_id ~ molecule, value.var = "value")
rna_prot_merge_fnl[, sample_no := (as.numeric(sub("MMA", "", rna_prot_merge_fnl$mma_id)))]
rna_prot_merge_fnl[, type := 
  ifelse(sample_no<151, "MMUT def.", 
    ifelse(sample_no>210, "unaffected", "unknown"))]
rna_prot_merge_fnl$type <- factor(rna_prot_merge_fnl$type, levels = c("MMUT def.", "unknown", "unaffected"))
rna_prot_merge_fnl[, type2 := ifelse(sample_no<151, "MMUT def.", "control")]
rna_prot_merge_fnl$type2 <- factor(rna_prot_merge_fnl$type2, levels = c("control", "MMUT def."))
rna_prot_merge_fnl[, type3 := ifelse(sample_no<151, "MMUT def.", "control")]

plot_genes <- data.table(rbind(arrange(cor_spearman_comp %>% top_n(4, wt = mean_cor), desc(mean_cor)), cor_spearman_comp %>% top_n(-4, wt = mean_cor)))$ENSG

plot_dt <- rna_prot_merge_fnl %>% 
	filter(ENSG %in% plot_genes)

HGNC_list_r <- data.table(read.table("Data/gene_sets/ENSG_HGNC_list_r.csv", header = TRUE))
setnames(HGNC_list_r, c("ENSG", "HGNC"))
HGNC_list <- HGNC_list_r[which(HGNC_list_r$ENSG %in% plot_dt$ENSG), ]
HGNC_list <- HGNC_list[match(plot_genes, HGNC_list$ENSG), ]
hgnc_order <- HGNC_list$HGNC
plot_dt <- merge(HGNC_list, plot_dt, by = "ENSG")
plot_dt$ENSG <- factor(plot_dt$ENSG, levels = plot_genes)
plot_dt$HGNC <- factor(plot_dt$HGNC, levels = hgnc_order)
# plot_dt1 <- plot_dt[order(unlist(sapply(plot_dt$ENSG, function(x) which(plot_genes == x)))),]

multiRnaProtPlot <- 
	ggplot(plot_dt, aes(x = rna_lvl, y = prot_lvl, color = type2)) +
	geom_point(size = 0.6, alpha = 0.6) +
 	geom_smooth(method = "lm", se = FALSE) +
  stat_cor(cor.coef.name = "rho", method = "spearman", label.y.npc = 0.15, r.accuracy = 0.01, size = 3, show.legend = FALSE) +
	scale_x_log10() +
	scale_y_log10() +
	annotation_logticks(sides = "bl", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
	facet_wrap(~HGNC, nrow = 2, scales = "free") +
	ggtitle("Transcript-protein expression: strongest correlations") +
	ylab("Protein level") +
	xlab("Transcript level") +
	scale_color_aaas() +
	theme_pubr() +
	rotate_x_text(angle = 90) +
	theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(size = 12), strip.background = element_rect(fill = "white"))

ggsave(paste(fig_path,"top_bottom_spearmanGenes.png", sep = ""), multiRnaProtPlot, device = png(), width = 10, height = 6)
dev.off()
ggsave(paste(fig_path_pdf,"SuppFig5/","top_bottom_spearmanGenes.pdf", sep = ""), multiRnaProtPlot, device = "pdf", width = 10, height = 6)
dev.off()


#######################################
# plot specific genes of interest

# MMUT Ensembl: ENSG00000146085
# OGDH Ensembl: ENSG00000105953
# GLUD1 Ensembl: ENSG00000148672
# SUCLA2: ENSG00000136143
# GLS: ENSG00000115419

plot_genes <- c("ENSG00000146085", "ENSG00000105953", "ENSG00000148672", "ENSG00000136143", "ENSG00000115419", "ENSG00000062485")

plot_dt_imp <- rna_prot_merge_fnl %>% 
	filter(ENSG %in% plot_genes)

HGNC_list_r <- data.table(read.table("Data/gene_sets/ENSG_HGNC_list_r.csv", header = TRUE))
setnames(HGNC_list_r, c("ENSG", "HGNC"))
HGNC_list <- HGNC_list_r[which(HGNC_list_r$ENSG %in% plot_dt_imp$ENSG), ]
HGNC_list <- HGNC_list[match(plot_genes, HGNC_list$ENSG), ]
hgnc_order <- HGNC_list$HGNC
plot_dt_imp <- merge(HGNC_list, plot_dt_imp, by = "ENSG")
plot_dt_imp$ENSG <- factor(plot_dt_imp$ENSG, levels = plot_genes)
plot_dt_imp$HGNC <- factor(plot_dt_imp$HGNC, levels = hgnc_order)

plot_dt_imp <- plot_dt_imp[!(HGNC == "SUCLA2" & rna_lvl < 3), ]

multiRnaProtPlot_important <- 
	ggplot(plot_dt_imp, aes(x = rna_lvl, y = prot_lvl, color = type2)) +
	geom_point(size = 0.6, alpha = 0.6) +
 	geom_smooth(method = "lm", se = FALSE) +
  	stat_cor(cor.coef.name = "rho", method = "spearman", label.y.npc = 0.15, r.accuracy = 0.01, size = 3, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), show.legend = FALSE) +
	scale_x_log10() +
	scale_y_log10() +
	annotation_logticks(sides = "bl", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
	facet_wrap(~HGNC, nrow = 2, scales = "free") +
	ggtitle("Transcript-protein expression") +
	ylab("Protein level") +
	xlab("Transcript level") +
	scale_color_aaas() +
	theme_pubr() +
	rotate_x_text(angle = 90) +
	theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(size = 12), strip.background = element_rect("white"))

ggsave(paste(fig_path,"ProtRNACorr_importantGenes.png", sep = ""), multiRnaProtPlot_important, device = png(), width = 8, height = 5)
dev.off()
ggsave(paste(fig_path_pdf,"SuppFig5/","ProtRNACorr_importantGenes.pdf", sep = ""), multiRnaProtPlot_important, device = "pdf", width = 8, height = 6)
dev.off()


multiAndimportant_plt <- 
plot_grid(multiRnaProtPlot, multiRnaProtPlot_important, align = "h", axis = "bt", rel_widths = c(1.26, 1))
ggsave(paste(fig_path_pdf,"SuppFig5/","ProtRNACorr_multiANDimportant.pdf", sep = ""), multiAndimportant_plt, device = "pdf", width = 18, height = 6)
dev.off()



#######################################
# check genes which showed a strong positive correlation in Quantitative Proteomics of the Cancer Cell Line Encyclopedia, https://doi.org/10.1016/j.cell.2019.12.023

# MET: ENSG00000105976
# WNT5A: ENSG00000114251
# IDH2: ENSG00000182054
# EGFR: ENSG00000146648


plot_genes1 <- c("ENSG00000105976", "ENSG00000114251", "ENSG00000182054", "ENSG00000146648")
# are genes present in merge table?
plot_genes1_mod <- unique(rna_prot_merge_fnl[which(rna_prot_merge_fnl$ENSG %in% plot_genes1), ]$ENSG)

plot_dt1 <- rna_prot_merge_fnl %>% 
	filter(ENSG %in% plot_genes1_mod)

HGNC_list1 <- HGNC_list_r[which(HGNC_list_r$ENSG %in% plot_dt1$ENSG), ]
HGNC_list1 <- HGNC_list1[match(plot_genes1_mod, HGNC_list1$ENSG), ]
hgnc_order1 <- HGNC_list1$HGNC
plot_dt1 <- merge(HGNC_list1, plot_dt1, by = "ENSG")
plot_dt1$ENSG <- factor(plot_dt1$ENSG, levels = plot_genes1_mod)
plot_dt1$HGNC <- factor(plot_dt1$HGNC, levels = hgnc_order1)
# plot_dt1 <- plot_dt[order(unlist(sapply(plot_dt$ENSG, function(x) which(plot_genes == x)))),]

multiRnaProtPlot_Gygi <- 
	ggplot(plot_dt1, aes(x = rna_lvl, y = prot_lvl, color = type2)) +
	geom_point(size = 0.3, alpha = 0.6) +
	# geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  	stat_cor(cor.coef.name = "rho", method = "spearman", label.y.npc = 0.15, r.accuracy = 0.01, size = 3, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
	scale_x_log10() +
	scale_y_log10() +
	facet_wrap(~HGNC, scales = "free") +
	theme_test() +
	rotate_x_text(angle = 90) +
	ggtitle("spearmanGenes_fromGygi") +
	scale_color_aaas()

ggsave(paste(fig_path,"spearmanGenes_fromGygi.png", sep = ""), multiRnaProtPlot_Gygi, device = png(), width = 150, height = 60, units = c("mm"), dpi = 600)
dev.off()


#######################################
# check TCA cycle associated genes

TCA_PROTEINS_TBL_MAN = "Data/gene_sets/tca_geneset_manual.csv"
tca_tbl = fread(TCA_PROTEINS_TBL_MAN, header = FALSE)[,V1]

tca_tbl_mod <- unique(rna_prot_merge_fnl[which(rna_prot_merge_fnl$ENSG %in% tca_tbl), ]$ENSG)

plot_dt2 <- rna_prot_merge_fnl %>% 
	filter(ENSG %in% tca_tbl_mod)

HGNC_list1 <- HGNC_list_r[which(HGNC_list_r$ENSG %in% plot_dt2$ENSG), ]
HGNC_list1 <- HGNC_list1[match(tca_tbl_mod, HGNC_list1$ENSG), ]
hgnc_order2 <- HGNC_list1$HGNC
plot_dt2 <- merge(HGNC_list1, plot_dt2, by = "ENSG")
plot_dt2$ENSG <- factor(plot_dt2$ENSG, levels = tca_tbl_mod)
plot_dt2$HGNC <- factor(plot_dt2$HGNC, levels = hgnc_order2)
# plot_dt2 <- plot_dt[order(unlist(sapply(plot_dt$ENSG, function(x) which(plot_genes == x)))),]

multiRnaProtPlot_TCA <- 
	ggplot(plot_dt2, aes(x = rna_lvl, y = prot_lvl, color = type3)) +
	geom_point(size = 0.3, alpha = 0.6) +
	# geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  	# stat_cor(cor.coef.name = "rho", method = "spearman", label.y.npc = 0.15, r.accuracy = 0.01, size = 3, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
	scale_x_log10() +
	scale_y_log10() +
	annotation_logticks(sides = "bl", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
	facet_wrap(~HGNC, nrow = 4, scales = "free") +
	theme_test() +
	rotate_x_text(angle = 90) +
	ggtitle("spearmanGenes_TCAcycle") +
	scale_color_aaas()

ggsave(paste(fig_path,"spearmanGenes_TCAcycle.png", sep = ""), multiRnaProtPlot_TCA, device = png(), width = 350, height = 150, units = c("mm"), dpi = 300)
dev.off()





# plot MMUT protein vs transcripti correlation
# MMUT Ensembl: ENSG00000146085

# prot/rna correlation color palette
corr_pal <- c("#9c1111" ,"#e86107")

MMUT_corr <- 
ggplot(rna_prot_merge_fnl[ENSG == "ENSG00000146085", ], aes(x = rna_lvl, y = prot_lvl, color = type3)) +
	geom_point(size = 1, alpha = 1) +
	geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  stat_cor(cor.coef.name = "rho", method = "spearman", label.y.npc = .2, r.accuracy = 0.01, size = 3, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
	scale_x_log10() +
	scale_y_log10() +
	xlab("Transcript level") +
	ylab("Protein level") +
	annotation_logticks(sides = "bl", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
	theme_pubr() +
	ggtitle("MMUT expression") +
	scale_color_manual(values = corr_pal) +
	theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(size = 12))

ggsave(paste(fig_path,"MMUTprotrnaCorr.png", sep = ""),MMUT_corr, device = png(), width = 4, height = 2)
dev.off()





#######################################
# protein-protein correlation of selected proteins

# MMUT Ensembl: ENSG00000146085
# ACO1: ENSG00000122729
# ACO2: ENSG00000100412


prot_exp_spec <- melt.data.table(prot_exp, id.vars = "ENSG")
prot_exp_spec <- prot_exp_spec[ENSG == "ENSG00000146085" | ENSG == "ENSG00000122729" | ENSG == "ENSG00000100412", ]
prot_exp_spec <- dcast(prot_exp_spec, variable ~ENSG)

prot_exp_spec[, sampleno := sub("MMA", "", variable)]
prot_exp_spec$sampleno <- as.numeric(prot_exp_spec$sampleno)
prot_exp_spec[, type := ifelse(sampleno <151, "MMUT def.", "control")]
setnames(prot_exp_spec, c("mmaid", "ACO2", "ACO1", "MMUT", "sampleno", "type"))


aco1_plt <- 
ggplot(prot_exp_spec, aes(x = MMUT, y = ACO1, color = type)) +
	geom_point(size = 0.6, alpha = 0.6) +
	geom_smooth(method = "lm", se = FALSE) +
  	stat_cor(cor.coef.name = "rho", method = "spearman", r.accuracy = 0.01, size = 3, aes(x = MMUT, y = ACO1, label = paste(..r.label.., ..p.label.., sep = "~`,`~")), inherit.aes = FALSE, show.legend = FALSE) +
	facet_wrap(type~.) +
	ylab("ACO1 protein level") +
	xlab("MMUT protein level") +
	scale_color_aaas() +
	theme_pubr() +
	theme(legend.position = "right", legend.title = element_blank(), strip.background = element_rect("white")) +
	rotate_x_text(angle = 90)


aco2_plt <- 
ggplot(prot_exp_spec, aes(x = MMUT, y = ACO2, color = type)) +
	geom_point(size = 0.6, alpha = 0.6) +
	geom_smooth(method = "lm", se = FALSE) +
  	stat_cor(cor.coef.name = "rho", method = "spearman", r.accuracy = 0.01, size = 3, aes(x = MMUT, y = ACO2, label = paste(..r.label.., ..p.label.., sep = "~`,`~")), inherit.aes = FALSE, show.legend = FALSE) +
	facet_wrap(type~.) +
	ylab("ACO2 protein level") +
	xlab("MMUT protein level") +
	scale_color_aaas() +
	theme_pubr() +
	theme(legend.position = "right", legend.title = element_blank(), strip.background = element_rect("white")) +
	rotate_x_text(angle = 90)


aco_plt <- 
ggarrange(aco2_plt, aco1_plt, nrow = 2, common.legend = TRUE, legend = "bottom")

ggsave(paste(fig_path,"ACOvsMMUT.png", sep = ""), aco_plt, device = png(), width = 5, height = 6, bg = "white")
dev.off()
ggsave(paste(fig_path_pdf,"SuppFig5/","ACOvsMMUT.pdf", sep = ""), aco_plt, device = "pdf", width = 5, height = 6, bg = "white")
dev.off()








#######################################
# test which samples are missing in the RNA dataset, i.e. did not meet quality control criteria

rna_samps <- as.numeric(gsub("MMA", "", colnames(rna_exp)[2:ncol(rna_exp)]))
seq_allsamps <- c(1:230)
seq_allsamps[!seq_allsamps %in% rna_samps]


# correlation of every protein sample with every transcript sample
# create vectors with correlation data

sample_p_cor_vector <- vector()
sample_s_cor_vector <- vector()

for(i in 1:221) {

cat(i, "\n")
sample_p_cor_vector <- c(sample_p_cor_vector, cor(unlist(prot_rna_merge[, i+1, with = FALSE]), unlist(prot_rna_merge[, i+222, with = FALSE]), use = "pairwise.complete.obs", method = "pearson"))
sample_s_cor_vector <- c(sample_s_cor_vector, cor(unlist(prot_rna_merge[, i+1, with = FALSE]), unlist(prot_rna_merge[, i+222, with = FALSE]), use = "pairwise.complete.obs", method = "spearman"))

}

cor_pearson_smpls <- data.table(sample_p_cor_vector)
cor_spearman_smpls <- data.table(sample_s_cor_vector)


# plot histograms

samples_corr_s_plot <- 
ggplot(cor_spearman_smpls, aes(x = sample_s_cor_vector)) +
	geom_histogram(bins = 16, color = "black", alpha = 0.7) +
	ggtitle("ProtRnaCorr_SampleLevel_spearman") +
	theme_bw()

samples_corr_p_plot <- 
ggplot(cor_pearson_smpls, aes(x = sample_p_cor_vector)) +
	geom_histogram(bins = 16, color = "black", alpha = 0.7) +
	ggtitle("ProtRnaCorr_SampleLevel_pearson") +
	theme_bw()


# add sample info to correlation data and make alternative histogram

cor_spearman_smpls[, mma_id := sub("\\.x", "", colnames(prot_rna_merge)[2:222])]
cor_spearman_smpls[, sample_no := (as.numeric(sub("MMA", "", cor_spearman_smpls$mma_id)))]
cor_spearman_smpls[, type :=
	ifelse(sample_no<151, "MMUT def.",
		ifelse(sample_no>210, "unaffected", "unknown"))]
cor_spearman_smpls$type <- factor(cor_spearman_smpls$type, levels = c("MMUT def.", "unknown", "unaffected"))
cor_spearman_smpls[, type2 :=
	ifelse(sample_no<151, "MMUT def.", "control")]
cor_spearman_smpls$type2 <- factor(cor_spearman_smpls$type2, levels = c("control", "MMUT def."))

mean(cor_spearman_smpls$sample_s_cor_vector)

samples_corr_s_type_plot <- 
ggplot(cor_spearman_smpls, aes(x = sample_s_cor_vector, color = type, fill = type)) +
	geom_histogram(bins = 16, alpha = 1) +
	geom_vline(xintercept = mean(cor_spearman_smpls$sample_s_cor_vector), color = "black", linetype = "dashed", alpha = 0.6) +
	ggtitle("Transcript-protein expression:", subtitle = "Sample-wise correlation") +
	ylab("Sample count") +
	xlab("Spearman correlation") +
	scale_color_npg() +
	scale_fill_npg() +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size =12), legend.position = "right", legend.title = element_blank())

ggsave(paste(fig_path,"ProtRNACorrSampleSpearman.png", sep = ""), samples_corr_s_type_plot, device = png(), width = 4, height = 3)
dev.off()

ggplot(cor_spearman_smpls, aes(x = sample_s_cor_vector, color = type2, fill = type2)) +
	geom_density(alpha = 0.5) +
	scale_fill_aaas() +
	scale_color_aaas()


samples_corr_s_type_plot2 <- 
ggplot(cor_spearman_smpls, aes(x = sample_s_cor_vector, color = type2, fill = type2)) +
	geom_histogram(bins = 16, alpha = 1) +
	geom_vline(xintercept = mean(cor_spearman_smpls$sample_s_cor_vector), color = "black", linetype = "dashed", alpha = 0.6) +
	ggtitle("Transcript-protein correlation") +
	ylab("Sample count") +
	xlab("Spearman correlation") +
	scale_color_aaas() +
	scale_fill_aaas() +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size =12), legend.position = c(0.25, 0.8), legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))

ggsave(paste(fig_path,"ProtRNACorrSampleSpearman2.png", sep = ""), samples_corr_s_type_plot2, device = png(), width = 5, height = 3)
dev.off()


# create combined plot of gene-wise and sample-wise correlation spearman
(prot_rna_corr_s_plot+theme(plot.title = element_blank())) + (samples_corr_s_type_plot2+theme(plot.title = element_blank())) + plot_annotation(title = "Transcript-protein correlation", theme = theme(plot.title = element_text(size = 12)))
ggsave(paste(fig_path,"ProtRNACorrSampleANDGeneSpearman.png", sep = ""), device = png(), width = 6, height = 3)
dev.off()
ggsave(paste(fig_path_pdf,"Fig4/","ProtRNACorrSampleANDGeneSpearman.pdf", sep = ""), device = "pdf", width = 6, height = 3)
dev.off()


cor_pearson_smpls[, type := cor_spearman_smpls$type]
cor_pearson_smpls[, type2 := cor_spearman_smpls$type2]


samples_corr_p_type_plot <- 
ggplot(cor_pearson_smpls, aes(x = sample_p_cor_vector, color = type, fill = type)) +
	geom_histogram(bins = 16, alpha = 1) +
	ggtitle("Transcript-protein expression:", subtitle = "Sample-wise correlation") +
	ylab("Sample count") +
	xlab("Pearson correlation") +
	scale_color_npg() +
	scale_fill_npg() +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size =12), legend.position = "right", legend.title = element_blank())

ggsave(paste(fig_path,"ProtRNACorrSamplePearson.png", sep = ""), samples_corr_p_type_plot, device = png(), width = 4, height = 3)


samples_corr_p_type_plot2 <- 
ggplot(cor_pearson_smpls, aes(x = sample_p_cor_vector, color = type2, fill = type2)) +
	geom_histogram(bins = 16, alpha = 1) +
	geom_vline(xintercept = mean(cor_pearson_smpls$sample_p_cor_vector), color = "black", linetype = "dashed", alpha = 0.6) +
	ggtitle("Transcript-protein expression:", subtitle = "Sample-wise correlation") +
	ylab("Sample count") +
	xlab("Pearson correlation") +
	scale_color_aaas() +
	scale_fill_aaas() +
	scale_x_continuous(breaks = c(-.65,-.32,0,.32,.65)) +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size =12), legend.position = c(0.8, 0.8), legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))



# create combined plot of gene-wise and sample-wise correlation pearson
(prot_rna_corr_p_plot+theme(plot.title = element_blank())) / (samples_corr_p_type_plot2+theme(plot.title = element_blank())) + plot_annotation(title = "Transcript-protein correlation", theme = theme(plot.title = element_text(size = 12)))
ggsave(paste(fig_path_pdf,"SuppFig5/","ProtRNACorrSampleANDGenePearson.pdf", sep = ""), device = "pdf", width = 4, height = 6)
dev.off()




# ggplot(cor_spearman_smpls, aes(x = sample_s_cor_vector, color = type, fill = type)) +
# 	geom_density(alpha = 0.7) +
# 	ggtitle("ProtRnaCorr_SampleLevel_spearman") +
# 	theme_bw()


# arrange plots

overall_corr_plots <- ggarrange(prot_rna_corr_p_plot, prot_rna_corr_s_plot, samples_corr_p_type_plot, samples_corr_s_type_plot)

ggsave(paste(fig_path,"corr_prot_rna.png", sep = ""), overall_corr_plots, device = png(), width = 190, height = 140, units = c("mm"), dpi = 600)
dev.off()







#####################################
# visualize correlation of gene subsets of interest
# corrr package

# MMUT Ensembl: ENSG00000146085
# OGDH Ensembl: ENSG00000105953
# GLUD1 Ensembl: ENSG00000148672
# gene_subset <- c("ENSG00000146085", "ENSG00000105953", "ENSG00000148672")

TCA_PROTEINS_TBL_MAN = "Data/gene_sets/tca_geneset_manual.csv"
tca_tbl = fread(TCA_PROTEINS_TBL_MAN, header = FALSE)[,V1]
class(tca_tbl)


# prepare proteomics dataset
prot_subset <- prot_exp[is.element(prot_exp$ENSG, tca_tbl), ]
prot_subset_t <- dcast(melt(prot_subset), variable ~ ENSG)
prot_subset_t[, variable := NULL]

HGNC_list_p <- data.table(read.table("Data/gene_sets/ENSG_HGNC_list_p.csv", header = TRUE))
setnames(HGNC_list_p, c("ENSG", "HGNC"))

new_names <- HGNC_list_p[which(HGNC_list_p$ENSG %in% colnames(prot_subset_t)), ]$HGNC
setnames(prot_subset_t, new_names)


# prepare transcript dataset
rna_subset <- rna_exp[is.element(rna_exp$ENSG, tca_tbl), ]
rna_subset_t <- dcast(melt(rna_subset), variable ~ ENSG)
rna_subset_t[, variable := NULL]

HGNC_list_r <- data.table(read.table("Data/gene_sets/ENSG_HGNC_list_r.csv", header = TRUE))
setnames(HGNC_list_r, c("ENSG", "HGNC"))

new_names_rna <- HGNC_list_r[which(HGNC_list_r$ENSG %in% colnames(rna_subset_t)), ]$HGNC
setnames(rna_subset_t, new_names_rna)



# calculate correlation matrix with corrr package
prot_corr_all <- correlate(prot_subset_t, method = "spearman")
prot_corr_MMUT <- correlate(prot_subset_t[1:150, ], method = "spearman")
prot_corr_rest <- correlate(prot_subset_t[151:230, ], method = "spearman")

# which row is the last MMUT sample in the RNA table? 143
rna_corr_all <- correlate(rna_subset_t, method = "spearman")
rna_corr_MMUT <- correlate(rna_subset_t[1:143, ], method = "spearman")
rna_corr_rest <- correlate(rna_subset_t[144:nrow(rna_subset_t), ], method = "spearman")


fashion(prot_corr_all)



###########################
# plot protein corrs

tca_tile_plot_all <- rplot(prot_corr_all, colors = c("skyblue1", "white", "indianred2"))
tca_tile_plot_all <- tca_tile_plot_all +
	ggtitle("ProtCorr_tca_subset_ALL") +
	rotate_x_text(angle = 90)

tca_tile_plot_MMUT <- rplot(prot_corr_MMUT, colors = c("skyblue1", "white", "indianred2"))
tca_tile_plot_MMUT <- tca_tile_plot_MMUT +
	ggtitle("ProtCorr_tca_subset_MMUT") +
	rotate_x_text(angle = 90)

tca_tile_plot_rest <- rplot(prot_corr_rest, colors = c("skyblue1", "white", "indianred2"))
tca_tile_plot_rest <- tca_tile_plot_rest +
	ggtitle("ProtCorr_tca_subset_UNKNandCTRL") +
	rotate_x_text(angle = 90)


tca_network_plot_all <- network_plot(prot_corr_all, min_cor = .3, colors = c("skyblue1", "white", "indianred2"))
tca_network_plot_all <- tca_network_plot_all +
	ggtitle("ProtCorr_tca_subset_ALL")

tca_network_plot_MMUT <- network_plot(prot_corr_MMUT, min_cor = .3, colors = c("skyblue1", "white", "indianred2"))
tca_network_plot_MMUT <- tca_network_plot_MMUT +
	ggtitle("ProtCorr_tca_subset_MMUT")

tca_network_plot_rest <- network_plot(prot_corr_rest, min_cor = .3, colors = c("skyblue1", "white", "indianred2"))
tca_network_plot_rest <- tca_network_plot_rest +
	ggtitle("ProtCorr_tca_subset_UNKNandCTRL")



tca_corr_plot_all <- ggarrange(tca_network_plot_all, tca_tile_plot_all)

ggsave(paste(fig_path,"tca_ProtNetworkANDtile_plot_ALL.png", sep = ""), tca_corr_plot_all, device = png(), width = 400, height = 200, units = c("mm"), dpi = 300)
dev.off()


tca_corr_plot_groups <- ggarrange(tca_network_plot_MMUT, tca_tile_plot_MMUT, tca_network_plot_rest, tca_tile_plot_rest)

ggsave(paste(fig_path,"tca_ProtNetworkANDtile_plot_MMUTvsUNKNandCTRL.png", sep = ""), tca_corr_plot_groups, device = png(), width = 400, height = 280, units = c("mm"), dpi = 300)
dev.off()





###########################
# plot RNA corrs

tca_tile_plot_all <- rplot(rna_corr_all, colors = c("skyblue1", "white", "indianred2"))
tca_tile_plot_all <- tca_tile_plot_all +
	ggtitle("RNACorr_tca_subset_ALL") +
	rotate_x_text(angle = 90)

tca_tile_plot_MMUT <- rplot(rna_corr_MMUT, colors = c("skyblue1", "white", "indianred2"))
tca_tile_plot_MMUT <- tca_tile_plot_MMUT +
	ggtitle("RNACorr_tca_subset_MMUT") +
	rotate_x_text(angle = 90)

tca_tile_plot_rest <- rplot(rna_corr_rest, colors = c("skyblue1", "white", "indianred2"))
tca_tile_plot_rest <- tca_tile_plot_rest +
	ggtitle("RNACorr_tca_subset_UNKNandCTRL") +
	rotate_x_text(angle = 90)


tca_network_plot_all <- network_plot(rna_corr_all, min_cor = .5, colors = c("skyblue1", "white", "indianred2"))
tca_network_plot_all <- tca_network_plot_all +
	ggtitle("RNACorr_tca_subset_ALL")

tca_network_plot_MMUT <- network_plot(rna_corr_MMUT, min_cor = .5, colors = c("skyblue1", "white", "indianred2"))
tca_network_plot_MMUT <- tca_network_plot_MMUT +
	ggtitle("RNACorr_tca_subset_MMUT")

tca_network_plot_rest <- network_plot(rna_corr_rest, min_cor = .5, colors = c("skyblue1", "white", "indianred2"))
tca_network_plot_rest <- tca_network_plot_rest +
	ggtitle("RNACorr_tca_subset_UNKNandCTRL")



tca_corr_plot_all <- ggarrange(tca_network_plot_all, tca_tile_plot_all)

ggsave(paste(fig_path,"tca_RNAnetworkANDtile_plot_ALL.png", sep = ""), tca_corr_plot_all, device = png(), width = 400, height = 200, units = c("mm"), dpi = 300)
dev.off()


tca_corr_plot_groups <- ggarrange(tca_network_plot_MMUT, tca_tile_plot_MMUT, tca_network_plot_rest, tca_tile_plot_rest)

ggsave(paste(fig_path,"tca_RNAnetworkANDtile_plot_MMUTvsUNKNandCTRL.png", sep = ""), tca_corr_plot_groups, device = png(), width = 400, height = 280, units = c("mm"), dpi = 300)
dev.off()