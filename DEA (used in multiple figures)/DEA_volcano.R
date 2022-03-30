##########################################################################
# volcano plot based on the DEA lmm results
##########################################################################

# required libraries

require(tidyverse)
require(data.table)
require(ggsci)
require(ggpubr)
library(ggrepel)
library(patchwork)

# colors
mypal <- pal_aaas("default", alpha = 1)(9)
mypal2 <- pal_npg("nrc", alpha = 1)(9)


##################################
# create paths
system("mkdir Figs/DEA/")
fig_path <- c("Figs/DEA/")
fig_path_pdf <- c("Figs/v8/pdf/")

##################################
# import DEA lmm output files

# rna_DEA1 <- fread("interimData/DEA/diff_exp_rnaseq_pathwayact_all_annotout.txt")
rna_DEA1 <- fread("interimData/DEAlmm/diff_exp_data_rnaseq_pathwayact_all.txt")
rna_DEA_cond <- fread("interimData/DEA/diff_exp_rnaseq_pathwayact_conditional_annotout.txt")

# prot_DEA1 <- fread("interimData/DEA/diff_exp_prot_pathwayact_all_annotout.txt")
prot_DEA1 <- fread("interimData/DEAlmm/diff_exp_data_prot_pathwayact_all.txt")
prot_DEA_cond <- fread("interimData/DEA/diff_exp_prot_pathwayact_conditional_annotout.txt")


# David's modifications on file

enzymes = fread("interimData/gene_sets/all_enzyme_ensg.tbl",header=FALSE)
prot_DEA1[,isenzyme:=is.element(ensembl,enzymes[,V4])]
rna_DEA1[,isenzyme:=is.element(ensembl,enzymes[,V4])]

pulld = fread("interimData/gene_sets/pulldown_mmm.txt",header=FALSE)
myg = pulld[V2<0.05,V4]
prot_DEA1[,ispulldown:=is.element(ensembl,myg)]
rna_DEA1[,ispulldown:=is.element(ensembl,myg)]

mitog = fread("interimData/gene_sets/mitoglist.txt",header=FALSE)[,V1]
tcaprot = fread("interimData/gene_sets/tca_proteins.txt", header=FALSE)[,V1]

prot_DEA1[,ismitog:=is.element(ensembl,mitog)]
prot_DEA1[,istca:=is.element(ensembl,tcaprot)]
prot_DEA1[,inpulldown:="not in pulldown"]
prot_DEA1[ispulldown==TRUE,inpulldown:="in pulldown"]
prot_DEA1[gname=="MMUT",inpulldown:="is MMUT"]
prot_DEA1[,p_theor:=rank(pval)/(dim(prot_DEA1)[1]+1)]
prot_DEA1[,pnorm_theor:=rank(pvalnorm)/(dim(prot_DEA1)[1]+1)]

prot_DEA=prot_DEA1[ismitog==TRUE & isenzyme==TRUE,]
prot_DEA[,p_theor:=rank(pval)/(dim(prot_DEA)[1]+1)]


rna_DEA1[,ismitog:=is.element(ensembl,mitog)]
rna_DEA1[,istca:=is.element(ensembl,tcaprot)]
rna_DEA1[,inpulldown:="not in pulldown"]
rna_DEA1[ispulldown==TRUE,inpulldown:="in pulldown"]
rna_DEA1[gname=="MMUT",inpulldown:="is MMUT"]
rna_DEA1[,p_theor:=rank(pval)/(dim(rna_DEA1)[1]+1)]
rna_DEA1[,pnorm_theor:=rank(pvalnorm)/(dim(rna_DEA1)[1]+1)]

rna_DEA=rna_DEA1[ismitog==TRUE & isenzyme==TRUE,]
rna_DEA[,p_theor:=rank(pval)/(dim(rna_DEA)[1]+1)]




### POSITIVE beta value --> RNA is up in samples with high pathway activity
### NEGATIVE beta value --> RNA is down in samples with high pathway activity
### samples with low pathway activity are mostly MMUT deficient samples
### betas and pval according to DL's fast_lmm model

# use automated enhanced volcano package
# EnhancedVolcano(rna_DEA,
#     lab = rna_DEA$gname,
#     x = 'betas',
#     y = 'pval',
#     pCutoff = 10e-2,
#     FCcutoff = 0.2,
#     xlim = c(-1, 1),
#     ylim = c(0,6),
#     title = "MMAomics volcano plot",
#     subtitle = "Differential expression",
#     xlab = bquote('Effect size (betas)'))


# create gene lists of interest
tca_genes_dt <- fread("Data/gene_sets/tca_geneset_manual_wHGNC.csv", header = TRUE)
tca_genes <- tca_genes_dt$hgnc_symbol
genes_of_interest <- c("MMUT", "OGDH", "GLUD1", "GLUD2", "PDK4", "DLD", "GOT1", "GOT2", "DLST")
genes_of_interest2_rna_up <- c("ALDH2", "PDK4", "SUCLA2")
# genes_of_interest2_rna_up <- c("PDK4")
genes_of_interest2_rna_down <- c("MMUT")
# genes_of_interest2_prot_up <- c("GLUD1", "ALDH2")
genes_of_interest2_prot_up <- c("GLUD1")
genes_of_interest2_prot_down <- c("MMUT", "OGDH")
# genes_label <- c("MMUT", "OGDH", "GLUD1", "ALDH2")
genes_label <- c("MMUT", "OGDH", "GLUD1", "ALDH2")
genes_label_rna <- c("MMUT", "PDK4", "ALDH2", "SUCLA2")


# prep RNA table
for(i in 1:length(tca_genes)) {
  rna_DEA[which(rna_DEA$gname %in% tca_genes[i]), label := as.character(tca_genes[i])]
}

for(i in 1:length(genes_of_interest)) {
  rna_DEA[which(rna_DEA$gname %in% genes_of_interest[i]), label1 := as.character(genes_of_interest[i])]
}

for(i in 1:length(genes_of_interest2_rna_up)) {
  rna_DEA[which(rna_DEA$gname %in% genes_of_interest2_rna_up[i]), label_rnaUP := as.character(genes_of_interest2_rna_up[i])]
}

for(i in 1:length(genes_of_interest2_rna_down)) {
  rna_DEA[which(rna_DEA$gname %in% genes_of_interest2_rna_down[i]), label_rnaDOWN := as.character(genes_of_interest2_rna_down[i])]
}

for(i in 1:length(genes_label_rna)) {
  rna_DEA[which(rna_DEA$gname %in% genes_label_rna[i]), label4 := as.character(genes_label_rna[i])]
}


rna_DEA[pval<0.001, label3 := gname]
rna_DEA[!is.na(label), ]
rna_DEA[!is.na(label1), ]
rna_DEA[!is.na(label3), ]


# prep protein table
for(i in 1:length(tca_genes)) {
  prot_DEA[which(prot_DEA$gname %in% tca_genes[i]), label := as.character(tca_genes[i])]
}

for(i in 1:length(genes_of_interest)) {
  prot_DEA[which(prot_DEA$gname %in% genes_of_interest[i]), label1 := as.character(genes_of_interest[i])]
}

for(i in 1:length(genes_of_interest2_prot_up)) {
  prot_DEA[which(prot_DEA$gname %in% genes_of_interest2_prot_up[i]), label_protUP := as.character(genes_of_interest2_prot_up[i])]
}

for(i in 1:length(genes_of_interest2_prot_down)) {
  prot_DEA[which(prot_DEA$gname %in% genes_of_interest2_prot_down[i]), label_protDOWN := as.character(genes_of_interest2_prot_down[i])]
}

for(i in 1:length(genes_label)) {
  prot_DEA[which(prot_DEA$gname %in% genes_label[i]), label4 := as.character(genes_label[i])]
}

prot_DEA[pval<0.001, label3 := gname]
prot_DEA[!is.na(label), ]
prot_DEA[!is.na(label1), ]
prot_DEA[!is.na(label3), ]




##################################
# plots volcanos for RNA and protein (limited genes of interest)


# rna_FDR_pulldown <- 
# ggplot(rna_DEA, aes(x = -log10(p_theor), y = -log10(pval), color = inpulldown)) +
#   geom_abline(intercept = 0, slope = 1) +
#   geom_abline(aes(intercept = 1, slope = 1, color = "FDR<0.1"), linetype = 2) +
#   geom_point() +
#   xlab("-log10(p-value) (theoretical)") +
#   ylab("-log10(p-value) (empirical)") +
#   ggtitle("DEA (PI activity) transcriptomics") +
#   geom_text_repel(data = rna_DEA[!is.na(label4), ], aes(label = gname), nudge_x = 0, nudge_y = 0.5, show.legend = FALSE) +
#   scale_color_manual(values = c("black", mypal[c(5,2,9)]), guide = guide_legend(override.aes = list(pch = c(NA, 16, 16, 16), linetype = c(2, 0,0,0)))) +
#   theme_pubr() +
#   theme(plot.title = element_text(size = 12), legend.position = c(0.26, 0.8), legend.title = element_blank(), legend.background = element_rect(fill = alpha('blue', 0)))

rna_FDR <- 
ggplot(rna_DEA, aes(x = -log10(p_theor), y = -log10(pval))) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(aes(intercept = 1, slope = 1, color = "FDR<0.1"), linetype = 2) +
  geom_point(color = "lightgrey") +
  geom_point(data = rna_DEA[!is.na(label4), ], color = "black", size = 3) +
  xlab("-log10(p-value) (theoretical)") +
  ylab("-log10(p-value) (empirical)") +
  ggtitle("DEA (PI activity) transcriptomics") +
  geom_text_repel(data = rna_DEA[!is.na(label4), ], aes(label = gname), nudge_x = 0, nudge_y = 0.5, show.legend = FALSE) +
  theme_pubr() +
  theme(plot.title = element_text(size = 12), legend.position = c(0.26, 0.9), legend.title = element_blank(), legend.background = element_rect(fill = alpha('blue', 0)))


# prot_FDR_pulldown <- 
# ggplot(prot_DEA, aes(x = -log10(p_theor), y = -log10(pval), color = inpulldown)) +
#   geom_abline(intercept = 0, slope = 1) +
#   geom_abline(aes(intercept = 1, slope = 1, color = "FDR<0.1"), linetype = 2) +
#   geom_point() +
#   xlab("-log10(p-value) (theoretical)") +
#   ylab("-log10(p-value) (empirical)") +
#   ggtitle("DEA (PI activity) proteomics") +
#   geom_text_repel(data = prot_DEA[!is.na(label4), ], aes(label = gname), nudge_x = 0, nudge_y = 0.5, show.legend = FALSE) +
#   scale_color_manual(values = c("black", mypal[c(5,2,9)]), guide = guide_legend(override.aes = list(pch = c(NA, 16, 16, 16), linetype = c(2, 0,0,0)))) +
#   theme_pubr() +
#   theme(plot.title = element_text(size = 12), legend.position = c(0.26, 0.8), legend.title = element_blank(), legend.background = element_rect(fill = alpha('blue', 0)))

prot_FDR <- 
ggplot(prot_DEA, aes(x = -log10(p_theor), y = -log10(pval))) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(aes(intercept = 1, slope = 1, color = "FDR<0.1"), linetype = 2) +
  geom_point(color = "lightgrey") +
  geom_point(data = prot_DEA[!is.na(label4), ], color = "black", size = 3) +
  xlab("-log10(p-value) (theoretical)") +
  ylab("-log10(p-value) (empirical)") +
  ggtitle("DEA (PI activity) proteomics") +
  geom_text_repel(data = prot_DEA[!is.na(label4), ], aes(label = gname), nudge_x = 0, nudge_y = 0.5, show.legend = FALSE) +
  theme_pubr() +
  theme(plot.title = element_text(size = 12), legend.position = c(0.26, 0.9), legend.title = element_blank(), legend.background = element_rect(fill = alpha('blue', 0)))




p0 <- 0.01


# dea_trans_pulldown <- 
# ggplot(rna_DEA, aes(x = -betas, y = -log10(pval), color = inpulldown)) +
#   geom_hline(yintercept = -log10(p0), linetype = "dashed", alpha = 0.6) +
#   geom_vline(xintercept = 0.2, linetype = "dashed", alpha = 0.6) +
#   geom_vline(xintercept = -0.2, linetype = "dashed", alpha = 0.6) +
#   geom_point(alpha = 0.2) +
#   geom_point(data = rna_DEA[!is.na(label_rnaUP), ], size = 3) +
#   geom_text_repel(data = rna_DEA[!is.na(label_rnaUP), ], aes(label = gname), nudge_x = 0.1, nudge_y = 0.5) +
#   geom_point(data = rna_DEA[!is.na(label_rnaDOWN), ], size = 3) +
#   geom_text_repel(data = rna_DEA[!is.na(label_rnaDOWN), ], aes(label = gname), nudge_x = 0, nudge_y = 0.5) +
#   xlab("Effect size") +
#   ylab("-log10(adjusted p value)") +
#   xlim(-0.45, 0.45) +
#   ylim(0, 6) +
#   scale_color_manual(values = mypal[c(5,2,9)]) +
#   theme_pubr() +
#   theme(plot.title = element_text(size = 12), legend.position = "none", legend.title = element_blank())

dea_trans <- 
ggplot(rna_DEA, aes(x = -betas, y = -log10(pval))) +
  geom_hline(yintercept = -log10(p0), linetype = "dashed", alpha = 0.6) +
  geom_vline(xintercept = 0.2, linetype = "dashed", alpha = 0.6) +
  geom_vline(xintercept = -0.2, linetype = "dashed", alpha = 0.6) +
  geom_point(color = "grey", alpha = 0.2) +
  geom_point(data = rna_DEA[!is.na(label4), ], size = 3) +
  geom_text_repel(data = rna_DEA[!is.na(label4), ], aes(label = gname), nudge_x = 0.1, nudge_y = 0.5) +
  xlab("Effect size") +
  ylab("-log10(adjusted p value)") +
  xlim(-0.45, 0.45) +
  ylim(0, 6) +
  scale_color_manual(values = mypal[c(5,2,9)]) +
  theme_pubr() +
  theme(plot.title = element_text(size = 12), legend.position = "none", legend.title = element_blank())


# dea_prot_pulldown <- 
# ggplot(prot_DEA, aes(x = -betas, y = -log10(pval), color = inpulldown)) +
#   geom_hline(yintercept = -log10(p0), linetype = "dashed", alpha = 0.6) +
#   geom_vline(xintercept = 0.2, linetype = "dashed", alpha = 0.6) +
#   geom_vline(xintercept = -0.2, linetype = "dashed", alpha = 0.6) +
#   geom_point(alpha = 0.2) +
#   geom_point(data = prot_DEA[!is.na(label_protUP), ], size = 3) +
#   geom_text_repel(data = prot_DEA[!is.na(label_protUP), ], aes(label = gname), nudge_x = -0.1, nudge_y = 0.5) +
#   geom_point(data = prot_DEA[!is.na(label_protDOWN), ], size = 3) +
#   geom_text_repel(data = prot_DEA[!is.na(label_protDOWN), ], aes(label = gname), nudge_x = 0, nudge_y = 0.5) +
#   xlab("Effect size") +
#   ylab("-log10(adjusted p value)") +
#   xlim(-0.45, 0.45) +
#   ylim(0, 6) +
#   scale_color_manual(values = mypal[c(5,2,9)]) +
#   theme_pubr() +
#   theme(plot.title = element_text(size = 12), legend.position = "none", legend.title = element_blank())

dea_prot <- 
ggplot(prot_DEA, aes(x = -betas, y = -log10(pval))) +
  geom_hline(yintercept = -log10(p0), linetype = "dashed", alpha = 0.6) +
  geom_vline(xintercept = 0.2, linetype = "dashed", alpha = 0.6) +
  geom_vline(xintercept = -0.2, linetype = "dashed", alpha = 0.6) +
  geom_point(color = "grey", alpha = 0.2) +
  geom_point(data = prot_DEA[!is.na(label4), ], size = 3) +
  geom_text_repel(data = prot_DEA[!is.na(label4), ], aes(label = gname), nudge_x = -0.1, nudge_y = 0.5) +
  xlab("Effect size") +
  ylab("-log10(adjusted p value)") +
  xlim(-0.45, 0.45) +
  ylim(0, 6) +
  scale_color_manual(values = mypal[c(5,2,9)]) +
  theme_pubr() +
  theme(plot.title = element_text(size = 12), legend.position = "none", legend.title = element_blank())




# dea_trans <- 
# ggplot(rna_DEA, aes(x = -betas, y = -log10(pval))) +
#   geom_hline(yintercept = -log10(p0), linetype = "dashed", alpha = 0.6) +
#   geom_vline(xintercept = 0.2, linetype = "dashed", alpha = 0.6) +
#   geom_vline(xintercept = -0.2, linetype = "dashed", alpha = 0.6) +
#   geom_point(alpha = 0.2) +
#   geom_point(data = rna_DEA[!is.na(label_rnaUP), ], inherit.aes = FALSE, aes(x = -betas, y = -log10(pval)), color = mypal2[1], size = 3) +
#   geom_text_repel(data = rna_DEA[!is.na(label_rnaUP), ], aes(label = gname), color = mypal2[1], nudge_x = 0.1, nudge_y = 0.5) +
#   geom_point(data = rna_DEA[!is.na(label_rnaDOWN), ], inherit.aes = FALSE, aes(x = -betas, y = -log10(pval)), color = mypal2[3], size = 3) +
#   geom_text_repel(data = rna_DEA[!is.na(label_rnaDOWN), ], aes(label = gname), color = mypal2[3], nudge_x = 0, nudge_y = 0.5) +
#   ggtitle("DEA (PI activity) transcriptomics") +
#   xlab("Effect size") +
#   ylab("-log10(adjusted p value)") +
#   xlim(-0.45, 0.45) +
#   ylim(0, 6) +
#   theme_pubr() +
#   theme(plot.title = element_text(size = 12))


# dea_prot <- 
# ggplot(prot_DEA, aes(x = -betas, y = -log10(pval))) +
#   geom_hline(yintercept = -log10(p0), linetype = "dashed", alpha = 0.6) +
#   geom_vline(xintercept = 0.2, linetype = "dashed", alpha = 0.6) +
#   geom_vline(xintercept = -0.2, linetype = "dashed", alpha = 0.6) +
#   geom_point(alpha = 0.2) +
#   geom_point(data = prot_DEA[!is.na(label_protUP), ], inherit.aes = FALSE, aes(x = -betas, y = -log10(pval)), color = mypal2[1], size = 3) +
#   geom_text_repel(data = prot_DEA[!is.na(label_protUP), ], aes(label = gname), color = mypal2[1], nudge_x = -0.1, nudge_y = 0.5) +
#   geom_point(data = prot_DEA[!is.na(label_protDOWN), ], inherit.aes = FALSE, aes(x = -betas, y = -log10(pval)), color = mypal2[3], size = 3) +
#   geom_text_repel(data = prot_DEA[!is.na(label_protDOWN), ], aes(label = gname), color = mypal2[3], nudge_x = 0, nudge_y = 0.5) +
#   ggtitle("DEA (PI activity) proteomics") +
#   xlab("Effect size") +
#   ylab("-log10(adjusted p value)") +
#   xlim(-0.45, 0.45) +
#   ylim(0, 6) +
#   theme_pubr() +
#   theme(plot.title = element_text(size = 12))


dea_plt_main <- 
ggarrange(rna_FDR, dea_trans, ncol = 1) + plot_spacer() + ggarrange(prot_FDR, dea_prot, ncol = 1) +  plot_layout(widths = c(1, 0.15, 1))


ggsave(paste(fig_path,"DEAVolcanoMainFDR.png", sep = ""),dea_plt_main, device = png(), width = 7, height = 7, bg = "white")
dev.off()
ggsave(paste(fig_path_pdf, "Fig3/","DEAVolcanoMainFDR.pdf", sep = ""),dea_plt_main, device = "pdf", width = 7, height = 7, bg = "white")
dev.off()


dea_plt_main_rna <- 
rna_FDR+theme(plot.title = element_blank()) + dea_trans + plot_annotation(title = "Differential expression analysis (PI activity) transcriptomics", theme = theme(plot.title = element_text(size = 12)))
ggsave(paste(fig_path_pdf, "Fig3/","DEAVolcanoMainFDR_rna.pdf", sep = ""),dea_plt_main_rna, device = "pdf", width = 6, height = 3, bg = "white")
dev.off()

dea_plt_main_prot <- 
prot_FDR+theme(plot.title = element_blank()) + dea_prot + plot_annotation(title = "Differential expression analysis (PI activity) proteomics", theme = theme(plot.title = element_text(size = 12)))
ggsave(paste(fig_path_pdf, "Fig3/","DEAVolcanoMainFDR_prot.pdf", sep = ""),dea_plt_main_prot, device = "pdf", width = 6, height = 3, bg = "white")
dev.off()



# alternative protein plot (different color)

# dea color palette
dea_pal <- c( "#084594","#4292C6")


dea_prot <-
ggplot(prot_DEA, aes(x = -betas, y = -log10(pval))) +
  geom_point(alpha = 0.2) +
  geom_point(data = prot_DEA[gname == "MMUT" | gname == "GLUD1" | gname == "OGDH", ], inherit.aes = FALSE, aes(x = -betas, y = -log10(pval)), color = dea_pal[2], size = 3) +
  geom_text_repel(data = prot_DEA[gname == "MMUT" | gname == "GLUD1" | gname == "OGDH", ], aes(label = gname), color = dea_pal[2], nudge_x = 0, nudge_y = 1) +
  geom_hline(yintercept = -log10(p0), linetype = "dashed", alpha = 0.6) +
  ggtitle("Differential protein expression") +
  xlab("Expression") +
  ylab("-log10(adjusted p value)") +
  theme_pubr() +
  theme(plot.title = element_text(size = 12))

ggsave(paste(fig_path,"DEAVolcanoProtOnly.png", sep = ""),dea_prot, device = png(), width = 4, height = 2)
dev.off()





##################################
# plots volcanos for RNA and protein
p0 <- 0.01

dea_plt_rna_p <- 
ggplot(rna_DEA, aes(x = -betas, y = -log(pval), label = label3)) +
  geom_point(alpha = 0.1) +
  geom_text_repel() +
  geom_point(data = rna_DEA[!is.na(label3), ], inherit.aes = FALSE, aes(x = -betas, y = -log(pval), color = mypal[2])) +
  guides(color = FALSE) +
  ggtitle("transcripts_DEA_lowPvals", subtitle = "positive effect size == low OHCblPlus") +
  theme_test()

dea_plt_rna <- 
ggplot(rna_DEA, aes(x = -betas, y = -log(pval), label = label1)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = -log(p0)) +
  geom_text(aes(0 , -log(p0), label = paste("pval=",p0, sep = ""), vjust = -1)) +
  geom_label_repel() +
  geom_point(data = rna_DEA[!is.na(label1), ], inherit.aes = FALSE, aes(x = -betas, y = -log(pval), color = mypal[2])) +
  guides(color = FALSE) +
  ggtitle("transcripts_DEA_interestingGenes", subtitle = "positive effect size == low OHCblPlus") +
  theme_test()


dea_plt_prot_p <- 
ggplot(prot_DEA, aes(x = -betas, y = -log(pval), label = label3)) +
  geom_point(alpha = 0.1) +
  geom_text_repel() +
  geom_point(data = prot_DEA[!is.na(label3), ], inherit.aes = FALSE, aes(x = -betas, y = -log(pval), color = mypal[2])) +
  guides(color = FALSE) +
  ggtitle("prot_DEA_lowPvals", subtitle = "positive effect size == low OHCblPlus") +
  theme_test()

dea_plt_prot <- 
ggplot(prot_DEA, aes(x = -betas, y = -log(pval), label = label1)) +
  geom_point(alpha = 0.1) +
  # geom_hline(yintercept = -log(p0)) +
  # geom_text(aes(0 , -log(p0), label = paste("pval=",p0, sep = ""), vjust = -1)) +
  geom_label_repel() +
  geom_point(data = prot_DEA[!is.na(label1), ], inherit.aes = FALSE, aes(x = -betas, y = -log(pval), color = mypal[2])) +
  guides(color = FALSE) +
  ggtitle("prot_DEA_interestingGenes", subtitle = "positive effect size == low OHCblPlus") +
  theme_test()

volcanos <- ggarrange(dea_plt_prot_p, dea_plt_prot, dea_plt_rna_p, dea_plt_rna)

ggsave(paste(fig_path,"DEA_volcanos.png", sep = ""), volcanos, width = 7, height = 8)




##################################
# same plots based on clinical severity score instead of OHCblPlus


rna_DEA_css <- fread("interimData/DEAlmm/diff_exp_data_rnaseq_wcss_all.txt")

prot_DEA_css <- fread("interimData/DEAlmm/diff_exp_data_prot_wcss_all.txt")

### POSITIVE beta value --> RNA is up in samples with high pathway activity
### NEGATIVE beta value --> RNA is down in samples with high pathway activity
### samples with low pathway activity are mostly MMUT deficient samples
### betas and pval according to DL's fast_lmm model


# create gene lists of interest
tca_genes_dt <- fread("Data/gene_sets/tca_geneset_manual_wHGNC.csv", header = TRUE)
tca_genes <- tca_genes_dt$hgnc_symbol
genes_of_interest <- c("MMUT", "OGDH", "GLUD1", "GLUD2", "PDK4", "DLD", "GOT1", "GOT2", "DLST")

# prep RNA table
for(i in 1:length(tca_genes)) {
  rna_DEA_css[which(rna_DEA_css$gname %in% tca_genes[i]), label := as.character(tca_genes[i])]
}

for(i in 1:length(genes_of_interest)) {
  rna_DEA_css[which(rna_DEA_css$gname %in% genes_of_interest[i]), label1 := as.character(genes_of_interest[i])]
}

rna_DEA_css[pval<0.001, label3 := gname]
rna_DEA_css[!is.na(label), ]
rna_DEA_css[!is.na(label1), ]
rna_DEA_css[!is.na(label3), ]


# prep protein table
for(i in 1:length(tca_genes)) {
  prot_DEA_css[which(prot_DEA_css$gname %in% tca_genes[i]), label := as.character(tca_genes[i])]
}

for(i in 1:length(genes_of_interest)) {
  prot_DEA_css[which(prot_DEA_css$gname %in% genes_of_interest[i]), label1 := as.character(genes_of_interest[i])]
}

prot_DEA_css[pval<0.005, label3 := gname]
prot_DEA_css[!is.na(label), ]
prot_DEA_css[!is.na(label1), ]
prot_DEA_css[!is.na(label3), ]


# plots volcanos for RNA and protein
p0 <- 0.01

dea_plt_rna_p <- 
ggplot(rna_DEA_css, aes(x = betas, y = -log(pval), label = label3)) +
  geom_point(alpha = 0.1) +
  geom_text_repel() +
  geom_point(data = rna_DEA_css[!is.na(label3), ], inherit.aes = FALSE, aes(x = betas, y = -log(pval), color = mypal[2])) +
  guides(color = FALSE) +
  ggtitle("RNA_DEA_lowPvals_cnCSS", subtitle = "positive effect size == high CSS") +
  theme_test()

dea_plt_rna <- 
ggplot(rna_DEA_css, aes(x = betas, y = -log(pval), label = label1)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = -log(p0)) +
  geom_text(aes(0 , -log(p0), label = paste("pval=",p0, sep = ""), vjust = -1)) +
  geom_label_repel() +
  geom_point(data = rna_DEA_css[!is.na(label1), ], inherit.aes = FALSE, aes(x = betas, y = -log(pval), color = mypal[2])) +
  guides(color = FALSE) +
  ggtitle("RNA_DEA_interestingGenes_cnCSS", subtitle = "positive effect size == high CSS") +
  theme_test()


dea_plt_prot_p <- 
ggplot(prot_DEA_css, aes(x = betas, y = -log(pval), label = label3)) +
  geom_point(alpha = 0.1) +
  geom_text_repel() +
  geom_point(data = prot_DEA_css[!is.na(label3), ], inherit.aes = FALSE, aes(x = betas, y = -log(pval), color = mypal[2])) +
  guides(color = FALSE) +
  ggtitle("prot_DEA_lowPvals_cnCSS", subtitle = "positive effect size == high CSS") +
  theme_test()

dea_plt_prot <- 
ggplot(prot_DEA_css, aes(x = betas, y = -log(pval), label = label1)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = -log(p0)) +
  # geom_text(aes(0 , -log(p0), label = paste("pval=",p0, sep = ""), vjust = -1)) +
  geom_label_repel() +
  geom_point(data = prot_DEA_css[!is.na(label1), ], inherit.aes = FALSE, aes(x = betas, y = -log(pval), color = mypal[2])) +
  guides(color = FALSE) +
  ggtitle("prot_DEA_interestingGenes_cnCSS", subtitle = "positive effect size == high CSS") +
  theme_test()

volcanos <- ggarrange(dea_plt_prot_p, dea_plt_prot, dea_plt_rna_p, dea_plt_rna)

ggsave(paste(fig_path,"DEA_volcanos_newcCSS.png", sep = ""), volcanos, width = 7, height = 8)


