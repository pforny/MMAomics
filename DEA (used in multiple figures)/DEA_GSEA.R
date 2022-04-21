

# libraries
library(fgsea)
library(data.table)
library(tidyverse)
library(pheatmap)
library(qusage)
library(ggsci)
library(ggpubr)
library(viridis)


# create figures path
system("mkdir Figs/DEA/")
system("mkdir Figs/DEA/GSEA/")
fig_path <- c("Figs/DEA/GSEA/")
fig_path_pdf <- c("Figs/v8/pdf/")



##################################################
# protein analysis


dea_prot <- fread("interimData/DEAlmm/diff_exp_data_prot_pathwayact_all.txt")


# dea_prot <- dea_prot[pval < 0.05, ]
setorder(dea_prot, -betas)
dea_prot[, order := 1:nrow(dea_prot)]

generanks_plt <- 
ggplot(dea_prot, aes(x = order, y = betas)) +
	geom_bar(stat = "identity") +
	xlab("gene rank") +
	ggtitle("Genes ordered according to effect size (lmm DEA)") +
	theme_test()


# pathways downloaded from http://www.gsea-msigdb.org/gsea/msigdb/index.jsp "MSigDB Collections", 28 Dec 2020
ranks <- dea_prot$betas
names(ranks) <- dea_prot$gname
pathways1 <- read.gmt("Data/GSEA/h.all.v7.2.symbols.gmt")
set.seed(42)
fgseaRes1 <- fgsea(pathways = pathways1, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes1[order(pval), ], 10)
enrich_oxphos_plt <- plotEnrichment(pathways1[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], ranks) + labs(title="HALLMARK_OXIDATIVE_PHOSPHORYLATION")

pathways2 <- read.gmt("Data/GSEA/c2.cp.kegg.v7.2.symbols.gmt")
set.seed(42)
fgseaRes2 <- fgsea(pathways = pathways2, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes2[order(pval), ], 20)

pathways3 <- read.gmt("Data/GSEA/c2.all.v7.2.symbols.gmt")
set.seed(42)
fgseaRes3 <- fgsea(pathways = pathways3, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes3[order(pval), ], 20)

setorder(fgseaRes1, -NES)
fgseaRes1$pathway <- factor(fgseaRes1$pathway, levels = fgseaRes1$pathway)

hallmarkpathways_plt <- 
ggplot(head(fgseaRes1[order(pval), ], 15), aes(x = NES, y = pathway, color = pval, size = size)) +
	geom_point() +
	geom_segment(aes(x = 0, xend = NES, y = pathway, yend = pathway), size = .5) +
	xlab("normalized enrichment score") +
	ggtitle("GSEA_prot_hallmrks") +
	theme_bw()


setorder(fgseaRes2, -NES)
fgseaRes2$pathway <- factor(fgseaRes2$pathway, levels = fgseaRes2$pathway)

keggpathways_plt <- 
ggplot(head(fgseaRes2[order(pval), ], 15), aes(x = NES, y = pathway, color = pval, size = size)) +
	geom_point() +
	geom_segment(aes(x = 0, xend = NES, y = pathway, yend = pathway), size = .5) +
	xlab("normalized enrichment score") +
	ggtitle("GSEA_prot_KEGG") +
	theme_bw()


gsea_prot_plt <- ggarrange(generanks_plt, enrich_oxphos_plt, hallmarkpathways_plt, keggpathways_plt, heights = c(1,2))

ggsave(paste(fig_path,"GSEA_prot.png", sep = ""), gsea_prot_plt, width = 16, height = 10)
dev.off()



# alternative polished plot (GESA on DEA results with pathway activity on all samples)

gsea_select <- head(fgseaRes1[order(pval), ], 8)
new_pathway_names <- c("Oxidative phosphorylation", "Hypoxia","E2F targets", "G2M checkpoint", "Inflammatory response", "Peroxisome", "Cholesterol homeostasis", "MYC targets")
data.table(gsea_select$pathway, new_pathway_names) # be careful here and doulbe check if manual names are correctly assigned
gsea_select[, pathway2 := new_pathway_names]
setorder(gsea_select, -NES)
gsea_select$pathway2 <- factor(gsea_select$pathway2, levels = gsea_select$pathway2)

gsea_PIact <- 
ggplot(gsea_select[order(pval), ], aes(x = NES, y = pathway2, color = pval, size = size)) +
	geom_point() +
	geom_segment(aes(x = 0, xend = NES, y = pathway2, yend = pathway2), size = .5) +
	xlab("Normalized enrichment score") +
	ylab("Gene set") +
	ggtitle("GSEA proteomics (PI activity)") +
	scale_color_viridis(name = "p value", option = "E") +
	guides(color = guide_colorbar(order = 1, barwidth = 15, barheight = 0.5), size = guide_legend(order = 0)) +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), legend.position = "bottom")

ggsave(paste(fig_path,"GSEA_prot_DEAwithPIact.png", sep = ""),gsea_PIact, width = 5.5, height = 3.7)
dev.off()





# alternative polished plot and analysis (GESA on DEA results with clinical score on all samples)

dea_protCSS <- fread("interimData/DEAlmm/diff_exp_data_prot_wcss_all.txt")


# dea_protCSS <- dea_protCSS[pval < 0.05, ]
setorder(dea_protCSS, -betas)
dea_protCSS[, order := 1:nrow(dea_protCSS)]

generanks_plt <- 
ggplot(dea_protCSS, aes(x = order, y = betas)) +
	geom_bar(stat = "identity") +
	xlab("gene rank") +
	ggtitle("Genes ordered according to effect size (lmm DEA)") +
	theme_test()


# pathways downloaded from http://www.gsea-msigdb.org/gsea/msigdb/index.jsp "MSigDB Collections", 28 Dec 2020
ranks <- dea_protCSS$betas
names(ranks) <- dea_protCSS$gname
pathways2 <- read.gmt("Data/GSEA/c2.cp.kegg.v7.2.symbols.gmt")
set.seed(42)
fgseaRes2CSS <- fgsea(pathways = pathways2, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes2CSS[order(pval), ], 20)


gsea_select2 <- head(fgseaRes2CSS[order(pval), ], 8)
new_pathway_names2 <- c("Glycolysis/Gluconeogeneiss", "Fructose/Mannose","TCA cycle","Ubiquitin mediated proteolysis", "Arginine/Proline", "Oxidative phosphorylation", "Galactose metabolism", "Parkinsons disease")
data.table(gsea_select2$pathway, new_pathway_names2) # be careful here and doulbe check if manual names are correctly assigned
gsea_select2[, pathway2 := new_pathway_names2]
setorder(gsea_select2, -NES)
gsea_select2$pathway2 <- factor(gsea_select2$pathway2, levels = gsea_select2$pathway2)


gsea_CSS <- 
ggplot(gsea_select2[order(pval), ], aes(x = NES, y = pathway2, color = pval, size = size)) +
	geom_point() +
	geom_segment(aes(x = 0, xend = NES, y = pathway2, yend = pathway2), size = .5) +
	xlab("Normalized enrichment score") +
	ylab("Gene set") +
	ggtitle("GSEA proteomics (clinical score)") +
	scale_color_viridis(name = "p value", option = "E") +
	guides(color = guide_colorbar(order = 1), size = guide_legend(order = 0)) +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), legend.position = "right")

ggsave(paste(fig_path,"GSEA_prot_DEAwithCSS.png", sep = ""), gsea_CSS, width = 5.5, height = 3.7)
dev.off()







##################################################
# transcript analysis

dea_rna <- fread("interimData/DEAlmm/diff_exp_data_rnaseq_pathwayact_all.txt")
dea_rna <- dea_rna[!duplicated(gname), ]



# dea_rna <- dea_rna[pval < 0.05, ]
setorder(dea_rna, -betas)
dea_rna[, order := 1:nrow(dea_rna)]

generanks_plt <- 
ggplot(dea_rna, aes(x = order, y = betas)) +
	geom_bar(stat = "identity") +
	xlab("gene rank") +
	ggtitle("Genes ordered according to effect size (lmm DEA)") +
	theme_test()


# pathways downloaded from http://www.gsea-msigdb.org/gsea/msigdb/index.jsp "MSigDB Collections", 28 Dec 2020
ranks <- dea_rna$betas
names(ranks) <- dea_rna$gname
pathways1 <- read.gmt("Data/GSEA/h.all.v7.2.symbols.gmt")
set.seed(42)
fgseaRes1RNA <- fgsea(pathways = pathways1, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes1RNA[order(pval), ], 30)
enrich_oxphos_plt <- plotEnrichment(pathways1[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], ranks) + labs(title="HALLMARK_OXIDATIVE_PHOSPHORYLATION")

pathways2 <- read.gmt("Data/GSEA/c2.cp.kegg.v7.2.symbols.gmt")
set.seed(42)
fgseaRes2 <- fgsea(pathways = pathways2, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes2[order(pval), ], 20)

pathways3 <- read.gmt("Data/GSEA/c5.go.mf.v7.2.symbols.gmt")
set.seed(42)
fgseaRes3 <- fgsea(pathways = pathways3, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes3[order(pval), ], 20)

setorder(fgseaRes1RNA, -NES)
fgseaRes1RNA$pathway <- factor(fgseaRes1RNA$pathway, levels = fgseaRes1RNA$pathway)

hallmarkpathways_plt <- 
ggplot(head(fgseaRes1RNA[order(pval), ], 15), aes(x = NES, y = pathway, color = pval, size = size)) +
	geom_point() +
	geom_segment(aes(x = 0, xend = NES, y = pathway, yend = pathway), size = .5) +
	xlab("normalized enrichment score") +
	ggtitle("GSEA_rna_hallmrks") +
	theme_bw()


setorder(fgseaRes2, -NES)
fgseaRes2$pathway <- factor(fgseaRes2$pathway, levels = fgseaRes2$pathway)

keggpathways_plt <- 
ggplot(head(fgseaRes2[order(pval), ], 15), aes(x = NES, y = pathway, color = pval, size = size)) +
	geom_point() +
	geom_segment(aes(x = 0, xend = NES, y = pathway, yend = pathway), size = .5) +
	xlab("normalized enrichment score") +
	ggtitle("GSEA_rna_KEGG") +
	theme_bw()


setorder(fgseaRes3, -NES)
fgseaRes3$pathway <- factor(fgseaRes3$pathway, levels = fgseaRes3$pathway)

ggplot(head(fgseaRes3[order(pval), ], 15), aes(x = NES, y = pathway, color = pval, size = size)) +
	geom_point() +
	geom_segment(aes(x = 0, xend = NES, y = pathway, yend = pathway), size = .5) +
	xlab("normalized enrichment score") +
	ggtitle("GSEA_rna_KEGG") +
	theme_bw()


gsea_rna_plt <- ggarrange(generanks_plt, enrich_oxphos_plt, hallmarkpathways_plt, keggpathways_plt, heights = c(1,2))

ggsave(paste(fig_path,"GSEA_rna.png", sep = ""), gsea_rna_plt, width = 16, height = 10)
dev.off()



# alternative polished plot (GESA on DEA results with pathway activity on all samples)

gsea_select3 <- head(fgseaRes1RNA[order(pval), ], 8)
new_pathway_names3 <- c("TNFalpha signalling via NFKB", "Epithelial mesenchymal transition", "Interferon gamma response", "E2F targets", "Hypoxia", "Myogenesis", "MYC targets", "Inflammatory response")
data.table(gsea_select3$pathway, new_pathway_names3) # be careful here and doulbe check if manual names are correctly assigned
gsea_select3[, pathway2 := new_pathway_names3]
setorder(gsea_select3, -NES)
gsea_select3$pathway2 <- factor(gsea_select3$pathway2, levels = gsea_select3$pathway2)

gsea_PIactRNA <- 
ggplot(gsea_select3[order(pval), ], aes(x = NES, y = pathway2, color = pval, size = size)) +
	geom_point() +
	geom_segment(aes(x = 0, xend = NES, y = pathway2, yend = pathway2), size = .5) +
	xlab("Normalized enrichment score") +
	ylab("Gene set") +
	ggtitle("GSEA transcriptomics (PI activity)") +
	guides(color = guide_colorbar(order = 1, barwidth = 15, barheight = 1), size = guide_legend(order = 0)) +
	scale_color_viridis(name = "p value", option = "E") +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), legend.position = "bottom")

ggsave(paste(fig_path,"GSEA_RNA_DEAwithPIact.png", sep = ""),gsea_PIactRNA, width = 5.5, height = 3.7)
dev.off()





# alternative polished plot and analysis (GESA on DEA results with clinical score on all samples)

dea_rnaCSS <- fread("interimData/DEAlmm/diff_exp_data_rnaseq_wcss_all.txt")
dea_rnaCSS <- dea_rnaCSS[!duplicated(gname), ]

# dea_rnaCSS <- dea_rnaCSS[pval < 0.05, ]
setorder(dea_rnaCSS, -betas)
dea_rnaCSS[, order := 1:nrow(dea_rnaCSS)]

generanks_plt <- 
ggplot(dea_rnaCSS, aes(x = order, y = betas)) +
	geom_bar(stat = "identity") +
	xlab("gene rank") +
	ggtitle("Genes ordered according to effect size (lmm DEA)") +
	theme_test()


# pathways downloaded from http://www.gsea-msigdb.org/gsea/msigdb/index.jsp "MSigDB Collections", 28 Dec 2020
ranks <- dea_rnaCSS$betas
names(ranks) <- dea_rnaCSS$gname
pathways2 <- read.gmt("Data/GSEA/c2.cp.kegg.v7.2.symbols.gmt")
set.seed(42)
fgseaRes2CSSRNA <- fgsea(pathways = pathways2, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes2CSSRNA[order(pval), ], 20)


gsea_select4 <- head(fgseaRes2CSSRNA[order(pval), ], 8)
new_pathway_names4 <- c("Fructose/Mannose", "Biosynthesis unsaturated FA","TCA cycle","Cytokine recepter interaction","Glycolysis/Gluconeogeneiss", "BCAA degradation", "Steroid biosynthesis","Starch/Sucrose")
data.table(gsea_select4$pathway, new_pathway_names4) # be careful here and doulbe check if manual names are correctly assigned
gsea_select4[, pathway2 := new_pathway_names4]
setorder(gsea_select4, -NES)
gsea_select4$pathway2 <- factor(gsea_select4$pathway2, levels = gsea_select4$pathway2)


gsea_CSS_RNA <- 
ggplot(gsea_select4[order(pval), ], aes(x = NES, y = pathway2, color = pval, size = size)) +
	geom_point() +
	geom_segment(aes(x = 0, xend = NES, y = pathway2, yend = pathway2), size = .5) +
	xlab("Normalized enrichment score") +
	ylab("Gene set") +
	ggtitle("GSEA transcriptomics (clinical score)") +
	scale_color_viridis(name = "p value", option = "E") +
	guides(color = guide_colorbar(order = 1), size = guide_legend(order = 0)) +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), legend.position = "right")

ggsave(paste(fig_path,"GSEA_RNA_DEAwithCSS.png", sep = ""), gsea_CSS_RNA, width = 5.5, height = 3.7)
dev.off()






gsea_combo_tbl <- rbind(gsea_select[, test := "Proteomics (PI activity)"], gsea_select2[, test := "Proteomics (clinical score)"], gsea_select3[, test := "Transcriptomics (PI activity)"], gsea_select4[, test := "Transcriptomics (clinical score)"])
colnames(gsea_combo_tbl)
setorder(gsea_combo_tbl, -NES)
gsea_combo_tbl$pathway2 <- factor(gsea_combo_tbl$pathway2, levels = unique(gsea_combo_tbl$pathway2), ordered=TRUE)


gsea_combinedProtRNA <- 
ggplot(gsea_combo_tbl, aes(x = NES, y = pathway2, color = pval, size = size)) +
	geom_point() +
	geom_segment(aes(x = 0, xend = NES, y = pathway2, yend = pathway2), size = .5) +
	facet_wrap(.~test, scales = "free") +
	xlab("Normalized enrichment score") +
	ylab("Gene set") +
	ggtitle("Gene set enrichment analysis") +
	scale_color_viridis(name = "p value", option = "E") +
	guides(color = guide_colorbar(order = 1), size = guide_legend(order = 0)) +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), legend.position = "right", strip.background = element_rect(fill="white"), axis.text.y = element_text(size = 10))


ggsave(paste(fig_path,"GSEA_RNAandProt_DEAcombined.png", sep = ""), gsea_combinedProtRNA, width = 10, height = 5)
dev.off()
ggsave(paste(fig_path_pdf,"Fig3/","GSEA_RNAandProt_DEAcombined.pdf", sep = ""), gsea_combinedProtRNA, device = "pdf", width = 10, height = 5)
dev.off()



gsea_Protcombo_tbl <- rbind(gsea_select[, test := "Proteomics (PI activity)"], gsea_select2[, test := "Proteomics (clinical score)"])
colnames(gsea_Protcombo_tbl)
setorder(gsea_Protcombo_tbl, -NES)
gsea_Protcombo_tbl$pathway2 <- factor(gsea_Protcombo_tbl$pathway2, levels = unique(gsea_Protcombo_tbl$pathway2), ordered=TRUE)


gsea_combinedProt <- 
ggplot(gsea_Protcombo_tbl, aes(x = NES, y = pathway2, color = pval, size = size)) +
	geom_point() +
	geom_segment(aes(x = 0, xend = NES, y = pathway2, yend = pathway2), size = .5) +
	facet_wrap(.~test, scales = "free", nrow = 2) +
	xlab("Normalized enrichment score") +
	ylab("Gene set") +
	ggtitle("Gene set enrichment analysis") +
	scale_color_viridis(name = "p value", option = "E") +
	guides(color = guide_colorbar(order = 1), size = guide_legend(order = 0)) +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), legend.position = "right", strip.background = element_rect(fill="white"), axis.text.y = element_text(size = 10))


ggsave(paste(fig_path,"GSEA_Prot_DEAcombined.png", sep = ""), gsea_combinedProt, width = 5.5, height = 5)
dev.off()













