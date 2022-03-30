# libraries
require(data.table)
require(ggplot2)
require(tidyverse)
require(ggpubr)
library(factoextra)
library(ggsci)
require(readxl)
library(DESeq2)
library(fgsea)
library(biomaRt)
library(pheatmap)
library(qusage)
library(viridis)






#colour_palette
##blue 
#82
#129
#163
myblue=rgb(60, 107, 131, maxColorValue = 255)
#orange
#231
#129
#48
myorange=rgb(231, 129, 48, maxColorValue = 255)
#brown
#180
#100
#36
mybrown=rgb(180, 100, 36, maxColorValue = 255)
#mygrey
#212
#212
#212
mygrey=rgb(180, 180, 180, maxColorValue = 255)


mypal1 <- pal_aaas("default", alpha = 1)(9)
mypal2 <- pal_npg("nrc", alpha = 1)(9)


### load some helper functions
source("csrproject/Code/normalizationFunctions.R")


# create figures path
system("mkdir Figs/mouse")
system("mkdir Figs/mouse/brainRNAseq")
fig_path <- c("Figs/mouse/brainRNAseq/")
system("mkdir interimData/mouse")
fig_path_pdf <- c("Figs/v8/pdf/")


######################################################################
## create human mouse gene name conversion table
######################################################################

# this table was created on 26 Aug 2021, 2 PM


# Mouse2Human <- function(MouseGenes){

#   human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#   mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#   genesMousetoHuman = getLDS(attributes = c("ensembl_gene_id","mgi_symbol"), 
#                              filters = "mgi_symbol", 
#                              values = MouseGenes , 
#                              mart = mouse, 
#                              attributesL = c("ensembl_gene_id", "hgnc_symbol"), 
#                              martL = human, 
#                              uniqueRows = TRUE)

#   colnames(genesMousetoHuman) <- c("Mouse.Gene_ID", "MGI", "Human.Gene_ID", "HGNC")

#   return(genesMousetoHuman) 

# }

# ## Get mouse genes
# mmusculus_genes <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"), useCache = FALSE)

# ## create the conversion table
# Mouse2HumanTable <- Mouse2Human(MouseGenes = mmusculus_genes$mgi_symbol)
# fwrite(Mouse2HumanTable, file = "Data/mouse/Mouse2HumanTable.csv")


######################################################################
## import data and clean
######################################################################


tbl <- data.table(read_excel("Data/mouse/result--ko-ki--over--ki-wt.xlsx"))
colnames(tbl)

tbl1 <- tbl[, c(2, 24:31)]
summary(duplicated(tbl1$gene_name))
tbl1[duplicated(gene_name), ]$gene_name
tbl2 <- tbl1[!duplicated(gene_name), ]
tbl2[tbl2$gene_name == "Mut", ]$gene_name <- "Mmut"
tbl3 <- melt.data.table(tbl2, id = "gene_name")
tbl4 <- dcast(tbl3, variable ~ gene_name, value.var = "value")
lbl <- c(rep("kiwt", 4), rep("koki", 4))
row.names(tbl4) <- paste(tbl4$variable, row.names(tbl4), lbl, sep = "_")

tbl4$variable <- NULL

tbl5 <- scalerRound(tbl4, 7)
tbl5 <- data.table(tbl5)

hist(tbl5[[1244]], breaks = 30)


#########
# run PCA

prot_exp_pca <- prcomp(tbl5)
plot(prot_exp_pca$x[,1], prot_exp_pca$x[,2])


df_out_prot <- as.data.frame(prot_exp_pca$x)
df_out_prot$no <- as.numeric(row.names(tbl5))

df_out_prot$type <- c(rep("kiwt", 4), rep("koki", 4))


#calculate percentages
perc_prot_pca <- vector("character", nrow(prot_exp_pca$x))

for(i in 1:nrow(prot_exp_pca$x)) {
	perc_prot_pca[i] <- round(((prot_exp_pca$sdev[i])^2) / sum(prot_exp_pca$sdev^2) * 100, 2)
	perc_prot_pca[i] <- paste(colnames(df_out_prot)[i], " (", paste(as.character(perc_prot_pca[i]), "%", ")", sep = ""), sep = "")
}



# scree plot with package factoextra
scree_plot <- fviz_eig(prot_exp_pca)

g_prot_pca2 <-
ggplot(df_out_prot, aes(x = PC1, y = PC2, color = type)) +
	geom_point() +
	ggtitle("RNAseq mouse brain") +
	xlab(perc_prot_pca[1]) +
	ylab(perc_prot_pca[2]) +
	scale_color_aaas() +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), legend.title = element_blank(), legend.position = "right")


mouse_pca <- 
ggarrange(scree_plot, g_prot_pca2, ncol = 1)

ggsave(paste(fig_path,"RNAmousePCA.png", sep = ""), mouse_pca, width = 4, height = 5)
dev.off()



################################################
# make simple comparisons

# compare Mmut
names <- colnames(tbl[, c(24:31)])
names_kiwt <- names[1:4]
names_koki <- names[5:8]
new_tbl <- copy(tbl3)

new_tbl[which(new_tbl$variable %in% names_kiwt), type := "kiwt"]
new_tbl[which(new_tbl$variable %in% names_koki), type := "koki"]


# hand-picked genes of interest
new_tbl[grep("pdk", new_tbl$gene_name, ignore.case = TRUE), ]
interest_genes <- c("Mmut", "Ogdh", "Glud1", "Pdk4")

ggplot(new_tbl[which(new_tbl$gene_name %in% interest_genes), ], aes(x = type, y = value)) +
	geom_boxplot() +
	facet_wrap(~ gene_name, scales = "free") +
	stat_compare_means(comparisons = list(c("kiwt", "koki")), method = "t.test", size = 3) +
	scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
	theme_pubr()

# test Ogdh without outlier
ggplot(new_tbl[gene_name == interest_genes[2] & value > 121, ], aes(x = type, y = value)) +
	geom_boxplot() +
	stat_compare_means(comparisons = list(c("kiwt", "koki")), method = "t.test", size = 3) +
	ggtitle(interest_genes[2]) +
	scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
	theme_pubr()

ggsave(paste(fig_path,"RNAmouse_Ogdh.png", sep = ""), width = 2, height = 3)
dev.off()


# TCA gene list (Florian)
flos_TCA <- fread("Data//gene_sets/List ENSG qPCR gene Florian.csv", header = FALSE)
setnames(flos_TCA, c("ENSG", "HGNC"))
flos_genes <- flos_TCA$HGNC


ggplot(filter(new_tbl, grepl(paste(paste("^", flos_genes, "$", sep = ""), collapse = "|"), gene_name, ignore.case = TRUE)), aes(x = type, y = value)) +
	geom_boxplot() +
	facet_wrap(~ gene_name, scales = "free") +
	stat_compare_means(comparisons = list(c("kiwt", "koki")), method = "t.test", size = 3) +
	scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
	theme_test()

ggsave(paste(fig_path,"RNAmouse_TCAgenes.png", sep = ""), width = 9, height = 9)
dev.off()




# GSEA below shows oxphos and hypoxia; hence check HIF 1 and 2.

unique(new_tbl[grep("hif", new_tbl$gene_name, ignore.case = TRUE), ]$gene_name)
interest_genes2 <- c("Hif1a", "Hif3a", "Hif1an")

ggplot(new_tbl[which(new_tbl$gene_name %in% interest_genes2), ], aes(x = type, y = value)) +
	geom_boxplot() +
	facet_wrap(~ gene_name, scales = "free") +
	stat_compare_means(comparisons = list(c("kiwt", "koki")), method = "t.test", size = 3) +
	scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
	theme_pubr()

ggsave(paste(fig_path,"RNAmouse_HIFgenes.png", sep = ""), width = 5, height = 3)
dev.off()





##############################################
##############################################
# apply DESeq2 package on RNA dataset

colnames(tbl)
deseq_tbl0 <- tbl %>% distinct(Identifier, .keep_all = TRUE)
cols <- colnames(tbl[, c(1, 24:31)])
deseq_tbl <- deseq_tbl0[, ..cols, with = FALSE]
setnames(deseq_tbl, c("Identifier", paste(colnames(deseq_tbl)[2:5], "_kiwt", sep = ""), paste(colnames(deseq_tbl)[6:9], "_koki", sep = "")))

# prep meta data tbl
mma_id <- colnames(deseq_tbl)[2:ncol(deseq_tbl)]
type <- c(rep("Mmut-ki/wt", 4), rep("Mmut-ko/ki", 4))
meta_tbl <- data.frame(type)
rownames(meta_tbl) <- colnames(deseq_tbl)[2:ncol(deseq_tbl)]

# make identifier row names and drop Identifier column
deseq_tbl_raw <- as.data.frame(deseq_tbl)
rownames(deseq_tbl_raw) <- deseq_tbl$Identifier
deseq_tbl_raw <- subset(deseq_tbl_raw, select = -c(Identifier))
colnames(deseq_tbl_raw)
all(rownames(meta_tbl) == colnames(deseq_tbl_raw)) #TRUE

# create DESeq2 object
dds_rna <- DESeqDataSetFromMatrix(countData = round(deseq_tbl_raw*1000), colData = meta_tbl, design = ~type)
dds_rna <- estimateSizeFactors(dds_rna)
sizeFactors(dds_rna)
normalized_rna_counts <- counts(dds_rna, normalized=TRUE)

##############################################
# quality control

# unsupervised clustering analysis
vsd_rna <- vst(dds_rna, blind=TRUE)
# extract the vst matrix from the object
vsd_mat_rna <- assay(vsd_rna)
# compute pairwise correlation values
vsd_cor_rna <- cor(vsd_mat_rna)
# view heatmap
plt <- pheatmap(vsd_cor_rna, col = viridis(4), annotation = meta_tbl, annotation_colors = list(type = c("Mmut-ki/wt" = mypal1[1], "Mmut-ko/ki" = mypal1[2])), show_rownames = FALSE, show_colnames = FALSE)
ggsave(paste0(fig_path_pdf,"SuppFig_mouseBrainTranscriptomics/RNA_DESeq2_heatmap.pdf"), plt[[4]], device = "pdf", width=5,height=3.5)
dev.off()

# plot PCA
deseq_pca_plt <- plotPCA(vsd_rna, intgroup = "type", ntop = 60)

deseq_pca_plt <- deseq_pca_plt + 
	scale_color_aaas() +
	theme_pubr() + 
	theme(legend.title = element_blank())

ggsave(paste0(fig_path_pdf,"SuppFig_mouseBrainTranscriptomics/RNA_DESeq2_PCA.pdf"), deseq_pca_plt, device = "pdf", width=3,height=3)
dev.off()

##############################################
# continuation of analysis

# run DESeq analysis
dds_rna <- DESeq(dds_rna)
res <- results(dds_rna, tidy = TRUE)
fwrite(res, file = "interimData/mouse/deseq-results-tidy-mouseBrain.csv")

# check dispersion
# calculating mean for each gene (each row)
mean_counts <- apply(deseq_tbl_raw, 1, mean)
# calculating variance for each gene (each row)
variance_counts <- apply(deseq_tbl_raw, 1, var)
dispersion_df <- data.frame(mean_counts, variance_counts)

dispersion_plt <- 
ggplot(dispersion_df, aes(x = mean_counts, y = variance_counts)) +
	geom_point() +
	scale_y_log10() +
	scale_x_log10() +
	xlab("Mean counts per gene") +
	ylab("Variance per gene")
ggsave(paste(fig_path,"RNA_dispersion_plot.png", sep = ""), dispersion_plt, width=4,height=4)

# plot dispersion estimates
png(paste(fig_path,"RNA_dispersion_estimates.png", sep = ""))
plotDispEsts(dds_rna)
dev.off()

# DESeq contrasts
rna_res <- results(dds_rna, contrast = c("type", "kiwt", "koki"), alpha = 0.05)
png(paste(fig_path,"RNA_DESeq_contrasts.png", sep = ""))
plotMA(rna_res, ylim=c(-2,2))
dev.off()
# LFC shrinkage
rna_res <- lfcShrink(dds_rna, contrast = c("type", "kiwt", "koki"), res = rna_res)
plotMA(rna_res, ylim=c(-2,2))

mcols(rna_res)
summary(rna_res)
rna_res1 <- data.frame(rna_res)
rna_res1 <- data.table(rna_res1)
nms <- rownames(rna_res)
rna_res1$ENSG <- nms
colnames(rna_res1)
rna_res1 <- rna_res1 %>% arrange(padj)
rna_res1 <- rna_res1 %>% mutate(threshold = padj < 0.05)
head(rna_res1)
rna_res1[ENSG == "Mut", ]




##############################################
# use DESeq2 results in GSEA
# DESeq2 does not need to be rerun as the results are saved in the interimData folder

res <- fread("interimData/mouse/deseq-results-tidy-mouseBrain.csv")

# convert mouse ENSM number to human gene symbol using the mouse2human gene table
m2h_tbl <- fread("Data/mouse/Mouse2HumanTable.csv")

human_smbls <- m2h_tbl[match(res$row, m2h_tbl$Mouse.Gene_ID), ]$HGNC
res$human_smbls <- human_smbls


# make unique and exclude NAs
res3 <- res %>% 
  dplyr::select(human_smbls, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(human_smbls) %>% 
  summarize(stat=mean(stat))
res3
summary(res3$human_smbls == "")
res3 <- res3[res3$human_smbls != "", ]
# -log(res3$stat) # when using p-value

# apply fgsea package
ranks <- deframe(res3)
head(ranks, 20)

# load pathways
# fls <- list.files("Data/GSEA/")

# for (i in 1:length(fls)) {
# 	print(paste0(i,"/",length(fls)))
# 	pathways1 <- read.gmt(paste0("Data/GSEA/", fls[i]))

# 	set.seed(42)
# 	fgseaRes1 <- fgsea(pathways = pathways1, 
# 	                  stats    = ranks,
# 	                  eps      = 0.0,
# 	                  minSize  = 15,
# 	                  maxSize  = 500)
# 	# head(fgseaRes1[order(pval), ], 10)
# 	# enrich_oxphos_plt <- plotEnrichment(pathways1[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], ranks) + labs(title="HALLMARK_OXIDATIVE_PHOSPHORYLATION")

# 	setorder(fgseaRes1, -NES)
# 	fgseaRes1$pathway <- factor(fgseaRes1$pathway, levels = fgseaRes1$pathway)

# 	pp <- ggplot(head(fgseaRes1[order(pval), ], 12), aes(x = NES, y = pathway, color = pval, size = size)) +
# 		geom_point() +
# 		geom_segment(aes(x = 0, xend = NES, y = pathway, yend = pathway), size = .5) +
# 		xlab("normalized enrichment score") +
# 		ylab("Gene set") +
# 		scale_color_viridis(name = "p value", option = "E") +
# 		guides(color = guide_colorbar(order = 1, barwidth = 1, barheight = 6), size = guide_legend(order = 0)) +
# 		ggtitle(gsub(".gmt", "", fls[i])) +
# 		theme_pubr() +
# 		theme(plot.title = element_text(size = 12), legend.position = "right",axis.text.y = element_text(size = 10))

# 	ggsave(paste0(fig_path,"GSEA_pval_",gsub(".gmt", "", fls[i]),".pdf"), pp, device = "pdf", width = 12, height = 3.5)
# }


# alternative polished plot (important to set.seed above, otherwise everytime results will slightly change and order of pathways might become incorrect)


pathways1 <- read.gmt("Data/GSEA/c2.cp.wikipathways.v7.4.symbols.gmt")

set.seed(42)
fgseaRes1 <- fgsea(pathways = pathways1, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)

gsea_select <- head(fgseaRes1[order(pval), ], 8)
gsea_select
new_pathway_names <- c("Cytoplasmic ribosomal proteins", "Non-alcoholic fatty liver disease", "ETC and oxidative phosphorylation", "Mitochondrial complex IV assembly", "Oxidative phosphorylation", "GABA receptor signaling", "Mitochondrial complex I assembly", "Hippomerlin signaling dysregulation")
data.table(gsea_select$pathway, new_pathway_names) # be careful here and doulbe check if manual names are correctly assigned
gsea_select[, pathway2 := new_pathway_names]
setorder(gsea_select, -NES)
gsea_select$pathway2 <- factor(gsea_select$pathway2, levels = gsea_select$pathway2)

gsea_mouserna <- ggplot(gsea_select[order(pval), ], aes(x = NES, y = pathway2, color = pval, size = size)) +
	geom_point() +
	geom_segment(aes(x = 0, xend = NES, y = pathway2, yend = pathway2), size = .5) +
	xlab("Normalized enrichment score") +
	ylab("Gene set") +
	scale_color_viridis(name = "p-value", option = "E") +
	guides(color = guide_colorbar(order = 1, barwidth = 1, barheight = 6), size = guide_legend(order = 0)) +
	ggtitle("Pathway enrichment analysis") +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), legend.position = "right",axis.text.y = element_text(size = 10))

ggsave(paste0(fig_path_pdf,"GSEA_stat_c2.cp.wikipathways.v7.4.pdf"), gsea_mouserna, device = "pdf", width = 6, height = 3.5)




# combine Mmut plot with GSEA plot

mmut_tbl <- copy(new_tbl)
mmut_tbl <- mmut_tbl[gene_name == "Mmut", ]
mean_ctrl <- mean(mmut_tbl[type == "kiwt", ]$value)
mmut_tbl[, norm_val := value/mean_ctrl]
mmut_tbl[, type2 := "Mmut-ki/wt"]
mmut_tbl[type == "koki", type2 := "Mmut-ko/ki"]


compare <- list(c("Mmut-ki/wt", "Mmut-ko/ki"))

mmut_exp_plt <- 
ggplot(mmut_tbl, aes(x = type2, y = norm_val, color = type2)) +
	geom_boxplot() +
	geom_jitter(shape =21, fill = "white", color = "black", width = 0.16) +
	stat_compare_means(comparisons = compare, method = "t.test", size = 3) +
	labs(title = "", x="", y="Relative Mmut expression")+
	scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
	scale_color_aaas() +
	guides(color=guide_legend(nrow=2,byrow=TRUE)) +
	theme_pubr() +
	theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), legend.position = "bottom", plot.title = element_text(size=12))



library(patchwork)

brain_mouse_trans <- 
mmut_exp_plt + gsea_mouserna + plot_layout(widths = c(1,2)) + plot_annotation(title = "Mouse brain transcriptomics", theme = theme(plot.title = element_text(size = 12)))

ggsave(paste0(fig_path_pdf,"Brain_mouse_transcriptomics_combo.pdf"), brain_mouse_trans, device = "pdf", width = 8, height = 4)



