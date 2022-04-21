# libraries
library(data.table)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggsci)
library(viridis)
library(ComplexHeatmap)

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

# create figures path
system("mkdir Figs/DEA/")
system("mkdir Figs/DEA/DESeq2/")
fig_path <- c("Figs/DEA/DESeq2/")
fig_path_pdf <- c("Figs/v8/pdf/")


mypal <- pal_aaas("default", alpha = 1)(9)
mypal2 <- pal_npg("nrc", alpha = 1)(9)


##############################################
##############################################
# apply DESeq2 package on RNA dataset
load("interimData/rnafiltered_mat.RData")
rna_exp0 <- data.table(ENSG = row.names(rnafiltered_mat), rnafiltered_mat)
rna_exp <- rna_exp0 %>% distinct(ENSG, .keep_all = TRUE)
hist(as.numeric(unlist(rna_exp[,3])), breaks = 200)
hist(as.numeric(unlist(rna_exp[4333,])))

# prep meta data tbl
mma_id <- colnames(rna_exp)[2:ncol(rna_exp)]
sample_no <- as.numeric(sub("MMA", "", mma_id))
type <- as.factor(ifelse(sample_no<151, "MMUT def.", ifelse(sample_no>210, "unaffected", "unknown")))
type1 <- as.factor(ifelse(sample_no<151, "MMUT def.", "control"))

meta_tbl <- data.frame(mma_id, sample_no, type, type1)
rownames(meta_tbl) <- colnames(rna_exp)[2:ncol(rna_exp)]

# make ENSG row names and drop ENSG column
rna_exp_raw <- as.data.frame(rna_exp)
rownames(rna_exp_raw) <- rna_exp$ENSG
rna_exp_raw <- subset(rna_exp_raw, select = -c(ENSG))
colnames(rna_exp_raw)

all(rownames(meta_tbl) == colnames(rna_exp_raw)) #TRUE


# create DESeq2 object
dds_rna <- DESeqDataSetFromMatrix(countData = round(rna_exp_raw), colData = meta_tbl, design = ~type1)
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
# ann_colors = list(type1 = c("MMUT def." = mypal[2], "control" = mypal[1]))
# plt <- pheatmap(vsd_cor_rna, annotation = select(meta_tbl, type1), annotation_color = ann_colors, show_colnames = FALSE, show_rownames = FALSE, border_color = "NA", color = inferno(100))

pdf(file = paste(fig_path_pdf, "SuppFig4/","RNA_DESeq2_heatmap.pdf", sep = ""), width = 6,height = 4.5)
Heatmap(vsd_cor_rna, col = cividis(100), name = "Correlation", column_title = "Count matrix (DESeq2 transcriptomics)", row_title = "Samples", show_row_names = FALSE, show_column_names = FALSE, top_annotation = HeatmapAnnotation(Group = meta_tbl$type1, col = list(Group = c("MMUT def." = mypal[2], "control" = mypal[1]))))
dev.off()


# plot PCA
plotPCA(vsd_rna, intgroup = "type")

install.packages("pheatmap")
##############################################
# continuation of analysis

# run DESeq analysis
dds_rna <- DESeq(dds_rna)

# check dispersion
# calculating mean for each gene (each row)
mean_counts <- apply(rna_exp_raw, 1, mean)
# calculating variance for each gene (each row)
variance_counts <- apply(rna_exp_raw, 1, var)
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
rna_res <- results(dds_rna, contrast = c("type1", "MMUT def.", "control"), alpha = 0.05)
png(paste(fig_path,"RNA_DESeq_contrasts.png", sep = ""))
plotMA(rna_res, ylim=c(-2,2))
dev.off()
# LFC shrinkage
rna_res <- lfcShrink(dds_rna, contrast = c("type1", "MMUT def.", "control"), res = rna_res)
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
rna_res1[ENSG == "ENSG00000146085", ]






##############################################
##############################################
# apply DESeq2 package on protein dataset
load("interimData/prot_mat.RData")
prot_exp0 <- data.table(ENSG = row.names(prot_mat), prot_mat)
prot_exp <- prot_exp0 %>% distinct(ENSG, .keep_all = TRUE)
hist(as.numeric(unlist(prot_exp[,3])), breaks = 200)
colnames(prot_exp)


# prep meta data tbl
mma_id <- colnames(prot_exp)[2:ncol(prot_exp)]
sample_no <- as.numeric(sub("MMA", "", mma_id))
type <- as.factor(ifelse(sample_no<151, "MMUT def.", ifelse(sample_no>210, "unaffected", "unknown")))
type1 <- as.factor(ifelse(sample_no<151, "MMUT def.", "control"))

meta_tbl <- data.frame(mma_id, sample_no, type, type1)
rownames(meta_tbl) <- colnames(prot_exp)[2:ncol(prot_exp)]

# make ENSG row names and drop ENSG column
prot_exp_raw <- as.data.frame(prot_exp)
rownames(prot_exp_raw) <- prot_exp$ENSG
prot_exp_raw <- subset(prot_exp_raw, select = -c(ENSG))
colnames(prot_exp_raw)

all(rownames(meta_tbl) == colnames(prot_exp_raw)) #TRUE


# create DESeq2 object
dds_prot <- DESeqDataSetFromMatrix(countData = round(prot_exp_raw), colData = meta_tbl, design = ~type1)
dds_prot <- estimateSizeFactors(dds_prot)
sizeFactors(dds_prot)
normalized_prot_counts <- counts(dds_prot, normalized=TRUE)

##############################################
# quality control

# unsupervised clustering analysis
vsd_prot <- vst(dds_prot, blind=TRUE)
# extract the vst matrix from the object
vsd_mat_prot <- assay(vsd_prot)
# compute pairwise correlation values
vsd_cor_prot <- cor(vsd_mat_prot)
# view heatmap
# ann_colors = list(type1 = c("MMUT def." = mypal[2], "control" = mypal[1]))
# plt <- pheatmap(vsd_cor_prot, main = "DESeq2 analsyis on proteomics", annotation = select(meta_tbl, type1), annotation_color = ann_colors, border_color = NA, show_colnames = FALSE, show_rownames = FALSE, scale = "none", color = cividis(1000, direction = -1))

pdf(file = paste(fig_path_pdf, "SuppFig4/","prot_DESeq2_heatmap.pdf", sep = ""), width = 6,height = 4.5)
Heatmap(vsd_cor_prot, col = cividis(100), name = "Correlation", column_title = "Count matrix (DESeq2 proteomics)", row_title = "Samples", show_row_names = FALSE, show_column_names = FALSE, top_annotation = HeatmapAnnotation(Group = meta_tbl$type1, col = list(Group = c("MMUT def." = mypal[2], "control" = mypal[1]))))
dev.off()

# plot PCA
plotPCA(vsd_prot, intgroup = "type1")


##############################################
# continuation of analysis

# run DESeq analysis
dds_prot <- DESeq(dds_prot)

# check dispersion
# calculating mean for each gene (each row)
mean_counts <- apply(prot_exp_raw, 1, mean)
# calculating variance for each gene (each row)
variance_counts <- apply(prot_exp_raw, 1, var)
dispersion_df <- data.frame(mean_counts, variance_counts)
dispersion_plt <-
ggplot(dispersion_df, aes(x = mean_counts, y = variance_counts)) +
	geom_point() +
	scale_y_log10() +
	scale_x_log10() +
	xlab("Mean counts per gene") +
	ylab("Variance per gene")
ggsave(paste(fig_path,"prot_dispersion_plot.png", sep = ""), dispersion_plt, width=4,height=4)

# plot dispersion estimates
png(paste(fig_path,"prot_dispersion_estimates.png", sep = ""))
plotDispEsts(dds_prot)
dev.off()

# DESeq contrasts
prot_res <- results(dds_prot, contrast = c("type1", "MMUT def.", "control"), alpha = 0.05)
png(paste(fig_path,"prot_DESeq_contrasts.png", sep = ""))
plotMA(prot_res, ylim=c(-2,2))
dev.off()
# LFC shrinkage
prot_res <- lfcShrink(dds_prot, contrast = c("type1", "MMUT def.", "control"), res = prot_res)
plotMA(prot_res, ylim=c(-2,2))

mcols(prot_res)
summary(prot_res)
prot_res1 <- data.frame(prot_res)
prot_res1 <- data.table(prot_res1)
colnames(prot_res1)
nms_prot <- rownames(prot_res)
prot_res1$ENSG <- nms_prot
prot_res1 <- prot_res1 %>% arrange(padj)
prot_res1 <- prot_res1 %>% mutate(threshold = padj < 0.05)
head(prot_res1)
prot_res1[ENSG == "ENSG00000146085", ]

ggplot(prot_res1, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
	geom_point() +
	theme_bw()
















