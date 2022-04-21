# libraries
library(OUTRIDER)
library(annotables)
library(data.table)
library(ggplot2)
library(ggpubr)
library("Hmisc")

library(GenomicFeatures)
library(AnnotationDbi)
library(dplyr)
library(patchwork)
library(ggsci)
library(htmltools)


# create figures path
system("mkdir Figs/Outrider/")
fig_path <- c("Figs/Outrider/")
fig_path_pdf <- c("Figs/v8/pdf/")

mypal <- pal_aaas("default", alpha = 1)(9)
mypal2 <- pal_npg("nrc", alpha = 1)(9)


##############################################
##############################################
# Load RNA dataset
load("interimData/rnafiltered_mat.RData")
rna_exp0 <- data.frame(ENSG = row.names(rnafiltered_mat), rnafiltered_mat)
counts0 <- rna_exp0 %>% distinct(ENSG, .keep_all = TRUE)
colnames(counts0)

# Outrider object:
rownames(counts0) <- counts0$ENSG
counts0 <- counts0[,-1]
counts <- counts0*10e3


# tmp metadata:
cond1= rep (c("MMA"), times=202)
cond2= rep (c("control"), times=19)
class1= c(cond1, cond2)
subtype1= rep(c("deficient"), times=143)
subtype2= rep(c("non_deficient"), times= 59)
subtype3= rep(c("control"), times=19)
class2= c(subtype1, subtype2, subtype3)

meta= data.frame(colnames(counts), class1, class2)
colnames(meta)= c("sampleID", "phenotype_1", "phenotype_2")

ods <- OutriderDataSet(countData= round(counts), colData = meta) 
ods <- estimateSizeFactors(ods)
# plot the sizeFactor across the dataset
ggplot(data=NULL, aes(y=sort(sizeFactors(ods)), x=1:ncol(ods))) + 
    geom_point() + 
    labs(x='Sample rank', y='Size factors', title="Size factor distribution") + 
    theme_bw()


## co-variation
options(repr.plot.width=6, repr.plot.height=5)
plotCountCorHeatmap(ods, colGroups=c("phenotype_1", "phenotype_2"), normalize=FALSE)


## gene/sample expression
ods <- plotCountGeneSampleHeatmap(ods, colGroups=c("phenotype_1", "phenotype_2"), normalized=FALSE, nRowCluster=3)


#Denoising autoencoder fitted to our data
plotCountCorHeatmap(ods_noconfounders, colGroups=c("phenotype_1", "phenotype_2"), normalize=TRUE, nRowCluster=3)





## Fitting the model
##Find q (takes a long time to compute (6 days):
ods_filtered_findparams <- findEncodingDim(ods,BPPARAM=MulticoreParam(1, progressbar=TRUE))
# visualize the hyper parameter optimization
# plotEncDimSearch(ods_filtered_findparams)


ods_model <- OUTRIDER(ods, q=38, verbose=TRUE, iterations=15)


## Exploration of expression outliers

tail(sort(aberrant(ods_model, by="sample")), n=10)
tail(sort(aberrant(ods_model, by="gene", zScoreCutoff=1)))


plotAberrantPerSample(ods_model, padjCutoff=0.05)



# add gene symbols
orig_rownames <- rownames(ods_model)
HGNC_list_r <- data.table(read.table("Data/gene_sets/ENSG_HGNC_list_r.csv", header = TRUE))
setnames(HGNC_list_r, c("ENSG", "HGNC"))
name_transfer_tbl <- data.table(orig_ENSG = orig_rownames)

HGNC_list_r[ENSG == name_transfer_tbl[7], ]$HGNC
HGNC_list_r[ENSG == "ISG153", ]

test_x <- ifelse(sum(HGNC_list_r$ENSG == name_transfer_tbl$orig_ENSG[7])>0,
	HGNC_list_r[ENSG == name_transfer_tbl$orig_ENSG[7], ]$HGNC,
	name_transfer_tbl$orig_ENSG[7])
name_transfer_tbl[orig_ENSG == name_transfer_tbl$orig_ENSG[7], new_names := test_x]

for (i in 1:length(name_transfer_tbl$orig_ENSG)) {
	print(i)
	test_x <- ifelse(sum(HGNC_list_r$ENSG == name_transfer_tbl$orig_ENSG[i])>0,
	HGNC_list_r[ENSG == name_transfer_tbl$orig_ENSG[i], ]$HGNC,
	name_transfer_tbl$orig_ENSG[i])
	name_transfer_tbl[orig_ENSG == name_transfer_tbl$orig_ENSG[i], new_names := test_x]
}

ods_model_gnames <- copy(ods_model)
rownames(ods_model_gnames) <- name_transfer_tbl$new_names

# saveRDS(ods_model_gnames, file = "interimData/Outrider/ods_model_gnames_q38.rds")

ods_model_gnames <- readRDS(file = "interimData/Outrider/ods_model_gnames_q44.rds")


res <- results(ods_model_gnames, all = TRUE)
# res <- results(ods_model_gnames, all = FALSE) # only significant results
#write.table(res, "results_ods_model_addmetadata_padj0.05.tsv", sep="\t", quote = FALSE, row.names = FALSE)

res_z <- results(ods_model_gnames, padjCutoff=0.1, zScoreCutoff=2)
#write.table(res_z, "results_ods_model_addmetadata_Zscore2.tsv", sep="\t", quote = FALSE, row.names = FALSE)



### Plot results
#Examples for specific sample or gene exploration: Do MUT and clinically relevant genes first!
res[geneID == "MMAB" & sampleID == "MMA152"]
res[geneID == "SUCLA2"]
res[which(res$sampleID == "MMA186"),1:5]


#############
#ACSF3:

acsf3_ranks <- plotExpressionRank(ods_model_gnames, "ACSF3", norm=TRUE, basePlot=TRUE) + theme_pubr() + theme(plot.title = element_text(size=12), legend.position = "none")

acsf3_ranks_html <- plotExpressionRank(ods_model_gnames, "ACSF3", norm=TRUE, basePlot=FALSE)
save_html(acsf3_ranks_html, file = paste0(fig_path, "ACSF3/acsf3_ranks_html.html"))

# volcano plot of sample MMA186 and MMA187:
MMA186=plotVolcano(ods_model_gnames, "MMA186", base=TRUE) + labs(title = "MMA186") + theme_pubr() + theme(plot.title = element_text(size=12))
MMA187=plotVolcano(ods_model_gnames, "MMA187", base=TRUE) + labs(title = "MMA187") + theme_pubr() + theme(plot.title = element_text(size=12))

save_html(plotVolcano(ods_model_gnames, "MMA186", base=FALSE), file = paste0(fig_path, "ACSF3/MMA186.html"))
save_html(plotVolcano(ods_model_gnames, "MMA187", base=FALSE), file = paste0(fig_path, "ACSF3/MMA187.html"))

acsf3_combo_plt <- 
acsf3_ranks + MMA186 + MMA187 + plot_layout(widths = c(2,1,1))

ggsave(filename=paste0(fig_path,"ACSF3_combo.pdf"), acsf3_combo_plt, device = "pdf", width = 7, height = 3)




#############
#SUCLA2:

sucla2_ranks <- plotExpressionRank(ods_model_gnames, "SUCLA2", norm=TRUE, basePlot=TRUE) + theme_pubr() + theme(plot.title = element_text(size=12), legend.position = "none")

sucla2_ranks_html <- plotExpressionRank(ods_model_gnames, "SUCLA2", norm=TRUE, basePlot=FALSE)
save_html(sucla2_ranks_html, file = paste0(fig_path, "SUCLA2/sucla2_ranks_html.html"))

# volcano plot of sample MMA164 and MMA192:
MMA164=plotVolcano(ods_model_gnames, "MMA164", base=TRUE) + labs(title = "MMA164") + theme_pubr() + theme(plot.title = element_text(size=12))
MMA192=plotVolcano(ods_model_gnames, "MMA192", base=TRUE) + labs(title = "MMA192") + theme_pubr() + theme(plot.title = element_text(size=12))

save_html(plotVolcano(ods_model_gnames, "MMA164", base=FALSE), file = paste0(fig_path, "SUCLA2/MMA164.html"))
save_html(plotVolcano(ods_model_gnames, "MMA192", base=FALSE), file = paste0(fig_path, "SUCLA2/MMA192.html"))

sucla2_combo_plt <- 
sucla2_ranks + MMA164 + MMA192 + plot_layout(widths = c(2,1,1))

ggsave(filename=paste0(fig_path,"SUCLA2_combo.pdf"), sucla2_combo_plt, device = "pdf", width = 7, height = 3)




#############
#MMAA:

mmaa_ranks <- plotExpressionRank(ods_model_gnames, "MMAA", norm=TRUE, basePlot=TRUE) + theme_pubr() + theme(plot.title = element_text(size=12), legend.position = "none")

mmaa_ranks_html <- plotExpressionRank(ods_model_gnames, "MMAA", norm=TRUE, basePlot=FALSE)
save_html(mmaa_ranks_html, file = paste0(fig_path, "MMAA/mmaa_ranks_html.html"))

# volcano plot of sample MMA206 and MMA209:
MMA206=plotVolcano(ods_model_gnames, "MMA206", base=TRUE) + labs(title = "MMA206") + theme_pubr() + theme(plot.title = element_text(size=12))
MMA209=plotVolcano(ods_model_gnames, "MMA209", base=TRUE) + labs(title = "MMA209") + theme_pubr() + theme(plot.title = element_text(size=12))

save_html(plotVolcano(ods_model_gnames, "MMA206", base=FALSE), file = paste0(fig_path, "MMAA/MMA206.html"))
save_html(plotVolcano(ods_model_gnames, "MMA209", base=FALSE), file = paste0(fig_path, "MMAA/MMA209.html"))

mmaa_combo_plt <- 
mmaa_ranks + MMA206 + MMA209 + plot_layout(widths = c(2,1,1))

ggsave(filename=paste0(fig_path,"MMAA_combo.pdf"), mmaa_combo_plt, device = "pdf", width = 7, height = 3)





#############
#MMAB:

mmab_ranks <- plotExpressionRank(ods_model_gnames, "MMAB", norm=TRUE, basePlot=TRUE) + theme_pubr() + theme(plot.title = element_text(size=12), legend.position = "none")

mmab_ranks_html <- plotExpressionRank(ods_model_gnames, "MMAB", norm=TRUE, basePlot=FALSE)
save_html(mmab_ranks_html, file = paste0(fig_path, "MMAB/mmab_ranks_html.html"))

# volcano plot of sample MMA152 and MMA210:
MMA152=plotVolcano(ods_model_gnames, "MMA152", base=TRUE) + labs(title = "MMA152") + theme_pubr() + theme(plot.title = element_text(size=12))

save_html(plotVolcano(ods_model_gnames, "MMA152", base=FALSE), file = paste0(fig_path, "MMAB/MMA152.html"))

mmab_combo_plt <- 
mmab_ranks + MMA152 + plot_spacer() + plot_layout(widths = c(2,1, 1))

ggsave(filename=paste0(fig_path,"MMAB_combo.pdf"), mmab_combo_plt, device = "pdf", width = 5, height = 3)



# arrange all outrider plots

outrider_all <- 
ggarrange(acsf3_combo_plt ,sucla2_combo_plt ,mmaa_combo_plt, mmab_combo_plt, labels = "auto")

ggsave(filename=paste0(fig_path,"outrider_complete_plt.pdf"), outrider_all, device = "pdf", width = 14, height = 6)




## other genes - put in any example gene
plotExpressionRank(ods_model_gnames, "NDN", norm=TRUE, basePlot=TRUE)

##QQ plot:
plotQQ(ods_model_gnames, "SUCLA2")


## Observed versus expected gene expression
plotExpectedVsObservedCounts(ods_model_gnames, "SUCLA2", basePlot=TRUE)



#NDN levels in SUCLA2 mutated samples
res[geneID == "SUCLA2"][,2]
group2=colnames(ods_model_gnames)[which(colnames(ods_model_gnames) %nin% c("MMA192", "MMA167", "MMA183", "MMA164"))]
group1=c("MMA192", "MMA167", "MMA183", "MMA164")
plotExpressionRank(ods_model_gnames, "NDN", norm=TRUE, groups = c(group1, group2),col = c("blue", "red"), groupColSet = "Paired", basePlot=TRUE)


