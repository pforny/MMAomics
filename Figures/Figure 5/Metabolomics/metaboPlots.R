# metabolmics data set analysis
library(data.table)
# library(readr)
# library(readxl)
# library(reshape2)
library(ggpubr)
library(extrafont)
library(gridExtra)
library(gtools)
# library(EnhancedVolcano)
# library(dplyr)
# library(plyr)
library(tidyverse)
library(ggsci)
library(Ipaper)
library(readxl)
library(ggrepel)
library(factoextra)
library(patchwork)


system("mkdir Figs/metabo")

fig_path <- c("Figs/metabo/")
fig_path_pdf <- c("Figs/v12/pdf/")
mypal1 <- pal_aaas("default", alpha = 1)(9)
mypal2 <- pal_npg("nrc", alpha = 1)(9)
  
### load some helper functions
source("csrproject/Code/normalizationFunctions.R")





id_tbl=fread("interimData/metabo_annot/metaboId.tbl")



# import file with calculations
tbl_calc <- data.table(read_excel("Data/metaboDat/intra/froese1 DIFF DATA.xlsx", sheet = 5))
setnames(tbl_calc, c("difference", "ionIdx", "index diff", "log2_fc","p_val","p_val_adj","pFDR","q_value","m_z","hmdb_id","kegg_id","other_id", "name","mz difference","formula","ion","annotation score"))


metabo_int <- fread("Data/metaboDat/intra/intracell_intensities1.csv")
metabo_int <- metabo_int %>%
  select(ionIdx, ionMz, contains("E1"), contains("E2"))

metabo_ids <- fread("Data/metaboDat/intra/intracell_annotation.csv")
metabo_ids <- metabo_ids %>%
  select(ionIdx, `label (bona fide)`)
setnames(metabo_ids, c("ionIdx", "label"))

# add name for each metabolite
mymatcher = match(metabo_int[,ionIdx],metabo_ids[,ionIdx])
metabo_names = metabo_ids[mymatcher, label]
metabo_int[,label := metabo_names]


# create meta data table
expid=colnames(metabo_int)[3:(ncol(metabo_int)-1)]
sampleid=sub("-\\d","",sub(".+E\\d ","",expid))
repid=sub(".+-","",expid)
batchid=sub("(E\\d+).+","\\1",sub(".+E","E",expid))
expcl=data.table(sampleid,repid,batchid,expid)
mymatcher=match(expcl[,sampleid],id_tbl[,sampleid])
mmaids=id_tbl[mymatcher,mmaid]
expcl[,mmaid:=mmaids]
expcl[,biorep:=sub("\\d+ / ","",expid)]


###################################################################
# create data table of metabolites


# collapse technical replicates into PSEUDO-biological replicates and melt table (and add experimental info)
biorepnames = unique(expcl[,biorep])
metabo_samplerep = do.call("cbind", lapply(biorepnames, function(x) {
  rowMeans(metabo_int %>% select(contains(x)))
}))
metabo_samplerep <- data.table(metabo_samplerep)
setnames(metabo_samplerep, biorepnames)
metabo_samplerep[, ionIdx := metabo_int[,"ionIdx"]]
metabo_samplerep[, label := metabo_int[, "label"]]
metabo_samplerep <- melt(metabo_samplerep, id.vars = c("ionIdx", "label"))
setnames(metabo_samplerep, c("ionIdx", "label", "biorep", "value"))
metabo_samplerep <- metabo_samplerep %>%
  mutate(br = biorep) %>%
  separate(br, c("expid", "rest"), sep = " ") %>%
  separate(rest, c("sampleid", "repid"), sep = "-") %>%
  mutate(si = sampleid) %>%
  separate(si, c("type", "cellno"), sep = 1)
metabo_samplerep <- metabo_samplerep[,-"ionIdx"]
metabo_samplerep_spread <- dcast(metabo_samplerep, biorep + expid + sampleid + repid + type + cellno ~ label, value.var = "value")
metabo_samplerep_spread[, type2 := type]
metabo_samplerep_spread$type <- factor(metabo_samplerep_spread$type, levels = c("C", "M"))
levels(metabo_samplerep_spread$type) <- c("control", "MMUT def.")

###################################################################
# create data table of metabolites as well as protein and RNA information


# collapse technical replicates into pseudo-biological replicates
biorepnames = unique(expcl[,biorep])
metabo_biorep = do.call("cbind", lapply(biorepnames, function(x) {
  rowMeans(metabo_int %>% select(contains(x)))
  }))
metabo_biorep <- data.table(metabo_biorep)
setnames(metabo_biorep, biorepnames)

# collapse pseudo-biological replicates into biological replicates
cellnames = unique(expcl[,sampleid])
cellnames_ = paste(cellnames, "-", sep="")
metabo_biorep = do.call("cbind", lapply(cellnames_, function(x) {
  apply((metabo_biorep %>% select(contains(x))), 1, median)
  })) # calculate median instead of mean to avoid influence of outliers
metabo_biorep <- data.table(metabo_biorep)
setnames(metabo_biorep, cellnames)
metabo_biorep[, label := metabo_int[, "label"]]
metabo_biorep <- melt(metabo_biorep, id.vars = c("label"))
setnames(metabo_biorep, c("label", "cellname", "value"))
metabo_biorep <- metabo_biorep %>%
  mutate(type = cellname) %>%
  separate(type, c("type", "cellno"), sep = 1)
metabo_biorep <- metabo_biorep[,-"cellno"]
metabo_biorep_spread <- dcast(metabo_biorep, cellname + type ~ label, value.var = "value")
mimatcher = match(metabo_biorep_spread[,cellname], id_tbl[,sampleid])
sample_names = id_tbl[mimatcher, mmaid]
metabo_biorep_spread[, mmaid := sample_names]
col_idx <- grep("mmaid", names(metabo_biorep_spread))
metabo_biorep_spread <- data.frame(metabo_biorep_spread)
metabo_biorep_spread <- metabo_biorep_spread[, c(col_idx, (1:ncol(metabo_biorep_spread))[-col_idx])]
metabo_biorep_spread <- data.table(metabo_biorep_spread)

metabo_biorep_spread$type <- factor(metabo_biorep_spread$type)
levels(metabo_biorep_spread$type) <- c("control", "MMUT def.")


##################################
# load protein matrix
load("interimData/prot_mat.RData")

prot_exp <- data.table(ENSG = row.names(prot_mat), prot_mat)
prot_exp[,PG.Qvalue:=NULL] 
prot_exp2 <- melt(prot_exp, id = c("ENSG"))

# there is a duplicate gene, to be removed:
prot_exp3 <- prot_exp2 %>% distinct(ENSG, variable, .keep_all = TRUE)
prot_exp4 <- dcast(prot_exp3, variable ~ ENSG, value.var = c("value"))
# row.names(prot_exp4) <- paste(prot_exp4$variable, row.names(prot_exp4), sep = "_")

rel_mmaids <- metabo_biorep_spread$mmaid
prot_exp5 <- prot_exp4[rel_mmaids,]
colnames(prot_exp5)[1] <- "mmaid"

# fetch HGNC gene ID
P_list <- data.table(read.table("Data/gene_sets/ENSG_HGNC_list_p.csv", header =  TRUE))
genes <- colnames(prot_exp5)[-1]

mymatch = match(genes, P_list[,ensembl_gene_id])
p_symbols = P_list[mymatch, hgnc_symbol]

p_symbols_1 <- paste("p", p_symbols, sep = "_")
colnames(prot_exp5)[-1] <- p_symbols_1

# remove duplicated columns
sum(duplicated(names(prot_exp5)))
prot_exp5[, which(duplicated(names(prot_exp5))) := NULL]


##################################
# load transcript matrix

load("interimData/rnafiltered_mat.RData")


rna_exp <- data.table(ENSG = row.names(rnafiltered_mat), rnafiltered_mat)
rna_exp2 <- melt(rna_exp, id = c("ENSG"))

length(unique(rna_exp2[variable == "MMA001",]$ENSG)) == length(rna_exp2[variable == "MMA001",]$ENSG)
# mma1_genes <- prot_exp2[variable == "MMA001",]$ENSG
# mma1_genes[duplicated(mma1_genes)]
# prot_exp2[ENSG == mma1_genes[duplicated(mma1_genes)],]

rna_exp3 <- dcast(rna_exp2, variable ~ ENSG, value.var = c("value"))
# row.names(rna_exp3) <- paste(rna_exp3$variable, row.names(rna_exp3), sep = "_")

rel_mmaids <- metabo_biorep_spread$mmaid
rna_exp4 <- rna_exp3[rel_mmaids,]
colnames(rna_exp4)[1] <- "mmaid"

# fetch HGNC gene ID
G_list <- data.table(read.table("Data/gene_sets/ENSG_HGNC_list_r.csv", header =  TRUE))
rna_genes <- colnames(rna_exp4)[-1]

mymatch = match(rna_genes, G_list[,ensembl_gene_id])
r_symbols = G_list[mymatch, hgnc_symbol]

r_symbols_1 <- paste("r", r_symbols, sep = "_")
colnames(rna_exp4)[-1] <- r_symbols_1

# remove duplicated columns
sum(duplicated(names(rna_exp4)))
rna_exp4[, which(duplicated(names(rna_exp4))) := NULL]


##################################
# merge metabolite and protein as well as transcript table
dt1 <- merge(metabo_biorep_spread, prot_exp5, by = "mmaid", all = TRUE)
dt2 <- merge(dt1, rna_exp4, by = "mmaid", all = TRUE)

dt2$type <- factor(dt2$type)
levels(dt2$type) <- c("control", "MMUT def.")


########################################################################
# PCA plot of metabolomics data
########################################################################
########################################################################
# pseudobiological replicates

pca_tbl <- metabo_samplerep_spread
row.names(pca_tbl) <- paste(pca_tbl$type2, pca_tbl$biorep, sep = "_")
colnames(metabo_samplerep_spread)[1:10]

pca_tbl$type2 <- NULL
pca_tbl$biorep <- NULL
pca_tbl$expid <- NULL
pca_tbl$sampleid <- NULL
pca_tbl$repid <- NULL
pca_tbl$type <- NULL
pca_tbl$cellno <- NULL


pca_tbl2 <- scalerRound(pca_tbl, 7)
pca_tbl2 <- data.table(pca_tbl2)

hist(pca_tbl[[500]], breaks = 30)
hist(pca_tbl2[[500]], breaks = 30)

# PCA plot
met_pca <- prcomp(pca_tbl2)
plot(met_pca$x[,1], met_pca$x[,2])

df_met_pca <- as.data.frame(met_pca$x)
df_met_pca$no <- row.names(pca_tbl)
df_met_pca$type <- sub("_.*", "", df_met_pca$no)

#calculate percentages
perc_met_pca <- vector("character", nrow(met_pca$x))

for(i in 1:nrow(met_pca$x)) {
  perc_met_pca[i] <- round(((met_pca$sdev[i])^2) / sum(met_pca$sdev^2) * 100, 2)
  perc_met_pca[i] <- paste(colnames(df_met_pca)[i], " (", paste(as.character(perc_met_pca[i]), "%", ")", sep = ""), sep = "")
}

# scree plot with package factoextra
scree_plot0 <- fviz_eig(met_pca)

scree_plot <- 
  scree_plot0 +
  theme_pubr() +
  theme(plot.title = element_text(size = 12))



met_pca_plt <-
ggplot(df_met_pca, aes(x = PC1, y = PC2, color = type)) +
  geom_point() +
  ggtitle("Metabolomics PCA") +
  xlab(perc_met_pca[1]) +
  ylab(perc_met_pca[2]) +
  scale_color_aaas() +
  theme_pubr() +
  theme(plot.title = element_text(size = 12), legend.title = element_blank(), legend.position = "right")


met_pca <- 
ggarrange(met_pca_plt, scree_plot, widths = c(1.2, 1))


ggsave(paste(fig_path, "met_pca.png", sep = ""), met_pca, device = png(), width = 7, height = 3)
dev.off()
ggsave(paste(fig_path_pdf, "SuppFig6/", "met_pca.pdf", sep = ""), met_pca, device = "pdf", width = 7, height = 3)
dev.off()






########################################################################
# volcano plot based on pre-calculated values
########################################################################


metabos_of_interest <- c("2-Methylcitric acid", "Pyruvate", "Oxoadipate", "Aminoadipate", "Hexose")
metabos_of_interest2 <- c("Glutamine", "Alanine")

tbl_calc[ , label1 := NULL]
tbl_calc[ , label2 := NULL]

for(i in 1:length(metabos_of_interest)) {
  tbl_calc[which(tbl_calc$name %in% metabos_of_interest[i]), label1 := as.character(metabos_of_interest[i])]
}
summary(!is.na(tbl_calc$label1))
for(i in 1:length(metabos_of_interest2)) {
  tbl_calc[which(tbl_calc$name %in% metabos_of_interest2[i]), label2 := as.character(metabos_of_interest2[i])]
}


metabo_volcano <- 
ggplot(tbl_calc, aes(x = log2_fc, y = -log10(p_val_adj))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, alpha = 0.6) +
  geom_point(alpha = 0.5) +
  geom_point(data = tbl_calc[!is.na(label1), ], inherit.aes = FALSE, aes(x = log2_fc, y = -log10(p_val_adj)), color = "darkorange", size = 3) +
  geom_text_repel(data = tbl_calc[!is.na(label1), ], aes(label = name), color = "darkorange", nudge_x = 1, nudge_y = 0) +
  geom_point(data = tbl_calc[!is.na(label2), ], inherit.aes = FALSE, aes(x = log2_fc, y = -log10(p_val_adj)), color = mypal2[3], size = 3) +
  geom_text_repel(data = tbl_calc[!is.na(label2), ], aes(label = name), color = mypal2[3], nudge_x = -0.7) +
  scale_x_continuous(limits = c(-1.8,2.1), breaks = c(-1,0,1)) +
  ggtitle("Differential metabolite expression") +
  xlab("log2(fold change)") +
  ylab("-log10(adjusted p value)") +
  theme_pubr() +
  theme(plot.title = element_text(size = 12))

ggsave(paste(fig_path,"Volcano_metabolites.png", sep = ""), metabo_volcano, device = png(), width = 4, height = 4)
dev.off()
ggsave(paste(fig_path_pdf,"Fig5/","Volcano_metabolites.pdf", sep = ""),metabo_volcano, device = "pdf", width = 4, height = 3)
dev.off()




########################################################################
# TCA cycle plots
########################################################################
# test which metabolites are in the dataset
colnames(metabo_samplerep_spread)[c(grep("Hexose", colnames(metabo_samplerep_spread), ignore.case=TRUE))]

# check all the possible IDs as reported in the annotation sheet
metabo_ids[c(grep("carnitine", metabo_ids$label, ignore.case=TRUE))]

# TCA metabolites
# tca_ions <- c("Pyruvate", "Alanine", "Phosphoenolpyruvate", "2-Methylcitric acid", "Acetyl citrate", "(Iso)Citrate", "2-oxoglutarate(2-)", "Succinate", "Fumarate", "Malate", "3-Propylmalate", "Aspartate", "Glutamate", "Glutamine", "Aminoadipate", "Oxoadipate")
tca_ions1 <- c("Fumarate","Malate", "Aspartate")
tca_ions2 <- c("Pyruvate", "Alanine", "2-Methylcitric acid", "(Iso)Citrate", "2-oxoglutarate(2-)", "Succinate", "Glutamate", "Glutamine", "Aminoadipate", "Oxoadipate")
tca_ions3 <- c("Pyruvate", "Alanine", "2-Methylcitric acid", "(Iso)Citrate", "2-oxoglutarate(2-)", "Succinate", "Glutamate", "Glutamine", "Aminoadipate", "Oxoadipate", "Fumarate","Malate", "Aspartate")
tca_mets <- metabo_samplerep_spread[,c("type2", "expid", tca_ions3), with = FALSE]
colnames(tca_mets)

# normalize to control
ctrl_means <- tca_mets[type2 == "C", lapply(.SD, mean)]

ff=lapply(c(1:length(tca_ions3)),function(i){
  print(i)
  tca_mets[, tca_ions3[i], with = FALSE]/as.numeric(ctrl_means[, tca_ions3[i], with = FALSE])
})
ff2=do.call("cbind",ff)
ff3=data.table(ff2)

tca_mets2 <- cbind(ff3, tca_mets[, "type2"])

tca_mets3 <- melt.data.table(tca_mets2, id.vars = c("type2"))


# exclude massive 2-oxoglutarate outlier
tca_mets3[variable == "2-oxoglutarate(2-)" & value >1.6, ]$value <- NA
tca_mets3[, type3 := "control"]
tca_mets3[type2 == "M", type3 := "MMUT def."]




# plot
compare <- list(c("control", "MMUT def."))

tca_mets_plt <- 
ggplot(tca_mets3[which(tca_mets3$variable %in% tca_ions2), ], aes(x = type3, y = value, color = type3))+
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width = 0.2,color = "black", fill = "white", shape = 21, show.legend = FALSE) +
  facet_wrap(~variable, scales = "free_y", nrow = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  stat_compare_means(comparisons = compare, size = 3, method = "wilcox.test", label = "p.format") +
  scale_fill_aaas() +
  scale_color_aaas() +
  theme_test() +
  ylab("Relative ion abundance normalized to control") +
  ggtitle("TCA metabolite levels") +
  theme(axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), strip.background = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color = "black"), axis.ticks.x = element_blank(), plot.title = element_text(size = 12))

ggsave(paste(fig_path,"TCA_mets_plots.png", sep = ""), tca_mets_plt, device = png(), width = 6, height = 5)
dev.off()
ggsave(paste0(fig_path_pdf,"Fig5/","TCA_mets_plots.pdf"), tca_mets_plt, device = "pdf", width = 6, height = 5)
dev.off()


tca_mets_plt2 <- 
ggplot(tca_mets3[which(tca_mets3$variable %in% tca_ions1), ], aes(x = type3, y = value, color = type3))+
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width = 0.2,color = "black", fill = "white", shape = 21, show.legend = FALSE) +
  facet_wrap(~variable, scales = "free_y", nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  stat_compare_means(comparisons = compare, size = 3, method = "wilcox.test", label = "p.format") +
  scale_fill_aaas() +
  scale_color_aaas() +
  theme_test() +
  ylab("Relative ion abundance\nnormalized to control") +
  ggtitle("TCA metabolite levels") +
  theme(axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), strip.background = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color = "black"), axis.ticks.x = element_blank(), plot.title = element_text(size = 12))

ggsave(paste(fig_path,"TCA_metsSelect_plots.png", sep = ""), tca_mets_plt2, device = png(), width = 3.6, height = 2.5)
dev.off()
ggsave(paste0(fig_path_pdf,"SuppFig_metabolomics/","TCA_metsSelect_plots.pdf"), tca_mets_plt2, device = "pdf", width = 3.6, height = 3)
dev.off()








# relationships of adjacent metabolites in the TCA cycle

pyr_mca <- 
ggplot(metabo_samplerep_spread[`2-Methylcitric acid`<25000,], aes(x = Pyruvate, y =  `2-Methylcitric acid`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.01, r.accuracy = 0.01, size = 3, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")))+
  scale_color_aaas() +
  theme_test() +
  theme(legend.position = "none")

pyr_citrate <- ggplot(metabo_samplerep_spread[`(Iso)Citrate`>1600000,], aes(x = Pyruvate, y =  `(Iso)Citrate`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  scale_color_aaas() +
  theme_test() +
  theme(legend.position = "none")

citrate_oxo <- ggplot(metabo_samplerep_spread[`2-oxoglutarate(2-)`<10000,], aes(x = `(Iso)Citrate`, y =  `2-oxoglutarate(2-)`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  scale_color_aaas() +
  theme_test() +
  theme(legend.position = "none")

oxo_succinate <- ggplot(metabo_samplerep_spread[`2-oxoglutarate(2-)`<10000,], aes(x = `2-oxoglutarate(2-)`, y =  Succinate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  scale_color_aaas() +
  theme_test() +
  theme(legend.position = "none")

succinate_fumarate <- ggplot(metabo_samplerep_spread[Fumarate<200000, ], aes(x = Succinate, y =  Fumarate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  scale_color_aaas() +
  theme_test() +
  theme(legend.position = "none")

fumarate_malate <- ggplot(metabo_samplerep_spread[Fumarate<200000, ], aes(x = Fumarate, y =  Malate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  scale_color_aaas() +
  theme_test() +
  theme(legend.position = "none")


oxo_glutamate <- ggplot(metabo_samplerep_spread[`2-oxoglutarate(2-)`<10000 & Glutamate>8e+06,], aes(x = `2-oxoglutarate(2-)`, y =  Glutamate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  scale_color_aaas() +
  theme_test() +
  theme(legend.position = "none")

oxo_glutamine <- ggplot(metabo_samplerep_spread[`2-oxoglutarate(2-)`<10000,], aes(x = `2-oxoglutarate(2-)`, y =  Glutamine, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  scale_color_aaas() +
  theme_test() +
  theme(legend.position = "none")

glutamate_aspartate <- ggplot(metabo_samplerep_spread, aes(x = Glutamate, y =  Aspartate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  scale_color_aaas() +
  theme_test() +
  theme(legend.position = "none")

glutamine_glutamate <- ggplot(metabo_samplerep_spread[Glutamate>8e+06, ], aes(x = Glutamine, y =  Glutamate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  scale_color_aaas() +
  theme_test() +
  theme(legend.position = "none")


TCA_sum_plt <- ggarrange(pyr_mca, pyr_citrate, citrate_oxo,oxo_succinate, succinate_fumarate, fumarate_malate, oxo_glutamate, oxo_glutamine, glutamate_aspartate, glutamine_glutamate)

ggsave(paste(fig_path,"TCA_relationship_plots.png", sep = ""), TCA_sum_plt, device = png(), width = 13, height = 7, bg = "white")
dev.off()










########################################################################
# top changes
########################################################################

# top metabolites
tbl_calc_sub <- data.table(tbl_calc, key = "p_val_adj")
top_metabos <- head(tbl_calc_sub, 26)$name
top_metabos1 <- c("type2", "expid", top_metabos[-13])

top_mets0 <- metabo_samplerep_spread[, ..top_metabos1]
top_mets <- melt(top_mets0, id.vars = c("type2", "expid"))

ggplot(top_mets, aes(x = type2, y = value, color = type2, fill = type2))+
  geom_boxplot(alpha = .5, width.errorbar = 0.2) +
  facet_wrap(~variable, scales = "free_y") +
  scale_fill_aaas() +
  scale_color_aaas() +
  theme_test() +
  ylab("Metabolite level") +
  theme(axis.title.x = element_blank(), legend.position = "none")

ggsave(paste(fig_path,"TOP_mets_plots.png", sep = ""), device = png(), width = 9, height = 6)
dev.off()








########################################################################
# metabolite-protein/transcript plots
########################################################################


# test if specific metabolite is present in dataset
# check in the dataset with one ID assigned to one metabolite, even if several IDs are possible
sum(unique(dt2[,label]) == "citric")
colnames(dt2)[c(grep("PDH", colnames(dt2)))]

# check all the possible IDs as reported in the annotation sheet
metabo_ids[c(grep("citr", metabo_ids$label))]


r_PDH_pyr <- ggplot(dt2, aes(x = r_PDHB, y = Pyruvate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

p_PDH_pyr <- ggplot(dt2, aes(x = p_PDHX, y = Pyruvate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

r_OGDH_Gln <- ggplot(dt2, aes(x = r_OGDH, y = Glutamine, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

p_OGDH_Gln <- ggplot(dt2, aes(x = p_OGDH/10e5, y = Glutamate/10e6, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

r_GLS_Gln <- ggplot(dt2, aes(x = r_GLS, y = Glutamine, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

p_GLS_Gln <- ggplot(dt2, aes(x = p_GLS, y = Glutamate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

r_GOT2_Glu <- ggplot(dt2, aes(x = r_GOT2, y = Glutamate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

p_GOT2_Glu <- ggplot(dt2, aes(x = p_GOT2, y = Glutamate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

r_GLUD1_Gln <- ggplot(dt2, aes(x = r_GLUD1, y = Glutamine, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

p_GLUD1_Gln <- ggplot(dt2, aes(x = p_GLUD1, y = Glutamine, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

r_GLUD1_Glu <- ggplot(dt2, aes(x = r_GLUD1, y = Glutamate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

p_GLUD1_Glu <- ggplot(dt2[p_GLUD1 < 8e5, ], aes(x = p_GLUD1, y = Glutamate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

r_GLUD1_ketoglt <- ggplot(dt2[dt2$`X2.oxoglutarate.2..`<8000,], aes(x = r_GLUD1, y = `X2.oxoglutarate.2..`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

p_GLUD1_ketoglt <- ggplot(dt2[dt2$`X2.oxoglutarate.2..`<8000,], aes(x = p_GLUD1, y = `X2.oxoglutarate.2..`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

r_OGDH_ketoglt <- ggplot(dt2[dt2$`X2.oxoglutarate.2..`<8000,], aes(x = `X2.oxoglutarate.2..`, y = r_OGDH, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

p_OGDH_ketoglt <- ggplot(dt2[dt2$`X2.oxoglutarate.2..`<8000,], aes(x = `X2.oxoglutarate.2..`, y = p_OGDH, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

r_MMUT_ketoglt <- ggplot(dt2[dt2$`X2.oxoglutarate.2..`<8000,], aes(x = `X2.oxoglutarate.2..`, y = r_MMUT, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()


metabolite_r_p_plts <- ggarrange(r_PDH_pyr, p_PDH_pyr, r_GLS_Gln, p_GLS_Gln, r_GOT2_Glu, p_GOT2_Glu, r_GLUD1_Gln, p_GLUD1_Gln, r_GLUD1_Glu, p_GLUD1_Glu, r_GLUD1_ketoglt, p_GLUD1_ketoglt, r_OGDH_ketoglt, p_OGDH_ketoglt)

ggsave(paste(fig_path,"metabolite_r_p_plots.png", sep = ""), metabolite_r_p_plts, device = png(), width = 14, height = 9, bg= "white")
dev.off()



# methylcitric acid plots

mca_rMMUT <- ggplot(dt2, aes(x = r_MMUT, y = `X2.Methylcitric.acid`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

mca_pMMUT <- ggplot(dt2, aes(x = p_MMUT, y = `X2.Methylcitric.acid`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

mca_rCS <- ggplot(dt2, aes(x = r_CS, y = `X2.Methylcitric.acid`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

mca_pCS <- ggplot(dt2, aes(x = p_CS, y = `X2.Methylcitric.acid`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

mca_rPHDX <- ggplot(dt2, aes(x = r_PDHX, y = `X2.Methylcitric.acid`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

mca_pPHDX <- ggplot(dt2, aes(x = p_PDHX, y = `X2.Methylcitric.acid`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

mca_rPDHA1 <- ggplot(dt2, aes(x = r_PDHA1, y = `X2.Methylcitric.acid`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

mca_pPDHA1 <- ggplot(dt2, aes(x = p_PDHA1, y = `X2.Methylcitric.acid`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

mca_rPDHB <- ggplot(dt2, aes(x = r_PDHB, y = `X2.Methylcitric.acid`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

mca_pPDHB <- ggplot(dt2, aes(x = p_PDHB, y = `X2.Methylcitric.acid`, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()


mca_plts <- ggarrange(mca_rMMUT, mca_pMMUT, mca_rCS, mca_pCS, mca_rPHDX, mca_pPHDX, mca_rPDHA1, mca_pPDHA1, mca_rPDHB, mca_pPDHB)

ggsave(paste(fig_path,"mca_plots.png", sep = ""), mca_plts, device = png(), width = 14, height = 7, bg = "white")
dev.off()






########################################################################
# collection of metabolite-transcript-peptide interactions
########################################################################


fumarate_malate_2 <- 
ggplot(metabo_samplerep_spread[Fumarate<200000, ], aes(x = Fumarate/10e4, y =  Malate/10e5, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 3, cor.coef.name = "rho", label.x.npc = "right", label.y.npc = 0.15, hjust = 1, show.legend = FALSE)+
  scale_color_aaas() +
  labs(y="Malate", x="Fumarate") +
  theme_test()


glutamine_glutamate_2 <- 
ggplot(metabo_samplerep_spread[Glutamate>8e+06, ], aes(x = Glutamine/10e6, y =  Glutamate/10e6, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 3, cor.coef.name = "rho", label.x.npc = "right", label.y.npc = 0.15, hjust = 1, show.legend = FALSE)+
  scale_color_aaas() +
  labs(y="Glutamate", x="Glutamine") +
  theme_test()


pyr_citrate_2 <- 
ggplot(metabo_samplerep_spread[`(Iso)Citrate`>1600000,], aes(x = Pyruvate/10e3, y =  `(Iso)Citrate`/10e5, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 3, cor.coef.name = "rho", label.x.npc = "left", label.y.npc = "top", hjust = 0, show.legend = FALSE)+
  scale_color_aaas() +
  labs(y="(Iso)Citrate", x="Pyruvate") +
  theme_test()


oxo_glutamate2 <- 
ggplot(metabo_samplerep_spread[`2-oxoglutarate(2-)`<10000 & Glutamate>8e+06,], aes(y = `2-oxoglutarate(2-)`/10e3, x =  Glutamate/10e6, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 3, cor.coef.name = "rho", label.x.npc = "right", label.y.npc = "top", hjust = 1, show.legend = FALSE)+
  scale_color_aaas() +
  labs(x="Glutamate", y="2-oxoglutarate") +
  theme_test()


glutamate_aspartate2 <- 
ggplot(metabo_samplerep_spread, aes(x = Glutamate/10e6, y =  Aspartate/10e4, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 3, cor.coef.name = "rho", label.x.npc = "left", label.y.npc = "top", hjust = 0, show.legend = FALSE)+
  scale_color_aaas() +
  labs(y="Aspartate", x="Glutamate") +
  theme_test()



r_MMUT_ketoglt_2 <- 
ggplot(dt2[dt2$`X2.oxoglutarate.2..`<8000,], aes(x = `X2.oxoglutarate.2..`/1000, y = r_MMUT, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 3, cor.coef.name = "rho", label.x.npc = "left", label.y.npc = 0.15, hjust = 0, show.legend = FALSE)+
  scale_color_aaas() +
  labs(y="MMUT transcript", x="2-oxoglutarate") +
  theme_test()

r_PDHB_pyr <- 
ggplot(dt2, aes(x = r_PDHB/10, y = Pyruvate/10e3, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 3, cor.coef.name = "rho", label.x.npc = "left", label.y.npc = "top", hjust = 0, show.legend = FALSE)+
  scale_color_aaas() +
  labs(y="Pyruvate", x="PDHB transcript") +
  theme_test()


p_PDHX_pyr <- 
ggplot(dt2, aes(x = p_PDHX/10e4, y = Pyruvate/10e3, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 3, cor.coef.name = "rho", label.x.npc = "left", label.y.npc = "top", hjust = 0, show.legend = FALSE)+
  scale_color_aaas() +
  labs(y="Pyruvate", x="PDHX protein") +
  theme_test()


mixed_mets_plt <- 
fumarate_malate_2 + glutamine_glutamate_2 + pyr_citrate_2 + oxo_glutamate2 + glutamate_aspartate2 + r_MMUT_ketoglt_2 + r_PDHB_pyr + p_PDHX_pyr + plot_layout(guides = "collect") + plot_annotation(title = "Transcript-protein-metabolite correlations", theme = theme(plot.title = element_text(size = 12))) & theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black"), legend.title = element_blank())

mixed_mets_plt2 <- 
fumarate_malate_2 + glutamine_glutamate_2 + pyr_citrate_2 + oxo_glutamate2 + glutamate_aspartate2 + r_MMUT_ketoglt_2 + r_PDHB_pyr + p_PDHX_pyr + plot_layout(guides = "collect", ncol=2) + plot_annotation(title = "Transcript-protein-metabolite correlations", theme = theme(plot.title = element_text(size = 12))) & theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black"), legend.title = element_blank(), legend.position = "bottom")

ggsave(paste0(fig_path_pdf,"Fig5_metabolomics/","MIXED_mets_plots2.pdf"), mixed_mets_plt2, device = "pdf", width = 5.2, height = 10)
dev.off()


mixed_mets_plt3 <- 
r_PDHB_pyr + p_PDHX_pyr + fumarate_malate_2 + pyr_citrate_2 + glutamine_glutamate_2 + oxo_glutamate2 + glutamate_aspartate2 + r_MMUT_ketoglt_2 + plot_layout(guides = "collect", ncol=2) + plot_annotation(title = "Transcript-protein-metabolite correlations", theme = theme(plot.title = element_text(size = 12))) & theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black"), legend.title = element_blank(), legend.position = "bottom")

ggsave(paste0(fig_path_pdf,"Fig5_metabolomics/","MIXED_mets_plots3.pdf"), mixed_mets_plt3, device = "pdf", width = 5.2, height = 10)
dev.off()





r_PDH_pyr <- ggplot(dt2, aes(x = r_PDHB, y = Pyruvate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()

p_PDH_pyr <- ggplot(dt2, aes(x = p_PDHX, y = Pyruvate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..p.label..)))+
  theme_bw()




ggplot(metabo_samplerep_spread[,], aes(x = "X2.Methylcitric.acid", y =  Pyruvate, color = type)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 3, cor.coef.name = "rho", label.x.npc = "left", label.y.npc = "bottom", hjust = 0, show.legend = FALSE)+
  scale_color_aaas() +
  labs(y="Pyruvate", x="2-Methylcitric acid") +
  theme_test()