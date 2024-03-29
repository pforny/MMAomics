

#### Circos plot of TCA genes

#libraries
#####
library(data.table)
library(dplyr)
library(tidyverse)
library(corrplot)
library(ggpubr)
library(ggplot2)
library(corrplot)
library(pheatmap)
library(gridExtra)
library(circlize)
library(tibble)
library(cowplot)
library(gridGraphics)
library(GeneNet)

library(ComplexHeatmap) #for legend
library(gridBase)
library(graph4lg) #to reorder the rows and columns of a symmetric matrix

# create figures path
system("mkdir Figs/TCA")
fig_path <- c("Figs/thesis/fig4/")
fig_path_pdf <- c("Figs/v8/pdf/")

################
# import data
################
TCA <- "Data/gene_sets/masterTCA_keggPlusmanual_circos.csv"

tca_tbl <- fread(TCA, header = TRUE)[,ENSG]

hgnc_r <- fread("Data/gene_sets/ENSG_HGNC_list_r.csv")
setnames(hgnc_r, c("ENSG", "HGNC"))
hgnc_p <- fread("Data/gene_sets/ENSG_HGNC_list_p.csv")
setnames(hgnc_p, c("ENSG", "HGNC"))

# load transcriptomics and protein data
load("interimData/rnafiltered_mat.RData")#14' 749 gene entries
load("interimData/prot_mat.RData") #4'406 gene entries


##### prep protein and transcriptome table
prot_exp <- data.table(ENSG = row.names(prot_mat), prot_mat)
summary(unique(prot_exp$ENSG) == prot_exp$ENSG)
prot_exp <- prot_exp %>% distinct(ENSG, .keep_all = TRUE)
prot_exp_samples <- prot_exp[, -"ENSG"]
prot_exp_prots <- t(prot_exp_samples)

rna_exp <- data.table(ENSG = row.names(rnafiltered_mat), rnafiltered_mat)
summary(unique(rna_exp$ENSG) == rna_exp$ENSG)
rna_exp_samples <- rna_exp[, -"ENSG"]
rna_exp_rnas <- t(rna_exp_samples)


##### merge protein and rna, CAVE: only for common genes and samples
rna_exp_common_genes <- rna_exp[which(rna_exp[, ENSG] %in% prot_exp[, ENSG]), ]
kept_samples_prot <- (which(colnames(prot_exp)[-1] %in% colnames(rna_exp)[-1])+1)
prot_exp_common_genes <- prot_exp[, ..kept_samples_prot]
prot_exp_common_genes[, ENSG := prot_exp[, ENSG]]

rna_exp_common_genes_MMUT <-rna_exp_common_genes[,1:144]
rna_exp_common_genes_rest <-rna_exp_common_genes[, 145:222]

prot_exp_common_genes<- prot_exp_common_genes %>% relocate(ENSG)
prot_exp_common_genes_MMUT <- prot_exp_common_genes[, 1:144]
prot_exp_common_genes_rest <- prot_exp_common_genes[, 145:222]

prot_exp_common_genes_rest$ENSG <- prot_exp_common_genes$ENSG
rna_exp_common_genes_rest$ENSG <- rna_exp_common_genes$ENSG

merged_rna_prot_MMUT <- merge(rna_exp_common_genes_MMUT, prot_exp_common_genes_MMUT,  by = "ENSG") # only genes shared in both
# column 1 -> ENSG
# columns 2:144 -> rna data (noted with .x in name)
# columns 145:287 -> protein data (noted with .y in name)

merged_rna_prot_rest <- merge(rna_exp_common_genes_rest, prot_exp_common_genes_rest, by ="ENSG")
# column 1 -> ENSG
# columns 2:79 -> rna
# columns 80:157 -> protein



##### calculate protein-rna correlations in the MMUT and REST group
spearman_cor_vector_MMUT <- vector()
spearman_cor_vector_rest <- vector()

for(i in 1:dim(merged_rna_prot_MMUT)[1]) {
  cat(i, "\n")
  spearman_cor_vector_MMUT <- c(spearman_cor_vector_MMUT, cor(unlist(merged_rna_prot_MMUT[i, 2:144]), unlist(merged_rna_prot_MMUT[i, 145:287]), use = "pairwise.complete.obs", method = "spearman"))
}

for(i in 1:dim(merged_rna_prot_rest)[1]){
  cat(i, "\n")
  spearman_cor_vector_rest <- c(spearman_cor_vector_rest, cor(unlist(merged_rna_prot_rest[i, 2:79]), unlist(merged_rna_prot_rest[i, 80:157]), use = "pairwise.complete.obs", method = "spearman"))
}


##### correlation data in data.table
cor_spearman_MMUT <- data.table(spearman_cor_vector_MMUT)
cor_spearman_MMUT[, ENSG := merged_rna_prot_MMUT[, ENSG]] 
cor_spearman_MMUT <- cor_spearman_MMUT[, c(2,1)]

cor_spearman_rest <- data.table(spearman_cor_vector_rest)
cor_spearman_rest[, ENSG := merged_rna_prot_rest[, ENSG]]
cor_spearman_rest <- cor_spearman_rest[,c(2,1)]# Use gene from TCA
spearman_MMUT_tca <- cor_spearman_MMUT[is.element(cor_spearman_MMUT$ENSG, tca_tbl),] 
spearman_rest_tca <- cor_spearman_rest[is.element(cor_spearman_rest$ENSG, tca_tbl),]

spearman_MMUT_tca <- as.data.frame(spearman_MMUT_tca)
spearman_MMUT_tca <- left_join(spearman_MMUT_tca, hgnc_r, by="ENSG")
spearman_MMUT_tca <- spearman_MMUT_tca[,c(3,2,1)]
spearman_MMUT_tca$ENSG <- NULL

spearman_rest_tca <- as.data.frame(spearman_rest_tca,)
spearman_rest_tca <- left_join(spearman_rest_tca, hgnc_r, by="ENSG")
spearman_rest_tca <-spearman_rest_tca[,c(3,2,1)]
spearman_rest_tca$ENSG<-NULL

df_spearman1 <- left_join(spearman_MMUT_tca, spearman_rest_tca, by="HGNC")
DF_spearman1 <- df_spearman1
rownames(df_spearman1) <- df_spearman1[,1]
df_spearman1[,1] <- NULL
colnames(df_spearman1)[1] <- "MMUT"
colnames(df_spearman1)[2] <- "rest"

df_spearman1 = as.matrix(df_spearman1)


##### import effect sizes from DEA lmm analysis
diffexp_prot <- read_table("interimData/DEAlmm/diff_exp_data_prot_pathwayact_all.txt", col_types = cols(betas = col_number()))
diffexp_rna <- read_table("interimData/DEAlmm/diff_exp_data_rnaseq_pathwayact_all.txt",col_types = cols(betas = col_number()))

diffexp_rna <- data.table(diffexp_rna)
diffexp_rna <- diffexp_rna %>% distinct(ensembl, .keep_all = TRUE)

diffexp_prot <- data.table(diffexp_prot)
diffexp_prot <- diffexp_prot %>% distinct(ensembl, .keep_all=TRUE)

diffexp_rna_common <- diffexp_rna[which(diffexp_rna[, ensembl] %in% diffexp_prot[,ensembl]),]
kept_samples_diffprot <- which(colnames(diffexp_prot)[-1] %in% colnames(diffexp_rna)[-1])+1
diffexp_prot_common <- diffexp_prot[, ..kept_samples_diffprot]

# filter for TCA genes
diffexp_rna_tca <- diffexp_rna_common[is.element(diffexp_rna_common$ensembl, tca_tbl), ]
diffexp_prot_tca <- diffexp_prot_common[is.element(diffexp_prot_common$ensembl, tca_tbl), ]

# ensure same order as in genome
genome1 <- DF_spearman1$HGNC

diffexp_rna_tca_df1 <- as.data.frame(diffexp_rna_tca)
de_rna_tca_df1 <- diffexp_rna_tca_df1[match(genome1, diffexp_rna_tca_df1$gname), ]
diffexp_prot_tca_df1 <- as.data.frame(diffexp_prot_tca)
de_prot_tca_df1 <- diffexp_prot_tca_df1[match(genome1, diffexp_prot_tca_df1$gname), ]


##### protein-protein correlation of TCA cycle and associated proteins for central chord diagram

# prepare protein dataset
prot_subset <- prot_exp[is.element(prot_exp$ENSG, tca_tbl), ]
prot_subset_t <- dcast(melt(prot_subset), variable ~ ENSG)
prot_subset_t[, variable := NULL]

new_names_prot <- hgnc_p[which(hgnc_p$ENSG %in% colnames(prot_subset_t)), ]$HGNC
setnames(prot_subset_t, new_names_prot)

prot_cor_all <- cor(prot_subset_t, use = "complete.obs", method = "spearman")
prot_cor_MMUT <- cor(prot_subset_t[1:150, ],  use = "complete.obs", method = "spearman")
prot_cor_rest <- cor(prot_subset_t[151:230, ],  use = "complete.obs", method = "spearman")

rownames <- row.names(prot_cor_MMUT)

# prepare transcript dataset
rna_subset <- rna_exp[is.element(rna_exp$ENSG, tca_tbl), ]
rna_subset_t <- dcast(melt(rna_subset), variable ~ ENSG)
rna_subset_t[, variable := NULL]
# HGNC_list_r <- data.table(read.table("Data/gene_sets/ENSG_HGNC_list_r.csv", header = TRUE))
# setnames(HGNC_list_r, c("ENSG", "HGNC"))
HGNC_list_r<- hgnc_r
new_names_rna <- HGNC_list_r[which(HGNC_list_r$ENSG %in% colnames(rna_subset_t)), ]$HGNC
setnames(rna_subset_t, new_names_rna)
rna_corr_all <- cor(rna_subset_t, method = "spearman")
rna_corr_MMUT <- cor(rna_subset_t[1:143, ], method = "spearman")
rna_corr_rest <- cor(rna_subset_t[144:nrow(rna_subset_t), ], method = "spearman")

### DF CHORD PLOT ###
prot_corr_eff <- prot_cor_rest - prot_cor_MMUT
prot_corr_eff1 <- prot_cor_MMUT - prot_cor_rest


##### generate df for plots

df_corr1 <- data.frame(sector = as.factor(c(rep(genome1, 2))),
                       x = c(rep(1, 37), rep(2, 37)),
                       y = c(df_spearman1[,"rest"], df_spearman1[,"MMUT"])) # y = c(df_spearman1[,"MMUT"], df_spearman1[,"rest"]))
df_dea1 <- data.frame(sector=as.factor(c(rep(genome1,2))),
                      x = c(rep(1,37), rep(2,37)), #1 encoding for rna, 2 encoding for protein
                      y2 = c(de_rna_tca_df1[,"betas"], de_prot_tca_df1[,"betas"]))

# change direction of bars
df_dea1 <- df_dea1 %>% 
  mutate( y2 = y2 *(-1))

# change character into factor and reorder levels to avoid alphabetical order of genes in plot
upper_boarder <- df_corr1$sector %>% unique() %>% length()
df_corr1$sector <- factor(df_corr1$sector, levels = c(as.character(unique(df_corr1$sector)[2:upper_boarder]), "MDH1")) # change order here!

gap_v <- c(rep(1, 27),12,rep(1, 9)) # gaps after sectors > 12 for gap in circos plot to label

### data diff. expression lvls rna and protein (most outer track)

# data complied in systems_plot_prot_rna_corr_PF.R 
# contains rna and protein lvls of all genes of interst for the circos plot
exp_lvl <- "Data/spearmanGenes_TCAcycle_expressionlevels_RnaProt_new.csv"
dt_exp_lvl <- fread(exp_lvl, header = TRUE, drop ="V1")

# gene_groups <- fread("Data/gene_sets/genes_in_groupsIso.csv", header = FALSE) #add the new genes here too; 37 entries
# gene_groups_orderd <- fread("Data/gene_sets/genes_in_groups_orderIso.csv", header = FALSE)

gene_groups <- fread("Data/gene_sets/genes_in_groups_01_02_21.csv", header = FALSE) #add the new genes here too; 37 entries
gene_groups_orderd <- fread("Data/gene_sets/genes_in_groups_order_01_03_21.csv", header = FALSE)

dt<-as.data.frame(dt_exp_lvl)
dt$HGNC<-factor(dt$HGNC, levels=c(as.character(unique(gene_groups$V1)[2:upper_boarder]),"IDH3G" ))
dt$type2<-factor(dt$type2)

# x = rna, y = protein
df = data.frame(sectors = dt$HGNC,
                x = dt$rna_lvl, y = dt$prot_lvl, type=dt$type2)

logav_exp_MMUT<-df %>% 
  group_by(sectors) %>% 
  filter(type =="MMUT") %>%
  summarise(logav_rna_mmut = mean(log2(x)), logav_prot_mmut = mean(log2(y)))

logav_exp_REST <- df%>%
  group_by(sectors) %>%
  filter(type =="REST") %>%
  summarise(logav_rna_rest = mean(log2(x)), logav_prot_rest = mean(log2(y)))


logav_exp <- left_join(logav_exp_MMUT, logav_exp_REST, by = "sectors")

fc_exp<-  logav_exp %>% 
  mutate(fc_rna = logav_rna_mmut- logav_rna_rest) %>% 
  mutate(fc_prot =logav_prot_mmut - logav_prot_rest ) %>% 
  as.matrix()

fc_exp_ <- data.frame(fc_exp[,-1], row.names=fc_exp[,1]) 
fc_exp_$fc_rna <- as.numeric(fc_exp_$fc_rna)
fc_exp_$fc_prot <- as.numeric(fc_exp_$fc_prot)
m_fc_exp<-data.matrix(fc_exp_[,5:6])

fc_exp_ <- data.frame(fc_exp[,-1], row.names=fc_exp[,1]) 
fc_exp_$fc_rna <- as.numeric(fc_exp_$fc_rna)
fc_exp_$fc_prot <- as.numeric(fc_exp_$fc_prot)
m_fc_exp<-data.matrix(fc_exp_[,5:6])


df_fc_exp <- data.frame(sector = as.factor(c(rep(fc_exp[,1], 2))),
                        x = c(rep(1, 37), rep(2, 37)),
                        y = c(m_fc_exp[, "fc_rna"], m_fc_exp[, "fc_prot"]))

df_fc_exp$sector <- factor(df_fc_exp$sector, levels = c(as.character(unique(gene_groups$V1)[4:upper_boarder2]),c("IDH3G", "OGDH", "DLD")))
######
#matrix zprot_corr_diff from  modify_chord_plt.R sollte direkt hier laufen
vMMUT<- prot_cor_MMUT[25,]
vCTRL <- prot_cor_rest[25,]

zvMMUT<- z.transform(vMMUT)
zvCTRL <- z.transform(vCTRL)

vdiff <- zvCTRL-zvMMUT
vdiff[25] <- 1

zprot_corr_diff <- as_tibble(prot_corr_eff) %>%
  mutate(across(, ~replace(.x, is.numeric(.x), 0))) %>%
  as.matrix()

zprot_corr_diff[25,] <- vdiff
zprot_corr_diff[,25] <- vdiff
row.names(zprot_corr_diff) <- rownames

# 
# 
# zprot_corr_diff <- as_tibble(prot_corr_diff) %>%
#   mutate(across(, ~replace(.x, is.numeric(.x), 0))) %>%
#   as.matrix()
# 
# zprot_corr_diff[25,] <- vdiff
# zprot_corr_diff[,25] <- vdiff
# row.names(zprot_corr_diff) <- rownames
# 


##
# prot_cor_rest
zvCTRLc<- zvCTRL
zvCTRLc[25] <- 1
zprot_cor_rest <- as_tibble(prot_cor_rest) %>%
  mutate(across(, ~replace(.x, is.numeric(.x), 0))) %>%
  as.matrix()
zprot_cor_rest[25,] <- zvCTRLc
zprot_cor_rest[,25] <- zvCTRLc
row.names(zprot_cor_rest) <- rownames

# MMUT chord data
zvMMUTc <- zvMMUT
zvMMUTc[25] <- 1
zprot_cor_MMUT <- as_tibble(prot_cor_MMUT) %>%
  mutate(across(, ~replace(.x, is.numeric(.x), 0))) %>%
  as.matrix()
zprot_cor_MMUT[25,] <- zvMMUTc
zprot_cor_MMUT[,25] <- zvMMUTc
row.names(zprot_cor_MMUT) <- rownames

# grouping of genes and order for plot
upper_boarder <- df_corr1$sector %>% unique() %>% length()
upper_boarder2 <- gene_groups$V1 %>% unique() %>% length()
df_corr1$sector<- factor(df_corr1$sector, levels=c(as.character(unique(gene_groups$V1)[4:upper_boarder2]),c("IDH3G", "OGDH", "DLD"))) #group_order$V1
# df_corr2$sector <-fct_relevel(df_corr2$sector, "SUCLG2")

df_dea1$sector<- factor(df_dea1$sector, levels=c(as.character(unique(gene_groups$V1)[4:upper_boarder2]),c("IDH3G", "OGDH", "DLD")))
# df_dea1$sector <-fct_relevel(df_dea1$sector, "DLST")

ordered <- gene_groups$V1
prot_corr_eff_ <- reorder_mat(prot_corr_eff, ordered)
zprot_corr_diff2<- reorder_mat(zprot_corr_diff, ordered)
zprot_cor_rest2<- reorder_mat(zprot_cor_rest, ordered)
zprot_cor_MMUT2 <- reorder_mat(zprot_cor_MMUT, ordered)



# zprot_cor_rest2[25,] <- zvCTRL
# zprot_cor_rest2[,25] <- zvCTRL
# row.names(zprot_cor_rest2) <- rownames
#legends
#####
col_fun <- colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("darkblue", "lightblue", "lightgrey", "orange", "darkred"))

lgd_links = Legend(at = c( -0.5, 0, 0.5,1), col_fun = col_fun, 
                   title = "Protein-Protein Relationship", title_position= "topleft", direction = 'horizontal')

lgd_points_spear = Legend(at = c("MMUT deficient", "CTRL"), type = "points", pch = 15, size = unit(4, "mm"), 
                          legend_gp = gpar(col =c("#238B45","#74C476")), title_position = "topleft",  background = "white", # red colors "#9c1111" ,"#e86107"
                          title = "Protein/RNA correlation")

lgd_points_DEA = Legend(at = c("RNA", "Protein"), type = "points", pch = 15, size = unit(4, "mm"),
                        legend_gp = gpar(col = c( "#084594","#4292C6")), title_position = "topleft", background = "white",
                        title = "Effect size (ES)")

lgd_groups = Legend(at = c("is MMUT","TCA", "anaplerotic", "anaplerotic glutamate", "regulatory", "cataplerotic", "pyruvate metabolism", "isoform"), type = "points",
                    legend_gp = gpar(col = c("#ba1116","#ed850e", "#0ddb55", "#10853b",  "#249db3","#3b45d9", "#5d0bb0","#807f7d")), title_position = "topleft",
                    background = "white",
                    title ="Pathway")

lgd_fc = Legend(at = c("RNA", "Protein"), type = "points", pch = 15, size = unit(4, "mm"),
                legend_gp = gpar(col =c("#9302bf","#d078eb")), title_position = "topleft", background = "white",
                title = "Fold change MMUT deficient vs CTRL")

lgd_list_vertical = packLegend(lgd_groups, lgd_fc, lgd_points_DEA,lgd_points_spear,lgd_links)



#### PLOTS


# rest
dev.new()
circos.clear()
chordDiagram(zprot_cor_rest2, grid.col = "grey", symmetric = TRUE, scale = FALSE, col = col_fun, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(zprot_cor_rest2))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  # circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.text( CELL_META$xcenter, CELL_META$cell.ylim[1]+0.33,
               CELL_META$sector.index, facing = "clockwise", niceFacing = T,
               adj = c(0.50, 0), cex = 0.8, col = "black",font =1)
}, bg.border = NA)
title("Protein-protein correlation in control")
circos.clear()
chord_CTRL <- recordPlot()
chord_CTRL_plot <- ggdraw(chord_CTRL)
dev.off()


# MMUT
dev.new()
circos.clear()
chordDiagram(zprot_cor_MMUT2, grid.col = "grey", symmetric = TRUE, scale = FALSE, col = col_fun, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(zprot_cor_MMUT2))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  # circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.text( CELL_META$xcenter, CELL_META$cell.ylim[1]+0.33,
               CELL_META$sector.index, facing = "clockwise", niceFacing = T,
               adj = c(0.50, 0), cex = 0.8, col = "black",font =1)
}, bg.border = NA)
title("Protein-protein correlation in MMUT def.")
circos.clear()
chord_MMUT <- recordPlot()
chord_MMUT_plt <- ggdraw(chord_MMUT)
dev.off()



# difference
dev.new()
circos.clear()
chordDiagram(zprot_corr_diff2, grid.col = "grey", symmetric = TRUE, scale = FALSE, col = col_fun, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(zprot_corr_diff2))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  # circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.text( CELL_META$xcenter, CELL_META$cell.ylim[1]+0.33,
               CELL_META$sector.index, facing = "clockwise", niceFacing = T,
               adj = c(0.50, 0), cex = 0.8, col = "black",font =1)
}, bg.border = NA)
title("Difference in protein-protein correlation")
circos.clear()
chord_diff <- recordPlot()
chord_diff_plot <- ggdraw(chord_diff)
dev.off()



chord_plts <- ggarrange(chord_CTRL_plot, chord_MMUT_plt, chord_diff_plot, nrow = c(1))

ggsave(paste(fig_path_pdf, "SuppFig5/","chordDgrms_iso.pdf", sep = ""), chord_plts, device = "pdf", width = 20, height = 7)
dev.off()





#########################################################
#Chord plots links between all genes (not just from MMUT)


#prot

#ctrl
dev.new()
circos.clear()
chordDiagram(prot_cor_rest, grid.col = "grey", symmetric = TRUE, scale = TRUE, col = col_fun, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(prot_cor_rest))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  # circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.text( CELL_META$xcenter, CELL_META$cell.ylim[1]+0.3,
               CELL_META$sector.index, facing = "clockwise", niceFacing = T,
               adj = c(0.3, 0), cex = 0.8, col = "black",font =1)
}, bg.border = NA)
title("All Protein-protein correlations (control)")
circos.clear()
chord_CTRLprot <- recordPlot()
chord_CTRLprot_plot <- ggdraw(chord_CTRLprot)
dev.off()

# MMUT
dev.new()
circos.clear()
chordDiagram(prot_cor_MMUT, grid.col = "grey", symmetric = TRUE, scale = TRUE, col = col_fun, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(prot_cor_MMUT))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  # circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.text( CELL_META$xcenter, CELL_META$cell.ylim[1]+0.3,
               CELL_META$sector.index, facing = "clockwise", niceFacing = T,
               adj = c(0.3, 0), cex = 0.8, col = "black",font =1)
}, bg.border = NA)
title("All Protein-protein correlations (MMUT def.)")
circos.clear()
chord_MMUTprot <- recordPlot()
chord_MMUTprot_plot <- ggdraw(chord_MMUTprot)
dev.off()

chord_plts_all_prot <- ggarrange(chord_CTRLprot_plot, chord_MMUTprot_plot, nrow = c(1))

ggsave(paste(fig_path_pdf, "Fig4/","chordDgrms_all_prot.pdf", sep = ""), chord_plts_all_prot, device = "pdf", width = 10, height = 5)
dev.off()




#RNA

#settings in RNA plots changed (org. like in prot. plots for plot with just rna chord plots!)

#ctrl
dev.new()
circos.clear()
chordDiagram(rna_corr_rest, grid.col = "grey", symmetric = TRUE, scale = FALSE, col = col_fun, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(rna_corr_rest))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  # circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.text( CELL_META$xcenter, CELL_META$cell.ylim[1]+0.25,
               CELL_META$sector.index, facing = "clockwise", niceFacing = T,
               adj = c(0.2, 0), cex = 0.8, col = "black",font =1)
}, bg.border = NA)
title("All RNA-RNA correlations (control)")
circos.clear()
chord_CTRLrna <- recordPlot()
chord_CTRLrna_plot <- ggdraw(chord_CTRLrna)
dev.off()

# MMUT
dev.new()
circos.clear()
chordDiagram(rna_corr_MMUT, grid.col = "grey", symmetric = TRUE, scale = FALSE, col = col_fun, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(rna_corr_MMUT))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  # circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.text( CELL_META$xcenter, CELL_META$cell.ylim[1]+0.25,
               CELL_META$sector.index, facing = "clockwise", niceFacing = T,
               adj = c(0.2, 0), cex = 0.8, col = "black",font =1)
}, bg.border = NA)
title("All RNA-RNA correlations (MMUT def.)")
circos.clear()
chord_MMUTrna <- recordPlot()
chord_MMUTrna_plot <- ggdraw(chord_MMUTrna)
dev.off()

chord_plts_all_rna <- ggarrange(chord_CTRLrna_plot, chord_MMUTrna_plot, nrow = c(1))

ggsave(paste(fig_path_pdf, "SuppFig5/","chordDgrms_all_rna.pdf", sep = ""), chord_plts_all_rna, device = "pdf", width = 10, height = 5)
dev.off()





