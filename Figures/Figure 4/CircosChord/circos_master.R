#### Circos plot of TCA genes

#libraries
#####
library(data.table)
library(dplyr)
library(tidyverse)
library(corrplot)
library(ggpubr)
library(ggplot2)
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
load("interimData/rnafiltered_mat.RData")
load("interimData/prot_mat.RData")

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
diffexp_prot <- read_table2("interimData/diff_exp_data_prot_pathwayact_all.txt", col_types = cols(betas = col_number()))
diffexp_rna <- read_table2("interimData/diff_exp_data_rnaseq_pathwayact_all.txt",col_types = cols(betas = col_number()))

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
prot_subset <- prot_exp[is.element(prot_exp$ENSG, tca_tbl), ] #only 37 genes of the 62 contained tca_tbl (based on masterTCA_keggPlusmanual.csv) are present in prot_exp data --> results in our list of genes for circos plot
prot_subset_t <- dcast(melt(prot_subset), variable ~ ENSG)
prot_subset_t[, variable := NULL]

new_names_prot <- hgnc_p[which(hgnc_p$ENSG %in% colnames(prot_subset_t)), ]$HGNC
setnames(prot_subset_t, new_names_prot)

prot_cor_all <- cor(prot_subset_t, use = "complete.obs", method = "spearman")
prot_cor_MMUT <- cor(prot_subset_t[1:150, ],  use = "complete.obs", method = "spearman")
prot_cor_rest <- cor(prot_subset_t[151:230, ],  use = "complete.obs", method = "spearman")

rownames <- row.names(prot_cor_MMUT)
# rna_corr_all <- cor(rna_subset_t, method = "spearman")
# rna_corr_MMUT <- cor(rna_subset_t[1:143, ], method = "spearman")
# rna_corr_rest <- cor(rna_subset_t[144:nrow(rna_subset_t), ], method = "spearman")

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


###### data diff. expression lvls rna and protein (most outer track)

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
# 
logav_exp_MMUT10<-df %>%
  group_by(sectors) %>%
  filter(type =="MMUT") %>%
  summarise(logav_rna_mmut = mean((x)), logav_prot_mmut = mean((y)))

# 
logav_exp_REST10 <- df%>%
  group_by(sectors) %>%
  filter(type =="REST") %>%
  summarise(logav_rna_rest = mean(x), logav_prot_rest = mean((y)))



# logav_exp <- left_join(logav_exp_MMUT, logav_exp_REST, by = "sectors")
logav_exp <- left_join(logav_exp_MMUT10, logav_exp_REST10, by = "sectors")

#CHANGE HERE TO CORRECT FOLD CHANGE FORMULA
fc_exp<-  logav_exp %>% 
  mutate(fc_rna = (logav_rna_mmut- logav_rna_rest)/ logav_rna_rest )%>% 
  # mutate(fc_rna = log2(fc_rna)) %>% 
  mutate(fc_prot =(logav_prot_mmut - logav_prot_rest)/ logav_prot_rest) %>% 
  # mutate(fc_prot = log2(fc_prot)) %>% 
  as.matrix()

fc_exp_ <- data.frame(fc_exp[,-1], row.names=fc_exp[,1]) 
fc_exp_$fc_rna <- as.numeric(fc_exp_$fc_rna)
fc_exp_$fc_prot <- as.numeric(fc_exp_$fc_prot)
m_fc_exp<-data.matrix(fc_exp_[,5:6])



df_fc_exp <- data.frame(sector = as.factor(c(rep(fc_exp[,1], 2))),
                        x = c(rep(1, 37), rep(2, 37)),
                        y = c(m_fc_exp[, "fc_rna"], m_fc_exp[, "fc_prot"]))

upper_boarder <- df_corr1$sector %>% unique() %>% length()
upper_boarder2 <- gene_groups$V1 %>% unique() %>% length()
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


#legends
#####
col_fun <- colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("blue", "lightblue", "#e3e3e3", "yellow","red"))

lgd_links = Legend(at = c( -0.5, 0, 0.5,1), col_fun = col_fun,
                   title = "Protein-Protein Relationship", title_position= "topleft", direction = 'horizontal')

lgd_points_spear = Legend(at = c( "CTRL","MMUT deficient"), type = "points", pch = 15, size = unit(4, "mm"), 
                          legend_gp = gpar(col =c("#9c1111" ,"#e86107")), title_position = "topleft",  background = "white", # green colors "#238B45","#74C476"
                          title = "Protein/RNA correlation")

lgd_points_DEA = Legend(at = c("RNA", "Protein"), type = "points", pch = 15, size = unit(4, "mm"),
                        legend_gp = gpar(col = c( "#084594","#4292C6")), title_position = "topleft", background = "white",
                        title = "Effect size (ES)")

# lgd_groups = Legend(at = c("is MMUT","TCA", "anaplerotic", "anaplerotic glutamate", "regulatory", "cataplerotic", "pyruvate metabolism", "isoform"), type = "points",
# legend_gp = gpar(col = c("#ba1116","#ed850e", "#0ddb55", "#10853b",  "#249db3","#3b45d9", "#5d0bb0","#807f7d")), title_position = "topleft",
# background = "white",
# title ="Pathway")

lgd_fc = Legend(at = c("RNA", "Protein"), type = "points", pch = 15, size = unit(4, "mm"),
                legend_gp = gpar(col =c("#9302bf","#d078eb")), title_position = "topleft", background = "white",
                title = "Fold change MMUT deficient vs CTRL")

lgd_list_vertical = packLegend(lgd_fc, lgd_points_DEA,lgd_points_spear,lgd_links)
# 
t_col <- function(color, percent = 65, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  ## Save the color
  invisible(t.col)
}

torange <- t_col("#ed850e")
tlightgreen <- t_col("#0ddb55")
tred <- t_col("#ba1116")
tgrey <- t_col("#807f7d")
tgreen <- t_col("#10853b")
tteal <- t_col("#249db3")
tblue <- t_col("#3b45d9")
tpurple <- t_col("#5d0bb0")
###############################################################################################################
#PLOT
###############################################################################################################

circos.clear()
dev.new()
circos.par(start.degree=84, gap.after = c(rep(2, 36), 12), track.height=0.18)
chordDiagram(zprot_corr_diff2, grid.col = 'grey', symmetric = TRUE, scale = TRUE, col = col_fun, annotationTrack = "grid", annotationTrackHeight = 0.025, preAllocateTracks = list(list(),list(), list())) # list determines number of tracks > now with 4 tracks (most inner not visible)
circle_size = unit(1, "snpc") # snpc unit gives you a square region
circos.clear()
circos.par( gap.after = c(rep(2, 36), 12), track.height=0.18)
pushViewport(viewport(x = 0, y = 0.5, width = circle_size+30, height = circle_size+30,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circos.par("points.overflow.warning" = FALSE, gap.after = gap_v) #,gap.degree=3
circos.initialize(sectors = df_corr1$sector, xlim =c(0,2)) #Dataframe statt matrix

# # circos.track(sectors=df_fc_exp$sector, y=1)
# circos.track(sectors = df_fc_exp$sector, x = df_fc_exp$x, y = df_fc_exp$y, panel.fun = function(x, y) { # inner track
#   xlim = CELL_META$xlim
#   ylim = CELL_META$ylim
#   value = y
#   circos.barplot(value, 1:2 - 0.5, col = c("#9302bf","#d078eb"), border = "transparent")
#   circos.axis(h = 0, labels = F, major.tick = F, col = "grey")
#   if(CELL_META$sector.numeric.index == 29) { # on first sector after gap on left side axis
#     circos.yaxis(side = "left", tick.length=0.1,labels.cex = 0.5, labels.niceFacing=T)}
# })

# outer track fold change
circos.track(sectors = df_fc_exp$sector, x = df_fc_exp$x, y = df_fc_exp$y, panel.fun = function(x, y) { # inner track
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  value = y
  circos.barplot(value, 1:2 - 0.5, col = c("#9302bf","#d078eb"), border = "transparent")
  circos.axis(h = 0, labels = F, major.tick = F, col = "grey")
  if(CELL_META$sector.numeric.index == 29) { # on first sector after gap on left side axis
    circos.yaxis(side = "left", tick.length=0.1,labels.cex = 0.5, labels.niceFacing=T)}
})

# middle track
circos.track(sectors = df_dea1$sector, x = df_dea1$x, y = df_dea1$y2, panel.fun = function(x, y) { 
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  value = y
  circos.barplot(value, 1:2 - 0.5, col = c("#084594","#4292C6"), border = "white")
  circos.axis(h = 0, labels = F, major.tick = F, col = "grey")
  if(CELL_META$sector.numeric.index == 29) { # on first sector after gap on left side
    circos.yaxis(side = "left",  tick.length=0.1, labels.cex = 0.5, labels.niceFacing=T)}
})

circos.track(sectors = df_corr1$sector, x = df_corr1$x, y = df_corr1$y, panel.fun = function(x, y) { # inner track
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  value = y
  circos.barplot(value, 1:2 - 0.5, col = c("#9c1111" ,"#e86107"), border = "transparent") # red colors "#9c1111" ,"#e86107" 
  circos.axis(h = 0, labels = F, major.tick = F, col = "grey")
  if(CELL_META$sector.numeric.index == 29) { # on first sector after gap on left side axis
    circos.yaxis(side = "left", tick.length=0.1,labels.cex = 0.5, labels.niceFacing=T)}
  # circos.text( CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(48),
  #              CELL_META$sector.index, facing = "bending.inside", niceFacing = T,
  #              adj = c(0, -0.5), cex = 0.75)
  circos.text( CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(45.7),
               CELL_META$sector.index, facing = "bending.inside", niceFacing = T,
               adj = c(0.50, 0), cex = 0.71, col = "grey30",font =2)
  
})

# labels in gap
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 28) { # gap sector
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 0.8,
                CELL_META$cell.xlim[2] + convert_x(7, "mm"), 1.95,
                col = "transparent", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(3, "mm"), 2.5,
                "FC", cex = 0.75, facing = "inside")
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -0.45, #l?nge Farbblock innen
                CELL_META$cell.xlim[2] + convert_x(7, "mm"), 0.73,
                # CELL_META$cell.xlim[2] + convert_x(13, "mm"), 0.9,
                col = "transparent", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(3, "mm"), 0.14,
                expression(rho), cex = 0.86, facing = "inside")
    # circos.text(CELL_META$cell.xlim[2] + convert_x(4.5, "mm"), 0.08,
    #             "", cex = 0.74, facing = "clockwise")
    # 
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 0.8,
                CELL_META$cell.xlim[2] + convert_x(7, "mm"), 1.95,
                col = "transparent", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(3, "mm"), 1.38,
                "ES", cex = 0.8, facing = "inside")
    # circos.text(CELL_META$cell.xlim[2] + convert_x(4.5, "mm"), 1.38,
    #             "", cex = 0.55, facing = "clockwise")
  }
}, bg.border = NA)


draw.sector(get.cell.meta.data("cell.start.degree", sector.index = head(gene_groups_orderd$V1[gene_groups_orderd$V2 == "TCA"], 1)),
            get.cell.meta.data("cell.end.degree", sector.index = tail(gene_groups_orderd$V1[gene_groups_orderd$V2 == "TCA"], 1)),
            rou1 = 1.054, #get.cell.meta.data("cell.top.radius", track.index = 2),
            rou2 = 1.045, #get.cell.meta.data("cell.bottom.radius", track.index = 3),
            col = torange,
            border = F)

draw.sector(get.cell.meta.data("cell.start.degree", sector.index = head(gene_groups$V1[gene_groups$V2 == "anaplerotic"], 1)),
            get.cell.meta.data("cell.end.degree", sector.index = tail(gene_groups$V1[gene_groups$V2 == "anaplerotic"], 1)),
            rou1 = 1.054, #get.cell.meta.data("cell.top.radius", track.index = 2),
            rou2 = 1.045, #get.cell.meta.data("cell.bottom.radius", track.index = 3),
            col = tlightgreen, #"#0ddb55"
            border = F)

draw.sector(get.cell.meta.data("cell.start.degree", sector.index = head(gene_groups$V1[gene_groups$V2 == "anaplerotic glutamate"], 1)),
            get.cell.meta.data("cell.end.degree", sector.index = tail(gene_groups$V1[gene_groups$V2 == "anaplerotic glutamate"], 1)),
            rou1 = 1.054, #get.cell.meta.data("cell.top.radius", track.index = 2),
            rou2 = 1.045, #get.cell.meta.data("cell.bottom.radius", track.index = 3),
            col = tgreen,
            border = F)

draw.sector(get.cell.meta.data("cell.start.degree", sector.index = head(gene_groups$V1[gene_groups$V2 == "regulatory"], 1)),
            get.cell.meta.data("cell.end.degree", sector.index = tail(gene_groups$V1[gene_groups$V2 == "regulatory"], 1)),
            rou1 = 1.054, #get.cell.meta.data("cell.top.radius", track.index = 2),
            rou2 = 1.045, #get.cell.meta.data("cell.bottom.radius", track.index = 3),
            col = tteal,
            border = F)

draw.sector(get.cell.meta.data("cell.start.degree", sector.index = head(gene_groups$V1[gene_groups$V2 == "cataplerotic"], 1)),
            get.cell.meta.data("cell.end.degree", sector.index = tail(gene_groups$V1[gene_groups$V2 == "cataplerotic"], 1)),
            rou1 = 1.054, #get.cell.meta.data("cell.top.radius", track.index = 2),
            rou2 = 1.045, #get.cell.meta.data("cell.bottom.radius", track.index = 3),
            col = tblue,
            border = F)

draw.sector(get.cell.meta.data("cell.start.degree", sector.index = head(gene_groups$V1[gene_groups$V2 == "pyruvate metabolism"], 1)),
            get.cell.meta.data("cell.end.degree", sector.index = tail(gene_groups$V1[gene_groups$V2 == "pyruvate metabolism"], 1)),
            rou1 = 1.054, #get.cell.meta.data("cell.top.radius", track.index = 2),
            rou2 = 1.045, #get.cell.meta.data("cell.bottom.radius", track.index = 3),
            col = tpurple,
            border = F)

# draw.sector(get.cell.meta.data("cell.start.degree", sector.index = head(gene_groups$V1[gene_groups$V2 == "is MMUT"], 1)),
#             get.cell.meta.data("cell.end.degree", sector.index = tail(gene_groups$V1[gene_groups$V2 == "is MMUT"], 1)),
#             rou1 = 1.054, #get.cell.meta.data("cell.top.radius", track.index = 2),
#             rou2 = 1.045, #get.cell.meta.data("cell.bottom.radius", track.index = 3),
#             col = tred,
#             border = F)

draw.sector(get.cell.meta.data("cell.start.degree", sector.index = head(gene_groups_orderd$V1[gene_groups_orderd$V2 == "Isoform"], 1)),
            get.cell.meta.data("cell.end.degree", sector.index = tail(gene_groups_orderd$V1[gene_groups_orderd$V2 == "Isoform"], 1)),
            rou1 = 1.054, #get.cell.meta.data("cell.top.radius", track.index = 2),
            rou2 = 1.045, #get.cell.meta.data("cell.bottom.radius", track.index = 3),
            col = tgrey,
            border = F)


#labels of groups
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { 
    circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(50),
                "TCA cycle", facing = "bending.inside", niceFacing = F,
                adj = c(0.50, 0), cex = 0.75, col = "#cf6e00",font =2)
  }})

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 13) { 
    circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(49.5),
                "anaplerotic", facing = "bending.inside", niceFacing = T,
                adj = c(0.50, 0), cex = 0.75, col = "#07a33d",font =2)
  }})

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 19) {
    circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(49.5),
                "anaplerotic glutamate", facing = "bending.inside", niceFacing = T,
                adj = c(0.35, 0), cex = 0.75, col = "#10853b",font =2)
  }})

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 22) {
    circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(50),
                "catapl.", facing = "bending.inside", niceFacing = T,
                adj = c(0.50, 0), cex = 0.75, col = "#3b45d9",font =2)
  }})

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 23) {
    circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(50),
                "Pyr metabolism.", facing = "bending.inside", niceFacing = T,
                adj = c(0.23, 0), cex = 0.75, col = "#5d0bb0",font =2)
  }})

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 26) {
    circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(50),
                "Isoform", facing = "bending.inside", niceFacing = T,
                adj = c(0.23, 0), cex = 0.75, col = "#807f7d" ,font =2)
  }})

# circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
#   if(CELL_META$sector.numeric.index == 29) {
#     circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(50),
#                 "MMUT", facing = "bending.inside", niceFacing = T,
#                 adj = c(0.5, 0), cex = 0.75, col = "#ba1116" ,font =2)
#   }})
par(xpd = NA)
# pushViewport(viewport(x = unit(0.15, "npc"), y = unit(-0.4, "npc"), just = c("left", "bottom")))

pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size,
                      just = c("center", "top")))
par(omi = gridOMI(), new = TRUE)
# pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size+1,
#                       just = c("center", "top")))
# draw(lgd_list_horizontal1, x= unit(0.7, "npc") - circle_size, just = c('right',"top"))
# draw(lgd_list_horizontal2, x= unit(1.05, "npc") - circle_size, just = c('left',"top"))
draw(lgd_list_horizontal1,x= unit(0.5, "npc"),  y = unit(0.65, "npc") - circle_size, just = c('center',"top"))
# draw(lgd_list_horizontal2, x= unit(0.32, "npc"), y = unit(0.52, "npc") - circle_size, just = c('center',"top"))
# draw(lgd_list_horizontal3, x= unit(0.23, "npc"), y = unit(0.445, "npc") - circle_size, just = c('center',"top"))
draw(lgd_list_horizontal2, x= unit(0.33, "npc"), y = unit(0.52, "npc") - circle_size, just = c('center',"top"))
draw(lgd_list_horizontal3, x= unit(0.24, "npc"), y = unit(0.445, "npc") - circle_size, just = c('center',"top"))


popViewport()

circos.clear()

dev.print(pdf, 'complete_circos_plot_Legend_horizontal.pdf')
