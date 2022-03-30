# Glutamine flux analysis
# Exp3
# measured at Kispi

# libraries
library(tidyr)
library(ggforce)
library(ggplot2)
library(dplyr)
library(scales)
library(ggsci)
library(data.table)
library(forcats)
library(ggpubr)
library(viridis)
library(readxl)
library(patchwork)
library(RColorBrewer)


system("mkdir Figs/flux")
system("mkdir Figs/flux/KispiExp3")

fig_path <- c("Figs/flux/KispiExp3/")
fig_path_pdf <- c("Figs/v12/pdf/")
mypal1 <- pal_aaas("default", alpha = 1)(9)
mypal2 <- pal_npg("nrc", alpha = 1)(9)


### source plotting function
source("Code/GlutamineFlux/plottingIsotop.R")



#### load meta data table

meta_tbl <- data.table(read_excel("Data/GlutamineFlux/Exp3_Kispi/211221_glutamineflux_samplelist_exp3.xlsx"))
meta_tbl[, MMAF := paste0("MMAF", MMAF_no)]
cell_lines <- unique(meta_tbl$cell_line)


#### load data table

tbl <- fread("Data/GlutamineFlux/Experiment3_new/Normalized_Data.csv")

colnames(tbl) <- gsub("\\d{3,}_", "", colnames(tbl))
tbl[, MMAF := regmatches(tbl$Sample_number, regexpr("MMAF\\d+", tbl$Sample_number))]
tbl <- tbl[55:114,] # only keep Exp3
### potentially exclude MMAF46_1 AND MMAF37_3 only in positive mode.
### potentially exclude MMAF22_2 AND MMAF17_1 only in negative mode.
pos <- grep("_pos", colnames(tbl))
neg <- grep("_neg", colnames(tbl))
tbl[Sample_number == "MMAF46_BH_1", (names(tbl)[pos]) := NA]

tbl[Sample_number == "MMAF46_BH_1", (names(tbl)[neg]) := NA]
# tbl[Sample_number == "MMAF46_BH_2", (names(tbl)[neg]) := NA]
# tbl[Sample_number == "MMAF46_BH_3", (names(tbl)[neg]) := NA]

tbl[Sample_number == "MMAF37_BD_3", (names(tbl)[pos]) := NA]
tbl[Sample_number == "MMAF22_AE_2", (names(tbl)[neg]) := NA]
tbl[Sample_number == "MMAF17_AC_1", (names(tbl)[neg]) := NA]

tbl[, protein := meta_tbl[match(tbl$MMAF, meta_tbl$MMAF), ]$Lowry_ugul]
tbl[, cell_line := meta_tbl[match(tbl$MMAF, meta_tbl$MMAF), ]$cell_line]
tbl[, cell_name := meta_tbl[match(tbl$MMAF, meta_tbl$MMAF), ]$cell_name]
tbl$cell_name <- factor(tbl$cell_name, levels = c("WT", "MMUT-KO", "DLST-KO", "OGDH-KO","WT no label","Blank"))
MMAF_names <- unique(tbl$MMAF)


# remove certain genotypes or clones # variable name has to remain tbl because of the customized plotting function later in the script
tbl_noOGDH <- tbl[cell_name != "OGDH-KO", ]
tbl_noOGDHblank <- tbl[cell_name != "OGDH-KO" & cell_name != "Blank", ]
tbl_noOGDHblankDLST <- tbl[cell_name != "OGDH-KO" & cell_name != "Blank" & cell_name != "DLST-KO", ]
tbl_noOGDHcontrols <- tbl[cell_name != "OGDH-KO" & cell_name != "Blank" & cell_name != "WT no label", ]
tbl_MMUTwtonly <- tbl[cell_name == "MMUT-KO" | cell_name == "WT", ]



# stat comparisons
compare <- list(c("WT", "MMUT-KO"), c("WT", "DLST-KO"))
compare2 <- list(c("WT", "MMUT-KO"))




#################################################
# test environment
#################################################


###



#################################################
# Pyruvate
#################################################

cols_gln <- colnames(tbl)[grep("pyru", ignore.case=TRUE, colnames(tbl))][c(1,3:5)]

tbl_noOGDHcontrols_pyr1 <- copy(tbl_noOGDHcontrols[, ..cols_gln])
tbl_noOGDHcontrols_pyr1[, sample_sum := rowSums(.SD, na.rm = FALSE), .SDcols=cols_gln]
incl_row <- tbl_noOGDHcontrols_pyr1$sample_sum>100000
tbl_noOGDHcontrols_pyr2 <- copy(tbl_noOGDHcontrols)
tbl_noOGDHcontrols_pyr2 <- tbl_noOGDHcontrols_pyr2[incl_row, ]

pyr_plt_poolsize <- plottingIsotop(tbl, cols_gln)
pyr_plt_poolsizeGroup <- plottingIsotopGrouped(tbl_noOGDHblankDLST, cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))
pyr_plt_poolsizeGroupTIC <- plottingTICgrouped(tbl_noOGDHcontrols_pyr2, cols_gln) + stat_compare_means(comparisons = compare, method = "wilcox.test", size=3, step.increase=0.2) +scale_y_continuous(expand=c(0.1,0))
pyr_plt_FractionsGroup <- plottingFractionsGrouped(tbl_noOGDHblankDLST, cols_gln)


#################################################
# Alanine
#################################################

cols_gln <- colnames(tbl)[grep("alanin", ignore.case=TRUE, colnames(tbl))][c(1,4:6)]

ala_plt_poolsize <- plottingIsotop(tbl, cols_gln)
ala_plt_poolsizeGroup <- plottingIsotopGrouped(tbl_noOGDHblankDLST, cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))
ala_plt_poolsizeGroupTIC <- plottingTICgrouped(tbl_noOGDHcontrols, cols_gln) + stat_compare_means(comparisons = compare, method = "wilcox.test", size=3, step.increase=0.2) +scale_y_continuous(expand=c(0.1,0))
ala_plt_FractionsGroup <- plottingFractionsGrouped(tbl_noOGDHblankDLST, cols_gln)


#################################################
# Glutamine
#################################################

cols_gln <- colnames(tbl)[grep("Glutamine", colnames(tbl))][1:6]

gln_plt_poolsize <- plottingIsotop(tbl, cols_gln)
gln_plt_poolsizeGroup <- plottingIsotopGrouped(tbl_noOGDHblankDLST, cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))
gln_plt_poolsizeGroupTIC <- plottingTICgrouped(tbl_noOGDHcontrols, cols_gln) + stat_compare_means(comparisons = compare, method = "wilcox.test", size=3, step.increase=0.2) +scale_y_continuous(expand=c(0.1,0))
gln_plt_FractionsGroup <- plottingFractionsGrouped(tbl_noOGDHblankDLST, cols_gln)
gln_ratio_plt <- plottingFractionsRatio(dataframeIsotops=tbl_MMUTwtonly, columns=cols_gln, numeratorIsotop="C5", denominatorIsotop="C4") + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) +scale_y_continuous(expand=c(0.1,0))
gln_labelling<- plottingMetaboliteLabelling(dataframeIsotops=tbl_MMUTwtonly, columns=cols_gln) + stat_compare_means(label.y = 0.9, comparisons = compare2, tip.length = 0.2, method = "wilcox.test", size=3) +scale_y_continuous(limits=c(0,1))
gln_labelling_box <- plottingMetaboliteLabellingBox(dataframeIsotops=tbl_MMUTwtonly, columns=cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))


#################################################
# Glutamic acid
#################################################

tbl_MMUTwtonly_glu <- copy(tbl_MMUTwtonly)
tbl_MMUTwtonly_glu[Sample_number == "MMAF33_BB_1" | Sample_number == "MMAF33_BB_2" | Sample_number == "MMAF33_BB_3", (names(tbl)[pos]) := NA]

cols_gln <- colnames(tbl)[grep("Glutamic", colnames(tbl))][1:6]

glu_plt_poolsize <- plottingIsotop(tbl, cols_gln)
glu_plt_poolsizeGroup <- plottingIsotopGrouped(tbl_noOGDHblankDLST, cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))
glu_plt_poolsizeGroupTIC <- plottingTICgrouped(tbl_noOGDHcontrols, cols_gln) + stat_compare_means(comparisons = compare, method = "wilcox.test", size=3, step.increase=0.2) +scale_y_continuous(expand=c(0.1,0))
glu_plt_FractionsGroup <- plottingFractionsGrouped(tbl_noOGDHblankDLST, cols_gln)
glu_ratio_plt <- plottingFractionsRatio(dataframeIsotops=tbl_MMUTwtonly, columns=cols_gln, numeratorIsotop="C5", denominatorIsotop="C4") + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) +scale_y_continuous(expand=c(0.1,0))
glu_labelling<- plottingMetaboliteLabelling(dataframeIsotops=tbl_MMUTwtonly_glu, columns=cols_gln) + stat_compare_means(label.y = 0.9, comparisons = compare2, tip.length = 0.2, method = "wilcox.test", size=3) +scale_y_continuous(limits=c(0,1))
glu_labelling_box <- plottingMetaboliteLabellingBox(dataframeIsotops=tbl_MMUTwtonly_glu, columns=cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))


#################################################
# Oxoglutarate
#################################################

tbl_noOGDHcontrols_oxo <- copy(tbl_noOGDHcontrols)
tbl_noOGDHcontrols_oxo[Sample_number == "MMAF45_BH_2", (names(tbl)[neg]) := NA]
tbl_MMUTwtonly_oxo <- copy(tbl_MMUTwtonly)
tbl_MMUTwtonly_oxo[Sample_number == "MMAF46_BH_2" | Sample_number == "MMAF43_BG_3", (names(tbl)[neg]) := NA]
tbl_MMUTwtonly_oxo2 <- copy(tbl_MMUTwtonly)
tbl_MMUTwtonly_oxo2[Sample_number == "MMAF33_BB_1" | Sample_number == "MMAF33_BB_2" | Sample_number == "", (names(tbl)[neg]) := NA]

cols_gln <- colnames(tbl)[grep("oxogluta", colnames(tbl), ignore.case=TRUE)][c(1,3:7)]

oxo_plt_poolsize <- plottingIsotop(tbl, cols_gln)
oxo_plt_poolsizeGroup <- plottingIsotopGrouped(tbl_noOGDHblankDLST, cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))
oxo_plt_poolsizeGroupTIC <- plottingTICgrouped(tbl_noOGDHcontrols_oxo, cols_gln) + stat_compare_means(comparisons = compare, method = "wilcox.test", size=3, step.increase=0.15) +scale_y_continuous(expand=c(0.1,0))
oxo_plt_FractionsGroup <- plottingFractionsGrouped(tbl_noOGDHblankDLST, cols_gln)
oxo_ratio_plt <- plottingFractionsRatio(dataframeIsotops=tbl_MMUTwtonly_oxo, columns=cols_gln, numeratorIsotop="C5", denominatorIsotop="C4") + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) +scale_y_continuous(expand=c(0.1,0))
oxo_labelling<- plottingMetaboliteLabelling(dataframeIsotops=tbl_MMUTwtonly_oxo2, columns=cols_gln) + stat_compare_means(label.y = 0.9, comparisons = compare2, tip.length = 0.2, method = "wilcox.test", size=3) +scale_y_continuous(limits=c(0,1))
oxo_labelling_box <- plottingMetaboliteLabellingBox(dataframeIsotops=tbl_MMUTwtonly_oxo2, columns=cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))


#################################################
# Citrate
#################################################

tbl_MMUTwtonly_cit1 <- copy(tbl_MMUTwtonly)
tbl_MMUTwtonly_cit1[Sample_number == "MMAF31_BA_2" |Sample_number == "MMAF33_BB_1" | Sample_number == "MMAF33_BB_2" | Sample_number == "MMAF33_BB_3"| Sample_number =="MMAF34_BB_1", (names(tbl)[neg]) := NA]
tbl_MMUTwtonly_cit2 <- copy(tbl_MMUTwtonly_cit1)
tbl_MMUTwtonly_cit2[Sample_number == "MMAF46_BH_2" |Sample_number == "MMAF46_BH_3" | Sample_number == "MMAF34_BB_2" | Sample_number == "MMAF34_BB_3", (names(tbl)[neg]) := NA]
tbl_noOGDHcontrols_cit <- copy(tbl_noOGDHcontrols)
tbl_noOGDHcontrols_cit[Sample_number == "MMAF31_BA_2" |Sample_number == "MMAF33_BB_1" | Sample_number == "MMAF33_BB_2" | Sample_number == "MMAF33_BB_3"| Sample_number =="MMAF34_BB_1" | Sample_number == "MMAF46_BH_2" |Sample_number == "MMAF46_BH_3" | Sample_number =="MMAF40_BE_1", (names(tbl)[neg]) := NA]

cols_gln <- colnames(tbl)[grep("Citr", colnames(tbl), ignore.case=FALSE)][c(1,3:8)]

cit_plt_poolsize <- plottingIsotop(tbl, cols_gln)
cit_plt_poolsizeGroup <- plottingIsotopGrouped(tbl_noOGDHblankDLST, cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))
cit_plt_poolsizeGroupTIC <- plottingTICgrouped(tbl_noOGDHcontrols, cols_gln) + stat_compare_means(comparisons = compare, method = "wilcox.test", size=3, step.increase=0.15) +scale_y_continuous(expand=c(0.1,0))
cit_plt_poolsizeGroupTIC2 <- plottingTICgrouped(tbl_MMUTwtonly_cit1, cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3, step.increase=0.15) +scale_y_continuous(expand=c(0.1,0))
cit_plt_FractionsGroup <- plottingFractionsGrouped(tbl_noOGDHblankDLST, cols_gln)
cit_ratio_pltC4C0 <- plottingFractionsRatio(dataframeIsotops=tbl_MMUTwtonly_cit1, columns=cols_gln, numeratorIsotop="C4", denominatorIsotop="C0") + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) +scale_y_continuous(expand=c(0.1,0))
cit_ratio_pltC5C4 <- plottingFractionsRatio(dataframeIsotops=tbl_MMUTwtonly_cit2, columns=cols_gln, numeratorIsotop="C5", denominatorIsotop="C4") + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) +scale_y_continuous(expand=c(0.1,0))
cit_ratio_pltC5C4_wDLST <- plottingFractionsRatio(dataframeIsotops=tbl_noOGDHcontrols_cit, columns=cols_gln, numeratorIsotop="C5", denominatorIsotop="C4") + stat_compare_means(comparisons = compare, method = "wilcox.test", size=3) +scale_y_continuous(expand=c(0.1,0))
cit_labelling<- plottingMetaboliteLabelling(dataframeIsotops=tbl_MMUTwtonly_cit1, columns=cols_gln) + stat_compare_means(label.y = 0.9, comparisons = compare2, tip.length = 0.2, method = "wilcox.test", size=3) +scale_y_continuous(limits=c(0,1))
cit_labelling_box <- plottingMetaboliteLabellingBox(dataframeIsotops=tbl_MMUTwtonly_cit1, columns=cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))


#################################################
# Aspartate
#################################################

tbl_noOGDHcontrols_asp <- copy(tbl_noOGDHcontrols)
tbl_noOGDHcontrols_asp[Sample_number == "MMAF46_BH_2" | Sample_number == "MMAF46_BH_3", (names(tbl)[neg]) := NA]
tbl_MMUTwtonly_asp <- copy(tbl_MMUTwtonly)
tbl_MMUTwtonly_asp[Sample_number == "MMAF46_BH_2" | Sample_number == "MMAF46_BH_3" | Sample_number == "MMAF33_BB_1" | Sample_number == "MMAF33_BB_2" | Sample_number == "MMAF33_BB_3", (names(tbl)[neg]) := NA]

cols_gln <- colnames(tbl)[grep("Aspar", colnames(tbl), ignore.case=FALSE)][c(2:6)]

asp_plt_poolsize <- plottingIsotop(tbl, cols_gln)
asp_plt_poolsizeGroup <- plottingIsotopGrouped(tbl_noOGDHblankDLST, cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))
asp_plt_poolsizeGroupTIC <- plottingTICgrouped(tbl_noOGDHcontrols_asp, cols_gln) + stat_compare_means(comparisons = compare, method = "wilcox.test", size=3, step.increase=0.15) +scale_y_continuous(expand=c(0.1,0))
asp_plt_FractionsGroup <- plottingFractionsGrouped(tbl_noOGDHblankDLST, cols_gln)
asp_ratio_plt <- plottingFractionsRatio(dataframeIsotops=tbl_MMUTwtonly_asp, columns=cols_gln, numeratorIsotop="C4", denominatorIsotop="C0") + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) +scale_y_continuous(expand=c(0.1,0))
asp_labelling<- plottingMetaboliteLabelling(dataframeIsotops=tbl_MMUTwtonly_asp, columns=cols_gln) + stat_compare_means(label.y = 0.9, comparisons = compare2, tip.length = 0.2, method = "wilcox.test", size=3) +scale_y_continuous(limits=c(0,1))
asp_labelling_box <- plottingMetaboliteLabellingBox(dataframeIsotops=tbl_MMUTwtonly_asp, columns=cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))


#################################################
# Fumarate
#################################################

tbl_noOGDHcontrols_fum <- copy(tbl_noOGDHcontrols)
tbl_noOGDHcontrols_fum[Sample_number == "MMAF46_BH_3" | Sample_number == "MMAF46_BH_2", (names(tbl)[neg]) := NA]
tbl_MMUTwtonly_fum <- copy(tbl_MMUTwtonly)
tbl_MMUTwtonly_fum[Sample_number == "MMAF35_BC_2", (names(tbl)[neg]) := NA]
tbl_MMUTwtonly_fum2 <- copy(tbl_MMUTwtonly)
tbl_MMUTwtonly_fum2[Sample_number == "MMAF46_BH_2" | Sample_number == "MMAF46_BH_3"| Sample_number == "MMAF33_BB_1" | Sample_number == "MMAF33_BB_2" | Sample_number == "MMAF35_BC_2", (names(tbl)[neg]) := NA]

cols_gln <- colnames(tbl)[grep("Fumaric", colnames(tbl), ignore.case=FALSE)][c(1:5)]

fum_plt_poolsize <- plottingIsotop(tbl, cols_gln)
fum_plt_poolsizeGroup <- plottingIsotopGrouped(tbl_noOGDHblankDLST, cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))
fum_plt_poolsizeGroupTIC <- plottingTICgrouped(tbl_noOGDHcontrols_fum, cols_gln) + stat_compare_means(comparisons = compare, method = "wilcox.test", size=3, step.increase=0.15) +scale_y_continuous(expand=c(0.1,0))
fum_plt_FractionsGroup <- plottingFractionsGrouped(tbl_noOGDHblankDLST, cols_gln)
fum_ratio_plt <- plottingFractionsRatio(dataframeIsotops=tbl_MMUTwtonly_fum, columns=cols_gln, numeratorIsotop="C4", denominatorIsotop="C0") + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) +scale_y_continuous(expand=c(0.1,0))
fum_labelling<- plottingMetaboliteLabelling(dataframeIsotops=tbl_MMUTwtonly_fum2, columns=cols_gln) + stat_compare_means(label.y = 0.9, comparisons = compare2, tip.length = 0.2, method = "wilcox.test", size=3) +scale_y_continuous(limits=c(0,1))
fum_labelling_box <- plottingMetaboliteLabellingBox(dataframeIsotops=tbl_MMUTwtonly_fum2, columns=cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))


#################################################
# Malate
#################################################

tbl_MMUTwtonly_mal <- copy(tbl_MMUTwtonly)
tbl_MMUTwtonly_mal[Sample_number == "MMAF46_BH_2" | Sample_number == "MMAF46_BH_3" | Sample_number == "MMAF33_BB_1"| Sample_number == "MMAF33_BB_2"| Sample_number == "MMAF33_BB_3", (names(tbl)[neg]) := NA]

cols_gln <- colnames(tbl)[grep("Mal", colnames(tbl), ignore.case=FALSE)][c(2,4:7)]

mal_plt_poolsize <- plottingIsotop(tbl, cols_gln)
mal_plt_poolsizeGroup <- plottingIsotopGrouped(tbl_noOGDHblankDLST, cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))
mal_plt_poolsizeGroupTIC <- plottingTICgrouped(tbl_noOGDHcontrols, cols_gln) + stat_compare_means(comparisons = compare, method = "wilcox.test", size=3, step.increase=0.15) +scale_y_continuous(expand=c(0.1,0))
mal_plt_FractionsGroup <- plottingFractionsGrouped(tbl_noOGDHblankDLST, cols_gln)
mal_ratio_plt <- plottingFractionsRatio(dataframeIsotops=tbl_MMUTwtonly_mal, columns=cols_gln, numeratorIsotop="C4", denominatorIsotop="C0") + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) +scale_y_continuous(expand=c(0.1,0))
mal_labelling<- plottingMetaboliteLabelling(dataframeIsotops=tbl_MMUTwtonly_mal, columns=cols_gln) + stat_compare_means(label.y = 0.9, comparisons = compare2, tip.length = 0.2, method = "wilcox.test", size=3) +scale_y_continuous(limits=c(0,1))
mal_labelling_box <- plottingMetaboliteLabellingBox(dataframeIsotops=tbl_MMUTwtonly_mal, columns=cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))


#################################################
# Succinate
#################################################

tbl_MMUTwtonly_suc <- copy(tbl_MMUTwtonly)
tbl_MMUTwtonly_suc[Sample_number == "MMAF33_BB_1"| Sample_number == "MMAF33_BB_2"| Sample_number == "MMAF33_BB_3", (names(tbl)[neg]) := NA]

cols_gln <- colnames(tbl)[grep("Succi", colnames(tbl), ignore.case=FALSE)][c(1,3:6)]

suc_plt_poolsize <- plottingIsotop(tbl, cols_gln)
suc_plt_poolsizeGroup <- plottingIsotopGrouped(tbl_noOGDHblankDLST, cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))
suc_plt_poolsizeGroupTIC <- plottingTICgrouped(tbl_noOGDHcontrols, cols_gln) + stat_compare_means(comparisons = compare, method = "wilcox.test", size=3, step.increase=0.15) +scale_y_continuous(expand=c(0.1,0))
suc_plt_FractionsGroup <- plottingFractionsGrouped(tbl_noOGDHblankDLST, cols_gln)
suc_ratio_plt <- plottingFractionsRatio(dataframeIsotops=tbl_MMUTwtonly, columns=cols_gln, numeratorIsotop="C4", denominatorIsotop="C0") + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) +scale_y_continuous(expand=c(0.1,0))
suc_labelling<- plottingMetaboliteLabelling(dataframeIsotops=tbl_MMUTwtonly_suc, columns=cols_gln) + stat_compare_means(label.y = 0.9, comparisons = compare2, tip.length = 0.2, method = "wilcox.test", size=3) +scale_y_continuous(limits=c(0,1))
suc_labelling_box <- plottingMetaboliteLabellingBox(dataframeIsotops=tbl_MMUTwtonly_suc, columns=cols_gln) + stat_compare_means(comparisons = compare2, method = "wilcox.test", size=3) + scale_y_continuous(expand=c(0.1,0))






#################################################
#################################################
# PLOT compilation - Fractions, TICs, Ratios
#################################################
#################################################



supp_TIC <- 
gln_plt_poolsizeGroupTIC + glu_plt_poolsizeGroupTIC + oxo_plt_poolsizeGroupTIC + pyr_plt_poolsizeGroupTIC + ala_plt_poolsizeGroupTIC +plot_layout(guides = 'collect', nrow=1) & theme(legend.position="bottom")
ggsave(supp_TIC, filename=paste0(fig_path,"Supp_TIC_GlnGluOxo.pdf"), device = "pdf", width = 8, height = 3)

fig_TIC <- 
cit_plt_poolsizeGroupTIC + fum_plt_poolsizeGroupTIC + mal_plt_poolsizeGroupTIC + suc_plt_poolsizeGroupTIC + asp_plt_poolsizeGroupTIC + plot_layout(guides = 'collect', nrow = 1) & theme(legend.position="bottom")
ggsave(fig_TIC, filename=paste0(fig_path,"Fig_TIC_CitAspFumMalSuc.pdf"), device = "pdf", width = 8, height = 3)



supp_fractionsGroup <- 
glu_plt_FractionsGroup + oxo_plt_FractionsGroup + plot_layout(guides = 'collect') & theme(legend.position="right")
ggsave(supp_fractionsGroup, filename=paste0(fig_path,"Supp_Fractions_GluOxo.pdf"), device = "pdf", width = 4, height = 2.5)

fig_fractionsGroup <- 
cit_plt_FractionsGroup + fum_plt_FractionsGroup + mal_plt_FractionsGroup + suc_plt_FractionsGroup + asp_plt_FractionsGroup + plot_layout(guides = 'collect', nrow = 1) & theme(legend.position="bottom") & guides(fill = guide_legend(nrow=1))
ggsave(fig_fractionsGroup, filename=paste0(fig_path,"Fig_Fractions_CitAspFumMalSuc.pdf"), device = "pdf", width = 8, height = 3)



supp_ratios <- 
cit_ratio_pltC5C4 + plot_layout(guides = 'collect', nrow = 1) & theme(legend.position="bottom")
ggsave(supp_ratios, filename=paste0(fig_path,"Supp_Ratios.pdf"), device = "pdf", width = 2, height = 3)

fig_ratios <- 
cit_ratio_pltC4C0 + fum_ratio_plt+ mal_ratio_plt + suc_ratio_plt + asp_ratio_plt + plot_layout(guides = 'collect', nrow = 1) & theme(legend.position="bottom")
ggsave(fig_ratios, filename=paste0(fig_path,"Fig_Ratios_CitAspFumMalSuc.pdf"), device = "pdf", width = 8, height = 3)





fig_labelling <- 
cit_labelling + fum_labelling+ mal_labelling + suc_labelling + asp_labelling + plot_layout(guides = 'collect', nrow = 1) & theme(legend.position="bottom")
ggsave(fig_labelling, filename=paste0(fig_path,"Fig_labelling.pdf"), device = "pdf", width = 8, height = 3)


fig_labelling_box <- 
cit_labelling_box + fum_labelling_box+ mal_labelling_box + suc_labelling_box + asp_labelling_box + plot_layout(guides = 'collect', nrow = 1) & theme(legend.position="bottom")
ggsave(fig_labelling_box, filename=paste0(fig_path,"Fig_labellingBox.pdf"), device = "pdf", width = 8, height = 3)



### citrate only plot

cit_plt_collection <- 
(cit_plt_poolsizeGroupTIC2 & theme(legend.position="none")) + cit_plt_FractionsGroup + cit_labelling_box + cit_ratio_pltC5C4 + plot_layout(nrow = 1)
ggsave(cit_plt_collection, filename=paste0(fig_path,"Citrate_plt_collection.pdf"), device = "pdf", width = 7, height = 3)



