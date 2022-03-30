# specific plots for MMUT and how biochemitry correaltes with protein and transcript levels


library(data.table)
# library(readr)
# library(readxl)
# library(reshape2)
library(ggpubr)
library(extrafont)
library(gridExtra)
library(gtools)
library(cowplot)
# library(EnhancedVolcano)
# library(dplyr)
# library(plyr)
# library(biomaRt)
library(tidyverse)
library(ggridges)
library(ggsci)
library(treemapify)
library(ggrepel)



system("mkdir Figs/MMUT")

fig_path <- c("Figs/MMUT/")
fig_path_pdf <- c("Figs/v12/pdf/")
mypal = c(pal_npg("nrc", alpha = 1)(9))
mypal1 <- pal_aaas("default", alpha = 1)(9)
mypal2 <- c(mypal1[2], mypal1[4], mypal1[1])
mypal3 <- c(mypal1[2], mypal1[6], mypal1[4], mypal1[1])


# MMUT Ensembl: ENSG00000146085

##################################
# load protein matrix
load("interimData/prot_mat.RData")

prot_exp <- data.table(ENSG = row.names(prot_mat), prot_mat)
prot_exp2 <- melt(prot_exp, id = c("ENSG"))
prot_exp2 <- prot_exp2 %>% distinct(ENSG, variable, .keep_all = TRUE)

# subset for MMUT
prot_exp3 <- prot_exp2[ENSG == "ENSG00000146085"]

setnames(prot_exp3, c("ENSG", "mma_id", "prot_exp"))
setorder(prot_exp3, "mma_id")

types <- c(rep("MMUT def.", 150), rep("unknown", 60), rep("unaffected", 20))
prot_exp3[, type := types]
prot_exp3[, sample_no := rep(1:230)]
prot_exp3$type <- factor(prot_exp3$type, levels = c("MMUT def.", "unknown", "unaffected"))


##################################
# load transcript matrix
load("interimData/rnafiltered_mat.RData")

rna_exp <- data.table(ENSG = row.names(rnafiltered_mat), rnafiltered_mat)
rna_exp1 <- melt(rna_exp, id = c("ENSG"))

# subset for MMUT
rna_exp2 <- rna_exp1[ENSG == "ENSG00000146085"]
# rna_exp_MMUT <- rna_exp1[ENSG == "ENSG00000146085"]
# # GLUD1 Ensembl: ENSG00000148672
# rna_exp_GLUD1 <- rna_exp1[ENSG == "ENSG00000148672"]
# # Ensembl: ENSG00000075624 (ACTB)
# rna_exp_ACTB <- rna_exp1[ENSG == "ENSG00000075624"]

# rna_exp_GLUD1$variable == rna_exp_MMUT$variable
# rna_exp_ACTB$variable == rna_exp_MMUT$variable

setnames(rna_exp2, c("ENSG", "mma_id", "rna_exp"))

rna_exp2[, sample_no := (as.numeric(sub("MMA", "", rna_exp2$mma_id)))]
rna_exp2[, type := 
  ifelse(sample_no<151, "MMUT def.", 
    ifelse(sample_no>210, "unaffected", "unknown"))]
rna_exp2$type <- factor(rna_exp2$type, levels = c("MMUT def.", "unknown", "unaffected"))
unique(rna_exp2[, type])

rna_exp2[, type2 := ifelse(sample_no<151, "MMUT def.", "control")]
rna_exp2$type2 <- factor(rna_exp2$type2, levels = c("MMUT def.", "control"))


##################################
# phenotype table with biochemistry and mut subtype
pheno_tbl0 <- fread("interimData/phenotypic_tables/meta_info_wcss.txt")
pheno_tbl0 <- data.table(pheno_tbl0)
pheno_tbl <- pheno_tbl0 %>%
  select("mma_id", "OHCblMinus", "OHCblPlus", "ratio", "SimultOHCblMinus", "SimultOHCblPlus", "mut_category", "init_mut_category", "AdoCblMinus", "AdoCblPlus", "SimultAdoCblMinus", "SimultAdoCblPlus", "wgs_bothtype", "wgs_qc", "rnaseq_qc", "prot_qc")

unique(pheno_tbl[,wgs_bothtype])
pheno_tbl[wgs_bothtype == "/", wgs_bothtype := ""]
pheno_tbl[wgs_bothtype == "splicing/missense", wgs_bothtype := "missense/splicing"]
pheno_tbl[wgs_bothtype == "truncating/missense", wgs_bothtype := "missense/truncating"]
pheno_tbl[wgs_bothtype == "splicing/truncating", wgs_bothtype := "truncating/splicing"]

pheno_tbl[wgs_bothtype == "deletion/deletion", wgs_bothtype2 := "del/del"]
pheno_tbl[wgs_bothtype == "missense/missense", wgs_bothtype2 := "mis/mis"]
pheno_tbl[wgs_bothtype == "missense/splicing", wgs_bothtype2 := "mis/spl"]
pheno_tbl[wgs_bothtype == "missense/truncating", wgs_bothtype2 := "mis/tru"]
pheno_tbl[wgs_bothtype == "splicing/splicing", wgs_bothtype2 := "spl/spl"]
pheno_tbl[wgs_bothtype == "truncating/splicing", wgs_bothtype2 := "tru/spl"]
pheno_tbl[wgs_bothtype == "truncating/truncating", wgs_bothtype2 := "tru/tru"]


pheno_tbl$AdoCblPlus <- as.numeric(pheno_tbl$AdoCblPlus)

summary(pheno_tbl$rnaseq_qc == "NA")

pheno_tbl0[wgs_bothtype2 == "spl/spl" | wgs_bothtype2 == "tru/spl" | wgs_bothtype2 == "mis/spl", c("mma_id", "wgs_mut1nt", "wgs_mut2nt")]
rna_exp2[mma_id == "MMA057" | mma_id == "MMA060" | mma_id == "MMA106" | mma_id == "MMA116" | mma_id == "MMA148", c("mma_id", "rna_exp")]

# merge tables for prot_pheno and rna_pheno

prot_pheno <- merge(prot_exp3, pheno_tbl, by = "mma_id", all = TRUE)

rna_pheno <- merge(rna_exp2, pheno_tbl, by = "mma_id", all = TRUE)
rna_pheno <- rna_pheno[!is.na(rna_pheno[,rna_exp])]

prot_pheno_val <- prot_pheno %>% select(mma_id, prot_exp)

prot_pheno_val1 <- prot_pheno_val[which(prot_pheno_val[, mma_id] %in% rna_pheno[, mma_id]),]
setnames(prot_pheno_val1, c("mma_id", "prot_exp"))

prot_rna_merge <- merge(prot_pheno_val1, rna_pheno, by = "mma_id")


# prot_pheno[sample_no < 151, type_MMUTmut := paste(type, mut_category)]
prot_pheno[sample_no < 151, type_MMUTmut := paste(mut_category)]
prot_pheno[sample_no > 150, type_MMUTmut := type]
# prot_pheno$type_MMUTmut <- factor(prot_pheno$type_MMUTmut, levels = c("MMUT mut0", "MMUT mut-", "UNKN", "CTRL"))
prot_pheno$type_MMUTmut <- factor(prot_pheno$type_MMUTmut, levels = c("mut0", "mut-", "unknown", "unaffected"))
prot_pheno$mut_category <- factor(prot_pheno$mut_category, levels = c("mut0", "mut-"))






#####################
# prepare table for lolliplot of MMUT gene

loli_tbl0 <- data.table(mutaa = c(pheno_tbl0$wgs_mut1aa[1:150], pheno_tbl0$wgs_mut2aa[1:150]),
  muttype = c(pheno_tbl0$wgs_mut1type[1:150], pheno_tbl0$wgs_mut2type[1:150]))

# splicing mutations cannot be plotted on the aa polypeptide --> exclude them
loli_tbl <- loli_tbl0[muttype != "splicing", ]

# count occurrence of unique mutations
loli_tbl <- as.data.table(table(loli_tbl$mutaa))
setnames(loli_tbl, c("aa_change", "freq"))

# put back mutation type
matcher <- match(loli_tbl$aa_change, loli_tbl0$mutaa)
loli_tbl[, type := loli_tbl0[matcher, ]$muttype]

# change short aa abbreviations and add position column
loli_tbl[, ref_aa := gsub("\\d+.+", "", gsub("p\\.\\(", "", aa_change))]
loli_tbl[, alt_aa := gsub("\\)", "", gsub("p\\.\\(\\D+\\d+", "", aa_change))]
loli_tbl[, pos_aa := gsub("\\D+", "", gsub("\\D+\\*\\d+\\)", "", gsub("p\\.\\(\\D+", "", aa_change)))]

loli_tbl$ref_aa
loli_tbl$alt_aa
loli_tbl$pos_aa

# import amino acid abbreviations file
# from https://www.ddbj.nig.ac.jp/ddbj/code-e.html 20 Apr 2021
aa_file <- read.csv("Data/aa_file.csv")
aa_file <- as.data.table(aa_file)


setnames(aa_file, c("ref_aa", "one_letter", "aa_name"))
setDT(loli_tbl)[, ref_aa_one := setDT(aa_file)[loli_tbl, one_letter, on = "ref_aa"]]
setnames(aa_file, c("alt_aa", "one_letter", "aa_name"))
setDT(loli_tbl)[, alt_aa_one := setDT(aa_file)[loli_tbl, one_letter, on = "alt_aa"]]
loli_tbl[, aa_change_short := ifelse(loli_tbl$type == "truncating", paste(loli_tbl$ref_aa_one, loli_tbl$pos_aa, "*", sep = ""), paste(loli_tbl$ref_aa_one, loli_tbl$pos_aa, loli_tbl$alt_aa_one, sep = "")) ]

print(loli_tbl[, c("aa_change_short", "aa_change")], nrows = nrow(loli_tbl))

loli_tbl[aa_change == "p.(Leu347del)", aa_change_short := "L347del"]
loli_tbl[aa_change == "p.(Leu347del)", type := "deletion"]

# adjust factor levels for mutation types
loli_tbl$type <- factor(loli_tbl$type, levels = c("missense", "truncating", "deletion"))
loli_tbl$pos_aa <- as.numeric(loli_tbl$pos_aa)



# transit peptide (mitochondrial leader sequence: 1-32
# chain (adenosyltransferase): 33-250
chain_tbl <- data.table(rbind(c("mito_peptide", 1, 32),c("adenosyltransferase", 33, 250)))
colnames(chain_tbl) <- c("feature", "start", "end")
chain_tbl$start <- as.numeric(chain_tbl$start)
chain_tbl$end <- as.numeric(chain_tbl$end)

# exon positions
# info from https://databases.lovd.nl/shared/refseq/MUT_NM_000255.3_table.html (21 July 2021, 3 PM)
# exon  c.startExon c.endExon g.startExon g.endExon lengthExon  lengthIntron
# 1 -265  -40 4812  5037  226 3597
# 2 -39 385 8635  9058  424 1023
# 3 386 753 10082 10449 368 1453
# 4 754 911 11903 12060 158 2323
# 5 912 1083  14384 14555 172 1870
# 6 1084  1332  16426 16674 249 2538
# 7 1333  1444  19213 19324 112 1030
# 8 1445  1560  20355 20470 116 2915
# 9 1561  1676  23386 23501 116 2667
# 10  1677  1808  26169 26300 132 1486
# 11  1809  1956  27787 27934 148 4582
# 12  1957  2124  32517 32684 168 3599
# 13  2125  *1368 36284 37780 1497    

exon_bound <- data.frame(ex = c(386/3, 754/3, 912/3, 1084/3, 1333/3, 1445/3, 1561/3, 1677/3, 1809/3, 1957/3, 2125/3))
exs <- exon_bound$ex

exon_labels <- data.frame(lbl = c(mean(c(-39, exs[1])), 
  mean(c(exs[1], exs[2])),
  mean(c(exs[2], exs[3])),
  mean(c(exs[3], exs[4])),
  mean(c(exs[4], exs[5])),
  mean(c(exs[5], exs[6])),
  mean(c(exs[6], exs[7])),
  mean(c(exs[7], exs[8])),
  mean(c(exs[8], exs[9])),
  mean(c(exs[9], exs[10])),
  mean(c(exs[10], exs[11])),
  mean(c(exs[11], 750))), exon = c(seq(2,13, by = 1)))

# total peptide length of MMUT: 1-750, resp. -13-750
chain_tbl <- data.table(cbind("MMUT peptide", "-13", "750"))
colnames(chain_tbl) <- c("feature", "start", "end")
chain_tbl$start <- as.numeric(chain_tbl$start)
chain_tbl$end <- as.numeric(chain_tbl$end)





##################################
##################################
# lolliplot of MMUT variants

loliplt <- 
ggplot(loli_tbl, aes(x = pos_aa, y = freq, color = type, label = aa_change_short)) +
  geom_point(size = 2) +
  geom_segment(aes(x = pos_aa, y = 0, xend = pos_aa, yend = freq)) +
  scale_color_uchicago() +
  geom_text_repel(nudge_y = 0.1, alpha = 1, size = 3, show.legend = FALSE, ylim = c(0,17), box.padding = 0.7, segment.curvature = -0.1, segment.angle = 20) +
  geom_rect(data = chain_tbl, aes(xmin = start, xmax = end, ymin = -3, ymax = -1), colour = "grey", fill = "grey", inherit.aes = FALSE) +
  geom_segment(data = exon_bound, aes(x = ex, xend = ex, y = -3, yend = -1), inherit.aes = FALSE, color = "white") + 
  geom_text(data = exon_labels, aes(x = lbl, y = -2, label = exon), inherit.aes = FALSE, color = "black") +
  annotate(geom = "text", x = max(loli_tbl[!is.na(loli_tbl$pos_aa), ]$pos_aa)+20, y = -2, label = "Exon", hjust = 0) +
  scale_x_continuous(limits = c(-20, 800), breaks = c(0, 200, 400, 600, 750)) +
  theme_pubr() +
  xlab("Amino acid number") +
  ylab("Allele count") +
  ggtitle(expression(paste(italic("MMUT"), ": distribution of damaging variants"))) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(), panel.grid.major.y = element_line(), legend.position = c(0.9, 0.85), legend.title = element_blank(), plot.title = element_text(size = 12))


ggsave(paste(fig_path,"LoliplotMMUT.png", sep = ""), loliplt, device = png(), width = 7, height = 3.5)
dev.off()
ggsave(paste(fig_path_pdf,"Fig1/","LoliplotMMUT.pdf", sep = ""), loliplt, device = "pdf", width = 7, height = 3.5)
dev.off()






##################################
##################################
# Frequency of mutation types


freq_tbl0 <- data.table(muttype = c(pheno_tbl0$wgs_mut1type[1:150], pheno_tbl0$wgs_mut2type[1:150]))

# count occurrence of mutation types
freq_tbl <- as.data.table(table(freq_tbl0$muttype))
setnames(freq_tbl, c("type", "freq"))
freq_tbl[, combo := paste(type, paste("(N=", freq, ")", sep = ""), sep = " ")]
combs <- freq_tbl$combo
freq_tbl$combo <- factor(freq_tbl$combo, levels = c(combs[2],combs[4],combs[1],combs[3]))


freqMMUT <- 
ggplot(freq_tbl, aes(area = freq, color = combo, fill = combo, label = type)) +
  geom_treemap() +
  geom_treemap_text(colour = "white", place = "centre", grow = FALSE) +
  scale_color_uchicago() +
  scale_fill_uchicago() +
  ggtitle(expression(paste("Frequency of ",italic("MMUT"), " variant types"))) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) + 
  theme(legend.position = "bottom", legend.title = element_blank(), axis.line = element_blank(), plot.title = element_text(size = 12))

ggsave(paste(fig_path,"VarFreqMMUT.png", sep = ""), freqMMUT, device = png(), width = 3, height = 4)
dev.off()
ggsave(paste(fig_path_pdf,"Fig1/","VarFreqMMUT.pdf", sep = ""), freqMMUT, device = "pdf", width = 3, height = 4)
dev.off()








#####################
# prepare table for lolliplot of ACSF3 gene

loliACSF3_tbl0 <- data.table(mutaa = c(pheno_tbl0[proposed_gene_wgs == "ACSF3", ]$wgs_mut1aa, pheno_tbl0[proposed_gene_wgs == "ACSF3", ]$wgs_mut2aa),
  muttype = c(pheno_tbl0[proposed_gene_wgs == "ACSF3", ]$wgs_mut1type, pheno_tbl0[proposed_gene_wgs == "ACSF3", ]$wgs_mut2type))

# splicing mutations cannot be plotted on the aa polypeptide --> exclude them
loliACSF3_tbl <- loliACSF3_tbl0[muttype != "splicing", ]

# count occurrence of unique mutations
loliACSF3_tbl <- as.data.table(table(loliACSF3_tbl$mutaa))
setnames(loliACSF3_tbl, c("aa_change", "freq"))

# put back mutation type
matcher <- match(loliACSF3_tbl$aa_change, loliACSF3_tbl0$mutaa)
loliACSF3_tbl[, type := loliACSF3_tbl0[matcher, ]$muttype]

# change short aa abbreviations and add position column
loliACSF3_tbl[, ref_aa := gsub("\\d+.+", "", gsub("p\\.\\(", "", aa_change))]
loliACSF3_tbl[, alt_aa := gsub("\\)", "", gsub("p\\.\\(\\D+\\d+", "", aa_change))]
loliACSF3_tbl[, pos_aa := gsub("\\D+", "", gsub("\\D+\\*\\d+\\)", "", gsub("p\\.\\(\\D+", "", aa_change)))]

loliACSF3_tbl$ref_aa
loliACSF3_tbl$alt_aa
loliACSF3_tbl$pos_aa



# import amino acid abbreviations file
# from https://www.ddbj.nig.ac.jp/ddbj/code-e.html 20 Apr 2021
aa_file <- read.csv("Data/aa_file.csv")
aa_file <- as.data.table(aa_file)


setnames(aa_file, c("ref_aa", "one_letter", "aa_name"))
setDT(loliACSF3_tbl)[, ref_aa_one := setDT(aa_file)[loliACSF3_tbl, one_letter, on = "ref_aa"]]
setnames(aa_file, c("alt_aa", "one_letter", "aa_name"))
setDT(loliACSF3_tbl)[, alt_aa_one := setDT(aa_file)[loliACSF3_tbl, one_letter, on = "alt_aa"]]
loliACSF3_tbl[, aa_change_short := ifelse(loliACSF3_tbl$type == "truncating", paste(loliACSF3_tbl$ref_aa_one, loliACSF3_tbl$pos_aa, "*", sep = ""), paste(loliACSF3_tbl$ref_aa_one, loliACSF3_tbl$pos_aa, loliACSF3_tbl$alt_aa_one, sep = "")) ]

print(loliACSF3_tbl[, c("aa_change_short", "aa_change")], nrows = nrow(loliACSF3_tbl))

loliACSF3_tbl[aa_change == "p.(Leu347del)", aa_change_short := "L347del"]
loliACSF3_tbl[aa_change == "p.(Leu347del)", type := "deletion"]

# adjust factor levels for mutation types
loliACSF3_tbl$type <- factor(loliACSF3_tbl$type, levels = c("missense", "truncating", "deletion"))
loliACSF3_tbl$pos_aa <- as.numeric(loliACSF3_tbl$pos_aa)


# exon positions
# info from https://databases.lovd.nl/shared/refseq/ACSF3_NM_174917.3_table.html (3 Aug 2021, 4:43 PM)
# exon  c.startExon c.endExon g.startExon g.endExon lengthExon  lengthIntron
# 1 -377  -194  5001  5184  184 4598
# 2 -193  -21 9783  9955  173 1898
# 3 -20 666 11854 12539 686 1256
# 4 667 822 13796 13951 156 9332
# 5 823 977 23284 23438 155 2092
# 6 978 1126  25531 25679 149 6313
# 7 1127  1239  31993 32105 113 12222
# 8 1240  1366  44328 44454 127 12004
# 9 1367  1501  56459 56593 135 536
# 10  1502  1613  57130 57241 112 8040
# 11  1614  *1556 65282 66955 1674   

exon_bound <- data.frame(ex = c(667/3, 823/3, 978/3, 1127/3, 1240/3, 1367/3, 1502/3, 1614/3))
exs <- exon_bound$ex

exon_labels <- data.frame(lbl = c(mean(c(-20, exs[1])), 
  mean(c(exs[1], exs[2])),
  mean(c(exs[2], exs[3])),
  mean(c(exs[3], exs[4])),
  mean(c(exs[4], exs[5])),
  mean(c(exs[5], exs[6])),
  mean(c(exs[6], exs[7])),
  mean(c(exs[7], exs[8])),
  mean(c(exs[8], 576))), exon = c(seq(3,11, by = 1)))

# total number of amino acids derived from here: https://www.ncbi.nlm.nih.gov/protein/NP_001120686.1 (3 Aug 2021, 4:49 PM)

# total peptide length of ACSF3: 1-750, resp. -13-750
chain_tbl <- data.table(cbind("ACSF3 peptide", "-20", "576"))
colnames(chain_tbl) <- c("feature", "start", "end")
chain_tbl$start <- as.numeric(chain_tbl$start)
chain_tbl$end <- as.numeric(chain_tbl$end)





##################################
##################################
# lolliplot of ACSF3 alleles

loliplt_acsf3 <- 
ggplot(loliACSF3_tbl, aes(x = pos_aa, y = freq, color = type, label = aa_change_short)) +
  geom_point(size = 2) +
  geom_segment(aes(x = pos_aa, y = 0, xend = pos_aa, yend = freq)) +
  scale_color_uchicago() +
  geom_text_repel(nudge_y = -0.1, alpha = 1, size = 3, show.legend = FALSE, ylim = c(0,20), box.padding = 0.5, segment.curvature = -0.1, segment.angle = 20) +
  geom_rect(data = chain_tbl, aes(xmin = start, xmax = end, ymin = -3, ymax = -1), colour = "grey", fill = "grey", inherit.aes = FALSE) +
  geom_segment(data = exon_bound, aes(x = ex, xend = ex, y = -3, yend = -1), inherit.aes = FALSE, color = "white") + 
  geom_text(data = exon_labels, aes(x = lbl, y = -2, label = exon), inherit.aes = FALSE, color = "black") +
  annotate(geom = "text", x = max(loliACSF3_tbl[!is.na(loliACSF3_tbl$pos_aa), ]$pos_aa)+20, y = -2, label = "Exon", hjust = 0) +
  scale_x_continuous(limits = c(-30, 600), breaks = c(0, 200, 400, 576)) +
  theme_pubr() +
  xlab("Amino acid number") +
  ylab("Allele count") +
  ggtitle(expression(paste(italic("ACSF3"), ": distribution of damaging variants"))) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(), panel.grid.major.y = element_line(), legend.position = c(0.1, 0.95), legend.title = element_blank(), plot.title = element_text(size = 12))


ggsave(paste(fig_path,"LoliplotACSF3.png", sep = ""), loliplt_acsf3, device = png(), width = 5, height = 3.5)
dev.off()
ggsave(paste(fig_path_pdf,"Fig1/","LoliplotACSF3.pdf", sep = ""), loliplt_acsf3, device = "pdf", width = 5, height = 3.5)
dev.off()










##################################
##################################
# "parts of whole" chart for identified genes in patients

freq_genes <- data.table(MMUT = 148,
  MMAA = 2,
  MMAB = 1,
  ACSF3 = 19,
  SUCLA2 = 4,
  TCN2 = 3,
  unknown = (210-(148+19+2+1+4+3)))

freq_genes_melt <- melt.data.table(freq_genes)
setnames(freq_genes_melt, c("gene", "frequency"))
freq_genes_melt[, combo := paste(gene, paste("(N=", frequency, ")", sep = ""), sep = " ")]
freq_genes_melt$combo <- factor(freq_genes_melt$combo, levels = c(freq_genes_melt$combo))


freq_genes_plt <-
ggplot(freq_genes_melt, aes(area = frequency, color = combo, fill = combo, label = gene)) +
  geom_treemap(alpha=1,size=0) +
  geom_treemap_text(colour = "white", place = "centre", grow = FALSE) +
  scale_fill_jama() +
  scale_color_jama()+
  ggtitle("Genes with identified biallelic variants") +
  theme_pubr() +
  theme(legend.position = "right", legend.title = element_blank(), axis.line = element_blank(), plot.title = element_text(size = 12))

ggsave(paste(fig_path,"IdentifiedGenes.png", sep = ""), freq_genes_plt, device = png(), width = 5, height = 3)
dev.off()
ggsave(paste(fig_path_pdf,"Fig1_diagnosis/","IdentifiedGenes.pdf", sep = ""), freq_genes_plt, device = "pdf", width = 5, height = 3)
dev.off()
  



##################################
##################################
# MMUT plots PROTEIN

# MMUT levels in all patients
pre_p_MMUT_all <-
ggplot(prot_pheno, aes(x = sample_no, y = prot_exp, color = type)) +
  geom_point(alpha = 1) +
  scale_y_log10() +
  annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  ggtitle("MMUT protein levels") +
  ylab("Protein level") +
  xlab("Sample number") +
  scale_color_manual(values = mypal2) +
  theme_pubr() +
  theme(legend.title = element_blank(), legend.position = "right", plot.title = element_text(size = 12))

ybox_p <- axis_canvas(pre_p_MMUT_all, axis = "y") +
  geom_boxplot(data = prot_pheno, aes(x = 0, y = prot_exp, fill = type, color = type, alpha = 0.7)) +
  scale_y_log10() +
  scale_color_manual(values = mypal2) +
  scale_fill_manual(values = mypal2)

plot_p_MMUT_all <- insert_yaxis_grob(pre_p_MMUT_all, ybox_p, position = "right", clip = "off")
plot_p_MMUT_all <- ggdraw(plot_p_MMUT_all)
plot_p_MMUT_all

ggsave(paste(fig_path,"MMUTProteinLevel.png", sep = ""), plot_p_MMUT_all, device = png(), width = 5, height = 3)
dev.off()



# MMUT levels in mut0 samples according to mutation type

# pre_p_MMUT_mut0 <-
# ggplot(prot_pheno[mut_category == "mut0" & wgs_bothtype != ""], aes(x = sample_no, y = prot_exp, color = wgs_bothtype)) +
#   geom_point(alpha = .7) +
#   scale_y_log10() +
#   ggtitle("p_MMUT_lvls_mut0") +
#   theme_bw()

pre_p_MMUT_muttype <-
ggplot(prot_pheno[wgs_bothtype != ""], aes(x = sample_no, y = prot_exp, color = wgs_bothtype)) +
  geom_point(alpha = 1) +
  scale_y_log10() +
  annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  ggtitle("p_MMUT_lvls_MMUT_mutations") +
  ylab("Protein level") +
  xlab("Sample number") +
  ggtitle("MMUT protein levels") +
  scale_color_jco() +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_pubr() +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(size = 12))

ybox_muttype <- axis_canvas(pre_p_MMUT_muttype, axis = "y")+
    geom_boxplot(data = prot_pheno[mut_category == "mut0" & wgs_bothtype != ""], aes(x = 0, y = prot_exp, fill = wgs_bothtype, color = wgs_bothtype, alpha = 0.7))+
    scale_y_log10() +
    scale_color_jco() +
    scale_fill_jco()

plot_p_MMUT_muttype <- insert_yaxis_grob(pre_p_MMUT_muttype, ybox_muttype, position = "right", clip = "off", width = grid::unit(0.5, "null"))
plot_p_MMUT_muttype <- ggdraw(plot_p_MMUT_muttype)
plot_p_MMUT_muttype



# MMUT levels in 1-150 MMUT-deficient samples according to mut subtype
pre_p_MMUT_mut <-
ggplot(prot_pheno[sample_no<151], aes(x = sample_no, y = prot_exp, color = mut_category)) +
  geom_point(alpha = 1) +
  scale_y_log10() +
  ggtitle("p_MMUT_lvls_MMUT_mut") +
  theme_bw()

ybox_mut <- axis_canvas(pre_p_MMUT_mut, axis = "y")+
    geom_boxplot(data = prot_pheno[sample_no<151], aes(x = 0, y = prot_exp, fill = mut_category, color = mut_category, alpha = 0.7))+
    scale_y_log10()

plot_p_MMUT_mut <- insert_yaxis_grob(pre_p_MMUT_mut, ybox_mut, position = "right")
plot_p_MMUT_mut <- ggdraw(plot_p_MMUT_mut)





##################################
##################################
# MMUT plots RNA

# MMUT levels in all patients
pre_r_MMUT_all <-
ggplot(rna_pheno, aes(x = sample_no, y = rna_exp, color = type)) +
  geom_point(alpha = 1) +
  scale_y_log10() +
  annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  ggtitle("MMUT transcript levels") +
  ylab("Transcript level") +
  xlab("Sample number") +
  scale_color_manual(values = mypal2) +
  theme_pubr() +
  theme(legend.title = element_blank(), legend.position = "right", plot.title = element_text(size = 12))

ybox_r <- axis_canvas(pre_r_MMUT_all, axis = "y") +
  geom_boxplot(data = rna_pheno, aes(x = 0, y = rna_exp, fill = type, color = type, alpha = 0.7)) +
  scale_y_log10() +
  scale_color_manual(values = mypal2) +
  scale_fill_manual(values = mypal2)

plot_r_MMUT_all <- insert_yaxis_grob(pre_r_MMUT_all, ybox_r, position = "right", clip = "off")
plot_r_MMUT_all <- ggdraw(plot_r_MMUT_all)

ggsave(paste(fig_path,"MMUTRNALevel.png", sep = ""), plot_r_MMUT_all, device = png(), width = 5, height = 3)
dev.off()



# MMUT levels in all patients (2 categories)

# expression color palette
exp_pal <- c("#9302bf","#d078eb")

# dea color palette
dea_pal <- c( "#084594","#4292C6")

# prot/rna correlation color palette
corr_pal <- c("#9c1111" ,"#e86107")



pre_r_MMUT_all <-
ggplot(rna_pheno, aes(x = sample_no, y = rna_exp, color = type2)) +
  geom_point(alpha = 1) +
  scale_y_log10() +
  annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  ggtitle("MMUT raw transcript levels") +
  ylab("Transcript level") +
  xlab("Sample number") +
  scale_color_manual(values = c(exp_pal[1], "darkgrey")) +
  theme_pubr() +
  theme(legend.title = element_blank(), legend.position = "none", plot.title = element_text(size = 12), legend.text = element_text(size = 12))

ybox_r <- axis_canvas(pre_r_MMUT_all, axis = "y") +
  geom_boxplot(data = rna_pheno, aes(x = 0, y = rna_exp, fill = type2, color = type2, alpha = 0.7)) +
  scale_y_log10() +
  scale_color_manual(values = c(exp_pal[1], "darkgrey")) +
  scale_fill_manual(values = c(exp_pal[1], "darkgrey"))

plot_r_MMUT_all2 <- insert_yaxis_grob(pre_r_MMUT_all, ybox_r, position = "right")
plot_r_MMUT_all2 <- ggdraw(plot_r_MMUT_all2)

ggsave(paste(fig_path,"MMUTRNALevel_type2.png", sep = ""), plot_r_MMUT_all, device = png(), width = 4, height = 2)
dev.off()




# MMUT levels in mut0 samples according to mutation type

# pre_r_MMUT_mut0 <-
# ggplot(rna_pheno[mut_category == "mut0" & wgs_bothtype != ""], aes(x = sample_no, y = rna_exp, color = wgs_bothtype)) +
#   geom_point(alpha = .7) +
#   scale_y_log10() +
#   ggtitle("r_MMUT_lvls_mut0") +
#   theme_bw()

pre_r_MMUT_muttype <-
ggplot(rna_pheno[wgs_bothtype != ""], aes(x = sample_no, y = rna_exp, color = wgs_bothtype)) +
  geom_point(alpha = 1) +
  scale_y_log10() +
  annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  ggtitle("r_MMUT_lvls_MMUT_mutations") +
  ylab("Transcript level") +
  xlab("Sample number") +
  ggtitle("MMUT transcript levels") +
  scale_color_jco() +
  theme_pubr() +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(size = 12))

r_ybox_muttype <- axis_canvas(pre_r_MMUT_muttype, axis = "y")+
    geom_boxplot(data = rna_pheno[mut_category == "mut0" & wgs_bothtype != ""], aes(x = 0, y = rna_exp, fill = wgs_bothtype, color = wgs_bothtype, alpha = 0.7), width = 10)+
    scale_y_log10() +
    scale_color_jco() +
    scale_fill_jco()


plot_r_MMUT_muttype <- insert_yaxis_grob(pre_r_MMUT_muttype, r_ybox_muttype, position = "right", clip = "off", width = grid::unit(0.5, "null"))
plot_r_MMUT_muttype <- ggdraw(plot_r_MMUT_muttype)
plot_r_MMUT_muttype


MMUTlvls <- 
ggarrange(plot_r_MMUT_muttype, plot_p_MMUT_muttype, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

ggsave(paste(fig_path,"MMUTlevelsMutationTypes.png", sep = ""), MMUTlvls, device = png(), width = 5.5, height = 7)
dev.off()

ggsave(paste(fig_path_pdf,"SuppFig2_diagnosis/","MMUTlevelsMutationTypes.pdf", sep = ""), MMUTlvls, device = "pdf", width = 5.5, height = 7)
dev.off()





# MMUT levels in 1-150 MMUT-deficient samples according to mut subtype
pre_r_MMUT_mut <-
ggplot(rna_pheno[sample_no<151], aes(x = sample_no, y = rna_exp, color = mut_category)) +
  geom_point(alpha = .7) +
  scale_y_log10() +
  ggtitle("r_MMUT_lvls_MMUT_mut") +
  theme_bw()

r_ybox_mut <- axis_canvas(pre_r_MMUT_mut, axis = "y")+
    geom_boxplot(data = rna_pheno[sample_no<151], aes(x = 0, y = rna_exp, fill = mut_category, color = mut_category, alpha = 0.7))+
    scale_y_log10()

plot_r_MMUT_mut <- insert_yaxis_grob(pre_r_MMUT_mut, r_ybox_mut, position = "right")
plot_r_MMUT_mut <- ggdraw(plot_r_MMUT_mut)






# combine transcript and protein levels of MMUT for Fig. 1
MMUT_levels_plot <- 
ggarrange(plot_r_MMUT_all, plot_p_MMUT_all, ncol = 2, common.legend = TRUE)
ggsave(paste(fig_path_pdf,"Fig1_diagnosis/","MMUT_levels_plot.pdf", sep = ""), MMUT_levels_plot, device = "pdf", width = 10, height = 3)
dev.off()





# combine transcript and protein plots regarding MMUT levels

MMUT_rna_prot <- ggarrange(plot_p_MMUT_all, plot_p_MMUT_muttype, plot_p_MMUT_mut, plot_r_MMUT_all, plot_r_MMUT_muttype, plot_r_MMUT_mut, widths = c(1, 1, 1))

ggsave(paste(fig_path,"MMUT_rna_prot_plot.png", sep = ""), MMUT_rna_prot, device = png(), width = 300, height = 120, units = c("mm"), dpi = 600)
dev.off()


##################################
# identify the MMUT-deficient samples with low MMUT transcript levels and compare to the protein level
lowRNA_ids <- rna_pheno[which(rna_pheno[, rna_exp] <3), mma_id]

prot_pheno_lowRNA <- prot_pheno[sample_no<151]
prot_pheno_lowRNA[, lowRNA := FALSE]
prot_pheno_lowRNA[which(prot_pheno_lowRNA$mma_id %in% lowRNA_ids), lowRNA := TRUE]
prot_pheno_lowRNA$lowRNA <- factor(prot_pheno_lowRNA$lowRNA, levels = c(TRUE, FALSE))

p_MMUT_lowRNA <-
ggplot(prot_pheno_lowRNA, aes(x = sample_no, y = prot_exp, color = lowRNA)) +
  geom_point(alpha = 0.7) +
  scale_y_log10() +
  ggtitle("p_MMUT_lvls_lowRNA") +
  theme_bw()

# correlation of MMUT transcript and protein
MMUT_rnaVSprot_plot <-
ggplot(prot_rna_merge[prot_exp>1500], aes(x = rna_exp, y = prot_exp)) +
  geom_point(alpha = 0.6) +
  scale_y_log10() +
  scale_x_log10() +
  # scale_y_sqrt() +
  # scale_x_sqrt() +
  # coord_trans(x="log10", y="log10") +
  geom_smooth(method = "lm", se = FALSE, color = "darkblue", size = 0.5) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  ggtitle("MMUT_rna_prot_ALL") +
  theme_bw()


MMUT_lowRNA <- ggarrange(plot_r_MMUT_all, p_MMUT_lowRNA, MMUT_rnaVSprot_plot, widths = c(1.3, 1, 1))

ggsave(paste(fig_path,"p_MMUT_lowRNA.png", sep = ""), MMUT_lowRNA, device = png(), width = 240, height = 150, units = c("mm"), dpi = 600)
dev.off()


##################################
##################################
# compare biochemistry to MMUT mutation type

pheno_tbl[, wgs_bothtype_ctrl := wgs_bothtype2]
pheno_tbl[, cell_no := as.numeric(sub("MMA", "", mma_id))]
pheno_tbl[cell_no > 150 & wgs_bothtype2 == "", wgs_bothtype_ctrl := "nonMMUT"]
pheno_tbl[, wgs_bothtype_ctrl]
pheno_tbl$wgs_bothtype_ctrl <- factor(pheno_tbl$wgs_bothtype_ctrl, levels = c("nonMMUT", "del/del", "mis/mis", "mis/spl", "mis/tru", "spl/spl", "tru/spl", "tru/tru", ""))


plt_pthwy_muttype_ridges <- 
ggplot(pheno_tbl[wgs_bothtype_ctrl != "" & wgs_bothtype_ctrl != "del/del"], aes(x = log(OHCblPlus), y = as.factor(wgs_bothtype_ctrl), fill = as.factor(wgs_bothtype_ctrl), color = as.factor(wgs_bothtype_ctrl))) +
  geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
  scale_color_aaas() +
  scale_fill_aaas() +
  theme_test()

# ggplot(pheno_tbl[wgs_bothtype2 != "" & wgs_bothtype2 != "del/del"], aes(x = log(OHCblPlus), y = as.factor(wgs_bothtype2), fill = as.factor(wgs_bothtype2), color = as.factor(wgs_bothtype2))) +
#   geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
#   scale_color_aaas() +
#   scale_fill_aaas() +
#   theme_test()


plt_pthwy_muttype_box <- 
ggplot(pheno_tbl[wgs_bothtype_ctrl != "" & wgs_bothtype_ctrl != "del/del"], aes(x = as.factor(wgs_bothtype_ctrl), y = log(OHCblPlus), fill = as.factor(wgs_bothtype_ctrl), color = as.factor(wgs_bothtype_ctrl))) +
  geom_boxplot(alpha = 0.6) +
  scale_color_aaas() +
  scale_fill_aaas() +
  theme_test()

# ggplot(prot_pheno[wgs_bothtype != "" & wgs_bothtype != "del/del"], aes(x = as.factor(wgs_bothtype), y = log(OHCblPlus), fill = as.factor(wgs_bothtype), color = as.factor(wgs_bothtype))) +
#   geom_boxplot(alpha = 0.6) +
#   scale_color_aaas() +
#   scale_fill_aaas() +
#   theme_test()

# add 0.01 to AdoCbl as 0 cannot be log

pheno_tbl[, AdoCblPlus_01 := AdoCblPlus]
pheno_tbl[, AdoCblMinus_01 := AdoCblMinus]
pheno_tbl[AdoCblPlus_01 == 0, ]$AdoCblPlus_01 <- 0.1
pheno_tbl$AdoCblPlus_01
pheno_tbl[AdoCblMinus_01 == 0, ]$AdoCblMinus_01 <- 0.1


# plt_AdoCblPlus_muttype_ridges <- 
# ggplot(pheno_tbl[wgs_bothtype_ctrl != "" & wgs_bothtype_ctrl != "del/del" & AdoCblPlus != 0], aes(x = log(AdoCblPlus), y = as.factor(wgs_bothtype_ctrl), fill = as.factor(wgs_bothtype_ctrl), color = as.factor(wgs_bothtype_ctrl))) +
#   geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
#   scale_color_aaas() +
#   scale_fill_aaas() +
#   theme_test()

# plt_AdoCblPlus_muttype_box <- 
# ggplot(pheno_tbl[wgs_bothtype_ctrl != "" & wgs_bothtype_ctrl != "del/del"], aes(x = as.factor(wgs_bothtype_ctrl), y = log(AdoCblPlus), fill = as.factor(wgs_bothtype_ctrl), color = as.factor(wgs_bothtype_ctrl))) +
#   geom_boxplot(alpha = 0.6) +
#   scale_color_aaas() +
#   scale_fill_aaas() +
#   theme_test()


plt_AdoCblPlus_muttype_ridges <- 
ggplot(pheno_tbl[wgs_bothtype_ctrl != "" & wgs_bothtype_ctrl != "del/del"], aes(x = log(AdoCblPlus_01), y = as.factor(wgs_bothtype_ctrl), fill = as.factor(wgs_bothtype_ctrl), color = as.factor(wgs_bothtype_ctrl))) +
  geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
  scale_color_aaas() +
  scale_fill_aaas() +
  theme_test()

plt_AdoCblPlus_muttype_box <- 
ggplot(pheno_tbl[wgs_bothtype_ctrl != "" & wgs_bothtype_ctrl != "del/del"], aes(x = as.factor(wgs_bothtype_ctrl), y = log(AdoCblPlus_01), fill = as.factor(wgs_bothtype_ctrl), color = as.factor(wgs_bothtype_ctrl))) +
  geom_boxplot(alpha = 0.6) +
  scale_color_aaas() +
  scale_fill_aaas() +
  theme_test()



plt_pthwy_muttype <-
ggarrange(plt_pthwy_muttype_ridges, plt_pthwy_muttype_box, plt_AdoCblPlus_muttype_ridges, plt_AdoCblPlus_muttype_box)

ggsave(paste(fig_path,"pathwayAdoCblPlus_muttype.png", sep = ""), plt_pthwy_muttype, device = png(), width = 10, height = 6)
dev.off()


prot_pheno[, wgs_bothtype]

##################################
##################################
# use biochemical information to compare to MMUT transcript and protein levels



# correlation of MMUT transcript and protein with BIOCHEMISTRY (MMUT-deficient only)

MMUT_PI_act_plot <-
ggplot(prot_pheno[AdoCblPlus != 0 & sample_no < 151], aes(x = AdoCblPlus, y = OHCblPlus)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm", color = "darkblue", size = 0.5) +
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~','~"))) +
  ggtitle("MMUT_PI_act_MMUT") +
  # label x and y 
  theme_bw()


r_MMUT_PI_plot <- 
ggplot(rna_pheno[sample_no<151], aes(x = OHCblPlus, y = rna_exp)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_point(aes(color = mut_category), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(label.y.npc = "bottom", size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  ggtitle("r_MMUT_PI_MMUT") +
  theme_bw()

r_MMUT_act_plot <- 
ggplot(rna_pheno[AdoCblPlus != 0 & sample_no < 151], aes(x = AdoCblPlus, y = rna_exp)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_point(aes(color = mut_category), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(label.y.npc = "bottom", size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  ggtitle("r_MMUT_act_MMUT") +
  theme_bw()


p_MMUT_PI_plot <- 
ggplot(prot_pheno[sample_no<151], aes(x = OHCblPlus, y = prot_exp)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_point(aes(color = mut_category), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(label.y.npc = "bottom", size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  ggtitle("p_MMUT_PI_MMUT") +
  theme_bw()

p_MMUT_act_plot <- 
ggplot(prot_pheno[AdoCblPlus != 0 & sample_no < 151], aes(x = AdoCblPlus, y = prot_exp)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_point(aes(color = mut_category), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(label.y.npc = "bottom", size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  ggtitle("p_MMUT_act_MMUT") +
  theme_bw()




MMUT_lvl_biochem <- ggarrange(r_MMUT_PI_plot, r_MMUT_act_plot, MMUT_PI_act_plot, p_MMUT_PI_plot, p_MMUT_act_plot)

ggsave(paste(fig_path,"MMUT_lvl_biochem_MMUT.png", sep = ""), MMUT_lvl_biochem, device = png(), width = 280, height = 140, units = c("mm"), dpi = 600)
dev.off()




# correlation of MMUT transcript and protein with BIOCHEMISTRY (ALL samples)

MMUT_PI_act_plot <-
ggplot(prot_pheno[AdoCblPlus != 0], aes(x = AdoCblPlus, y = OHCblPlus)) +
  geom_point(aes(color = type, shape = type), alpha = 1) +
  geom_smooth(method = "lm", color = "black", size = 0.5) +
  stat_cor(size = 3, cor.coef.name = "rho", aes(label = paste(..r.label.., ..p.label.., sep = "~','~"))) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks(sides = "lb", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  ggtitle("MMUT_PI_act_ALL") +
  ylab("PI activity\n[pmol/mg protein/16 h]") +
  xlab("MMUT enzyme activity\n[pmol/min/mg protein]") +
  scale_color_manual(values=mypal2) +
  theme_pubr() +
  theme(plot.title = element_blank(), legend.position = "right", legend.title = element_blank())

ggsave(paste(fig_path,"MMUTenzymeActivityVSPIactivity.png", sep = ""), MMUT_PI_act_plot, device = png(), width = 5, height = 3)
dev.off()
ggsave(paste(fig_path_pdf,"SuppFig2_diagnosis/","MMUTenzymeActivityVSPIactivity.pdf", sep = ""), MMUT_PI_act_plot, device = "pdf", width = 5, height = 3)
dev.off()



r_MMUT_PI_plot <- 
ggplot(rna_pheno, aes(x = OHCblPlus, y = rna_exp)) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(sides = "lb", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  geom_point(aes(color = type, shape = type), alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(label.y.npc = "bottom", label.x.npc = 0.55, size = 3, cor.coef.name = "rho", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  scale_color_manual(values=mypal2) +
  ggtitle("r_MMUT_PI_ALL") +
  ylab("MMUT transcript level") +
  xlab("PI activity\n[pmol/mg protein/16 h]") +
  theme_pubr() +
  theme(plot.title = element_blank(), legend.position = "right", legend.title = element_blank())

r_MMUT_act_plot <- 
ggplot(rna_pheno[AdoCblPlus != 0], aes(x = AdoCblPlus, y = rna_exp)) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(sides = "lb", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  geom_point(aes(color = type, shape = type), alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(label.y.npc = "bottom", label.x.npc = .55, size = 3, cor.coef.name = "rho",aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  scale_color_manual(values=mypal2) +
  ggtitle("r_MMUT_act_ALL") +
  ylab("MMUT transcript level") +
  xlab("MMUT enzyme activity\n[pmol/min/mg protein]") +
  theme_pubr() +
  theme(plot.title = element_blank(), legend.position = "right", legend.title = element_blank())


p_MMUT_PI_plot <- 
ggplot(prot_pheno, aes(x = OHCblPlus, y = prot_exp)) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(sides = "lb", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  geom_point(aes(color = type, shape = type), alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(label.y.npc = "bottom", label.x.npc = .55, size = 3,cor.coef.name = "rho", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  scale_color_manual(values=mypal2) +
  ggtitle("p_MMUT_PI_ALL") +
  ylab("MMUT protein level") +
  xlab("PI activity\n[pmol/mg protein/16 h]") +
  theme_pubr() +
  theme(plot.title = element_blank(), legend.position = "right", legend.title = element_blank())

p_MMUT_act_plot <- 
ggplot(prot_pheno[AdoCblPlus != 0], aes(x = AdoCblPlus, y = prot_exp)) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(sides = "lb", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  geom_point(aes(color = type, shape = type), alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(label.y.npc = "bottom", label.x.npc = .55, size = 3, cor.coef.name = "rho",aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  scale_color_manual(values=mypal2) +
  ggtitle("p_MMUT_act_ALL") +
  ylab("MMUT protein level") +
  xlab("MMUT enzyme activity\n[pmol/min/mg protein]") +
  theme_pubr() +
  theme(plot.title = element_blank(), legend.position = "right", legend.title = element_blank())




MMUT_lvl_biochem <- ggarrange(r_MMUT_PI_plot, r_MMUT_act_plot, MMUT_PI_act_plot, p_MMUT_PI_plot, p_MMUT_act_plot)

ggsave(paste(fig_path,"MMUT_lvl_biochem_ALL.png", sep = ""), MMUT_lvl_biochem, device = png(), width = 280, height = 140, units = c("mm"), dpi = 600)
dev.off()



MMUT_lvl_biochem2 <- ggarrange(r_MMUT_PI_plot, p_MMUT_PI_plot, r_MMUT_act_plot, p_MMUT_act_plot, common.legend = TRUE, legend = "bottom")

ggsave(paste(fig_path,"MMUTenzymeActivityVSPIactivity2.png", sep = ""), MMUT_lvl_biochem2, device = png(), width = 7.1, height = 6, bg = "white")
dev.off()
ggsave(paste(fig_path_pdf,"SuppFig2_diagnosis/","MMUTenzymeActivityVSPIactivity2.pdf", sep = ""), MMUT_lvl_biochem2, device = "pdf", width = 7.1, height = 6, bg = "white")
dev.off()

# correlation of MMUT transcript and protein with BIOCHEMISTRY (nonMMUT only)


MMUT_PI_act_plot <-
ggplot(prot_pheno[AdoCblPlus != 0 & sample_no > 150], aes(x = AdoCblPlus, y = OHCblPlus)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_point(aes(color = type), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", size = 0.5) +
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~','~"))) +
  ggtitle("MMUT_PI_act_nonMMUT") +
  # label x and y 
  theme_bw()


r_MMUT_PI_plot <- 
ggplot(rna_pheno[sample_no > 150], aes(x = OHCblPlus, y = rna_exp)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_point(aes(color = type), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(label.y.npc = "bottom", size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  ggtitle("r_MMUT_PI_nonMMUT") +
  theme_bw()

r_MMUT_act_plot <- 
ggplot(rna_pheno[AdoCblPlus != 0 & sample_no > 150], aes(x = AdoCblPlus, y = rna_exp)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_point(aes(color = type), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(label.y.npc = "bottom", size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  ggtitle("r_MMUT_act_nonMMUT") +
  theme_bw()


p_MMUT_PI_plot <- 
ggplot(prot_pheno[sample_no > 150], aes(x = OHCblPlus, y = prot_exp)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_point(aes(color = type), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(label.y.npc = "bottom", size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  ggtitle("p_MMUT_PI_nonMMUT") +
  theme_bw()

p_MMUT_act_plot <- 
ggplot(prot_pheno[AdoCblPlus != 0 & sample_no > 150], aes(x = AdoCblPlus, y = prot_exp)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_point(aes(color = type), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor(label.y.npc = "bottom", size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  ggtitle("p_MMUT_act_nonMMUT") +
  theme_bw()




MMUT_lvl_biochem <- ggarrange(r_MMUT_PI_plot, r_MMUT_act_plot, MMUT_PI_act_plot, p_MMUT_PI_plot, p_MMUT_act_plot)

ggsave(paste(fig_path,"MMUT_lvl_biochem_nonMMUT.png", sep = ""), MMUT_lvl_biochem, device = png(), width = 280, height = 140, units = c("mm"), dpi = 600, bg = "white")
dev.off()



##################################
##################################
# pure BIOCHEMICAL plots (w/o transcript and protein info)

# build table with new types (MMUT divided into mut- and mut0)

act_MMUT_tbl <- prot_pheno %>%
  select("AdoCblMinus", "AdoCblPlus", "type")
act_MMUT_tbl <- act_MMUT_tbl[!is.na(act_MMUT_tbl$AdoCblPlus), ]
act_MMUT_tbl_melt <- melt(act_MMUT_tbl, id.vars = "type")
setnames(act_MMUT_tbl_melt, c("type", "AdoCbl", "MMUT_activity"))
act_MMUT_tbl_melt[AdoCbl == "AdoCblMinus", AdoCbl := "-"]
act_MMUT_tbl_melt[AdoCbl == "AdoCblPlus", AdoCbl := "+"]

act_MMUT_tbl_mut <- prot_pheno %>%
  select("AdoCblMinus", "AdoCblPlus", "type_MMUTmut")
act_MMUT_tbl_mut <- act_MMUT_tbl_mut[!is.na(act_MMUT_tbl_mut$AdoCblPlus), ]
act_MMUT_tbl_mut_melt <- melt(act_MMUT_tbl_mut, id.vars = "type_MMUTmut")
setnames(act_MMUT_tbl_mut_melt, c("type_MMUTmut", "AdoCbl", "MMUT_activity"))
act_MMUT_tbl_mut_melt[AdoCbl == "AdoCblMinus", AdoCbl := "-"]
act_MMUT_tbl_mut_melt[AdoCbl == "AdoCblPlus", AdoCbl := "+"]


# build table for PI activity plots

act_PI_tbl <- prot_pheno %>%
  select("OHCblMinus", "OHCblPlus", "type_MMUTmut")
act_PI_tbl <- act_PI_tbl[!is.na(act_PI_tbl$OHCblPlus), ]
act_PI_tbl_melt <- melt(act_PI_tbl, id.vars = "type_MMUTmut")
setnames(act_PI_tbl_melt, c("type_MMUTmut", "OHCbl", "Prop_incorp"))
act_PI_tbl_melt[OHCbl == "OHCblMinus", OHCbl := "-"]
act_PI_tbl_melt[OHCbl == "OHCblPlus", OHCbl := "+"]


# plot MMUT activity according to different types

act_MMUT_mut_plot_ALL <- 
ggplot(act_MMUT_tbl_mut_melt, aes(x = AdoCbl, y = MMUT_activity, color = type_MMUTmut, fill = type_MMUTmut)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10() +
  annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  facet_grid(~type_MMUTmut) +
  ggtitle("act_MMUT_mut") +
  ylab("MMUT enzyme activity\n[pmol/min/mg protein]") +
  scale_color_manual(values = mypal3) +
  scale_fill_manual(values = mypal3) +
  theme_pubr() +
  theme(legend.position = "none", plot.title = element_blank(), strip.background = element_rect(fill = "white"))

ggsave(paste(fig_path,"MMUTenzymeActivityPlusAndMinusAdo.png", sep = ""), act_MMUT_mut_plot_ALL, device = png(), width = 5, height = 3)
dev.off()
ggsave(paste(fig_path_pdf,"SuppFig2_diagnosis/","MMUTenzymeActivityPlusAndMinusAdo.pdf", sep = ""), act_MMUT_mut_plot_ALL, device = "pdf", width = 5, height = 3)
dev.off()


# plot MMUT activity (AdoCblPlus only) according to different types

act_MMUT_mut_plot_ALL_AdoCblPlus <- 
ggplot(act_MMUT_tbl_mut_melt[AdoCbl == "+" & MMUT_activity > 0], aes(x = type_MMUTmut, y = MMUT_activity, color = type_MMUTmut, fill = type_MMUTmut)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(alpha = 0.5, width = 0.18) +
  scale_y_log10() +
  ggtitle("act_MMUT_mut_Ado+") +
  theme_bw() +
  scale_color_npg() +
  scale_fill_npg()

act_MMUT_plot_ALL_AdoCblPlus <- 
ggplot(act_MMUT_tbl_melt[AdoCbl == "+" & MMUT_activity > 0], aes(x = type, y = MMUT_activity, color = type, fill = type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(alpha = 1, width = 0.18, shape =21, color = "white") +
  scale_y_log10() +
  annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  stat_compare_means(comparisons = list(c("MMUT def.", "unknown"), c("MMUT def.", "unaffected"), c("unknown", "unaffected")), method = "wilcox.test", size = 3, step.increase = 0.15) +
  ggtitle("act_MMUT_Ado+") +
  ylab("MMUT enzyme activity\n[pmol/min/mg protein]") +
  xlab("") +
  scale_color_manual(values = mypal2) +
  scale_fill_manual(values = mypal2) +
  theme_pubr() +
  theme(legend.title = element_blank(), plot.title = element_blank(), legend.position = c(0.7, 0.2), axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(paste(fig_path,"MMUTenzymeActivity.png", sep = ""), act_MMUT_plot_ALL_AdoCblPlus, device = png(), width = 2.7, height = 3)
dev.off()
ggsave(paste(fig_path_pdf,"Fig1_diagnosis/","MMUTenzymeActivity.pdf", sep = ""), act_MMUT_plot_ALL_AdoCblPlus, device = "pdf", width = 2.7, height = 3)
dev.off()




# plot PI activity according to different types

act_PI_plot_ALL <- 
ggplot(act_PI_tbl_melt, aes(x = OHCbl, y = Prop_incorp, color = type_MMUTmut, fill = type_MMUTmut)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10() +
  annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  facet_grid(~type_MMUTmut) +
  ggtitle("act_PI") +
  ylab("PI activity\n[pmol/mg protein/16 h]") +
  scale_color_manual(values = mypal3) +
  scale_fill_manual(values = mypal3) +
  theme_pubr() +
  theme(legend.position = "none", plot.title = element_blank(), strip.background = element_rect(fill = "white"))

ggsave(paste(fig_path,"PIActivity.png", sep = ""), act_PI_plot_ALL, device = png(), width = 5, height = 3)
dev.off()
ggsave(paste(fig_path_pdf,"SuppFig2_diagnosis/","PIActivity.pdf", sep = ""), act_PI_plot_ALL, device = "pdf", width = 5, height = 3)
dev.off()



# PI ratio plot as above, this time with density

pre_mut0vsMinus_plot <-
ggplot(prot_pheno[sample_no<151], aes(x = log(ratio), color = mut_category, fill = mut_category)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = log(1.5)) +
  ggtitle("mut0 vs mut-, PI ratio") +
  theme_bw() +
  scale_color_npg() +
  scale_fill_npg()  
xbox_p <- axis_canvas(pre_mut0vsMinus_plot, axis = "x") +
  geom_boxplot(data = prot_pheno[sample_no<151], aes(x = log(ratio), color = mut_category, fill = mut_category, alpha = 0.6)) +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  scale_color_npg() +
  scale_fill_npg()

mut0vsMinus_plot <- insert_xaxis_grob(pre_mut0vsMinus_plot, xbox_p, position = "top")
mut0vsMinus_plot <- ggdraw(mut0vsMinus_plot)



act_MMUTandPI_plot <- ggarrange(act_PI_plot_ALL, act_MMUT_mut_plot_ALL_AdoCblPlus, mut0vsMinus_plot, act_MMUT_mut_plot_ALL, act_MMUT_plot_ALL_AdoCblPlus, nrow = 2, ncol = 3, widths = c(1.5, 1, 1.5), legend = "none")

ggsave(paste(fig_path,"act_MMUTandPI_ALL.png", sep = ""), act_MMUTandPI_plot, device = png(), width = 250, height = 120, units = c("mm"), dpi = 600)
dev.off()















