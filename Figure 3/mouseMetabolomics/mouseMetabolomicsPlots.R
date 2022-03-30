# script to generate mouse metaobloimcs plots
# run PAEM script first
# next run ApplyPAEM script next

# libraries
require(ggplot2)
require(openxlsx)
require(stringr)
require(data.table)
require(ggpubr)
library(matrixTests)
library(ggsci)
library(patchwork)


# extract the list of datasets and annotate the modes manually
fls <- list.files("Data/mouse/mouseMetabolomics/diffreports_kokikiwt/")
# pos_fls <- c(2,3,5,7,9,11)#koki vs wt
# neg_fls <- c(1,4,6,8,10,12)# koki vs wt
pos_fls <- c(2,4,6,8,10,12)#koki vs kiwt
neg_fls <- c(1,3,5,7,9,11)# koki vs kiwt



###########################################################
###########################################################
# plots


###########################################################
# create plot based on TCA cycle p-values in tstatpositive, negative mode

datalist = list()

for (q in neg_fls) {
  cat(q, "\n")
  tbl <- read.csv(paste0("Results/mouseMetabolomics/enrichedPathways/enrichedPathwaysTstatPos_negMode_",gsub("diffreport_MMA_","",fls[q])))
  tbl$tissue <- gsub("diffreport_MMA_","",fls[q])
  datalist[[q]] <- tbl
}

combolist = do.call(rbind, datalist)
combolist <- data.table(combolist)
combolist$tissue2 <- gsub("_.*","", combolist$tissue)


TCA_enrich_plt <- 
ggplot(combolist[Pathway == "Citrate cycle (TCA cycle)", ], aes(x=-log10(pvalue), y=reorder(tissue2, pvalue))) +
  geom_bar(stat = "identity", fill = "grey", width = 0.8) +
  geom_vline(xintercept = -log10(0.05), linetype = 2) +
  ggtitle("TCA cycle enrichment in mouse") +
  xlab("-log10(p-value)") +
  ylab("Mouse tisse") +
  theme_pubr() +
  theme(plot.title = element_text(size=12), axis.title.y = element_blank(), axis.text.y = element_text(size=10))

ggsave("Figs/mouseMetabolomics/TCAenrichmentInTissues.pdf", TCA_enrich_plt, device = "pdf", width = 4.5, height = 2.5)
dev.off()

ggsave("Figs/mouseMetabolomics/TCAenrichmentInTissues.png", TCA_enrich_plt, device = png(), width = 4.5, height = 2.5)
dev.off()


###########################################################
# plot methylmalonate as sanity check

# table for negative mode
datalist0 = list()
for (q in neg_fls) {
  cat(q, "\n")
  tbl <- read.csv(paste0("Results/mouseMetabolomics/tstatpos/",gsub(".csv","",gsub("diffreport_MMA_","",fls[q])), "_annotation_accuracy_0.01.csv"))
  tbl$tissue <- gsub("_.*","",gsub("diffreport_MMA_","",fls[q]))
  datalist0[[q]] <- tbl
}

summary(datalist0)
datalist2 <- datalist0[c(1,3,5,7,9,11)] # if using negative mode
summary(datalist2)
old_names <- colnames(datalist2[[1]])
new_names <- c(old_names[1:15], paste0("kiwt_",seq(1,15)), paste0("koki_",seq(1,15)), old_names[46:48])
dt_neg <- do.call(rbind, lapply(datalist2, setNames, new_names))
dt_neg$mode <- "negative"


# table for positive mode
datalist0 = list()
for (q in pos_fls) {
  cat(q, "\n")
  tbl <- read.csv(paste0("Results/mouseMetabolomics/tstatpos/",gsub(".csv","",gsub("diffreport_MMA_","",fls[q])), "_annotation_accuracy_0.01.csv"))
  tbl$tissue <- gsub("_.*","",gsub("diffreport_MMA_","",fls[q]))
  datalist0[[q]] <- tbl
}

summary(datalist0)
datalist2 <- datalist0[c(2, 4, 6, 8, 10,12)] # if using positive mode
summary(datalist2)
old_names <- colnames(datalist2[[1]])
new_names <- c(old_names[1:15], paste0("kiwt_",seq(1,15)), paste0("koki_",seq(1,15)), old_names[46:48])
dt_pos <- do.call(rbind, lapply(datalist2, setNames, new_names))
dt_pos$mode <- "positive"


dt <- rbind(dt_neg, dt_pos)
colnames(dt)


dt[, 16:45][dt[,16:45] == 0] <- NA
dt1 <- data.table(dt)
unique(dt1$tissue)

# calculate row means per type (wt vs mutant)
dt1[, mean_kiwt := rowMeans(dt1[,16:30], na.rm = TRUE)]
dt1[, mean_koki := rowMeans(dt1[,31:45], na.rm = TRUE)]

dt2 <- dt1[!is.na(dt1$mean_koki), ]
dt2 <- dt2[!is.na(dt2$mean_kiwt), ]


# calculate stats
koki_names <- names(dt2[, 31:45])
wt_names <- names(dt2[, 16:30])
stats_tbl <- data.table(row_t_welch(dt2[, koki_names, with = FALSE], dt2[, wt_names, with = FALSE]))
# compare pvalues
data.table(stats_tbl$pvalue, dt2$pvalue)

# calculate fold change and compare
dt2[, fold2 := mean_koki/mean_kiwt]
data.table(dt2$fold, dt2$fold2)


# normalize to kiwt
ff=lapply(c(1:dim(dt2)[1]),function(i){
  print(i)
  dt2[i, c(wt_names, koki_names), with = FALSE]/dt2[i, ]$mean_kiwt
})

ff2=do.call("rbind",ff)
ff3=data.table(ff2)
norm_wt_names <- paste0("norm_",wt_names)
norm_koki_names <- paste0("norm_",koki_names)
colnames(ff3)
setnames(ff3, c(norm_wt_names, norm_koki_names))

dt3 <- cbind(dt2, ff3)
dt3_neg <- dt3[mode == "negative", ]
dt3_pos <- dt3[mode == "positive", ]


# labels for plots
wt_lbl <- "Mmut-ki/wt"
koki_lbl <- "Mmut-ko/ki"


# identify MMA (KEGG ID C02170) in negative mode
dtMMA <- dt3_neg[grep("C02170", dt3_neg$KEGGID), ]
dtMMA <- dtMMA[rtmed <60, ]
colnames(dtMMA)
dtMMA1 <- dtMMA[, c(norm_wt_names, norm_koki_names, "tissue"), with = FALSE]
dtMMA2 <- melt.data.table(dtMMA1, id.vars = c("tissue"))
dtMMA2[grep("norm_kiwt", dtMMA2$variable), type := wt_lbl]
dtMMA2[grep("norm_koki", dtMMA2$variable), type := koki_lbl]
dtMMA2$type=factor(dtMMA2$type,levels=c(wt_lbl, koki_lbl))


# comparisons and y limits:
compare <- list(c(wt_lbl, koki_lbl))
max_val <- max(dtMMA2$value, na.rm = TRUE)
min_val <- min(dtMMA2$value, na.rm = TRUE)


# plot
mma_plt <- 
ggplot(dtMMA2, aes(x=type, y=value, color = type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2,color = "black", fill = "white", shape = 21, show.legend = FALSE) +
  scale_y_log10(limits = c(min_val,2*max_val)) +
  stat_compare_means(comparisons = compare, size = 3, method = "t.test", label = "p.format") +
  facet_grid(.~tissue, scales = "free_y") +
  labs(title="Methylmalonic acid levels",x ="", y = "Ion abundance relative to control",fill = '') +
  # annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  scale_color_aaas() +
  scale_fill_aaas() +
  theme_pubr() +
  theme(strip.background = element_rect("white"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(size=12))

ggsave("Figs/mouseMetabolomics/MMAInmouseTissues.pdf", mma_plt, device = "pdf", width = 5, height = 3.5)
dev.off()

ggsave("Figs/mouseMetabolomics/MMAInmouseTissues.png", mma_plt, device = png(), width = 5, height = 3.5)
dev.off()



# identify propionylcarnitine (KEGG ID C02170) in positive mode
dtpropcar <- dt3_pos[grep("C03017", dt3_pos$KEGGID), ]
# select ions with similar retention time
dtpropcar <- dtpropcar[rtmed > 400 & rtmed < 500, ]
dtpropcar <- dtpropcar[1:6,]
colnames(dtpropcar)
dtpropcar1 <- dtpropcar[, c(norm_wt_names, norm_koki_names, "tissue"), with = FALSE]
dtpropcar2 <- melt.data.table(dtpropcar1, id.vars = c("tissue"))
dtpropcar2[grep("norm_kiwt", dtpropcar2$variable), type := wt_lbl]
dtpropcar2[grep("norm_koki", dtpropcar2$variable), type := koki_lbl]
dtpropcar2$type=factor(dtpropcar2$type,levels=c(wt_lbl, koki_lbl))


# comparisons and y limits:
compare <- list(c(wt_lbl, koki_lbl))
max_val <- max(dtpropcar2$value, na.rm = TRUE)
min_val <- min(dtpropcar2$value, na.rm = TRUE)


# plot
propcar_plt <- 
ggplot(dtpropcar2, aes(x=type, y=value, color = type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2,color = "black", fill = "white", shape = 21, show.legend = FALSE) +
  ylim(0, 1.1 * max_val) +
  stat_compare_means(comparisons = compare, size = 3, method = "t.test", label = "p.format") +
  facet_grid(.~tissue, scales = "free_y") +
  labs(title="Propionylcarntine levels",x ="", y = "Ion abundance relative to control",fill = '') +
  # annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  scale_color_aaas() +
  scale_fill_aaas() +
  theme_pubr() +
  theme(strip.background = element_rect("white"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(size=12))

ggsave("Figs/mouseMetabolomics/PropCarInmouseTissues.pdf", propcar_plt, device = "pdf", width = 5, height = 3.5)
dev.off()

ggsave("Figs/mouseMetabolomics/PropCarInmouseTissues.png", propcar_plt, device = png(), width = 5, height = 3.5)
dev.off()





# combine MMA and TCA enrichment plot


mouse_metab_plt <- 
mma_plt + TCA_enrich_plt + plot_layout(widths = c(2, 1.3))

ggsave("Figs/v8/pdf/Fig3/mouse_metabolomics_combo.pdf", mouse_metab_plt, device = "pdf", width = 7.3, height = 3.5)
dev.off()


# identify methylmalonate semialdehyde (KEGG ID C06002) in positive mode??




