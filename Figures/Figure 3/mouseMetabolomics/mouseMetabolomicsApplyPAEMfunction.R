# script to apply PAEM function developed by Sarah on multiple datasets
# script to run in PF directory

# libraries
require(ggplot2)
require(openxlsx)
require(stringr)
require(data.table)
require(ggpubr)
library(matrixTests)
library(ggsci)


# extract the list of datasets and annotate the modes manually
fls <- list.files("Data/mouseMetabolomics/diffreports_kokikiwt/")
# pos_fls <- c(2,3,5,7,9,11)#koki vs wt
# neg_fls <- c(1,4,6,8,10,12)# koki vs wt
pos_fls <- c(2,4,6,8,10,12)#koki vs kiwt
neg_fls <- c(1,3,5,7,9,11)# koki vs kiwt
file_path <- c("Data/mouseMetabolomics/diffreports_kokikiwt/")



# run for loops using the modified PAEM function

# sort metabolites according to pvals

for (q in pos_fls) {
  cat(q, "\n")
  pathlist <- PAEM(file = paste0(file_path,fls[q]), mode = "positive", accuracy = 0.01, metabolitesDatabase='Data/mouseMetabolomics/helperFiles/mref_kegg.xlsx',pathwayDatabase='Data/mouseMetabolomics/helperFiles/hsa_HumanSpecificMetab.csv',exportAnnotation=TRUE,exportPath='Results/mouseMetabolomics/pval/',sortList='pvalue')
  pos_sign=pathlist$pvalue<0.05
  df_sign=data.frame(x = -log10(pathlist$pvalue[pos_sign]),y = pathlist$Pathway[pos_sign])
  barplot_sign = ggplot(df_sign, aes(x,reorder(y, x))) + 
  geom_bar(stat="identity", width=0.5)+
  labs(title = paste0("KEGG pathways: ",gsub("diffreport_MMA_","",gsub(".csv","",fls[q]))),x = "-log10(p-value)")+
  theme_minimal() +
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.title.x = element_text(size=10))
  ggsave(filename=paste0('Figs/mouseMetabolomics/','TopKEGG_pval_',gsub("diffreport_MMA_","",gsub(".csv","",fls[q])),'.pdf'), barplot_sign, device = "pdf", width = 7, height = 4)
}

for (q in neg_fls) {
  cat(q, "\n")
  pathlist <- PAEM(file = paste0(file_path,fls[q]), mode = "negative", accuracy = 0.01, metabolitesDatabase='Data/mouseMetabolomics/helperFiles/mref_kegg.xlsx',pathwayDatabase='Data/mouseMetabolomics/helperFiles/hsa_HumanSpecificMetab.csv',exportAnnotation=TRUE,exportPath='Results/mouseMetabolomics/pval/',sortList='pvalue')
  pos_sign=pathlist$pvalue<0.05
  df_sign=data.frame(x = -log10(pathlist$pvalue[pos_sign]),y = pathlist$Pathway[pos_sign])
  barplot_sign = ggplot(df_sign, aes(x,reorder(y, x))) + 
  geom_bar(stat="identity", width=0.5)+
  labs(title = paste0("KEGG pathways: ",gsub("diffreport_MMA_","",gsub(".csv","",fls[q]))),x = "-log10(p-value)")+
  theme_minimal() +
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.title.x = element_text(size=10))
  ggsave(filename=paste0('Figs/mouseMetabolomics/','TopKEGG_pval_',gsub("diffreport_MMA_","",gsub(".csv","",fls[q])),'.pdf'), barplot_sign, device = "pdf", width = 7, height = 4)
}


# sort metabolites according to fold

for (q in pos_fls) {
  cat(q, "\n")
  pathlist <- PAEM(file = paste0(file_path,fls[q]), mode = "positive", accuracy = 0.01, metabolitesDatabase='Data/mouseMetabolomics/helperFiles/mref_kegg.xlsx',pathwayDatabase='Data/mouseMetabolomics/helperFiles/hsa_HumanSpecificMetab.csv',exportAnnotation=TRUE,exportPath='Results/mouseMetabolomics/fold/',sortList='fold')
  pos_sign=pathlist$pvalue<0.05
  df_sign=data.frame(x = -log10(pathlist$pvalue[pos_sign]),y = pathlist$Pathway[pos_sign])
  barplot_sign = ggplot(df_sign, aes(x,reorder(y, x))) + 
  geom_bar(stat="identity", width=0.5)+
  labs(title = paste0("KEGG pathways: ",gsub("diffreport_MMA_","",gsub(".csv","",fls[q]))),x = "-log10(p-value)")+
  theme_minimal() +
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.title.x = element_text(size=10))
  ggsave(filename=paste0('Figs/mouseMetabolomics/','TopKEGG_fold_',gsub("diffreport_MMA_","",gsub(".csv","",fls[q])),'.pdf'), barplot_sign, device = "pdf", width = 7, height = 4)
}

for (q in neg_fls) {
  cat(q, "\n")
  pathlist <- PAEM(file = paste0(file_path,fls[q]), mode = "negative", accuracy = 0.01, metabolitesDatabase='Data/mouseMetabolomics/helperFiles/mref_kegg.xlsx',pathwayDatabase='Data/mouseMetabolomics/helperFiles/hsa_HumanSpecificMetab.csv',exportAnnotation=TRUE,exportPath='Results/mouseMetabolomics/fold/',sortList='fold')
  pos_sign=pathlist$pvalue<0.05
  df_sign=data.frame(x = -log10(pathlist$pvalue[pos_sign]),y = pathlist$Pathway[pos_sign])
  barplot_sign = ggplot(df_sign, aes(x,reorder(y, x))) + 
  geom_bar(stat="identity", width=0.5)+
  labs(title = paste0("KEGG pathways: ",gsub("diffreport_MMA_","",gsub(".csv","",fls[q]))),x = "-log10(p-value)")+
  theme_minimal() +
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.title.x = element_text(size=10))
  ggsave(filename=paste0('Figs/mouseMetabolomics/','TopKEGG_fold_',gsub("diffreport_MMA_","",gsub(".csv","",fls[q])),'.pdf'), barplot_sign, device = "pdf", width = 7, height = 4)
}


# sort metabolites according to tstat-positive

for (q in pos_fls) {
  cat(q, "\n")
  pathlist <- PAEM(file = paste0(file_path,fls[q]), mode = "positive", accuracy = 0.01, metabolitesDatabase='Data/mouseMetabolomics/helperFiles/mref_kegg.xlsx',pathwayDatabase='Data/mouseMetabolomics/helperFiles/hsa_HumanSpecificMetab.csv',exportAnnotation=TRUE,exportPath='Results/mouseMetabolomics/tstatpos/',sortList='tstat-positive')
  write.csv(pathlist, paste0("Results/mouseMetabolomics/enrichedPathways/enrichedPathwaysTstatPos_posMode_",gsub("diffreport_MMA_","",fls[q])))
  pos_sign=pathlist$pvalue<0.05
  df_sign=data.frame(x = -log10(pathlist$pvalue[pos_sign]),y = pathlist$Pathway[pos_sign])
  barplot_sign = ggplot(df_sign, aes(x,reorder(y, x))) + 
  geom_bar(stat="identity", width=0.5)+
  labs(title = paste0("KEGG pathways: ",gsub("diffreport_MMA_","",gsub(".csv","",fls[q]))),x = "-log10(p-value)")+
  theme_minimal() +
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.title.x = element_text(size=10))
  ggsave(filename=paste0('Figs/mouseMetabolomics/','TopKEGG_tstatPos_',gsub("diffreport_MMA_","",gsub(".csv","",fls[q])),'.pdf'), barplot_sign, device = "pdf", width = 7, height = 4)
}

for (q in neg_fls) {
  cat(q, "\n")
  pathlist <- PAEM(file = paste0(file_path,fls[q]), mode = "negative", accuracy = 0.01, metabolitesDatabase='Data/mouseMetabolomics/helperFiles/mref_kegg.xlsx',pathwayDatabase='Data/mouseMetabolomics/helperFiles/hsa_HumanSpecificMetab.csv',exportAnnotation=TRUE,exportPath='Results/mouseMetabolomics/tstatpos/',sortList='tstat-positive')
  write.csv(pathlist, paste0("Results/mouseMetabolomics/enrichedPathways/enrichedPathwaysTstatPos_negMode_",gsub("diffreport_MMA_","",fls[q])))
  pos_sign=pathlist$pvalue<0.05
  df_sign=data.frame(x = -log10(pathlist$pvalue[pos_sign]),y = pathlist$Pathway[pos_sign])
  barplot_sign = ggplot(df_sign, aes(x,reorder(y, x))) + 
  geom_bar(stat="identity", width=0.5)+
  labs(title = paste0("KEGG pathways: ",gsub("diffreport_MMA_","",gsub(".csv","",fls[q]))),x = "-log10(p-value)")+
  theme_minimal() +
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.title.x = element_text(size=10))
  ggsave(filename=paste0('Figs/mouseMetabolomics/','TopKEGG_tstatPos_',gsub("diffreport_MMA_","",gsub(".csv","",fls[q])),'.pdf'), barplot_sign, device = "pdf", width = 7, height = 4)
}


# sort metabolites according to tstat-negative

for (q in pos_fls) {
  cat(q, "\n")
  pathlist <- PAEM(file = paste0(file_path,fls[q]), mode = "positive", accuracy = 0.01, metabolitesDatabase='Data/mouseMetabolomics/helperFiles/mref_kegg.xlsx',pathwayDatabase='Data/mouseMetabolomics/helperFiles/hsa_HumanSpecificMetab.csv',exportAnnotation=TRUE,exportPath='Results/mouseMetabolomics/tstatneg/',sortList='tstat-negative')
  write.csv(pathlist, paste0("Results/mouseMetabolomics/enrichedPathways/enrichedPathwaysTstatNeg_posMode_",gsub("diffreport_MMA_","",fls[q])))
  pos_sign=pathlist$pvalue<0.05
  df_sign=data.frame(x = -log10(pathlist$pvalue[pos_sign]),y = pathlist$Pathway[pos_sign])
  barplot_sign = ggplot(df_sign, aes(x,reorder(y, x))) + 
  geom_bar(stat="identity", width=0.5)+
  labs(title = paste0("KEGG pathways: ",gsub("diffreport_MMA_","",gsub(".csv","",fls[q]))),x = "-log10(p-value)")+
  theme_minimal() +
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.title.x = element_text(size=10))
  ggsave(filename=paste0('Figs/mouseMetabolomics/','TopKEGG_tstatNeg_',gsub("diffreport_MMA_","",gsub(".csv","",fls[q])),'.pdf'), barplot_sign, device = "pdf", width = 7, height = 4)
}
  
for (q in neg_fls) {
  cat(q, "\n")
  pathlist <- PAEM(file = paste0(file_path,fls[q]), mode = "negative", accuracy = 0.01, metabolitesDatabase='Data/mouseMetabolomics/helperFiles/mref_kegg.xlsx',pathwayDatabase='Data/mouseMetabolomics/helperFiles/hsa_HumanSpecificMetab.csv',exportAnnotation=TRUE,exportPath='Results/mouseMetabolomics/tstatneg/',sortList='tstat-negative')
  write.csv(pathlist, paste0("Results/mouseMetabolomics/enrichedPathways/enrichedPathwaysTstatNeg_negMode_",gsub("diffreport_MMA_","",fls[q])))
  pos_sign=pathlist$pvalue<0.05
  df_sign=data.frame(x = -log10(pathlist$pvalue[pos_sign]),y = pathlist$Pathway[pos_sign])
  barplot_sign = ggplot(df_sign, aes(x,reorder(y, x))) + 
  geom_bar(stat="identity", width=0.5)+
  labs(title = paste0("KEGG pathways: ",gsub("diffreport_MMA_","",gsub(".csv","",fls[q]))),x = "-log10(p-value)")+
  theme_minimal() +
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.title.x = element_text(size=10))
  ggsave(filename=paste0('Figs/mouseMetabolomics/','TopKEGG_tstatNeg_',gsub("diffreport_MMA_","",gsub(".csv","",fls[q])),'.pdf'), barplot_sign, device = "pdf", width = 7, height = 4)
}


# sort metabolites according to tstat-absolute

for (q in pos_fls) {
  cat(q, "\n")
  pathlist <- PAEM(file = paste0(file_path,fls[q]), mode = "positive", accuracy = 0.01, metabolitesDatabase='Data/mouseMetabolomics/helperFiles/mref_kegg.xlsx',pathwayDatabase='Data/mouseMetabolomics/helperFiles/hsa_HumanSpecificMetab.csv',exportAnnotation=TRUE,exportPath='Results/mouseMetabolomics/tstatabs/',sortList='tstat-absolute')
  pos_sign=pathlist$pvalue<0.05
  df_sign=data.frame(x = -log10(pathlist$pvalue[pos_sign]),y = pathlist$Pathway[pos_sign])
  barplot_sign = ggplot(df_sign, aes(x,reorder(y, x))) + 
  geom_bar(stat="identity", width=0.5)+
  labs(title = paste0("KEGG pathways: ",gsub("diffreport_MMA_","",gsub(".csv","",fls[q]))),x = "-log10(p-value)")+
  theme_minimal() +
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.title.x = element_text(size=10))
  ggsave(filename=paste0('Figs/mouseMetabolomics/','TopKEGG_tstatAbs_',gsub("diffreport_MMA_","",gsub(".csv","",fls[q])),'.pdf'), barplot_sign, device = "pdf", width = 7, height = 4)
}

for (q in neg_fls) {
  cat(q, "\n")
  pathlist <- PAEM(file = paste0(file_path,fls[q]), mode = "negative", accuracy = 0.01, metabolitesDatabase='Data/mouseMetabolomics/helperFiles/mref_kegg.xlsx',pathwayDatabase='Data/mouseMetabolomics/helperFiles/hsa_HumanSpecificMetab.csv',exportAnnotation=TRUE,exportPath='Results/mouseMetabolomics/tstatabs/',sortList='tstat-absolute')
  pos_sign=pathlist$pvalue<0.05
  df_sign=data.frame(x = -log10(pathlist$pvalue[pos_sign]),y = pathlist$Pathway[pos_sign])
  barplot_sign = ggplot(df_sign, aes(x,reorder(y, x))) + 
  geom_bar(stat="identity", width=0.5)+
  labs(title = paste0("KEGG pathways: ",gsub("diffreport_MMA_","",gsub(".csv","",fls[q]))),x = "-log10(p-value)")+
  theme_minimal() +
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.title.x = element_text(size=10))
  ggsave(filename=paste0('Figs/mouseMetabolomics/','TopKEGG_tstatAbs_',gsub("diffreport_MMA_","",gsub(".csv","",fls[q])),'.pdf'), barplot_sign, device = "pdf", width = 7, height = 4)
}






# create plot based on TCA cycle p-values in tstatnegative, positive mode

datalist = list()

for (q in pos_fls) {
  cat(q, "\n")
  tbl <- read.csv(paste0("Results/mouseMetabolomics/enrichedPathways/enrichedPathwaysTstatNeg_posMode",gsub("diffreport_MMA_","",fls[q])))
  tbl$tissue <- gsub("diffreport_MMA_","",fls[q])
  datalist[[q]] <- tbl
}

combolist = do.call(rbind, datalist)
combolist <- data.table(combolist)
combolist$tissue2 <- gsub("_.*","", combolist$tissue)
combolist[tissue2 == "Brain"]$tissue2 <- "brain"


ggplot(combolist[Pathway == "Citrate cycle (TCA cycle)", ], aes(x=-log10(pvalue), y=reorder(tissue2, pvalue))) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = -log10(0.05), linetype = 2) +
  ggtitle("TCA cycle pathway enrichment in mouse tissues") +
  xlab("-log10(p value)") +
  ylab("Mouse tisse") +
  theme_pubr() +
  theme(plot.title = element_text(size=12))



# show methylmalonate as sanity check

datalist0 = list()

for (q in neg_fls) {
  cat(q, "\n")
  tbl <- read.csv(paste0("Results/mouseMetabolomics/tstatneg/",gsub(".csv","",gsub("diffreport_MMA_","",fls[q])), "_annotation_accuracy_0.01.csv"))
  tbl$tissue <- gsub("_.*","",gsub("diffreport_MMA_","",fls[q]))
  datalist0[[q]] <- tbl
}

summary(datalist0)
# datalist2 <- datalist0[-c(1, 4, 6, 8, 10)] # if using positive mode
datalist2 <- datalist0[-c(2,3,5,7,9,11)] # if using negative mode
summary(datalist2)
old_names <- colnames(datalist2[[1]])
new_names <- c(old_names[1:15], paste0("koki_",seq(1,15)), paste0("wt_",seq(1,15)), old_names[46:48])
dt <- do.call(rbind, lapply(datalist2, setNames, new_names))
colnames(dt)


dt[, 16:45][dt[,16:45] == 0] <- NA
dt1 <- data.table(dt)
unique(dt1$tissue)

# calculate row means per type (wt vs mutant)
dt1[, mean_koki := rowMeans(dt1[,16:30], na.rm = TRUE)]
dt1[, mean_wt := rowMeans(dt1[,31:45], na.rm = TRUE)]

dt2 <- dt1[!is.na(dt1$mean_koki), ]
dt2 <- dt2[!is.na(dt2$mean_wt), ]


# calculate stats
wt_names <- names(dt2[, 31:45])
koki_names <- names(dt2[, 16:30])
stats_tbl <- data.table(row_t_welch(dt2[, koki_names, with = FALSE], dt2[, wt_names, with = FALSE]))
# compare pvalues
data.table(stats_tbl$pvalue, dt2$pvalue)

# calculate fold change and compare
dt2[, fold2 := mean_koki/mean_wt]
data.table(dt2$fold, dt2$fold2)



dt2[1, c(wt_names, koki_names), with = FALSE]/dt2[1, ]$mean_wt

# normalize to wt
ff=lapply(c(1:dim(dt2)[1]),function(i){
  print(i)
  dt2[i, c(wt_names, koki_names), with = FALSE]/dt2[i, ]$mean_wt
})

ff2=do.call("rbind",ff)
ff3=data.table(ff2)
norm_wt_names <- paste0("norm_",wt_names)
norm_koki_names <- paste0("norm_",koki_names)
colnames(ff3)
setnames(ff3, c(norm_wt_names, norm_koki_names))

dt3 <- cbind(dt2, ff3)

# identify MMA (KEGG ID C02170)
wt_lbl <- "wildtype"
koki_lbl <- "MMUT-ko/ki"
dtMMA <- dt3[grep("C02170", dt3$KEGGID), ]
colnames(dtMMA)
dtMMA1 <- dtMMA[, c(norm_wt_names, norm_koki_names, "tissue"), with = FALSE]
dtMMA2 <- melt.data.table(dtMMA1, id.vars = c("tissue"))
dtMMA2[grep("norm_wt", dtMMA2$variable), type := wt_lbl]
dtMMA2[grep("norm_koki", dtMMA2$variable), type := koki_lbl]
dtMMA2$type=factor(dtMMA2$type,levels=c(wt_lbl, koki_lbl))

# comparisons:
compare <- list(c(wt_lbl, koki_lbl))
max_val <- max(dtMMA2$value, na.rm = TRUE)
min_val <- min(dtMMA2$value, na.rm = TRUE)


# plot

ggplot(dtMMA2, aes(x=type, y=value, color = type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, aes(fill=type),color = "white", shape = 21, show.legend = FALSE) +
  scale_y_log10(limits = c(min_val,2*max_val)) +
  stat_compare_means(comparisons = compare, size = 3, method = "t.test", label = "p.format") +
  facet_grid(.~tissue, scales = "free_y") +
  labs(title="Methylmalonic acid levels",x ="", y = "Ion abundance relative to wildtype",fill = '') +
  # annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  scale_color_aaas() +
  scale_fill_aaas() +
  theme_pubr() +
  theme(strip.background = element_rect("white"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(size=12))


