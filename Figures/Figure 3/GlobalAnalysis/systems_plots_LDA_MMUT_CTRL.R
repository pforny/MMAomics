#LDA of MMUT vs CTRL, including density plots
###############################################################


library(data.table)
library(dplyr)
library(tidyverse)
library(corrplot)
library(ggpubr)
library(ggplot2)
library(corrplot)
library(pheatmap)
library(tibble)
library(caret)
library(clipr)
library(ggrepel)
library(ggsci)
library(patchwork)
# LDA
library(MASS)



mypal = c(pal_aaas("default", alpha = 1)(7))

fig_path <- c("Figs/systems_plots/")
fig_path_pdf <- c("Figs/v8/pdf//")




################
# import data
################
TCA <- "Data/masterTCA_keggPlusmanual_circos.csv" #

tca_tbl <- fread(TCA, header = TRUE)[,ENSG]

hgnc_r <- fread("Data/ENSG_HGNC_list_r.csv")
setnames(hgnc_r, c("ENSG", "HGNC"))
hgnc_p <- fread("Data/ENSG_HGNC_list_p.csv")
setnames(hgnc_p, c("ENSG", "HGNC"))

# load transcriptomics and protein data
load("interimData/rnafiltered_mat.RData")
load("interimData/prot_mat.RData")



###############################################################
###############################################################
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
rna_exp_common_genes_rest <-rna_exp_common_genes[, c(1,145:222)]

prot_exp_common_genes<- prot_exp_common_genes %>% relocate(ENSG)
prot_exp_common_genes_MMUT <- prot_exp_common_genes[, 1:144]
prot_exp_common_genes_rest <- prot_exp_common_genes[, c(1,145:222)]

# prot_exp_common_genes_rest$ENSG <- prot_exp_common_genes$ENSG
# rna_exp_common_genes_rest$ENSG <- rna_exp_common_genes$ENSG

merged_rna_prot_MMUT <- merge(rna_exp_common_genes_MMUT, prot_exp_common_genes_MMUT,  by = "ENSG") # only genes shared in both
# column 1 -> ENSG
# columns 2:144 -> rna data (noted with .x in name)
# columns 145:287 -> protein data (noted with .y in name)

merged_rna_prot_rest <- merge(rna_exp_common_genes_rest, prot_exp_common_genes_rest, by ="ENSG")
# column 1 -> ENSG
# columns 2:79 -> rna
# columns 80:157 -> protein

# rna exp ALL GENES separated
rna_expression_MMUT <-rna_exp[, 1:144]
rna_expression_REST <-rna_exp[, c(1,145:222)]

# protein exp ALL GENES separated
prot_expression_MMUT <- prot_exp[, 1:151]
prot_expression_REST <- prot_exp[, c(1,152:231)]


######################################
######################################
# RNA DATA prep
######################################
######################################

# rna_exp_common_genes_MMUT$category <- "MMUT"
# rna_exp_common_genes_rest$category<- "CTRL"
rna_exp_MMUT<- as.data.frame(t(rna_exp_common_genes_MMUT))
# rna_exp_MMUT<- as.data.frame(t(rna_expression_MMUT))
names(rna_exp_MMUT)<- rna_exp_MMUT[1,]
rna_exp_MMUT<- rna_exp_MMUT[-1,]
rna_exp_MMUT <- rownames_to_column(rna_exp_MMUT, "patient")

rna_exp_CTRL<- as.data.frame(t(rna_exp_common_genes_rest))
# rna_exp_CTRL<- as.data.frame(t(rna_expression_REST))
names(rna_exp_CTRL) <- rna_exp_CTRL[1,]
rna_exp_CTRL<- rna_exp_CTRL[-1,]
rna_exp_CTRL <- rownames_to_column(rna_exp_CTRL, "patient")

rna_exp_MMUT$cat <- "MMUT def."
rna_exp_CTRL$cat <- "control"

merged_RNA <- bind_rows(rna_exp_MMUT, rna_exp_CTRL)

#######

#prep data for LDA
test_data <- as.data.table(merged_RNA)
# test_data <- bind_cols(test_data[,1:3001], merged_RNA[,4320])
# names(test_data)[4320]<-paste("cat")

d<-test_data
d<-d[,-1]
d$cat<-as.factor(d$cat) #column 'cat' needs to be a factor!
d<-as.data.frame(d)
ncol(d)-1
d[,1:(ncol(d)-1)] <- sapply(d[,1:(ncol(d)-1)], as.numeric)

set.seed(123)
training.samples <- d$cat %>%
  createDataPartition(p = .5, list = FALSE)
train.data <- d[training.samples, ]
test.data <- d[-training.samples, ]

# Estimate preprocessing parameters
#normalization
preproc.param <- train.data %>%
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)

# Fit  model
model <- lda(cat~., data = train.transformed)

#all coefficents driving the separation
coeffi.dfRNA <- as.data.frame(model$scaling)

coeffi.df_rna <- cbind( "ENSG" = rownames(coeffi.dfRNA), coeffi.dfRNA)
coeffi.df_rna <- left_join(coeffi.df_rna, hgnc_r) %>% 
  mutate(absolute.LD1 = abs(.$LD1))


# Make predictions
predictions <- model %>% predict(test.transformed)
#predicted classes
# head(predictions$class, 10)
# # Predicted probabilities class membership
# head(predictions$posterior, 6) 
# # Linear discriminants
head(predictions$x, 5)

# Model accuracy
mean(predictions$class==test.transformed$cat) 
#accuracy  0.7142857

#plot LDA
lda.data <- cbind(train.transformed, predict(model)$x)
ggplot(lda.data, aes(LD1, 0)) +
  geom_point(aes(color = cat))+
  theme_light()






######################################
######################################
# PROTEIN DATA prep
######################################
######################################

summary(prot_exp_common_genes_MMUT$ENSG == prot_exp_common_genes_rest$ENSG)

prot_exp_MMUT<- as.data.frame(t(prot_exp_common_genes_MMUT))
colnames(prot_exp_common_genes_MMUT)
# prot_exp_MMUT<- as.data.frame(t(prot_expression_MMUT))
names(prot_exp_MMUT)<- prot_exp_MMUT[1,]
prot_exp_MMUT<- prot_exp_MMUT[-1,]
prot_exp_MMUT <- rownames_to_column(prot_exp_MMUT, "patient")


prot_exp_CTRL<- as.data.frame(t(prot_exp_common_genes_rest))
colnames(prot_exp_common_genes_rest)
# prot_exp_CTRL<- as.data.frame(t(prot_expression_REST))
names(prot_exp_CTRL) <- prot_exp_CTRL[1,]
prot_exp_CTRL<- prot_exp_CTRL[-1,]
prot_exp_CTRL <- rownames_to_column(prot_exp_CTRL, "patient")

prot_exp_MMUT$cat <- "MMUT def."
prot_exp_CTRL$cat <- "control"

colnames(prot_exp_MMUT)[1:6]
colnames(prot_exp_CTRL)[1:6]

merged_prot <- bind_rows(prot_exp_MMUT, prot_exp_CTRL)


#####


test_data_prot <- as.data.table(merged_prot)
# test_data <- bind_cols(test_data[,1:3001], merged_RNA[,4320])
# names(test_data_prot)[4320]<-paste("cat")

prot<-test_data_prot
prot<-prot[,-1]
prot$cat<-as.factor(prot$cat) #column 'cat' needs to be a factor!
prot<-as.data.frame(prot)
prot[,1:(ncol(prot)-1)] <- sapply(prot[,1:(ncol(prot)-1)], as.numeric)


set.seed(123)
training.samples_prot <- prot$cat %>%
  createDataPartition(p = 0.5, list = FALSE)
train.data_prot <- prot[training.samples_prot, ]
test.data_prot <- prot[-training.samples_prot, ]

# Estimate preprocessing parameters
#normalization
preproc.param.prot <- train.data_prot %>%
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed.prot <- preproc.param.prot %>% predict(train.data_prot)
test.transformed.prot <- preproc.param.prot %>% predict(test.data_prot)

# Fit  model
# model_prot <- lda(cat~., data = train.transformed.prot)
model_prot<- lda(cat~., data = train.transformed.prot)

#all coefficents driving the separation
coeffi.df <- as.data.frame(model_prot$scaling)

coeffi.df_ <- cbind( "ENSG" = rownames(coeffi.df), coeffi.df)
coeffi.df_ <- left_join(coeffi.df_, hgnc_r) %>% 
  mutate(absolute.LD1 = abs(.$LD1))
# write.csv(coeffi.df_, "LDA_MMUT_CTRL_linear_components.csv")

top20drivers <- coeffi.df %>% 
  top_n(., 20)
df.top20 <- cbind("ENSG" = rownames(top20drivers), top20drivers)
df.top20 <- left_join(df.top20, hgnc_r) 

bottom20divers <- coeffi.df %>% 
  top_n(., -20)
df.bottom20 <- cbind("ENSG" = rownames(bottom20divers), bottom20divers)
df.bottom20 <- left_join(df.bottom20, hgnc_r)

# Make predictions
predictions_prot <- model_prot %>% predict(test.transformed.prot[,-4406])
#predicted classes
# head(predictions_prot$class, 10)
# # Predicted probabilities class membership
# head(predictions_prot$posterior, 6) 
# # Linear discriminants
# head(predictions_prot$x, 5) 

# Model accuracy
mean(predictions_prot$class==test.transformed.prot[,4406]) 
#accuracy 0.6363636


#plot LDA
lda.data.prot <- cbind(train.transformed.prot, predict(model_prot)$x)
colnames(lda.data.prot) <- make.unique(names(lda.data.prot))
ggplot(lda.data.prot, aes(LD1, 1)) +
  geom_point(aes(color = cat))+
  theme_light()
# ggsave("LAD_prot_1D_88.jpg")





######
#density plot rna
plt_rna1 <- 
ggplot(lda.data, aes(x=LD1, fill = cat))+
  geom_density(color = 0, alpha = 0.6) +
  scale_fill_aaas() +
  ylab("Density") +
  xlab("Linear discriminant 1 (LD1)") +
  ggtitle("Transcript-driven separation") +
  theme_pubr() +
  theme(plot.title = element_text(size = 12), legend.position = c(0.2, 0.9), legend.title = element_blank())



#__________________________________________________________
#Plot LD1 of genes ranked with TCA genes labeled
dftest_prot <-coeffi.df_[order(-coeffi.df_$absolute.LD1),]
dftest_prot$nr <- seq(1,length(coeffi.df_[,1]), by =1)
dftest_prot <- as.data.table(dftest_prot)

dftest_rna <-coeffi.df_rna[order(-coeffi.df_rna$absolute.LD1),]
dftest_rna$nr <- seq(1,length(coeffi.df_rna[,1]), by =1)
dftest_rna <- as.data.table(dftest_rna)
colnames(dftest_rna)


ggplot(dftest_prot, aes(x = nr, y = absolute.LD1)) +
  geom_bar(stat = "identity")+
  geom_label_repel(data = dftest_prot[HGNC == "MMUT" | HGNC == "OGDH" | HGNC == "GLUD1", ], aes(label = HGNC), nudge_y = 1)

dftest_rna_interest <- dftest_rna[HGNC == "MMUT" | HGNC == "OGDH" | HGNC == "SUCLA2" | HGNC == "PDHB", ]

plt_rna2 <- 
ggplot(dftest_rna, aes(x = nr, y = absolute.LD1)) +
  geom_bar(stat = "identity", fill = "grey", alpha = 0.8) +
  geom_point(data = dftest_rna_interest, color = "white", fill = "black", shape = 21, size = 3) +
  geom_label_repel(data = dftest_rna_interest, aes(label = HGNC), color = "black", nudge_x = 400) +
  xlab("Transcript rank") +
  ylab("Coefficient of LD1") +
  ggtitle("Driving genes of separation") +
  theme_pubr() +
  theme(plot.title = element_text(size = 12))
plt_rna2



TCA <- "Data/masterTCA_keggPlusmanual_circos.csv"
tca_tbl <- fread(TCA, header = TRUE)[,ENSG]
tca_tbl <- as.data.table(tca_tbl)
colnames(tca_tbl) <- "ENSG"

hgnc_r <- fread("Data/gene_sets/ENSG_HGNC_list_r.csv")
setnames(hgnc_r, c("ENSG", "HGNC"))
tca_genes_tbl <- inner_join(tca_tbl, hgnc_r, by="ENSG")
tcagenes_new <- tca_genes_tbl$HGNC


df.plot<-full_join(dftest_rna, tca_genes_tbl, by = "ENSG")

ggplot(df.plot, aes(x=nr, y= absolute.LD1, label = HGNC.y)) +
  geom_bar(stat = "identity", width = )+
  geom_text_repel(min.segment.length = 1, max.overlaps = 50)#, arrow = arrow(length = unit(0.015, "npc")))


LDA_plt <- 
ggarrange(plt_rna1, plt_rna2, widths = c(1, 1))

ggsave(paste(fig_path, "LDA_transcriptomics.png", sep = ""), LDA_plt, device = png(), width = 9, height = 3, bg = "white")
dev.off()
ggsave(paste(fig_path_pdf,"Fig3/", "LDA_transcriptomics.pdf", sep = ""), LDA_plt, device = "pdf", width = 8, height = 3, bg = "white")
dev.off()



LDA_plt2 <- 
plt_rna1/plt_rna2

ggsave(paste(fig_path, "LDA_transcriptomics2.png", sep = ""), LDA_plt2, device = png(), width = 4, height = 6, bg = "white")
dev.off()




