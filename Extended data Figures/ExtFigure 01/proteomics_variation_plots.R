# these plots illustrate the quality of the data of the proteomics dataset.

library(data.table)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(patchwork)


fig_path <- c("Figs/QC/")
fig_path_pdf <- c("Figs/v12/pdf/")
mypal1 <- pal_aaas("default", alpha = 1)(9)
mypal2 <- c(mypal1[2], mypal1[4], mypal1[1], mypal1[9])



# function to caculate coefficient of variation
cv_na <- function(x) {
  return(sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
}


##################################
# load protein matrix
load("interimData/prot_mat.RData")

prot_exp <- data.table(ENSG = row.names(prot_mat), prot_mat)
prot_exp <- prot_exp %>% distinct(ENSG, .keep_all = TRUE)
colnames(prot_exp)
prot_exp <- prot_exp[, ENSG := NULL]

#wenguang: tried different cut-off, such as 2^16, 100000 and 2^17, 2^16 gives the most reasonable results, in terms of # proteins identified per run.
intensity_cutoff <- 2^15.5
num_presence <- apply(prot_exp, 1, function(x) length(which(x > intensity_cutoff)))

# overlap in detected proteins among samples

ovrlp_prot <- data.table(
  length(which(num_presence > 1 )),
  length(which(num_presence > 0.25*230 )),
  length(which(num_presence > 0.50*230 )),
  length(which(num_presence > 0.75*230 )),
  length(which(num_presence ==   1*230 ))
)
ovrlp_prot <- melt(ovrlp_prot)
setnames(ovrlp_prot, c("threshold", "number_of_proteins"))
ovrlp_prot[, threshold := c("1 or more", ">25%", ">50%", ">75%", "100%")]
ovrlp_prot$threshold <- factor(ovrlp_prot$threshold, levels = c("1 or more", ">25%", ">50%", ">75%", "100%"))

# how many proteins are in more than 75% of samples?
ovrlp_prot[threshold == ">75%", ]

ovrlp_prot_plot <- 
ggplot(ovrlp_prot, aes(x = threshold, y = number_of_proteins)) +
  geom_col() +
  xlab("Overlap of proteins quantified") +
  ylab("Number of proteins") +
  theme_pubr() +
  theme(plot.title = element_text(size=12))

ggsave(paste(fig_path,"OverlapOfProteinsQuantified.png", sep = ""), ovrlp_prot_plot, device = png(), width = 4, height = 3)


# comparison of variation coeffcient among sample types

var_tbl <- prot_exp[which(num_presence > 0.50*230), ]

coef_tbl <- data.table(
  apply(var_tbl[, 1:150, with = FALSE], 1, function(x) cv_na(x[which(x > intensity_cutoff)])),
  apply(var_tbl[, 151:210, with = FALSE], 1, function(x) cv_na(x[which(x > intensity_cutoff)])),
  apply(var_tbl[, 211:230, with = FALSE], 1, function(x) cv_na(x[which(x > intensity_cutoff)])),
  apply(var_tbl[, 1:230, with = FALSE], 1, function(x) cv_na(x[which(x > intensity_cutoff)]))
)
setnames(coef_tbl, c("MMUT def.", "unknown", "unaffected", "all"))
coef_tbl <- melt(coef_tbl)
setnames(coef_tbl, c("type", "coefficient_of_variation"))

coef_var_prot_plot <- 
ggplot(coef_tbl, aes(x = type, y = coefficient_of_variation, color = type, fill = type)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(color = "black", width = 0.1, alpha = 0.4) +
  ylab("Coefficient of variation") +
  ylim(0, 1.5) +
  scale_color_manual(values=mypal2) +
  scale_fill_manual(values=mypal2) +
  theme_pubr() +
  theme(axis.title.x = element_blank(), legend.position = "none")

ggsave(paste(fig_path,"CoefficientOfVariationProteomics.png", sep = ""), coef_var_prot_plot, device = png(), width = 5, height = 3)
dev.off()

proteomics_bas <- 
ovrlp_prot_plot / coef_var_prot_plot + plot_annotation(title = "Proteomics QC", theme = theme(plot.title = element_text(size = 12)))

ggsave(paste(fig_path,"proteomics_basics.png", sep = ""), proteomics_bas, device = png(), width = 4, height = 5)
dev.off()
ggsave(paste(fig_path_pdf, "SuppFig1_qc/","proteomics_basics.pdf", sep = ""), proteomics_bas, device = "pdf", width = 4.5, height = 5)
dev.off()



