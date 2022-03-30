# quality control plots for genomics

library(data.table)
# library(readr)
# library(readxl)
# library(reshape2)
library(ggpubr)
library(extrafont)
library(gridExtra)
library(gtools)
require(cowplot)
# library(EnhancedVolcano)
# library(dplyr)
# library(plyr)
# library(biomaRt)
library(tidyverse)
library(ggridges)
library(ggsci)
library(patchwork)


system("mkdir Figs/QC")

fig_path <- c("Figs/QC/")
fig_path_pdf <- c("Figs/v8/pdf/")


##################################
# genomics QC

gen_qc <- data.table(read.table("Data/QC/WGSsummaryStats_MMAproject.tsv", header = TRUE))
colnames(gen_qc)
gen_qc[, sample_no := (as.numeric(sub("MMA", "", gen_qc$SAMPLE_NAME)))]
gen_qc[, type := 
  ifelse(sample_no<151, "MMUT def.", 
    ifelse(sample_no>210, "unknown", "unaffected"))]
gen_qc$type <- factor(gen_qc$type, levels = c("MMUT def.", "unknown", "unaffected"))
unique(gen_qc[, type])

mean_cov <- gen_qc$MEAN_COVERAGE

cov_sec <- c("≥5x", "≥15x", "≥30x", "≥40x", "≥50x")

gen_qc_cov <- gen_qc %>% select(type, MEAN_COVERAGE)
gen_qc_cov[, coverage := ifelse(MEAN_COVERAGE >= 50, cov_sec[5], ifelse(MEAN_COVERAGE >= 40, cov_sec[4], ifelse(MEAN_COVERAGE >= 30, cov_sec[3], ifelse(MEAN_COVERAGE >= 15, cov_sec[2], ifelse(MEAN_COVERAGE >= 5, cov_sec[1])))))]
gen_qc_cov[, cov50 := ifelse(MEAN_COVERAGE >= 50, 1, 0)]
gen_qc_cov[, cov40 := ifelse(MEAN_COVERAGE >= 40, 1, 0)]
gen_qc_cov[, cov30 := ifelse(MEAN_COVERAGE >= 30, 1, 0)]
gen_qc_cov[, cov15 := ifelse(MEAN_COVERAGE >= 15, 1, 0)]
gen_qc_cov[, cov5 := ifelse(MEAN_COVERAGE >= 5, 1, 0)]


coverage_per_type00 <- gen_qc_cov[, -c("MEAN_COVERAGE", "coverage")]
coverage_per_type0 <- melt.data.table(coverage_per_type00, id.vars = c("type"))
coverage_per_type <- data.table(table(coverage_per_type0[value != 0, ]))

sapply(coverage_per_type, class)
coverage_per_type$type <- factor(coverage_per_type$type, levels = c("MMUT def.", "unknown", "unaffected"))
coverage_per_type$variable <- factor(coverage_per_type$variable, levels = c("cov5", "cov15", "cov30", "cov40", "cov50"))
levels(coverage_per_type$variable) <- cov_sec


mean_cov_hist <-
ggplot(gen_qc_cov, aes(x = MEAN_COVERAGE)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.7) +
  ggtitle("Histogram of mean coverage") +
  theme_test() +
  theme(axis.title.x=element_blank())

mean_cov_dens  <-
ggplot(gen_qc_cov, aes(x = MEAN_COVERAGE)) +
  geom_density(fill = "black", alpha = 0.4) +
  ggtitle("Histogram_mean_coverage") +
  theme_bw() +
  theme(axis.title.x = element_blank())


# mean_cov_tbl <- data.table(
#   length(which(mean_cov < 15)),
#   sum(mean_cov >= 15 & mean_cov < 25),
#   sum(mean_cov >= 25 & mean_cov < 30),
#   sum(mean_cov >= 30 & mean_cov < 40),
#   sum(mean_cov >= 40 & mean_cov < 45),
#   sum(mean_cov >= 45)
# )
# mean_cov_tbl <- melt(mean_cov_tbl)
# setnames(mean_cov_tbl, c("mean_coverage", "Count"))
# mean_cov_tbl[, mean_coverage := c("0-15x", "15-25x", "25-30x", "30-40x", "40-45x", ">45x")]
# mean_cov_tbl$mean_coverage <- factor(mean_cov_tbl$mean_coverage, levels = c("0-15x", "15-25x", "25-30x", "30-40x", "40-45x", ">45x"))

# mean_cov_colplot <- 
# ggplot(mean_cov_tbl, aes(x = mean_coverage, y = Count)) +
#   geom_col() +
#   ggtitle("Mean_coverage") +
#   theme_bw() +
#   theme(axis.title.x = element_blank())


mean_cov_tbl1 <- data.table(
  length(which(mean_cov >= 5)),
  length(which(mean_cov >= 15)),
  length(which(mean_cov >= 30)),
  length(which(mean_cov >= 40)),
  length(which(mean_cov >= 50))
)
mean_cov_tbl1 <- melt(mean_cov_tbl1)
setnames(mean_cov_tbl1, c("mean_coverage", "Count"))
mean_cov_tbl1[, mean_coverage := c("at least 5x", "at least 15x", "at least 30x", "at least 40x", "at least 50x")]
mean_cov_tbl1$mean_coverage <- factor(mean_cov_tbl1$mean_coverage, levels = c("at least 5x", "at least 15x", "at least 30x", "at least 40x", "at least 50x"))
levels(mean_cov_tbl1$mean_coverage) <- c("≥5x", "≥15x", "≥30x", "≥40x", "≥50x")


mean_cov_colplot1 <- 
ggplot(mean_cov_tbl1, aes(x = mean_coverage, y = Count)) +
  geom_col() +
  geom_text(aes(label = Count), vjust = -0.5) +
  ylim(0,250) +
  ylab("Sample count") +
  xlab("Mean genomic coverage") +
  theme_pubr() +
  theme(plot.title = element_text(size = 12))

ggsave(paste(fig_path,"MeanGenomicCoverage.png", sep = ""), mean_cov_colplot1, device = png(), width = 4, height = 3)
dev.off()



mean_cov_colplotType <- 
ggplot(coverage_per_type, aes(x = variable, y = N, fill = type)) +
  geom_col() +
  geom_text(data = mean_cov_tbl1, aes(x = mean_coverage, y = Count, label = Count), vjust = -0.5, inherit.aes = FALSE) +
  ylim(0,250) +
  ylab("Sample count") +
  xlab("Mean genomic coverage") +
  scale_fill_npg() +
  theme_pubr() +
  theme(plot.title = element_text(size = 12), legend.position = "right",  legend.title = element_blank())

ggsave(paste(fig_path,"MeanGenomicCoveragePerType.png", sep = ""), mean_cov_colplotType, device = png(), width = 5, height = 3)
dev.off()





gen_qc_frac <- gen_qc %>% select(SAMPLE_NAME, "PCT_1X", "PCT_5X", "PCT_10X", "PCT_15X", "PCT_20X", "PCT_25X", "PCT_30X", "PCT_40X", "PCT_50X", "PCT_60X", "PCT_70X", "PCT_80X", "PCT_90X", "PCT_100X")
# setnames(gen_qc_frac, c("SAMPLE_NAME", "1x", "5x", "10x", "15x", "20x", "25x", "30x", "40x", "50x", "60x", "70x", "80x", "90x", "100x"))
setnames(gen_qc_frac, c("SAMPLE_NAME", 1, seq(5,30, by=5), seq(40, 100, by=10)))
gen_qc_frac_melt <- melt(gen_qc_frac, id.vars = c("SAMPLE_NAME"))
setnames(gen_qc_frac_melt, c("SAMPLE_NAME", "Coverage", "FractionOfGenome"))
gen_qc_frac_melt$Coverage <- as.numeric(levels(gen_qc_frac_melt$Coverage))[gen_qc_frac_melt$Coverage]


mean_cov_persample_plot <-
ggplot(gen_qc_frac_melt, aes(x = Coverage, y = FractionOfGenome, group = SAMPLE_NAME)) +
  geom_line(alpha = 0.16) +
  xlab("Sequencing coverage per sample (x times)") +
  ylab("Faction of genome") +
  theme_pubr() +
  theme(legend.position = "none", plot.title = element_text(size = 12))

ggsave(paste(fig_path,"SequencingCoveragePerSample.png", sep = ""), mean_cov_persample_plot, device = png(), width = 5, height = 3)
dev.off()


gen_qc_reads <- gen_qc %>% select(SAMPLE_NAME, PF_READS, PF_READS_ALIGNED, PF_HQ_ALIGNED_READS)
gen_qc_reads[, SAMPLE_NAME := (as.numeric(sub("MMA", "", gen_qc_reads$SAMPLE_NAME)))]
setnames(gen_qc_reads, c("SAMPLE_NAME", "Reads", "Aligned r.", "HQ aligned r."))
gen_qc_reads_melt <- melt(gen_qc_reads, id.vars = c("SAMPLE_NAME"))
setnames(gen_qc_reads_melt, c("SAMPLE_NAME", "read_type", "reads"))

# what is the mean of high quality aligned reads?
mean(gen_qc_reads_melt[read_type == "HQ aligned reads", ]$reads)


readno_plot <- 
ggplot(gen_qc_reads_melt, aes(x = SAMPLE_NAME, y = log(reads))) +
  geom_point(alpha = 0.2) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~read_type) +
  ylim(20,21.2) +
  ggtitle("NumberOfReads") +
  theme_bw()

readno_violin <- 
ggplot(gen_qc_reads_melt, aes(x = read_type, y = reads)) +
  geom_violin(alpha = 1, fill = "grey", color = "grey") +
  geom_boxplot(color = "black", width = 0.1, alpha = 0.4) +
  ylab("Number of reads") +
  scale_color_npg() +
  scale_fill_npg() +
  theme_pubr() +
  theme(legend.position = "none", axis.title.x = element_blank())

ggsave(paste(fig_path,"NumberOfReads.png", sep = ""), readno_violin, device = png(), width = 5, height = 3)
dev.off()


readno_ridges <- 
  ggplot(gen_qc_reads_melt, aes(x = log(reads), y = read_type,  fill = read_type)) +
  geom_density_ridges(alpha = 0.6, bandwidth = 4) +
  theme_bw() +
  theme(legend.position = "none")
  


genome_qc_plt <- ggarrange(readno_violin, mean_cov_persample_plot, mean_cov_colplot1, mean_cov_hist)

ggsave(paste(fig_path,"Genome_QC.png", sep = ""), genome_qc_plt, device = png(), width = 220, height = 140, units = c("mm"), dpi = 600)
dev.off()



# figure compilation for paper
genomic_qc <- 
readno_violin / mean_cov_persample_plot + plot_annotation(title = "Genomics QC", theme = theme(plot.title = element_text(size = 12)))
ggsave(paste(fig_path_pdf, "SuppFig1/", "GenomicsQC.pdf", sep = ""), device = "pdf", width = 4.5, height = 5)
dev.off()


#################################################
# special flattened violin plot with raw data
# code from https://www.data-to-viz.com/caveat/boxplot.html

# "%||%" <- function(a, b) {
#   if (!is.null(a)) a else b
# }

# geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
#                         position = "dodge", trim = TRUE, scale = "area",
#                         show.legend = NA, inherit.aes = TRUE, ...) {
#   layer(
#     data = data,
#     mapping = mapping,
#     stat = stat,
#     geom = GeomFlatViolin,
#     position = position,
#     show.legend = show.legend,
#     inherit.aes = inherit.aes,
#     params = list(
#       trim = trim,
#       scale = scale,
#       ...
#     )
#   )
# }

# GeomFlatViolin <-
#   ggproto("GeomFlatViolin", Geom,
#           setup_data = function(data, params) {
#             data$width <- data$width %||%
#               params$width %||% (resolution(data$x, FALSE) * 0.9)

#             # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
#             data %>%
#               group_by(group) %>%
#               mutate(ymin = min(y),
#                      ymax = max(y),
#                      xmin = x,
#                      xmax = x + width / 2)

#           },

#           draw_group = function(data, panel_scales, coord) {
#             # Find the points for the line to go all the way around
#             data <- transform(data, xminv = x,
#                               xmaxv = x + violinwidth * (xmax - x))

#             # Make sure it's sorted properly to draw the outline
#             newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
#                              plyr::arrange(transform(data, x = xmaxv), -y))

#             # Close the polygon: set first and last point the same
#             # Needed for coord_polar and such
#             newdata <- rbind(newdata, newdata[1,])

#             ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
#           },

#           draw_key = draw_key_polygon,

#           default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
#                             alpha = NA, linetype = "solid"),

#           required_aes = c("x", "y")
# )


# ggplot(gen_qc_reads_melt, aes(x = read_type, y = log(reads), fill = read_type)) + 
#     geom_flat_violin(legend = "none") +
#     coord_flip()










##################################
# transcriptomics QC



rna_qc_tbl0 <- fread("Data/QC/MMA_PHRT_R_baseQmetrix.csv")
setnames(rna_qc_tbl0, c("mma_id", "run_id", "no1", "no2", "score_mean", "score_sd", "depth"))

rna_qc_tbl <- rna_qc_tbl0 %>% group_by(mma_id, run_id) %>%
  mutate(score_mean2 = mean(score_mean), score_sd2 = mean(score_sd))

rna_qc_tbl <- as.data.table(rna_qc_tbl)

head(rna_qc_tbl, 20)

# check maximum number of n-plicates of mma-ids
txtbl <- rna_qc_tbl %>% group_by(mma_id) %>% summarize(n())
max(txtbl[2])

# compute weighed mean and sd depending on depth
mma_ids <- unique(rna_qc_tbl0$mma_id)


for (i in 1:length(mma_ids)) {
  cat(i, "\n")
  work_tbl <- rna_qc_tbl[mma_id == mma_ids[i], ]
  dims <- dim(work_tbl)[1]
  sum_depth <- sum(unique(work_tbl$depth))
  if(dims == 2) {
    mean_tmp <- work_tbl$score_mean2[1]/sum_depth*work_tbl$depth[1]
    sd_tmp <- work_tbl$score_sd2[1]/sum_depth*work_tbl$depth[1]
  } 
  if(dims == 4) {
    mean_tmp <- work_tbl$score_mean2[1]/sum_depth*work_tbl$depth[1] + work_tbl$score_mean2[3]/sum_depth*work_tbl$depth[3]
    sd_tmp <- work_tbl$score_sd2[1]/sum_depth*work_tbl$depth[1] + work_tbl$score_sd2[3]/sum_depth*work_tbl$depth[3]
  }
  if(dims == 6) {
    mean_tmp <- work_tbl$score_mean2[1]/sum_depth*work_tbl$depth[1] + work_tbl$score_mean2[3]/sum_depth*work_tbl$depth[3] + work_tbl$score_mean2[5]/sum_depth*work_tbl$depth[5]
    sd_tmp <- work_tbl$score_sd2[1]/sum_depth*work_tbl$depth[1] + work_tbl$score_sd2[3]/sum_depth*work_tbl$depth[3] + work_tbl$score_sd2[5]/sum_depth*work_tbl$depth[5] 
  }
  if(dims == 8) {
    mean_tmp <- work_tbl$score_mean2[1]/sum_depth*work_tbl$depth[1] + work_tbl$score_mean2[3]/sum_depth*work_tbl$depth[3] + work_tbl$score_mean2[5]/sum_depth*work_tbl$depth[5] + work_tbl$score_mean2[7]/sum_depth*work_tbl$depth[7]
    sd_tmp <- work_tbl$score_sd2[1]/sum_depth*work_tbl$depth[1] + work_tbl$score_sd2[3]/sum_depth*work_tbl$depth[3] + work_tbl$score_sd2[5]/sum_depth*work_tbl$depth[5] + work_tbl$score_sd2[7]/sum_depth*work_tbl$depth[7]
  }
  rna_qc_tbl[mma_id == mma_ids[i], score_mean3 := mean_tmp]
  rna_qc_tbl[mma_id == mma_ids[i], score_sd3 := sd_tmp]
}

rna_qc_tbl

rna_qc_tbl_uniq <- rna_qc_tbl[match(unique(rna_qc_tbl$mma_id), rna_qc_tbl$mma_id), ]

library(forcats)

ggplot(rna_qc_tbl_uniq, aes(x = fct_reorder(mma_id, -score_mean3), y = score_mean3)) +
  geom_point(size = 0.1) +
  geom_errorbar(aes(ymin = score_mean3-score_sd3, ymax = score_mean3+score_sd3), size = 0.1) +
  ylim(0, 50) +
  ylab("Mean base quality score") +
  xlab("Ranked samples") +
  ggtitle("PHRED trascriptomic quality score ") +
  theme_pubr() +
  theme(plot.title = element_text(size=12), axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(paste(fig_path,"PHREDscoreTranscriptomics.png", sep = ""), device = png(), width = 4, height = 3)
dev.off()






# intermediary file received from Cédric to display PHRED score
mean_base_tbl <- data.table(read.table("Data/QC/MMA_PHRT_R_meanBaseQperCycle.tsv", header=TRUE, sep=";"))
colnames(mean_base_tbl)
unique(mean_base_tbl$cycleNb)

mean_base_tbl[, new_cycle := "1-2"]
mean_base_tbl[cycleNb > 2, new_cycle := "3-5"]
mean_base_tbl[cycleNb > 5, new_cycle := "6-10"]
mean_base_tbl[cycleNb > 10, new_cycle := "11-20"]
mean_base_tbl[cycleNb > 20, new_cycle := "21-50"]
mean_base_tbl[cycleNb > 50, new_cycle := ">50"]

mean_base_tbl$new_cycle <- factor(mean_base_tbl$new_cycle, levels = c("1-2", "3-5", "6-10", "11-20", "21-50", ">50"))

# what is the median Phred score at 3-5 cycles?
median(mean_base_tbl[cycleNb == 4, ]$meanBaseQ)

qc_transcript <- 
ggplot(mean_base_tbl, aes(x = new_cycle, y = meanBaseQ)) + 
  geom_boxplot(outlier.size = 0.1, outlier.shape = 1) + 
  ylab("Phred quality score") +
  xlab("Number of cycles") +
  ggtitle("Transcriptomics QC") +
  theme_pubr() +
  theme(plot.title = element_text(size=12)) 

ggsave(paste(fig_path,"PHREDscoreTranscriptomics2.png", sep = ""),qc_transcript, device = png(), width = 4.5, height = 3)
dev.off()
ggsave(paste(fig_path_pdf, "SuppFig1/","PHREDscoreTranscriptomics2.pdf", sep = ""), qc_transcript, device = "pdf", width = 4.5, height = 2.5)
dev.off()




