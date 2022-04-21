# install old version of MOFA (manually)
# library(remotes)
# install_github("pmbio/mofa")
# also install mofadata

# libraries
require(data.table)
require(ggplot2)
require(tidyverse)
library(MOFA)
library(org.Hs.eg.db)
library(ggpubr)
library(patchwork)
library(ggrepel)


# load human gene database
hs <- org.Hs.eg.db

# create figures path
system("mkdir Figs/systems_plots")

fig_path <- c("Figs/systems_plots/")
fig_path_pdf <- c("Figs/v8/pdf/")



# slightly modified plotting function from MOFA package

plotEnrichmentDetailed <- function(object, factor, feature.sets, fsea.results, 
                                   adjust = TRUE, alpha = 0.1, max.genes = 5, max.pathways = 10, text_size = 3) {
  
  # Sanity checks
  stopifnot(length(factor)==1) 
  if(is.numeric(factor)) factor <- factorNames(object)[factor]
  if(!factor %in% colnames(fsea.results$pval)) 
    stop(paste0("No feature set enrichment calculated for factor ", factor, ".\n Use runEnrichmentAnalysis first."))
  
  # Fetch and prepare data  
  
  # foo
  foo <- melt(fsea.results$feature.statistics, factorsAsStrings=TRUE)
  colnames(foo) <- c("feature","factor","feature.statistic")
  foo <- foo[foo$factor==factor,]
  foo$feature <- as.character(foo$feature)
  
  # bar
  bar <- melt(feature.sets)
  bar <- bar[bar$value==1,c(1,2)]
  colnames(bar) <- c("pathway","feature")
  bar$pathway <- as.character(bar$pathway)
  bar$feature <- as.character(bar$feature)
  
  # baz
  if (adjust) {
    baz <- melt(fsea.results$pval.adj)
  } else {
    baz <- melt(fsea.results$pval)
  }
  colnames(baz) <- c("pathway","factor","pvalue")
  baz$pathway <- as.character(baz$pathway)
  baz <- baz[baz$factor==factor,]
  
  # Filter out pathways by p-values
  baz <- baz[baz$pvalue<=alpha,,drop=FALSE]
  if(nrow(baz)==0) {
    stop("No siginificant pathways at the specified alpha threshold. \n
         For an overview use plotEnrichmentHeatmap() or plotEnrichmentBars().")
  }
  if (nrow(baz) > max.pathways)
    baz <- head(baz[order(baz$pvalue),],n=max.pathways)
  
  # order pathways according to significance
  baz$pathway <- factor(baz$pathway, levels = baz$pathway[order(baz$pvalue, decreasing = TRUE)])
  
  # Merge
  foobar <- merge(foo, bar, by="feature")
  tmp <- merge(foobar, baz, by=c("pathway","factor"))
  
  # Select the top N features with the largest feature.statistic (per pathway)
  tmp_filt <- top_n(group_by(tmp, pathway), n=max.genes, abs(feature.statistic))
  
  # Add number of features and p-value per pathway
  pathways <- unique(tmp_filt$pathway)
  
  # Add Ngenes and p-values to the pathway name
  df <- data.frame(pathway=pathways, nfeatures=rowSums(feature.sets)[pathways])
  df <- merge(df, baz, by="pathway")
  df$pathway_long_name <- sprintf("%s\n (Ngenes = %d) \n (p-val = %0.2g)",df$pathway, df$nfeatures, df$pvalue)
  tmp <- merge(tmp, df[,c("pathway","pathway_long_name")], by="pathway")
  tmp_filt <- merge(tmp_filt, df[,c("pathway","pathway_long_name")], by="pathway")
  
  # sort pathways by p-value
  order_pathways <- df$pathway_long_name[order(df$pvalue,decreasing=TRUE) ]
  tmp$pathway_long_name <- factor(tmp$pathway_long_name, levels=order_pathways)
  tmp_filt$pathway_long_name <- factor(tmp_filt$pathway_long_name, levels=order_pathways)
  
  p <- ggplot(tmp, aes_string(x="pathway_long_name", y="feature.statistic")) +
    geom_point(size=2, color="lightgrey") +
    geom_point(aes_string(x="pathway_long_name", y="feature.statistic"), size=3, color="black", data=tmp_filt) +
    geom_text_repel(aes_string(x="pathway_long_name", y="feature.statistic", label="feature"), size=text_size, color="black", force=1, data=tmp_filt) +
    labs(x="", y="Feature statistic", title="") +
    coord_flip() +
    theme(
      axis.line = element_line(color="black"),
      axis.text.y = element_text(size=rel(1.2), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(1.2), vjust=0.5, color='black'),
      axis.title.y=element_blank(),
      legend.position='none',
      panel.background = element_blank()
    )
  return(p)
}





# use_python("/home/user/python", required = TRUE)
# use_python("/usr/local/bin/python3.9", required = TRUE)


mod <- loadModel("./Data/MOFA/output_MOFAobject")

factorNames(mod) #get or set factor names
MOFA::featureNames(mod) #get or set feature names
MOFA::sampleNames(mod) #get or set sample names
viewNames(mod) #get or set view names
getDimensions(mod) #get dimensions (number of samples, features, etc.)
getFactors(mod) #get model factors
getWeights(mod) #get model weights
getTrainData(mod) #get training data
# impute(mod)
# getImputedData(mod) #get imputed data


#############################################
# Disentangle sources of variation

plotVarianceExplained(mod) #plot the variance explained by each factor in each view. This is the key plot of MOFA and should always be done before inspecting factors or weights.
calculateVarianceExplained(mod) #calculate and return the variance explained by each factor in each view.

class(MOFA::featureNames(mod))
summary(MOFA::featureNames(mod))



#############################################
# Inspect loadings

plotTopWeights(mod, view = "gene", factor = "LF7") #plot the top loadings for a given factor and view
plotWeights(mod, view = "protein", factor = "LF1") #plot all loadings for a given factor and view
plotWeightsHeatmap(mod, view = "protein") #plot a heatmap of the loadings from multiple factors in a given view


#############################################
# Inspect factors

plotFactorCor(mod) #correlation plot between factors. Ideally, they should be uncorrelated
plotFactorScatter(mod, factors = c("LF3", "LF2")) #scatterplot between two factors, this is similar to doing a PCA plot
plotFactorScatters(mod) #pairwise combination of scatterplots between multiple factors
plotFactorBeeswarm(mod, factors = "LF1") #beeswarm plot for a single factor








#############################################
# Inspect the data

plotDataOverview(mod) #plot overview of the input data, including the number of samples, views, features, and the missing assays.
plotDataHeatmap(mod, view = "protein", factor = "LF7") #heatmap of the training data using only top features for a given factor. This is very useful to map the factors and features back to the original data
plotDataScatter(mod, view = "protein", factor = "LF7") #scatterplot of the data using only top features for a given factor


#############################################
# Feature set enrichment analysis

# generate feature list (gene list)

data("reactomeGS", package = "MOFAdata")

fsea.results_prot <- runEnrichmentAnalysis(mod, view = "protein", feature.sets = reactomeGS)
fsea.results_rna <- runEnrichmentAnalysis(mod, view = "rna", feature.sets = reactomeGS)
fsea.results_gene <- runEnrichmentAnalysis(mod, view = "gene", feature.sets = reactomeGS)

#do feature set enrichment analysis. Takes a bit amount of options, check the example on the vignette
plotEnrichment(mod, fsea.results = fsea.results_prot, factor = "LF7") #plot the top most enriched feature sets per factor
plotEnrichmentDetailed(object = mod, fsea.results = fsea.results_prot, factor = "LF7", feature.sets = reactomeGS, max.pathways = 6, max.genes = 3, text_size = 3) #plot a more detailed output of the top most enriched feature sets per factor
plotEnrichmentBars(fsea.results_prot, alpha = 0.01) #plot the number of enriched feature sets per factor as a barplot


plotEnrichmentBars(fsea.results_prot, alpha = 0.01)

ggsave(paste(fig_path, "MOFA_enrichmentBars.png", sep = ""), device = png(), width = 4, height = 4)
dev.off()


margs = c(0.5, 0.5, 0.5, 0.5)
margin = theme(plot.margin = unit(margs, "cm"))

gene_mofa_detailed <- 
plotEnrichmentDetailed(object = mod, fsea.results = fsea.results_gene, factor = "LF8", feature.sets = reactomeGS, max.pathways = 5, max.genes = 5, text_size = 3.5, alpha = 0.1) + ggtitle("Genomics: feature statistics per set") + theme(axis.text.y = element_text(size = 10), plot.title = element_text(size = 12)) + margin

ggsave(paste(fig_path, "MofaFeatureDetailedGene.png", sep = ""), gene_mofa_detailed, device = png(), width = 12, height = 4)
# ggsave(paste(fig_path_pdf, "Fig3/", "MofaFeatureDetailedGene.pdf", sep = ""), gene_mofa_detailed, device = "pdf", width = 12, height = 4)
dev.off()


rna_mofa_detailed <- 
plotEnrichmentDetailed(object = mod, fsea.results = fsea.results_rna, factor = "LF8", feature.sets = reactomeGS, max.pathways = 5, max.genes = 5, text_size = 3.5, alpha = 0.1) + ggtitle("Transcriptomics: feature statistics per set") + theme(axis.text.y = element_text(size = 10), plot.title = element_text(size = 12)) + margin

ggsave(paste(fig_path, "MofaFeatureDetailedRNA.png", sep = ""), rna_mofa_detailed, device = png(), width = 12, height = 4)
# ggsave(paste(fig_path_pdf, "Fig3/", "MofaFeatureDetailedRNA.pdf", sep = ""), rna_mofa_detailed, device = "pdf", width = 12, height = 4)
dev.off()

prot_mofa_detailed <- 
plotEnrichmentDetailed(object = mod, fsea.results = fsea.results_prot, factor = "LF7", feature.sets = reactomeGS, max.pathways = 6, max.genes = 4, text_size = 3.5, alpha = 0.1) + ggtitle("Proteomics: feature statistics per set") + theme(axis.text.y = element_text(size = 10), plot.title = element_text(size = 12)) + margin

ggsave(paste(fig_path_pdf, "Fig3/", "MofaFeatureDetailedProt.pdf", sep = ""), prot_mofa_detailed, device = "pdf", width = 8, height = 5)
dev.off()




gene_mofa_gsea <- 
plotEnrichment(mod, fsea.results = fsea.results_gene, factor = "LF8", max.pathways = 6)+ggtitle("Genomics")

rna_mofa_gsea <- 
plotEnrichment(mod, fsea.results = fsea.results_rna, factor = "LF8", max.pathways = 6)+ggtitle("Transcriptomics")

prot_mofa_gsea <- 
plotEnrichment(mod, fsea.results = fsea.results_prot, factor = "LF7", max.pathways = 6)+ggtitle("Proteomics")



patchwo <- gene_mofa_gsea / rna_mofa_gsea / prot_mofa_gsea
mofa_feature_set <- patchwo + plot_annotation(title = "Multi-omics factor analysis: feature set enrichment test", theme = theme(plot.title = element_text(hjust = 1)))

ggsave(paste(fig_path_pdf, "Fig3/", "MofaFeatureSet.pdf", sep = ""), mofa_feature_set, device = "pdf", width = 11, height = 6)














