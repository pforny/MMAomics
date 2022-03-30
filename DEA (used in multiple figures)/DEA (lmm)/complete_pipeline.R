#required R packages: Matrix, ggplot2,data.table


source("Code/setupPaths.R")

######################################################
## gets: git@github.com:dlampart/csrproject.git
##
## ## to be extended to automatize as much as possible the retrieval  of external data.
######################################################
source("Code/getExternals.R")

######################################################
#input: Data/20190301_PHRT-5_MMA_Protein-Report_short.csv
#input: Data/gencode.v31lift37.basic.annotation.gff3.gz
#input: Data/HUMAN_9606_idmapping_selected.tab
#input: Data/MMA_PHRT_R_gene_rpkm.bed
######################################################

source("Code/prepare_exp_mats.R")

######################################################
#input: Data/phenotypic_tables/MMA_master_pheno_tbl_1_6_19.csv
#input: Data/phenotypic_tables/MMA_master_pheno_tbl_28_8_19.csvxs
######################################################

source("Code/preprocess_meta_tbl.R")

######################################################
#input: Data/gene_sets/all_enzyme_ids2.tbl
#input: Data/gene_sets/all_enzyme_ids.tbl
#input: Data/gene_sets/mitochondrion_uniprot.txt
#input: Data/gene_sets/Pull_down_MUT_MCEE_MMAB_pooled_edited.csv
######################################################

source("Code/prepare_genesets.R")

######################################################
######################################################

source("Code/diff_expression_analysis_lmm.R")

######################################################
######################################################

source("Code/prepare_diff_expplot.R")

######################################################
######################################################
