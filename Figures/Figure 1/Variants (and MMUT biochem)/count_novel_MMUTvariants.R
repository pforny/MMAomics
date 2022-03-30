# This script is used to compare the detected MMUT variants in the samples 1-150 to the ClinVar database to identify the "novel" variants.
# Potentially other databases can also be used to compare.

# libraries
require(data.table)
require(tidyverse)


system("mkdir Figs/variants")

#####################################################################
# import the clinvar table
# clinvar table was downloaded on 22 Oct 2020 from https://www.ncbi.nlm.nih.gov/clinvar/?term=MMUT%5Bgene%5D
clinvar_tbl <- fread("Data/variants/clinvar_result.csv")

clinvars1 <- clinvar_tbl$Name

# extract mutations (w/o RefSeq etc.)
clinvars <- str_extract(clinvars1, "c\\.(.*)$")
clinvars <- str_split(clinvars," ", simplify = TRUE)[,1]
clinvars <- clinvars[!is.na(clinvars)]


#####################################################################
# import the variants detected in our cohort
# Master table v3 from 23 Oct 2020
variant_tbl <- fread("Data/variants/Master Table MMAomics_v3.csv", skip = 2)
# variant_tbl <- variant_tbl[-(231:nrow(variant_tbl)),]
variant_tbl <- variant_tbl[,-(2:31)]

vars1 <- variant_tbl$"Mutation 1 (nt)"
vars2 <- variant_tbl$"Mutation 2 (nt)"

# limit to the first 150 (MMUT) patients and merge the two vectors and make the individual mutations unique.

varsmix <- c(vars1[1:150], vars2[1:150])
varsmix <- unique(varsmix)
vars <- str_extract(varsmix, "c\\.(.*)$")
vars <- vars[!is.na(vars)]


#####################################################################
# check which variants (vars) are present in the ClinVar (clinvars) list --> the reverse value will indicate how many variants are novel.

notnvl_vars <- vars %in% clinvars
sum(notnvl_vars) # sum of TRUE, indicating the variants which are NOT novel

summary(notnvl_vars)
summary(notnvl_vars)["FALSE"]

# result:
# > summary(notnvl_vars)
#    Mode   FALSE    TRUE 
# logical      54      60


# which variants are novel?
nvl_vars <- vars[!notnvl_vars]
write.table(nvl_vars,file = "Figs/variants/novel_variants_ClinVarComparison.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
