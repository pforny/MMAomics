require(data.table)
pheno_tbl0=fread("less Data/phenotypic_tables/MMA_master_pheno_tbl_1_6_19.csv | perl -nle 'if($.>1){print}'")

# alternatively the csv file provided in the extended data of the manuscript can be used.

pheno_tbl=pheno_tbl0[,1:46,with=FALSE]
pheno_tbl=pheno_tbl[!is.na(ID),]

new_names=c("phenoid","mma_id","forny_id","case","gender","n_passage","date_collection","date_freezing","origin_contry","onset_age","onset_age_upper_limit","asymptomatic","metabolic_acidosis","hyper_ammonia","ketosis","hypotonia","hypothermia","dyspnea","lethargy","convulsions","feeding_problems","vomiting","dehydration","dev_delay","irratibility","autism","concurrent_infection","other_clinical","OHCblMinus","OHCblPlus","ratio","mut_category","AdoCblMinus","AdoCblPlus","sanger_mut1nt","sanger_mut1aa","sanger_mut2nt","sanger_mut2aa","sanger_comment","proposed_gene_wgs","wgs_mut1nt","wgs_mut1aa","wgs_mut2nt","wgs_mut2aa","wgs_zygosity","wgs_comment")
setnames(pheno_tbl,new_names)

checking_tbl=data.table(name_new=names(pheno_tbl),name_orig=names(pheno_tbl0)[1:46])## looks good
system("mkdir interimData/phenotypic_tables")
write.table(pheno_tbl,file="interimData/phenotypic_tables/meta_info_1_6_19.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
##################################################
##################################################
mmut_id="ENSG00000146085"

load("interimData/rnaseq_mat.RData")
load("interimData/prot_mat.RData")
mut_rnaseq=rnaseq_mat[rownames(rnaseq_mat)==mmut_id,]
mut_prot=prot_mat[rownames(prot_mat)==mmut_id,]
ee=match(pheno_tbl[,mma_id], names(mut_prot))
pheno_tbl[,prot_mut_level:=NULL]
pheno_tbl[,prot_mut_level:=1]
prot_mut2=mut_prot[ee[!is.na(ee)]]
pheno_tbl[!is.na(ee),prot_mut_level:=unlist(prot_mut2)]
pheno_tbl[is.na(ee),prot_mut_level:=NA]


ee2=match(pheno_tbl[,mma_id], names(mut_rnaseq))
pheno_tbl[,rnaseq_mut_level:=NULL]
pheno_tbl[,rnaseq_mut_level:=1]
rnseq_mut2=mut_rnaseq[ee2[!is.na(ee2)]]
pheno_tbl[!is.na(ee2),rnaseq_mut_level:=rnseq_mut2]
pheno_tbl[is.na(ee2),rnaseq_mut_level:=NA]

write.table(pheno_tbl, file="interimData/phenotypic_tables/meta_info_combined.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
##################################################
##################################################

require(data.table)
pheno_tbl0=fread("less Data/phenotypic_tables/MMA_master_pheno_tbl_28_8_19.csv | perl -nle 'if($.>2){print}'")

pheno_tbl=pheno_tbl0[,1:45,with=FALSE]
pheno_tbl=pheno_tbl[`Cell Line Number`!="",]

new_names=c("mma_id","forny_id","case","gender","n_passage","date_collection","date_freezing","sample_origin_contry","patient_origin_contry","consanguinity","onset_age","onset_age_upper_limit","OHCblMinus","OHCblPlus","ratio","SimultOHCblMinus","SimultOHCblPlus","mut_category","init_mut_category","AdoCblMinus","AdoCblPlus","SimultAdoCblMinus","SimultAdoCblPlus","wgs_qc","rnaseq_qc","prot_qc","sanger_mut1nt","sanger_mut1aa","sanger_mut2nt","sanger_mut2aa","sanger_comment","proposed_gene_wgs","wgs_mut1nt","wgs_mut1aa","wgs_mut1type","wgs_mut2nt","wgs_mut2aa","wgs_mut2type","wgs_bothtype","wgs_zygosity","knownornew","foundbypf","causefound","sanger_confirmed","rna_expression")
setnames(pheno_tbl,new_names)

checking_tbl=data.table(name_new=names(pheno_tbl),name_orig=names(pheno_tbl0)[1:45])## looks good
write.table(pheno_tbl,file="interimData/phenotypic_tables/meta_info_28_8_19.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
##################################################
mmut_id="ENSG00000146085"

load("interimData/rnafiltered_mat.RData")
load("interimData/prot_mat.RData")
mut_rnaseq=rnaseq_mat[rownames(rnaseq_mat)==mmut_id,]
mut_prot=prot_mat[rownames(prot_mat)==mmut_id,]
ee=match(pheno_tbl[,mma_id], names(mut_prot))
pheno_tbl[,prot_mut_level:=NULL]
pheno_tbl[,prot_mut_level:=1]
prot_mut2=mut_prot[ee[!is.na(ee)]]
pheno_tbl[!is.na(ee),prot_mut_level:=unlist(prot_mut2)]
pheno_tbl[is.na(ee),prot_mut_level:=NA]


ee2=match(pheno_tbl[,mma_id], names(mut_rnaseq))
pheno_tbl[,rnaseq_mut_level:=NULL]
pheno_tbl[,rnaseq_mut_level:=1]
rnseq_mut2=mut_rnaseq[ee2[!is.na(ee2)]]
pheno_tbl[!is.na(ee2),rnaseq_mut_level:=rnseq_mut2]
pheno_tbl[is.na(ee2),rnaseq_mut_level:=NA]
pheno_tbl[,rna_expression:=NULL]
write.table(pheno_tbl, file="interimData/phenotypic_tables/meta_info_combined28_8_19.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


##################################################
##################################################
require(data.table)
pheno_tbl0=fread("less Data/phenotypic_tables/MMA_master_pheno_tbl_2_12_20.csv | perl -nle 'if($.>2){print}'")

pheno_tbl=pheno_tbl0[,1:49,with=FALSE]
pheno_tbl=pheno_tbl[`Cell Line Number`!="",]

new_names=c("mma_id","forny_id","case","gender","n_passage","date_collection","date_freezing","sample_origin_contry", "sample_origin_town","patient_origin_contry","consanguinity","onset_age","onset_age_upper_limit","OHCblMinus","OHCblPlus","ratio","SimultOHCblMinus","SimultOHCblPlus","mut_category","init_mut_category","AdoCblMinus","AdoCblPlus","SimultAdoCblMinus","SimultAdoCblPlus","wgs_qc","rnaseq_qc","prot_qc","sanger_mut1nt","sanger_mut1aa","sanger_mut2nt","sanger_mut2aa","sanger_comment","proposed_gene_wgs","incomplete_gene_wgs","wgs_mut1nt","wgs_mut1aa","wgs_mut1type","wgs_mut2nt","wgs_mut2aa","wgs_mut2type","wgs_bothtype","wgs_zygosity","knownornew","foundbyCNV","foundbypf","causefound","sanger_confirmed","prot_mut_level","rna_expression")
setnames(pheno_tbl,new_names)

checking_tbl=data.table(name_new=names(pheno_tbl),name_orig=names(pheno_tbl0)[1:49])## looks good
write.table(pheno_tbl,file="interimData/phenotypic_tables/meta_info_2_12_20.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

##################################################
##################################################
require(data.table)
require(readxl)

pheno_tbl0=data.table(read_excel("Data/phenotypic_tables/MMA_master_pheno_tbl_15_7_21_CF_PF.xlsx", skip = 3))

pheno_tbl=pheno_tbl0[,1:ncol(pheno_tbl0),with=FALSE]
pheno_tbl=pheno_tbl[`Cell Line Number`!="",]

new_names=c("mma_id","forny_id","case","gender","n_passage","date_collection","date_freezing","sample_origin_contry", "sample_origin_town","patient_origin_contry","consanguinity","onset_age","onset_age_upper_limit","OHCblMinus","OHCblPlus","ratio","SimultOHCblMinus","SimultOHCblPlus","mut_category","init_mut_category","AdoCblMinus","AdoCblPlus","SimultAdoCblMinus","SimultAdoCblPlus","wgs_qc","rnaseq_qc","prot_qc","sanger_mut1nt","sanger_mut1aa","sanger_mut2nt","sanger_mut2aa","sanger_comment","proposed_gene_wgs","incomplete_gene_wgs","ncmut1_NCBIrefSeq", "mut1_mutalyzer", "wgs_mut1nt","wgs_mut1aa","wgs_mut1aa_short","wgs_mut1type","ncmut2_NCBIrefSeq", "mut2_mutalyzer","wgs_mut2nt","wgs_mut2aa","wgs_mut2aa_short", "wgs_mut2type", "comment_wgs", "wgs_bothtype","wgs_zygosity","knownornew","foundbyCNV","foundbypf","causefound","sanger_confirmed","prot_mut_level","rna_expression")
setnames(pheno_tbl,new_names)

checking_tbl=data.table(name_new=names(pheno_tbl),name_orig=names(pheno_tbl0)[1:ncol(pheno_tbl0)])## looks good
write.table(pheno_tbl,file="interimData/phenotypic_tables/meta_info_15_7_21.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

##################################################
# add MMUT protein and RNA levels
mmut_id="ENSG00000146085"

load("interimData/rnafiltered_mat.RData")
load("interimData/prot_mat.RData")
mut_rnaseq=rnafiltered_mat[rownames(rnafiltered_mat)==mmut_id,]
mut_prot=prot_mat[rownames(prot_mat)==mmut_id,]
ee=match(pheno_tbl[,mma_id], names(mut_prot))
pheno_tbl[,prot_mut_level:=NULL]
pheno_tbl[,prot_mut_level:=1]
prot_mut2=mut_prot[ee[!is.na(ee)]]
pheno_tbl[!is.na(ee),prot_mut_level:=unlist(prot_mut2)]
pheno_tbl[is.na(ee),prot_mut_level:=NA]


ee2=match(pheno_tbl[,mma_id], names(mut_rnaseq))
pheno_tbl[,rna_expression:=NULL]
pheno_tbl[,rnaseq_mut_level:=NULL]
pheno_tbl[,rnaseq_mut_level:=1]
rnseq_mut2=mut_rnaseq[ee2[!is.na(ee2)]]
pheno_tbl[!is.na(ee2),rnaseq_mut_level:=rnseq_mut2]
pheno_tbl[is.na(ee2),rnaseq_mut_level:=NA]


##################################################
# add clinical information
clin_info <- fread("Data/phenotypic_tables/Clinical_information_v3_2_12_20.csv", skip = 5, drop = c(1, 3:27))
pheno_tbl <- merge(pheno_tbl, clin_info, by = "mma_id")

write.table(pheno_tbl, file="interimData/phenotypic_tables/meta_info_combined2_12_20.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

##################################################
## adding severity score
##################################################
scoreTbl0=fread("Data/phenotypic_tables/Final\ Clinical\ Severity\ Score_v2.csv")
scoreTbl1 = scoreTbl0[,list(mma_id=MMAid, category, age_of_onset_d, acCSS,acnCSS, cCSS,cnCSS ,neuro_abnormalitites_chr,global_dev_delay_chr,hypotonia_chr,basal_ganglia_abnormality_chr,nCSSonlyNeuro=neuro_abnormalities_nCSS,kidney_impairment,haemat_abnormalitites,failure_to_thrive_or_tube_feeding)]

scoreTbl1[,cssModif:=cCSS] #alternative scores: acCSS, acnCSS, cCSS, cnCSS; a, acute; c, chronic; n, neuro
scoreTbl1[cssModif>2,cssModif:=2]

setkey(scoreTbl1, mma_id)
setkey(pheno_tbl, mma_id)
pheno_tbl_css = merge(pheno_tbl,scoreTbl1)

write.table(pheno_tbl_css, file="interimData/phenotypic_tables/meta_info_wcss.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

