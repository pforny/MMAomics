# libraries
require(data.table)
require(ggplot2)
require(ggsci)
require(tidyverse)
library(ggridges)
library(lubridate)
require(ggpubr)
library(ggrepel)
library(pheatmap)
library(waffle)
library(viridis)
library(ggplotify)


tbl=fread("interimData/phenotypic_tables/meta_info_wcss.txt")

#colour_palette
##blue 
#82
#129
#163
myblue=rgb(60, 107, 131, maxColorValue = 255)
#orange
#231
#129
#48
myorange=rgb(231, 129, 48, maxColorValue = 255)
#brown
#180
#100
#36
mybrown=rgb(180, 100, 36, maxColorValue = 255)
#mygrey
#212
#212
#212
mygrey=rgb(180, 180, 180, maxColorValue = 255)


system("mkdir Figs/clinical")
fig_path <- c("Figs/clinical/")
fig_path_pdf <- c("Figs/v12/pdf/")
mypal1 <- pal_aaas("default", alpha = 1)(9)
mypal2 <- c(mypal1[2], mypal1[4], mypal1[1])



##########################################
## cohort characterization
##########################################

tbl_coll <- tbl
tbl_coll[1:150, type := "MMUT def."]
tbl_coll[151:210, type := "unknown"]
tbl_coll[211:230, type := "unaffected"]
tbl_coll$type <- factor(tbl_coll$type, levels = c("MMUT def.", "unknown", "unaffected"))


tbl[,proposed_gene_ismut:="NOT MMUT"]
tbl[proposed_gene_wgs=="MMUT",proposed_gene_ismut:="MMUT"]
tbl[, proposed_gene_ismut]
tbl[,log(OHCblPlus)]



sample_year_plt <- 
ggplot(tbl_coll, aes(x = date_collection, fill = type)) +
	geom_bar() +
	xlab("Year at collection") +
	ylab("Number of samples") +
	scale_fill_manual(values = mypal2) +
	theme_pubr() +
	theme(legend.position = "right", legend.title = element_blank())


tbl_coll20 <- tbl_coll %>% select(type)
tbl_coll2 <- as.data.table(table(tbl_coll20))
setorder(tbl_coll2, -N)
setnames(tbl_coll2, c("type", "N"))
tbl_coll2$type <- factor(tbl_coll2$type, levels = c("MMUT def.", "unknown", "unaffected"))


sample_waffle_plt <- 
ggplot(tbl_coll2, aes(values = N, fill = type)) + 
	geom_waffle(flip = TRUE, color = "white") +
  	scale_fill_manual(values = mypal2) +
  	ylab("Number of samples") +
 	theme_pubr() +
	scale_y_continuous(labels = function(x) x * 10, expand = c(0,0)) +
	theme(legend.title = element_blank(), legend.position = "right", axis.line.y = element_blank()) +
	guides(x = "none", fill = guide_legend(reverse = TRUE))

sample_plt <- 
ggarrange(sample_year_plt, sample_waffle_plt, nrow = 2, common.legend = TRUE)

ggsave(paste(fig_path, "sample_plt.png", sep = ""), sample_plt, width = 3.5, height = 5.5, bg = "white")
dev.off()
ggsave(paste(fig_path_pdf, "SuppFig1_qc/", "sample_plt.pdf", sep = ""), sample_plt, device = "pdf", width = 3.5, height = 5.5, bg = "white")
dev.off()



##########################################
## show that OHCblPlus can act as a proxy for mut.
##########################################

summary(lm(tbl[,log10(OHCblPlus)]~tbl[,proposed_gene_ismut=="MMUT"]))

OHCblPlus_proxy_plot <- 
ggplot(data=tbl,aes(x=log10(OHCblPlus),fill=proposed_gene_ismut)) +
	geom_histogram() +
	facet_grid(proposed_gene_ismut~.) +
	guides(fill=FALSE) +
	ggtitle("OHCblPlus as MMUT proxy") +
	theme_test() +
	scale_fill_aaas()
ggsave(paste(fig_path,"OHCblPlusProxyForMMUT.png", sep = ""),width=4,height=3)
#ggplot(data=tbl,aes(x=log10(OHCblPlus),fill=proposed_gene_ismut,colour=proposed_gene_ismut,group=proposed_gene_ismut))+geom_histogram(alpha=0.5,position="identity")# +geom_histogram()#+scale_fill_manual(values=c(myblue,myorange))+theme_bw()+facet_grid(proposed_gene_ismut~.)+ guides(fill=FALSE)

ggplot(data=tbl,aes(x=log10(OHCblPlus),fill=proposed_gene_ismut,group=proposed_gene_ismut)) +
	geom_histogram(alpha=0.65,position="identity",bins=20) +
	theme_test() +
	scale_fill_aaas() +
	theme(legend.position = c(0.75, 0.8)) +
	guides(fill=guide_legend(title=NULL))
summary(lm(tbl[,log10(OHCblPlus)]~tbl[,proposed_gene_ismut=="MMUT"]))
ggsave(paste(fig_path,"OHCblPlusProxyForMMUT_v2.png", sep = ""),width=4,height=3)
# +geom_histogram()#+scale_fill_manual(values=c(myblue,myorange))+theme_bw()+facet_grid(proposed_gene_ismut~.)+ guides(fill=FALSE)

##########################################
## shows that OHCblPlus predictable from clinical symptoms.
##########################################
tbl1=tbl[,has_any_other_clinical_symptom:=!is.element(other_clinical_findings_atypical,c("","-"))]
tbl1[is.na(metabolic_acidosis),metabolic_acidosis:=0]
tbl1[is.na(acidosis), acidosis:=0]
tbl1[is.na(any_delay), any_delay:=0]
tbl1[is.na(global_dev_delay_chronic), global_dev_delay_chronic:=0]
tbl1[is.na(hyperammonemia), hyperammonemia:=0]
tbl1[is.na(abnormal_muscle_tone), abnormal_muscle_tone:=0]
tbl1[is.na(any_neurological_abnormalities_chronic), any_neurological_abnormalities_chronic:=0]
tbl1[is.na(any_neurological_abnormalities_acuteANDchronic), any_neurological_abnormalities_acuteANDchronic:=0]

tbl1[,clinically_atypical:=(acidosis==0 & hyperammonemia==0 & (has_any_other_clinical_symptom + any_delay)>0)]

tbl1[,`is_clinically_atypical`:="clinically typical"]
tbl1[clinically_atypical==TRUE,`is_clinically_atypical`:="clinically atypical"]
tbl1$is_clinically_atypical <- factor(tbl1$is_clinically_atypical, levels = c("clinically typical", "clinically atypical"))
summary(lm(tbl1[,log10(OHCblPlus)] ~tbl1[,is_clinically_atypical]))
## results of different combinations
# any_delay, acidosis, has_any_other_clinical_symptom		R2 = .25, p 1.3e-15
# any_delay, metabolic_acidosis, has_any_other_clinical_symptom		R2 = .22, p 1.7e-13
# global_dev_delay_chronic, acidosis, has_any_other_clinical_symptom		R2 = .15, p 1.8e-9
# any_delay, acidosis		R2 = .12, p 5.7e-8

atypical_plot <- 
ggplot(tbl1,aes(x=log10(OHCblPlus),group=is_clinically_atypical,fill=is_clinically_atypical)) +
	scale_fill_aaas() +
	geom_histogram(aes(y=..density..)) +
	theme_test() +
	facet_grid(is_clinically_atypical~.) +
	ggtitle("Typical have low OHCblPlus") +
	guides(fill=FALSE)
ggsave(paste(fig_path,"OHCblPlusPredictableFromClinical.png", sep = ""),width=4,height=3)

ggplot(tbl1,aes(x=log10(OHCblPlus),group=is_clinically_atypical,fill=is_clinically_atypical)) +
	scale_fill_aaas() +
	geom_histogram(alpha=0.65,position="identity",bins=20) +
	theme_test() +
	theme(legend.position = c(0.75, 0.8)) +
	guides(fill=guide_legend(title=NULL))
ggsave(paste(fig_path,"OHCblPlusPredictableFromClinical_v2.png", sep = ""),width=4,height=3)

##########################################
tbl2=tbl[proposed_gene_wgs=="MMUT",]
tbl2[,AdoCblPlusNumeric:=as.numeric(AdoCblPlus)]
myresids=summary(lm(tbl2[,log10(OHCblPlus)]~tbl2[,log10(OHCblMinus)]))$residuals
tbl2[,pwresid:=myresids]
tbl3=tbl2[!is.na(AdoCblPlusNumeric),]
##########################################
## outliers w.r.t. OHCblPlus activity tend to have residual enzymatic activity.
##########################################
tbl3[,AdoCblPlusBelow100:=AdoCblPlusNumeric<100]
tbl3[AdoCblPlusBelow100==TRUE, AdoCblPlusBelow100_txt:="AdoCblPlus<100"]
tbl3[AdoCblPlusBelow100==FALSE, AdoCblPlusBelow100_txt:="AdoCblPlus>100"]
summary(lm(tbl3[,log10(OHCblPlus)] ~tbl3[,AdoCblPlusBelow100]))

OHCblPlus_outliers_plot <-
ggplot(data=tbl3,aes(x=log(OHCblPlus),fill=AdoCblPlusBelow100_txt)) +
	geom_histogram() +
	facet_grid(AdoCblPlusBelow100_txt~.) +
	scale_fill_aaas() +
	theme_test() +
	guides(fill=FALSE)
ggsave(paste(fig_path,"OHCBlOutliers.png", sep = ""), OHCblPlus_outliers_plot,width=4,height=3)



##########################################
## low expression (DRIVEN BY NMD!) seemingly makes it unlikely that disease is late onset.(only mut cases considered). Association not strong.
##########################################
tbl4=tbl2[!is.na(onset_age),]
tbl4[,onset_category:="between"]
tbl4[onset_age>120,onset_category:="above 4 months"]
tbl4[onset_age<4,onset_category:="below 4 days"]
summary(lm(tbl4[,onset_category=="above 4 months"] ~tbl4[,onset_category=="below 4 days"]))

ggplot(data=tbl4[onset_category!="between" ,],aes(x=rnaseq_mut_level,fill=onset_category)) +
	geom_histogram(bins=20) +
	facet_grid(onset_category~.) +
	scale_fill_aaas() +
	theme_test() +
	guides(fill=FALSE)

NMD_drives_AoO_plot <- 
ggplot(data=tbl4[onset_category!="between" ,],aes(x=rnaseq_mut_level,fill=onset_category, color=onset_category)) +
	geom_density(alpha = 0.6) +
	scale_fill_aaas() +
	scale_color_aaas() +
	ggtitle("NMD drives age of onset") +
	theme_test()
ggsave(paste(fig_path,"NMD_driving_earlyOnset_onlyMMUT.png", sep = ""),NMD_drives_AoO_plot, width=4,height=2)

##########################################
## shows NMD on rnaseq.
##########################################
#tbl2 is table with MMUT only
tbl2[,mut2_nonsense:="not_nonsense"]
tbl2[grepl("\\*",wgs_mut2aa) | grepl("X",wgs_mut2aa) | grepl("Ter",wgs_mut2aa),mut2_nonsense:="nonsense"]
ggplot(data=tbl2[wgs_zygosity=="hom",],aes(x=rnaseq_mut_level,fill=mut2_nonsense)) +
	geom_histogram(bins=20) +
	facet_grid(mut2_nonsense~.) +
	scale_fill_aaas() +
	theme_test() +
	guides(fill=FALSE)

tbl2[, mut_nonsense:=NULL]
tbl2[, mut_nonsense:="no nonsense"]
tbl2[wgs_mut1type == "truncating" | wgs_mut2type == "truncating", mut_nonsense:="1 nonsense"]
tbl2[wgs_bothtype=="truncating/truncating", mut_nonsense:="2 nonsense"]
tbl2$mut_nonsense <- factor(tbl2$mut_nonsense, levels = c("2 nonsense", "1 nonsense", "no nonsense"))
# ggplot(data=tbl2,aes(x=rnaseq_mut_level,fill=mut_nonsense)) +
# 	geom_histogram(bins=20) +
# 	facet_grid(mut_nonsense~.) +
# 	scale_fill_aaas() +
# 	theme_test() +
# 	guides(fill=FALSE)
ggplot(data=tbl2,aes(x=rnaseq_mut_level,fill=mut_nonsense)) +
	geom_histogram(bins = 30, position = "identity", alpha = 0.6) +
	scale_fill_aaas() +
	theme_test()

NMD_mut_type_plt <-
ggplot(data=tbl2,aes(x=rnaseq_mut_level,fill=mut_nonsense, color = mut_nonsense)) +
	geom_density(alpha = 0.6) +
	scale_fill_aaas() +
	scale_colour_aaas() +
	ggtitle("NMD by mutation type") +
	theme_test()
ggsave(paste(fig_path,"NMD_mutation_type_onlyMMUT.png", sep = ""),NMD_mut_type_plt, width=4,height=2)


########################################
# compile panel of the above plots

clin1_plots <- ggarrange(OHCblPlus_proxy_plot, NMD_drives_AoO_plot, atypical_plot, NMD_mut_type_plt, widths = c(1, 1.3))
ggsave(paste(fig_path,"clinical1_biochem.png", sep = ""),clin1_plots,width=8,height=6)




##########################################
## p-values of clinical associations 
##########################################
#ee=fread("interimData/phenotypic_tables/ClinicalInformationTable.csv",header=TRUE,sep=",")
ee00=fread("interimData/phenotypic_tables/meta_info_wcss.txt",header=TRUE)
# matcher=match(tbl2[,mma_id],ee[,`Cell Line Number`])
# ee2=ee[matcher,]

# descriptive text variables or not meaningful (e.g. due to redundancy) variables to correlate to; they are excluded from the ee table
exclude_vars <- c(
	"mma_id",
	"forny_id",
	"MMA_unknown_units",
	"MMA2_urine_after_treat",
	"ammonia_ugDL",
	"glucose_mgDL",
	"onset_age_upper_limit",
	"init_mut_category",
	"wgs_qc",
	"rnaseq_qc",
	"prot_qc",
	"sanger_mut1nt",
	"sanger_mut1aa",
	"sanger_mut2nt",
	"sanger_mut2aa",
	"sanger_comment",
	"incomplete_gene_wgs",
	"knownornew",
	"foundbyCNV",
	"foundbypf",
	"causefound",
	"sanger_confirmed",
	"category",
	"age_of_onset_d",
	"acCSS",
	"acnCSS",
	"cnCSS",
	"other_clinical_findings",
	"other_clinical_findings_atypical",
	"comments_diagnostics",
	"comments_treatment",
	"comment_responsiveness",
	"age_at_follow_up_treatment",
	"antibiotic_treatment_type",
	"age_at_follow_up",
	"weight",
	"weight_centile",
	"length",
	"length_centile",
	"head_circumference",
	"head_circum_centile",
	"global_dev_delay_chronic",
	"comments_milestones",
	"comments_MRI",
	"any_neurological_abnormalities_acuteANDchronic",
	"any_neurological_abnormalities_acuteANDchronic_text",
	"neurological_summary",
	"comments_cardiomyopathy",
	"impaired_kidney_fct_text",
	"hemat_abnormality_text",
	"comments_overall",
	"failure_to_thrive_OR_tube_feeding",
	"kidney_impairment",
	"haemat_abnormalitites",
	"neuro_abnormalitites_chr",
	"global_dev_delay")

ee00 <- data.frame(ee00)
ee1 <- ee00 %>%	select(-all_of(exclude_vars))

# set all clinical symptoms (rated as present or absent) to 0, if not reported to be present

# correct consanguinity variable
ee1[c(59, 64, 154), ]$consanguinity <- 1
ee1$consanguinity <- as.numeric(ee1$consanguinity)

pres_abs_vars <- c("consanguinity","acidosis", "metabolic_acidosis", "metabolic_ketoacidosis","ketosis","hyperammonemia","abnormal_muscle_tone","musc_hypotonia","musc_hypertonia","hypothermia","fct_respiratory_abnormality" , "dyspnea","tachypnea","hyperventilation","reduced_consciousness","lethargy","drowsiness","somnolence","coma","seizures","general_tonic_clonic_seizure","any_GI_problem","feeding_problem","failure_to_thrive","vomiting","dehydration","any_delay","motor_delay","behavioral_abnormality","irritability","autistic_behaviour","concurrent_infection","urine_ketones","dialysis", "hemodialysis","peritoneal_dialysis","insulin","diet","carnitine","cobalamin","bicarb","glucose_IV", "responsive_to_acute_treatment","cobalamin_responsive","antibiotic_treatment","protein_restriction","tube_feeding_day","tube_feeding_night","tube_feeding_overall", "growth_hormone","weight_centile_quant","length_centile_quant","head_circumfernce_quant","language_delay", "any_neurological_abnormalities_chronic","cardiomyopathy","impaired_kidney_fct", "hemat_abnormality","anemia","neutropenia","pancytopenia","skin_abnormalities","hearing_impairment","osteoporosis","failure_to_thrive_chronic", "global_dev_delay_chr","hypotonia_chr","basal_ganglia_abnormality_chr","failure_to_thrive_or_tube_feeding")

ee1[, pres_abs_vars][is.na(ee1[pres_abs_vars])] <- 0 # replace the NA's with 0

# delete some variables as only very few (<3) patients positive
ee1 <- subset(ee1, select = -c(cardiomyopathy))
ee1 <- subset(ee1, select = -c(drowsiness))
ee1 <- subset(ee1, select = -c(motor_delay))
ee1 <- subset(ee1, select = -c(autistic_behaviour))
ee1 <- subset(ee1, select = -c(hemodialysis))
ee1 <- subset(ee1, select = -c(pancytopenia))
ee1 <- subset(ee1, select = -c(growth_hormone))


# continuous variables; they need to be changed to numeric or checked whether they are numeric
cont_vars <- c(
	"n_passage",
	"onset_age",
	"OHCblMinus",
	"OHCblPlus",
	"ratio",
	"SimultOHCblMinus",
	"SimultOHCblPlus",
	"AdoCblMinus",
	"AdoCblPlus",
	"SimultAdoCblMinus",
	"SimultAdoCblPlus",
	"prot_mut_level",
	"rnaseq_mut_level",
	"MMA_urine",
	"MMA_plasma",
	"ammonia_umolL",
	"pH",
	"base_excess",
	"pCO2",
	"glucose_mmolL",
	"MMA_urine_after_treat",
	"cobalamin_dose_PO",
	"cobalamin_dose_IM",
	"carnitine_dose",
	"natural_protein_amount",
	"synthetic_protein_amount",
	"total_protein_amount",
	"meals_per_day",
	"calories_per_day",
	"IQ",
	"uric_acid",
	"creatinine",
	"urea",
	"GFR",
	"cCSS",
	"cssModif",
	"nCSSonlyNeuro")

summary(ee1[, cont_vars])  # check if there are NAs
sapply(ee1[, cont_vars], class) # all values are numeric or integer
ee1[, cont_vars] <- apply(ee1[, cont_vars], 2, function(x) as.numeric(as.character(x)))

# time/date variables
date_vars <- c("date_collection", "date_freezing")
ee1[, date_vars]
sapply(ee1[, date_vars], class)
ee1$date_collection
ee1$date_freezing <- dmy(ee1$date_freezing)
ee1$date_freezing <- substring(ee1$date_freezing,1,4)
ee1$date_freezing <- as.integer(ee1$date_freezing)
ee1$date_freezing
sapply(ee1[, date_vars], class)


# special variables
special_vars <- c("case", "gender", "sample_origin_contry", "sample_origin_town", "patient_origin_contry", "mut_category", "proposed_gene_wgs", "wgs_mut1nt", "wgs_mut1aa", "wgs_mut1type", "wgs_mut2nt","wgs_mut2aa","wgs_mut2type","wgs_bothtype","wgs_zygosity")

head(ee1[, special_vars])
sapply(ee1[, special_vars], class)
# case only two levels 1 and 0, leave as is
ee1$case
# gender only two levels F and M, change to 1(=F) and 0(=M)
ee1 <- data.table(ee1)
ee1[gender == "F", ]$gender <- 1
ee1[gender == "M", ]$gender <- 0
ee1$gender <- as.numeric(ee1$gender)
# mut_category variable needs to be changed to a numeric variable, mut0(=0), mut-(=1)
ee1[mut_category == ""]$mut_category <- NA
ee1[mut_category == "mut0"]$mut_category <- 0
ee1[mut_category == "mut-"]$mut_category <- 1
ee1$mut_category <- as.numeric(ee1$mut_category)
# change wgs_zygosity variable to hom(=0), het(=1)
ee1[wgs_zygosity == ""]$wgs_zygosity <- NA
ee1[wgs_zygosity == "hom"]$wgs_zygosity <- 0
ee1[wgs_zygosity == "het"]$wgs_zygosity <- 1
ee1$wgs_zygosity <- as.numeric(ee1$wgs_zygosity)
# change IQ variable so that two patients with IQ 20 are excluded and not driving the correlation
ee1[IQ == 20]$IQ <- NA
ee1$IQ

# some variables are factor variables, which must be excluded from the analysis of numeric variables
factor_vars <- c("sample_origin_contry", "sample_origin_town", "patient_origin_contry", "proposed_gene_wgs", "wgs_mut1nt", "wgs_mut1aa", "wgs_mut1type", "wgs_mut2nt", "wgs_mut2aa", "wgs_mut2type", "wgs_bothtype", "ncmut1_NCBIrefSeq", "mut1_mutalyzer", "wgs_mut1aa_short", "ncmut2_NCBIrefSeq", "mut2_mutalyzer", "wgs_mut2aa_short", "comment_wgs")

ee1 <- data.frame(ee1)
ee1 <- ee1 %>%	select(-all_of(factor_vars))

ee2 <- data.table(ee1)


# save list of clinical/biochemical variables
vars_list <- colnames(ee2)
write.table(vars_list, file = "Results/listOfVariables_numeric.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

vars_list_exp <- data.table(vars_list, var_category = c(
	"clinical",
	"clinical",
	"biochemical",
	"biochemical",
	"biochemical",
	"clinical",
	"clinical",
	"biochemical",
	"biochemical",
	"biochemical",
	"biochemical",
	"biochemical",
	"biochemical",
	"biochemical",
	"biochemical",
	"biochemical",
	"biochemical",
	"multi-omic",
	"multi-omic",
	"multi-omic",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical/biochemical",
	"clinical/biochemical",
	"clinical/biochemical",
	"clinical/biochemical",
	"clinical/biochemical",
	"clinical/biochemical",
	"clinical/biochemical",
	"clinical/biochemical",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/biochemical",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical/treatment",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical/biochemical",
	"clinical/biochemical",
	"clinical/biochemical",
	"clinical/biochemical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical",
	"clinical"
	)
)

print(vars_list_exp, nrow = 119)
write.csv(vars_list_exp, file = "Results/listOfVariables_categorized.csv", row.names=FALSE, col.names = TRUE, quote=FALSE)

vars_list_exp_no <- data.table(table(vars_list_exp$var_category))
setnames(vars_list_exp_no, c("category", "N"))
vars_list_exp_no$category <- factor(vars_list_exp_no$category , levels = c("clinical", "clinical/treatment", "clinical/biochemical", "biochemical", "multi-omic"))
levels(vars_list_exp_no$category) <- c("Symptoms", "Treatment", "Clinical chemistry", "Biochemical", "MMUT multi-omic")
write.csv(vars_list_exp_no, file = "Results/listOfVariables_categorized_numbers.csv", row.names = FALSE, col.names = TRUE, quote = FALSE)
vars_list_exp_no[, N_alt := paste("N=", vars_list_exp_no$N, sep = "")]

vars_plt <- 
ggplot(vars_list_exp_no, aes(x = category, y = N)) +
	geom_col() +
	geom_text(aes(x = category, y = N+4, label = N_alt), size = 3) +
	ylab("Count") +
	ggtitle("Phenotype variable category") +
	theme_pubr() +
	rotate_x_text(angle = 45) +
	theme(axis.title.x = element_blank(), plot.title = element_text(size = 12))

ggsave(paste(fig_path,"PhenotypeVarialbes_categories.png", sep = ""),vars_plt,width = 2.8, height = 3)
dev.off()
ggsave(paste(fig_path_pdf,"Fig2/","PhenotypeVarialbes_categories.pdf", sep = ""), vars_plt, device = "pdf", width = 2.8, height = 3)
dev.off()

ggplot(data.table(table(sapply(ee2, class))), aes(x = V1, y = N)) +
	geom_col() +
	ylab("Count") +
	ggtitle("Phenotype variable type") +
	theme_pubr() +
	rotate_x_text(angle = 45) +
	theme(axis.title.x = element_blank(), plot.title = element_text(size = 12))

ggsave(paste(fig_path,"PhenotypeVarialbes_types.png", sep = ""),width = 2.8, height = 3)
dev.off()



### add 0.01 to AdoCbl activity as zero cannot be log10
tbl[AdoCblPlus==0, ]$AdoCblPlus <- 0.01
tbl[AdoCblMinus==0, ]$AdoCblMinus <- 0.01
tbl[, AdoCblPlus]



ff=lapply(c(1:dim(ee2)[2]),function(i){
	print(i)
	ee3=as.numeric(unlist(ee2[,i,with=FALSE]))### remove NA's and empty
	ee4=ee3
	myna=is.na(ee4)
	if(sum(!myna)<2 || sd(ee3,na.rm=TRUE)==0){
		return(NULL)
	}
	xx=scale(ee4[!myna])

	yy=tbl[!myna,scale(log10(MMA_urine))]
	tt=summary(lm(yy~xx))
	pvalmmaurine=tt$coef[2,4]
	ESmmaurine=tt$coef[2,1]
	yy=tbl[!myna,scale(log10(ammonia_umolL))]
	tt=summary(lm(yy~xx))
	pvalammonia=tt$coef[2,4]
	ESammonia=tt$coef[2,1]
	yy=tbl[!myna,scale(log10(pH))]
	tt=summary(lm(yy~xx))
	pvalph=tt$coef[2,4]
	ESph=tt$coef[2,1]
	yy=tbl[!myna,scale(log10(onset_age))]
	tt=summary(lm(yy~xx))
	pvalonset=tt$coef[2,4]
	ESonset=tt$coef[2,1]
	yy=tbl[!myna,scale(log10(prot_mut_level))]
	tt=summary(lm(yy~xx))
	pvalprotmut=tt$coef[2,4]
	ESprotmut=tt$coef[2,1]
	yy=tbl[!myna,scale(log10(rnaseq_mut_level))]
	tt=summary(lm(yy~xx))
	pvalrnamut=tt$coef[2,4]
	ESrnamut=tt$coef[2,1]
	yy=tbl[!myna,scale(log10(AdoCblMinus))]
	tt=summary(lm(yy~xx))
	pvaladominus=tt$coef[2,4]
	ESadominus=tt$coef[2,1]
	yy=tbl[!myna,scale(log10(AdoCblPlus))]
	tt=summary(lm(yy~xx))
	pvaladoplus=tt$coef[2,4]
	ESadoplus=tt$coef[2,1]
	yy=tbl[!myna,scale(log10(OHCblMinus))]
	tt=summary(lm(yy~xx))
	pvalminus=tt$coef[2,4]
	ESminus=tt$coef[2,1]
	yy=tbl[!myna,scale(log10(ratio))]
	tt=summary(lm(yy~xx))
	pvalrat=tt$coef[2,4]
	ESrat=tt$coef[2,1]
	yy=tbl[!myna,scale(log10(OHCblPlus))]
	tt=summary(lm(yy~xx))
	pvalplus=tt$coef[2,4]
	ESplus=tt$coef[2,1]
	yy=tbl[!myna,scale(cCSS)]
	tt=summary(lm(yy~xx))
	pvalccss=tt$coef[2,4]
	ESccss=tt$coef[2,1]
	data.table(
		MMA_urine=pvalmmaurine,
		ammonia_umolL=pvalammonia,
		pH=pvalph,
		onset_age=pvalonset,
		prot_mut_level=pvalprotmut,
		rnaseq_mut_level=pvalrnamut,
		AdoCblMinus=pvaladominus,
		AdoCblPlus=pvaladoplus,
		OHCblMinus=pvalminus,
		ratio=pvalrat, 
		OHCblPlus=pvalplus,
		cCSS=pvalccss,
		ES_MMA_urine=ESmmaurine,
		ES_ammonia_umolL=ESammonia,
		ES_pH=ESph,
		ES_onset_age=ESonset,
		ES_prot_mut_level=ESprotmut,
		ES_rnaseq_mut_level=ESrnamut,
		ES_AdoCblMinus=ESadominus,
		ES_AdoCblPlus=ESadoplus,
		ES_OHCblMinus=ESminus,
		ES_ratio=ESrat, 
		ES_OHCblPlus=ESplus,
		ES_cCSS=ESccss,
		myname=names(ee2)[i]
		)
})


new_phens <- c("Case/control",
"Gender",
"Cell passage no.",
"Cell collection date",
"Cell freezing date",
"Consanguinity",
"Age at onset",
"PI activity -",
"PI activity +",
"PI ratio",
"PI activity - (rel. day control)",
"PI activity + (rel. day control)",
"mut-type",
"MMUT activity -",
"MMUT activity +",
"MMUT activity - (rel. day control)",
"MMUT activity + (rel. day control)",
"Zygosity",
"MMUT protein expression",
"MMUT transcript expression",
"Acidosis",
"Metabolic acidosis",
"Metabolic ketoacidosis",
"Ketosis",
"Hyperammonemia",
"Abnormal muscle tone",
"Muscular hypotonia",
"Muscular hypertonia",
"Hypothermia",
"Respiratory Abnormality",
"Dyspnea",
"Tachypnea",
"Hyperventilation",
"Reduced consciousness",
"Lethargy",
"Somnolence",
"Coma",
"Seizures",
"Generalized tonic clonic seizures",
"Gastrointestinal problem",
"Feeding difficulties",
"Failure to thrive (on acute onset)",
"Vomiting",
"Dehydration",
"Developmental delay",
"Behavioral abnormality",
"Irritability",
"Cnocurrent infection",
"MMA (urine)",
"MMA (plasma)",
"Ammonia [umol/L]",
"pH",
"Base excess",
"pCO2",
"Glucose (plasma) [mmol/L]",
"Ketones (urine)",
"Dialysis",
"Peritoneal dialysis",
"Insulin (treatment)",
"Diet (treatment)",
"Carnitine (treatment)",
"Cobalamin (treatment)",
"Bicarbonate (treatment)",
"Glucose IV (treatment)",
"Responsive to acute treatment",
"Responsive to cobalamin treatment",
"MMA (urine) after treatment",
"Cobalamin dose PO",
"Cobalamin dose IM",
"Carnitine dose",
"Antibiotic treatment",
"Protein restriction",
"Protein (natural) amount",
"Protein (synthetic) amount",
"Protein (total) amount",
"Meals per day",
"Calories per day",
"Tube feeding (day-time)",
"Tube feeding (night-time)",
"Tube feeding (overall)",
"Weight centile",
"Length centile",
"Head circumference",
"Speech delay",
"IQ",
"Neurological abnormalities (chronic)",
"Impaired kidney function",
"Uric acid (plasma)",
"Creatinine (plasma)",
"Urea",
"(e)GFR",
"Hematological abnormalities",
"Anemia",
"Neutropenia",
"Skin abnormalities",
"Hearing impairment",
"Osteoporosis",
"Failure to thrive (chronic)",
"Clinical severity score",
"Global developmental delay",
"Muscular hypotonia (chronic)",
"Basal ganglia abnormalities",
"Clinical severity score (neurological)",
"Failure to thrive or tube feeding",
"Clinical severity score (collapsed)")



ff2=do.call("rbind",ff)
ff3=ff2[1:dim(ff2)[1],]
pval_tbl <- ff3[,c(1:12, 25)]
pval_tbl$myname <- new_phens
pval_tbl_mel <- melt(pval_tbl, id.vars = c("myname"))
setnames(pval_tbl_mel, c("myname", "molecular_variable", "pval"))
ES_tbl <- ff3[, 13:25]
ES_tbl$myname <- new_phens
ES_tbl_mel <- melt(ES_tbl, id.vars = c("myname"))
setnames(ES_tbl_mel, c("myname", "molecular_variable", "ES"))

ESpval_mer <- data.table(pval_tbl_mel, ES_tbl_mel$ES)
setnames(ESpval_mer, c("myname", "molecular_variable", "pval", "ES"))
ESpval_mer[, pval_adj := p.adjust(pval, method = "fdr")]
ESpval_mer$molecular_variable <- factor(ESpval_mer$molecular_variable, levels = c("cCSS", "OHCblPlus", "OHCblMinus", "ratio", "AdoCblPlus", "AdoCblMinus", "rnaseq_mut_level", "prot_mut_level", "onset_age","MMA_urine", "pH", "ammonia_umolL"))
ESpval_mer2 <- copy(ESpval_mer)

levels(ESpval_mer2$molecular_variable) <- c("Clinical score", "PI activity +", "PI activity -", "PI ratio", "MMUT activity +", "MMUT activity -", "MMUT transcript", "MMUT protein", "Age at onset", "MMA (urine)", "pH", "Ammonia")

ESpval_mer2[, pval_sig := ifelse(pval_adj>0.05, 0, 1)]
ESpval_sigs <- data.table(with(ESpval_mer2, table(molecular_variable, pval_sig)))
ESpval_sigs1 <- ESpval_sigs[pval_sig == 1, ]
ESpval_sigs1$molecular_variable <- factor(ESpval_sigs1$molecular_variable, levels = c("Clinical score", "PI activity +", "PI activity -", "PI ratio", "MMUT activity +", "MMUT activity -", "MMUT transcript", "MMUT protein", "Age at onset", "MMA (urine)", "pH", "Ammonia"))



# ff33 <- melt(ff3, id.vars = c("myname"))
# setnames(ff33, c("myname", "molecular_variable", "pval"))
# ff33[, pval_adj := p.adjust(pval, method = "fdr")]
# ff33$molecular_variable <- factor(ff33$molecular_variable, levels = c("cCSS", "OHCblPlus", "OHCblMinus", "ratio", "AdoCblPlus", "AdoCblMinus", "rnaseq_mut_level", "prot_mut_level", "onset_age","MMA_urine", "pH", "ammonia_umolL"))

# plot p value distribution of clinical and biochemical variables
ggplot(data = ESpval_mer2, aes(x = pval_adj)) +
	geom_histogram(bins = 20) +
	facet_wrap(molecular_variable~.) +
	xlab("p value (adjusted)") +
	ylab("Count") +
	scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9)) +
	ggtitle("Distribution of pvals for all clinical variables") +
	theme_pubr() +
	guides(fill = FALSE)

ggsave(paste(fig_path,"pvals_for_clinical_variables_ALL.png", sep = ""),width = 6, height = 6)
dev.off()

summary(ESpval_mer2$ES)


#≤

multi_reg_clinvars <- 
ggplot(ESpval_mer2[ES > -1.2 & ES < 1.2, ], aes(x = pval_adj, y = ES, color = ES)) +
	geom_point(alpha = 0.6) +
	geom_density(aes(x = pval_adj), inherit.aes = FALSE) +
	ylim(-1.5, 2) +
	geom_vline(xintercept = 0.05, linetype = "dashed", size = 0.5, alpha = 0.6) +
  	geom_text(data = ESpval_sigs1, aes(label = paste("p<0.05: N=",N, sep = "")), x = Inf, y = Inf, hjust=1, vjust=2, inherit.aes = FALSE, size = 3) +
	facet_wrap(molecular_variable~.) +
	scale_color_viridis(name = "ES", option = "inferno") +
	ggtitle("Multiple regression analysis") +
	xlab("p value (adjusted)") +
	ylab("Effect size (ES)") +
	scale_x_continuous(breaks = c(0, 0.4, 0.8)) +
	theme_pubr() +
	theme(legend.position = "right", plot.title = element_text(size = 12), strip.background = element_rect(color = "white"))

ggsave(paste(fig_path,"pvals_for_clinical_variables_ALL_withES.png", sep = ""),multi_reg_clinvars,width = 8, height = 6)
dev.off()
ggsave(paste(fig_path_pdf,"Fig2/","pvals_for_clinical_variables_ALL_withES.pdf", sep = ""),multi_reg_clinvars, device = "pdf", width = 7, height = 5)
dev.off()



CS_plt <- 
ggplot(data = ESpval_mer2[molecular_variable == "Clinical score", ], aes(x = pval_adj)) +
	geom_histogram(bins = 20) +
	xlab("p value (adjusted)") +
	ylab("Count") +
	ylim(0, 50) +
	annotate(geom = "text", x = 0.5, y = 25, label = paste("Sign. p values: N=", ESpval_sigs[molecular_variable == "Clinical score" & pval_sig == 1, ]$N, sep = "")) +
	scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9)) +
	ggtitle("Clinical score") +
	theme_pubr() +
	theme(plot.title = element_text(size = 12))

PI_plt <- 
ggplot(data = ESpval_mer2[molecular_variable == "PI activity +", ], aes(x = pval_adj)) +
	geom_histogram(bins = 20) +
	xlab("p value (adjusted)") +
	ylab("Count") +
	ylim(0, 50) +
	annotate(geom = "text", x = 0.5, y = 25, label = paste("Sign. p values: N=", ESpval_sigs[molecular_variable == "PI activity +" & pval_sig == 1, ]$N, sep = "")) +
	scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9)) +
	ggtitle("PI activity +") +
	theme_pubr() +
	theme(plot.title = element_text(size = 12))


ggsave(paste(fig_path,"pvals_for_ClinicalScoreAndPIact.png", sep = ""), ggarrange(CS_plt, PI_plt), width = 5, height = 3)
dev.off()


# only top3 (or 3 selected) variables
ggplot(data = ESpval_mer[molecular_variable == "cCSS" | molecular_variable == "OHCblPlus" | molecular_variable == "ratio", ], aes(x = pval_adj, fill = molecular_variable)) +
	geom_histogram(bins = 20) +
	facet_grid(molecular_variable~.) +
	ggtitle("Distribution of pvals for all clinical variables") +
	theme_test() +
	scale_fill_aaas() +
	guides(fill = FALSE)

ggsave(paste(fig_path,"pvals_for_clinical_variables_top3.png", sep = ""),width = 6,height = 6)


# plot the effect size of the clinial/biochemical parameters against their significance pval_adj - volcano plot

class(ESpval_mer$molecular_variable)

ESpval_mer_mod <- ESpval_mer[molecular_variable == "OHCblPlus", ]
ESpval_mer_mod$myname

ESpval_mer_mod <- ESpval_mer_mod[myname != "PI activity +" & myname != "PI activity + (rel. day control)" & myname != "PI activity -" & myname != "PI activity - (rel. day control)" & myname != "MMUT activity +" & myname != "MMUT activity -" & myname != "MMUT activity + (rel. day control)" & myname != "MMUT activity - (rel. day control)", ]
setorder(ESpval_mer_mod, pval)
ESpval_mer_mod

marked_vars <- c("Clinical severity score", "mut-type", "Acidosis", "Base excess", "MMA (plasma)", "pH", "Hyperventilation", "Impaired kidney function", "Protein restriction", "pCO2", "(e)GFR", "Age at onset")
marked_tbl <- ESpval_mer_mod %>% filter(myname %in% marked_vars)


clinical_vars_volcano_plt <- 
ggplot(ESpval_mer_mod, aes(x = ES, y = -log10(pval_adj), label = myname)) +
  geom_point(alpha = 0.4) +
  geom_point(data = marked_tbl, color = "red", size = 2) +
  geom_text_repel(data = marked_tbl, nudge_y = 0.4) +
  xlab("Effect size") +
  ylab("-log10(adjusted p value)") +
  xlim(-1,1)+
  guides(color = FALSE) +
  ggtitle("Correlation of PI activity to phenotypic variables") +
  theme_pubr() +
  theme(plot.title = element_text(size = 12))

ggsave(paste(fig_path,"pvalsANDeffectsizes_for_clinical_variables_OHCblPlus.png", sep = ""),clinical_vars_volcano_plt, width = 4, height = 3)
dev.off()
ggsave(paste(fig_path_pdf,"Fig2/","pvalsANDeffectsizes_for_clinical_variables_OHCblPlus.pdf", sep = ""),clinical_vars_volcano_plt, device = "pdf", width = 4, height = 3)
dev.off()

ESpval_mer_mod[myname == "MMA_urine"]
ESpval_mer_mod[myname == "MMA_plasma"]


# correlation matrix of clinical and biochemical variables
ee2x <- ee2
colnames(ee2x)
colnames(ee2x) <- new_phens
sapply(ee2x, class)
ee2x <- apply(ee2x, 2, function(x) as.numeric(as.character(x)))
apply(ee2x, 2, class)
ee2x <- data.frame(ee2x)

# exclude columns where sd is 0 (exclude NAs)
ncol(ee2x)
colnames(ee2x)
sapply(ee2x, sd, na.rm = TRUE) != 0
ee2xx <- Filter(function(x) sd(x, na.rm = TRUE) !=0, ee2x)
ncol(ee2xx)
colnames(ee2xx) == colnames(ee2x)

# generate correlation matrix
ee2xx_cormat <- cor(ee2xx, use = "pairwise.complete.obs")
head(ee2xx_cormat)
colnames(ee2xx_cormat) <- new_phens
rownames(ee2xx_cormat) <- new_phens



# plot correlation matrix
clinical_cor_htmp <-
pheatmap(ee2xx_cormat, color = inferno(100), border_color = "NA")

clinical_cor_htmp_gg <- as.ggplot(clinical_cor_htmp) + 
	ggtitle("Correlation matrix of clinical variables") +
	theme(plot.title = element_text(size = 22))



ggsave(paste(fig_path,"phenotype_variables_corrMatrix.png", sep = ""), clinical_cor_htmp_gg, width = 15, height = 15, bg = "white")
ggsave(paste(fig_path_pdf,"Fig2/","phenotype_variables_corrMatrix.pdf", sep = ""), clinical_cor_htmp_gg, device = "pdf", width = 15, height = 15, bg = "white")
dev.off()


# save_pheatmap_png <- function(x, filename, width = 1200, height = 1000, res = 150) {
#   png(filename, width = width, height = height, res = res)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }
 
# save_pheatmap_png(clinical_cor_htmp, paste(fig_path, "phenotype_variables_corrMatrix.png", sep = ""), width = 2100, height = 2100, res = 150)


#############################################################
# IQ variable analysis
tbl$IQ
summary(tbl$IQ)
ggplot(tbl, aes(x = IQ, fill = proposed_gene_wgs)) +
	geom_histogram() +
	theme_bw()

ggsave(paste(fig_path,"IQ_distribution.png", sep = ""),width = 6, height = 2)

ggplot(tbl[IQ>20, ], aes(x = IQ, y = OHCblPlus)) +
	geom_point() +
	geom_smooth(method = "lm") +
	ggtitle("No relationship between IQ and OHCblPlus", subtitle = "Excluding two outliers at IQ == 20") +
	theme_bw()

ggsave(paste(fig_path,"IQ_OHCblPlus.png", sep = ""),width = 6, height = 3)


#############################################################
# linear regression: analysis of factor variables (above only numeric variables are analysed)

ee = fread("interimData/phenotypic_tables/meta_info_wcss.txt",header=TRUE)

factor_vars <- c("sample_origin_contry", "sample_origin_town", "patient_origin_contry", "proposed_gene_wgs", "wgs_mut1nt", "wgs_mut1aa", "wgs_mut1type", "wgs_mut2nt", "wgs_mut2aa", "wgs_mut2type", "wgs_bothtype")

ee_fac <- ee[, ..factor_vars]

# make wgs_bothtype variable homogenous
ee_fac[wgs_bothtype == "/", wgs_bothtype := NA]
ee_fac[wgs_bothtype == "", wgs_bothtype := NA]
ee_fac[wgs_bothtype == "splicing/missense", wgs_bothtype := "missense/splicing"]
ee_fac[wgs_bothtype == "truncating/missense", wgs_bothtype := "missense/truncating"]
ee_fac[wgs_bothtype == "splicing/truncating", wgs_bothtype := "truncating/splicing"]
unique(ee_fac$wgs_bothtype)

lapply(c(1:dim(ee_fac)[2]), function(i){
	ee_fac2 = as.factor(unlist(ee_fac[, i, with = FALSE]))
	myna = is.na(ee_fac2)
	xx_fac = ee_fac2[!myna]
	yy_fac = tbl[!myna, scale(log10(OHCblPlus))]
	summary(lm(yy_fac~xx_fac))
})





#############################################################
# plot to visualize how well OHCblPlus correlates with cCSS
plt_pathwayCorrModcCSS <- 
ggplot(data = tbl, aes(x = log(OHCblPlus), y = as.factor(cssModif), fill = as.factor(cssModif), color = as.factor(cssModif))) +
	geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
	scale_fill_aaas() +
	scale_color_aaas() +
	theme_test()

# ggplot(data = tbl, aes(x = log(OHCblPlus), fill = as.factor(cssModif), color = as.factor(cssModif))) +
# 	geom_density(alpha = 0.6) +
# 	scale_fill_aaas() +
# 	scale_color_aaas() +
# 	theme_test()

plt_pathwayCorrcCSS <- 
ggplot(data = tbl, aes(x = log(OHCblPlus), y = as.factor(cCSS), fill = as.factor(cCSS), color = as.factor(cCSS))) +
	geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
	scale_fill_aaas() +
	scale_color_aaas() +
	theme_test()

plt_pathwayCSS <-
ggarrange(plt_pathwayCorrModcCSS, plt_pathwayCorrcCSS, nrow = 2)

ggsave(paste(fig_path,"pathwayCorrcCSS.png", sep = ""), plt_pathwayCSS, width = 6, height = 6)


# same plot with histograms
plt_pathwayCorrModcCSS <- 
ggplot(data = tbl, aes(x = log(OHCblPlus), fill = as.factor(cssModif), color = as.factor(cssModif))) +
	geom_histogram() +
	scale_fill_aaas() +
	scale_color_aaas() +
	theme_test()

plt_pathwayCorrcCSS <- 
ggplot(data = tbl, aes(x = log(OHCblPlus), fill = as.factor(cCSS), color = as.factor(cCSS))) +
	geom_histogram() +
	scale_fill_aaas() +
	scale_color_aaas() +
	theme_test()

plt_pathwayCSS <-
ggarrange(plt_pathwayCorrModcCSS, plt_pathwayCorrcCSS, nrow = 2)

ggsave(paste(fig_path,"pathwayCorrcCSS_hist.png", sep = ""), plt_pathwayCSS, width = 4, height = 6)

##########################################
## p-values of clinical associations (cCSS, OHCblPlus, ratio only) David's method
# ##########################################
# #ee=fread("interimData/phenotypic_tables/ClinicalInformationTable.csv",header=TRUE,sep=",")
# ee0=fread("interimData/phenotypic_tables/meta_info_wcss.txt",header=TRUE)
# # matcher=match(tbl2[,mma_id],ee[,`Cell Line Number`])
# # ee2=ee[matcher,]

# # set all clinical symptoms (rated as present or absent) to 0, if not reported to be present
# pres_abs_vars <- c("acidosis", "metabolic_acidosis", "metabolic_ketoacidosis", "ketosis", "hyperammonemia", "abnormal_muscle_tone", "musc_hypotonia", "musc_hypertonia", "hypothermia", "fct_respiratory_abnormality", "dyspnea", "tachypnea", "hyperventilation", "reduced_consciousness", "lethargy", "drowsiness", "somnolence", "coma", "seizures", "general_tonic_clonic_seizure", "any_GI_problem", "feeding_problem", "failure_to_thrive", "vomiting", "dehydration", "any_delay", "motor_delay", "behavioral_abnormality", "irritability", "autistic_behaviour", "concurrent_infection", "dialysis", "hemodialysis", "peritoneal_dialysis", "insulin", "diet", "carnitine", "cobalamin", "bicarb", "glucose_IV", "responsive_to_acute_treatment", "cobalamin_responsive", "protein_restriction", "global_dev_delay_chr", "language_delay", "any_neurological_abnormalities_chronic", "any_neurological_abnormalities_acuteANDchronic", "impaireded_kidney_fct", "hemat_abnormality", "anemia", "neutropenia", "pancytopenia", "skin_abnormalities", "hearing_impairment", "osteoporosis", "failure_to_thrive_chronic")
# ee1 <- data.frame(ee0)
# ee1[, pres_abs_vars][is.na(ee1[pres_abs_vars])] <- 0 # replace the NA's with 0
# ee2 <- data.table(ee1)

# ff=lapply(c(1:dim(ee2)[2]),function(i){
# 	print(i)
# 	ee3=as.numeric(unlist(ee2[,i,with=FALSE]))### remove NA's and empty
# 	ee4=ee3
# 	myna=is.na(ee4)
# 	if(sum(!myna)<2 || sd(ee3,na.rm=TRUE)==0){
# 		return(NULL)
# 	}
# 	xx=scale(ee4[!myna])
# #	yy=tbl2[!myna,log10(OHCblPlus)]
# 	yy=tbl[!myna,scale(log10(ratio))]
# 	tt=summary(lm(yy~xx))
# 	pvalrat=tt$coef[2,4]
# 	zrat=tt$coef[2,3]
# 	sloperat=tt$coef[2,1]
# 	yy=tbl[!myna,scale(log10(OHCblPlus))]
# 	tt=summary(lm(yy~xx))
# 	pvalplus=tt$coef[2,4]
# 	zplus=tt$coef[2,3]
# 	slopeact=tt$coef[2,1]
# 	yy=tbl[!myna,scale(cCSS)]
# 	tt=summary(lm(yy~xx))
# 	pvalccss=tt$coef[2,4]
# 	zccss=tt$coef[2,3]
# 	slopeccss=tt$coef[2,1]
# 	data.table(pval_ratio=pvalrat, z_ratio=zrat, slope_ratio = sloperat, 
# 		pval_OHCblPlus=pvalplus, z_OHCblPlus=zplus, slope_OHCblPlus = slopeact, 
# 		pval_cCSS=pvalccss, z_cCSS=zccss, slope_cCSS=slopeccss,
# 		myname=names(ee2)[i])
# })
# ff2=do.call("rbind",ff)
# ff3=ff2[3:dim(ff2)[1],]
# ff4=ff3[,list(myname, pval=pval_ratio, pval_adj=p.adjust(pval_ratio, method = "fdr"), slope=slope_ratio, molecular_variable="PI_ratio")]
# ff5=ff3[,list(myname, pval=pval_OHCblPlus, pval_adj=p.adjust(pval_OHCblPlus, method = "fdr"), slope=slope_OHCblPlus, molecular_variable="OHCblPlus")]
# ff6=ff3[,list(myname, pval=pval_cCSS, pval_adj=p.adjust(pval_cCSS, method = "fdr"), slope=slope_cCSS, molecular_variable="cCSS")]
# ff7=rbind(ff4,ff5,ff6)

# ggplot(data = ff7,aes(x = pval_adj, fill = molecular_variable)) +
# 	geom_histogram(bins = 20) +
# 	facet_grid(molecular_variable~.) +
# 	scale_fill_aaas() +
# 	ggtitle("Distribution of pvals for all clinical variables") +
# 	theme_test() +
# 	guides(fill = FALSE)

# ggsave(paste(fig_path,"pvals_for_clinical_variables_top3.png", sep = ""),width=4,height=4)

# mean(ff4[,pval<0.0125])
# ff4[pval<0.01,]
# ff5[pval<0.0125,]
# ff6[pval_adj<0.0125,]
# # mynames	pval 	molecular_variable
# # Age of Onset (d)	8.553449e-05	Ratio
# # Vit. B12-responsive? 0=no, 1=yes	4.090529e-05	Ratio
# # Age (months)	7.096545e-04	Ratio
# # no. of meals/day	1.853398e-07	Ratio
# # age (months)	2.472534e-04	Ratio
# # length (cm)	9.022156e-03	Ratio
# #	myname	pval	molecular_variable
# #	pH	0.0013084031	Activity\n(saturated)
# #	Vit. B12-responsive? 0=no, 1=yes	0.0003255492	Activity\n(saturated)
# fdr_act=0.01/mean(ff5[,pval<0.01])
# fdr_ratio=0.01/mean(ff4[,pval<0.01])

# #                                    myname         pval molecular_variable
# #  1:                      Age of Onset (d) 8.553449e-05              Ratio
# #  2:      Vit. B12-responsive? 0=no, 1=yes 4.090529e-05              Ratio
# #  3:                          Age (months) 7.096545e-04              Ratio
# #  4:                      no. of meals/day 1.853398e-07              Ratio
# #  5:                          age (months) 2.472534e-04              Ratio
# #  6:                           length (cm) 9.022156e-03              Ratio
#  7:                                    pH 1.308403e-03           Activity
#  8:      Vit. B12-responsive? 0=no, 1=yes 3.255492e-04           Activity
#  9:                            HP:0001263 8.669584e-03           Activity
# 10: impaired kidney function? 0=no, 1=yes 6.770888e-03           Activity





##################
# x-y diagrams for continous variables which are significantly associatd with OHCblPlus

# prep tables

ee0x = fread("interimData/phenotypic_tables/meta_info_wcss.txt",header=TRUE)
ee0x <- data.frame(ee0x)

# set all clinical symptoms (rated as present or absent) to 0, if not reported to be present
pres_abs_vars <- c("consanguinity","acidosis", "metabolic_acidosis", "metabolic_ketoacidosis","ketosis","hyperammonemia","abnormal_muscle_tone","musc_hypotonia","musc_hypertonia","hypothermia","fct_respiratory_abnormality" , "dyspnea","tachypnea","hyperventilation","reduced_consciousness","lethargy","drowsiness","somnolence","coma","seizures","general_tonic_clonic_seizure","any_GI_problem","feeding_problem","failure_to_thrive","vomiting","dehydration","any_delay","motor_delay","behavioral_abnormality","irritability","autistic_behaviour","concurrent_infection","urine_ketones","dialysis", "hemodialysis","peritoneal_dialysis","insulin","diet","carnitine","cobalamin","bicarb","glucose_IV", "responsive_to_acute_treatment","cobalamin_responsive","antibiotic_treatment","protein_restriction","tube_feeding_day","tube_feeding_night","tube_feeding_overall", "growth_hormone","weight_centile_quant","length_centile_quant","head_circumfernce_quant","language_delay", "any_neurological_abnormalities_chronic","cardiomyopathy","impaired_kidney_fct", "hemat_abnormality","anemia","neutropenia","pancytopenia","skin_abnormalities","hearing_impairment","osteoporosis","failure_to_thrive_chronic", "global_dev_delay_chr","hypotonia_chr","basal_ganglia_abnormality_chr","failure_to_thrive_or_tube_feeding")

ee0x[, pres_abs_vars][is.na(ee0x[pres_abs_vars])] <- 0 # replace the NA's with 0




tbl_cont <- ee0x %>% select(c("proposed_gene_wgs", "mma_id", "OHCblPlus", "onset_age", "MMA_urine", "MMA_plasma", "GFR", "pH", "base_excess", "pCO2", "cCSS"))
tbl_cont <- data.table(tbl_cont)
tbl_cont[, onset_age := log(tbl_cont$onset_age)]
tbl_cont[, MMA_urine := log(tbl_cont$MMA_urine)]
tbl_cont <- melt(tbl_cont, id.vars = c("proposed_gene_wgs", "mma_id", "OHCblPlus"))
tbl_cont[, sample_no := as.numeric(gsub("MMA", "", mma_id))]
tbl_cont[, type := ifelse(sample_no<151, "MMUT def.", ifelse(sample_no>210, "unaffected", "unknown"))]
tbl_cont[, type1 := ifelse(sample_no<151, "MMUT def.", "control")]
tbl_cont$type1 <- factor(tbl_cont$type1, levels = c("control", "MMUT def."))
levels(tbl_cont$variable) <- c("Age at onset", "MMA (urine)", "MMA (plasma)", "(e)GFR", "pH", "Base excess", "pCO2", "CSS")


PIwithContvars <- 
ggplot(tbl_cont[proposed_gene_wgs != "", ], aes(x = value, y = OHCblPlus)) +
	geom_point(alpha = 1, size = 1, aes(color = type1)) +
	geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
	stat_cor(label.y.npc = 0.1, size = 3, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), cor.coef.name = c("rho")) +
	scale_y_log10() +
  	annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
	facet_wrap(~variable, scales = "free_x", nrow = 2) +
	ggtitle("Association of clinical variables (continuous) with PI activity +", subtitle = "Samples with MMA gene defect (MMUT, MMAA, MMAB, ACSF3, SUCLA2, SUCLG1, TCN2)") +
  	ylab("PI activity\n[pmol/mg protein/16 h]") +
	scale_color_aaas() +
	theme_pubr() +
	theme(legend.position = "bottom", legend.title = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size = 12), plot.subtitle = element_blank(), strip.background = element_rect("white"))

ggsave(paste(fig_path,"clinical_variables_cont.png", sep = ""), PIwithContvars, device = png(), width = 9,height = 5)
dev.off() 
ggsave(paste(fig_path_pdf, "SuppFig3/","clinical_variables_cont.pdf", sep = ""), device = "pdf", PIwithContvars, width = 7.5,height = 5)
dev.off() 


##################
# correlation of PI activity with Clinical score

ggplot(tbl_cont[proposed_gene_wgs != "" & variable == "CSS", ], aes(x = value, y = OHCblPlus)) +
	geom_point(alpha = 1, size = 1, aes(color = type1)) +
	geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
	stat_cor(label.y.npc = 0.7, label.x.npc = 0.5, size = 3, cor.coef.name = c("rho"), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
	scale_y_log10() +
  	annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
	xlab("Clinical score") +
  	ylab("PI activity\n[pmol/mg protein/16 h]") +
	scale_color_aaas() +
	theme_pubr() +
	theme(legend.position = c(0.8, 0.9), legend.title = element_blank(), plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10))

ggsave(paste(fig_path,"PIvsClinicalScore.png", sep = ""), width = 3.5, height = 2.5)
dev.off() 



##################
# boxplots for 1-0 variables which are significantly associatd with OHCblPlus

tbl_disc <- ee0x %>% select(c("proposed_gene_wgs", "mma_id", "OHCblPlus", "acidosis", "hyperventilation", "any_delay", "diet", "bicarb", "glucose_IV", "cobalamin_responsive", "protein_restriction", "tube_feeding_overall", "impaired_kidney_fct"))
tbl_disc <- data.table(tbl_disc)
tbl_disc <- melt(tbl_disc, id.vars = c("proposed_gene_wgs", "mma_id", "OHCblPlus"))
tbl_disc[, sample_no := as.numeric(gsub("MMA", "", mma_id))]
tbl_disc[, type := ifelse(sample_no<151, "MMUT def.", ifelse(sample_no>210, "unaffected", "unknown"))]
tbl_disc[, type1 := ifelse(sample_no<151, "MMUT def.", "control")]
tbl_disc$type <- factor(tbl_disc$type, levels = c("MMUT def.", "unaffected", "unknown"))
tbl_disc$type1 <- factor(tbl_disc$type1, levels = c("control", "MMUT def."))
tbl_disc[, value2 := ifelse(value == 1, "Present", "Absent")]
levels(tbl_disc$variable) <- c("Acidosis", "Hyperventilation", "Developmental delay", "Diet (treatment)", "Bicarbonate (treatment)", "Glucose IV (treatment)", "Responsive to cobalamin", "Protein restriction", "Tube feeding", "Impaired kidney function")



# ggplot(tbl_disc, aes(x = as.factor(value), y = log(OHCblPlus))) +
# 	#geom_boxplot(alpha = 0.6, aes(color = as.factor(value))) +
# 	geom_violin(alpha = 0.6) +
# 	geom_jitter(alpha = 0.6, size = 1, width = 0.2, aes(color = type)) +
# 	facet_wrap(~variable) +
# 	theme_test() +
# 	scale_color_aaas() +
# 	ggtitle("Association of clinical symptoms (discrete) with OHCblPlus")

# ggsave(paste(fig_path,"clinical_variables_disc.png", sep = ""),width=8,height=6)

max(tbl_disc[proposed_gene_wgs!="",]$OHCblPlus)

PIwithDiscvars <- 
ggplot(tbl_disc[proposed_gene_wgs!="",], aes(x = as.factor(value2), y = OHCblPlus)) +
	geom_boxplot(alpha = 1, outlier.shape = NA, aes(color = as.factor(value))) +
	geom_jitter(width = 0.16, size = 1, alpha = 1, pch = 21, fill = "white", aes(color = as.factor(value))) +
	stat_compare_means(comparisons = list(c("Present","Absent")), method = "t.test", size = 3) +
	scale_y_log10(limits = c(50,50000)) +
  	annotation_logticks(sides = "l", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
	facet_wrap(~variable) +
	scale_color_jama() +
	ggtitle("Association of clinical variables (discrete) with PI activity +", subtitle = "Samples with MMA gene defect (MMUT, MMAA, MMAB, ACSF3, SUCLA2, SUCLG1, TCN2)") +
	guides(color = FALSE) +
  	ylab("PI activity\n[pmol/mg protein/16 h]") +
	xlab("Presence of clinical variable") +
	theme_pubr() +
	theme(plot.title = element_text(size = 12), plot.subtitle = element_blank(), strip.background = element_rect("white"))


ggsave(paste(fig_path,"clinical_variables_disc_box.png", sep = ""), PIwithDiscvars, device = png(), width = 9, height = 7)
dev.off()
ggsave(paste(fig_path_pdf, "SuppFig3/","clinical_variables_disc_box.pdf", sep = ""), PIwithDiscvars, device = "pdf", width = 8, height = 7.5)
dev.off()


##################
# plots for some variables to correlate with cCSS
# cCSS (chronic clinical score) consists of: onset_score, neuro_abnormality, kidney impairment, haematological abnormality, failure to thrive

#≥

tbl_css0 <- ee0x %>% select(c("proposed_gene_wgs", "mma_id", "cCSS", "acidosis", "hyperammonemia", "musc_hypotonia", "fct_respiratory_abnormality", "dehydration", "dehydration", "dialysis", "diet", "carnitine", "bicarb", "glucose_IV", "protein_restriction"))
tbl_css0 <- data.table(tbl_css0)
tbl_css <- melt(tbl_css0, id.vars = c("proposed_gene_wgs", "mma_id", "cCSS"))
tbl_css <- data.table(tbl_css)
tbl_css[, cCSS_mod := ifelse(cCSS>2,2,cCSS)]
tbl_css[, cCSS_mod2 := cCSS_mod]
tbl_css$cCSS_mod2 <- factor(tbl_css$cCSS_mod2, levels = c("0", "1", "2"))
levels(tbl_css$cCSS_mod2) <- c("0", "1", ">=2")
tbl_css[, value2 := ifelse(value == 1, "Present", "Absent")]
levels(tbl_css$variable) <- c("Acidosis", "Hyperammonemia", "Muscular hypotonia", "Respiratory abnormality", "Dehydration", "Dialysis (treatment)", "Diet (treatment)", "Carnitine (treatment)", "Bicarbonate (treatment)", "Glucose IV (treatment)", "Protein restriction")


ggplot(tbl_css[proposed_gene_wgs!="",], aes(x = as.factor(value), y = cCSS_mod)) +
	geom_boxplot(alpha = 0.6, outlier.shape = NA, aes(color = as.factor(value))) +
	geom_jitter(width = 0.16, size = 1, alpha = 0.3, aes(color = as.factor(value))) +
	facet_wrap(~variable) +
	theme_test() +
	scale_color_aaas() +
	ggtitle("Association of clinical variables (discrete) with clinical severity score", subtitle = "Samples with MMA gene defect (MMUT, ACSF3, SUCLA2, SUCLG1, TCN2)") +
	guides(color = FALSE) +
	theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10))

ggsave(paste(fig_path,"clinical_variables_disc_cCSS_box.png", sep = ""),width=6.5,height=6)


CSSwithDiscvars <- 
ggplot(tbl_css[proposed_gene_wgs!="",], aes(x = cCSS_mod2, fill = as.factor(value2))) +
	geom_bar(alpha = 0.8, position = "fill") +
	facet_wrap(~variable) +
	theme_pubr() +
	scale_fill_jama() +
	labs(fill = "Presence of clinical variable") +
	ggtitle("Association of clinical variables (discrete) with clinical severity score (CSS)", subtitle = "Samples with MMA gene defect (MMUT, MMAA, MMAB, ACSF3, SUCLA2, SUCLG1, TCN2)") +
	xlab("Clinical score") +
	ylab("Relative abundance of clinical variable") +
	theme(legend.position = "bottom", plot.title = element_text(size = 12), plot.subtitle = element_blank(), strip.background = element_rect("white"))

ggsave(paste(fig_path,"clinical_variables_disc_cCSS.png", sep = ""), CSSwithDiscvars, width = 8, height = 7)
dev.off()
ggsave(paste(fig_path_pdf, "SuppFig3/", "clinical_variables_disc_cCSS.pdf", sep = ""), CSSwithDiscvars, device = "pdf", width = 7.5, height = 7.5)
dev.off()




tbl_cssAoO <- ee0x %>% select(c("proposed_gene_wgs", "mma_id", "cCSS", "age_of_onset_d"))
tbl_cssAoO <- data.table(tbl_cssAoO)
tbl_cssAoO[, cCSS_mod := ifelse(cCSS>2,2,cCSS)]
tbl_cssAoO[, cCSS_mod2 := cCSS_mod]
tbl_cssAoO$cCSS_mod2 <- factor(tbl_cssAoO$cCSS_mod2, levels = c("0", "1", "2"))
levels(tbl_cssAoO$cCSS_mod2) <- c("0", "1", ">=2")

tbl_cssAoO <- data.table(tbl_cssAoO)
tbl_cssAoO[, onset := ifelse(age_of_onset_d>28, "Late (>28 days)", "Early (<=28 days)")]



ccss_onset_plt0 <-
ggplot(tbl_cssAoO, aes(x = log(age_of_onset_d), fill = as.factor(cCSS_mod), color = as.factor(cCSS_mod))) +
	geom_density(alpha = 0.6) +
	scale_fill_aaas() +
	scale_colour_aaas() +
	ggtitle("Late onset has low cCSS", subtitle = "all samples") +
	theme_test()

ccss_onset_plt1 <-
ggplot(tbl_cssAoO, aes(x = age_of_onset_d, y = cCSS_mod2, fill = cCSS_mod2, color = cCSS_mod2)) +
	geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7, color = "black") +
	scale_x_log10(breaks = c(1, 28, 100, 10000), labels = c(1, 28, 100, 10000)) +
  	annotation_logticks(sides = "b", short = unit(0.5,"mm"), mid = unit(0.5,"mm"), long = unit(1,"mm")) +
  	geom_vline(xintercept = 28, linetype = "dashed", alpha = 0.6) +
	scale_fill_rickandmorty() +
	coord_cartesian(clip = "off") +
	ylab("Clinical severity score") +
	xlab("Age at onset [days]") + 
	theme_pubr() +
	theme(legend.position = "none",  plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10))

ccss_onset_plt2 <-
ggplot(tbl_cssAoO[onset != "", ], aes(x = cCSS_mod2, fill = onset)) +
	geom_bar(alpha = 0.8, position = "fill") +
	theme_pubr() +
	coord_flip() +
	scale_fill_uchicago() +
	ylab("Proportion") +
	xlab("Clinical severity score") +
	labs(fill = "Age at onset") +
	theme(legend.position = "bottom",  plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10))

ccss_onset_plt <- 
ggarrange(ccss_onset_plt1, ccss_onset_plt2, ncol = 1)

ggsave(paste(fig_path,"onset_cCSS.png", sep = ""),ccss_onset_plt, width=5,height=5)
dev.off()
ggsave(paste(fig_path_pdf,"SuppFig3/","onset_cCSS.pdf", sep = ""),ccss_onset_plt, device = "pdf", width = 5, height = 5)
dev.off()



##########################################
##########################################
##########################################
##########################################
##########################################