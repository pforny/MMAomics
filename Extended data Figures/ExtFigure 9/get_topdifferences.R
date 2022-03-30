require(Matrix)
require(GSA)
require(data.table)
require(ggplot2)
require(ggsci)
require(ggpubr)



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


#source("../tools/csrproject-master/Code/fast_lmm.R")

source("csrproject/Code/normalizationFunctions.R")
source("csrproject/Code/fast_lmm.R")
source("csrproject/Code/fast_lmm_group.R")


fig_path <- c("Figs/DEA/TopDifferences/")
fig_path_pdf <- c("Figs/v12/pdf/SuppFig9_metabolomics/")
mypal1 <- pal_aaas("default", alpha = 1)(9)
mypal2 <- pal_npg("nrc", alpha = 1)(9)
mypal3 <- c(mypal1[2], mypal1[4], mypal1[1])


#phenoFresh=fread("Data/mma_tbl.csv")
#PROT=TRUE
PROT=TRUE
CONDITIONAL=FALSE



pheno0=fread("interimData/phenotypic_tables/meta_info_combined28_8_19.txt")
#pheno0=pheno0[case==1,]
pheno0[is.na(OHCblPlus),OHCblPlus:=10000]
pheno0[,AdoCblPlus:=as.numeric(AdoCblPlus)]
pheno0[is.na(AdoCblPlus) & case==0,AdoCblPlus:=3000]
pheno0[is.na(AdoCblPlus) & case==1,AdoCblPlus:=50]

if(CONDITIONAL){
	pheno0=pheno0[is.element(as.numeric(sub("MMA","",mma_id)),1:150),]
	cond_str="_conditional"
}else{
	cond_str=""
}

if(PROT){
	ww=load("interimData/prot_mat.RData")
	exp_mat=prot_mat
	type_used_str="prot"
	}else{
	ww=load("interimData/rnafiltered_mat.RData")
	exp_mat=rnafiltered_mat
	type_used_str="rnaseq"
}

mmaids=intersect(pheno0[,mma_id],colnames(exp_mat))
exp_mat_used=exp_mat[,match(mmaids,colnames(exp_mat))]
pheno=pheno0[match(mmaids,mma_id),]
mapping_tbl=fread("interimData/protein_coding.bed", header=TRUE, sep="\t")
intersected_ensg=intersect(rownames(exp_mat_used),mapping_tbl[,ensembl])
exp_mat_final=exp_mat_used[match(intersected_ensg,rownames(exp_mat_used)),]
mapping_final=mapping_tbl[match(intersected_ensg,ensembl),]


pheno[,normedpw:=qnorm(rank(OHCblPlus)/(1+nrow(pheno)))]
myy=pheno[,normedpw]
exp_sc=scalerRound(exp_mat_final,7)



ee3=fread(paste("interimData/DEA/diff_exp_",type_used_str,"_pathwayact_all_annotout.txt",sep=""))
subtbl=ee3[isenzyme==TRUE & ismitog==TRUE & chirank<20 & gname!="MMUT",]

#ee=fast_lmm(my_y=myy,my_X=t(exp_sc))

get_res=function(ensembl_list,getFullDataOut=FALSE){
	myels=is.element(rownames(exp_sc),ensembl_list)
	myindex_l=list()
	myindex_l[[1]]=myels	
	ff=fast_lmm_group(my_y=myy,my_X=t(exp_sc),indexList=myindex_l)
	mymod=summary(lm(myy~t(exp_sc)[,myindex_l[[1]]]-1))
	RR=t(exp_sc)[,myindex_l[[1]]]%*%ff$all_betas_l[[1]]
	if(getFullDataOut){
		outl=list()
		outl[["X"]]=t(exp_sc)[,myindex_l[[1]]]
		outl[["beta"]]=ff$all_betas_l
		return(outl)
	}
	rownames(RR)=colnames(exp_sc)
	orderedNames=rownames(RR)[order(RR)]
	RR[is.element(rownames(RR),tail(orderedNames,10))]
	RR[is.element(rownames(RR),head(orderedNames,10))]
	mymatcher=match(orderedNames,pheno[,mma_id])
	pheno_m=pheno[mymatcher,]
	plot(pheno_m[,(OHCblPlus>3000)+rnorm(230)/10])
	tblName=data.table(mma_id=orderedNames,pw=pheno_m[,OHCblPlus],proposed_g=pheno_m[,proposed_gene_wgs],iscase=pheno_m[,case], mymut_category=pheno_m[,mut_category])
	tblName[,active_pw:=(pw>3000)+rnorm(230)/10]
	tblName[,ordered_samples:=c(1:230)]
	tblName[,type:="is mut"]
	tblName[,ismut:=proposed_g=="MMUT"]
	tblName[ismut==FALSE,type:="not mut"]
	tblName[iscase==0,type:="control"]
	return(tblName)
}

table_printing=function(tt){
	subtbl=head(tt[proposed_g=="MMUT",],10)
	ff=data.table(mma_id=subtbl[,mma_id],type="first 10 with mut")
	subtbl2=tail(tt[proposed_g!="MMUT",],10)
	ff2=data.table(mma_id=subtbl2[,mma_id],type="last 10 without mut")
	subtbl3=tail(tt[iscase!=1,],10)
	ff3=data.table(mma_id=subtbl3[,mma_id],type="last 10 controls")
	ff4=rbind(rbind(ff,ff2),ff3)
	return(ff4)
}

ensembl_list=subtbl[,ensembl]


tt=get_res(ensembl_list)
tt[,ismut:=proposed_g=="MMUT"]

ggplot(tt,aes(group=type,fill=type,x=ordered_samples))+geom_histogram(alpha=0.8)
ggsave(paste0(fig_path,"predictions_4top_prots.pdf"))
outtbl=table_printing(tt)

write.table(outtbl,file="interimData/outl_samples_4top_prots.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
ensembl_list2=subtbl[gname=="OGDH",ensembl]
tt2=get_res(ensembl_list2)
ggplot(tt2,aes(group=type,fill=type,x=ordered_samples))+geom_histogram(alpha=0.8)
ggsave(paste0(fig_path,"predictions_ogdh_prots.pdf"))
outtbl2=table_printing(tt2)
write.table(outtbl2,file="interimData/outl_samples_OGDH.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)



ensembl_list_ogdhgl=c("ENSG00000148672","ENSG00000105953") # GLUD1, OGDH
#subtbl[,ensembl]


tt0=get_res(ensembl_list_ogdhgl,getFullDataOut=TRUE)
tt=get_res(ensembl_list_ogdhgl,getFullDataOut=FALSE)

tt[,ismut:=proposed_g=="MMUT"]
tt[, type2 := "unaffected"]
tt[type == "is mut", type2 := "MMUT def."]
tt[type == "not mut", type2 := "unknown"]
tt$type2 <- factor(tt$type2, levels = c("MMUT def.", "unknown", "unaffected"))
tt[, type3 := "control"]
tt[type == "is mut", type3 := "MMUT def."]
tt[type == "not mut", type3 := "NA"]
tt$type3 <- factor(tt$type3)

model_ogdh_glud1 <- 
ggplot(tt,aes(group=type2,fill=type2,x=ordered_samples))+
	geom_density(alpha=0.6, color = NA)+
	labs(y = "Density", x="Sample rank for OGDH/GLUD1 model") +
	scale_fill_manual(values = mypal3)+
	theme_test() +
	theme(legend.title = element_blank(), legend.position = "right",  axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))

model_ogdh_glud1_2 <- 
ggplot(tt[type3 != "NA", ],aes(group=type3,fill=type3,x=ordered_samples))+
	geom_density(alpha=0.6, color = NA)+
	labs(y = "Density", x="Sample rank for OGDH/GLUD1 model") +
	scale_fill_aaas()+
	theme_pubr() +
	theme(legend.title = element_blank(), legend.position = "right",  axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))

ggsave(paste0(fig_path,"modelrank_ogdhg_glud1_prots.pdf"), model_ogdh_glud1)
ggsave(paste0(fig_path_pdf,"modelrank_ogdhg_glud1_prots.pdf"), model_ogdh_glud1, width = 4, height = 2.5)
ggsave(paste0(fig_path_pdf,"modelrank_ogdhg_glud1_prots_2.pdf"), model_ogdh_glud1_2)

##### 
##### add
##### 
# updated_gnames_match=match(tt[,mma_id],phenoFresh[,`Cell Line Number`])
# newgenes=phenoFresh[updated_gnames_match,Gene]
# tt[,proposed_g_old:=proposed_g]
# tt[,proposed_g_updated:=newgenes]
# tt[,proposed_g:=NULL]

tt[,proposed_g_updated:=proposed_g]
#####
ggplot(tt,aes(group=type,fill=type,x=ordered_samples))+geom_histogram(alpha=0.8)

ggsave(paste0(fig_path,"predictions_ogdhg_glud1_prots.pdf"))
outtbl=table_printing(tt)
tt2=tt[1:230,]
tt2[,proposed_g_old:=NULL]
tt2[,ismut:=NULL]
write.table(tt2,file="interimData/all_samples_ogdhg_glud1_prots.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
write.table(outtbl,file="interimData/outl_samples__ogdhg_glud1_prots.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

tt3=tt2[is.element(proposed_g_updated, c("","ACSF3","MMUT")),]
ggplot(tt3,aes(group=proposed_g_updated,fill=proposed_g_updated,x=ordered_samples))+geom_density(alpha=0.7)+xlab("sample rank for (OGDH/GLUD1 model)") + theme_bw() + scale_fill_manual(values=c(myblue,myorange, mygrey))



ensembl_list_fullenz=c("ENSG00000105953","ENSG00000148672","ENSG00000119689","ENSG00000091140","ENSG00000125166")

tt0_enz=get_res(ensembl_list_fullenz,getFullDataOut=TRUE)
tt_enz=get_res(ensembl_list_fullenz,getFullDataOut=FALSE)

tt_enz[,ismut:=proposed_g=="MUT"]

ggplot(tt_enz,aes(group=type,fill=type,x=ordered_samples))+geom_histogram(alpha=0.8)
ggsave(paste0(fig_path,"predictions_fullenz_prots.pdf"))
outtbl_enz=table_printing(tt_enz)

write.table(outtbl_enz,file="interimData/outl_samples_fullenz_prots.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


# aa=tt_enz[,list(mma_id,ordered_samples)]
# aa2=tt[,list(mma_id,ordered_samples)]
# setkey(aa,mma_id)
# setkey(aa2,mma_id)
# ee=merge(aa,aa2)
# mytbl=fread("outl_samples_4top_prots.txt")
# aa3=mytbl[,list(mma_id,ordered_samples)]
# setkey(aa3,mma_id)
# ee2=merge(ee,aa3)


