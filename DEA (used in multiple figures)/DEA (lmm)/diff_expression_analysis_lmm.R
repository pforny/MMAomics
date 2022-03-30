require(Matrix)
require(data.table)
require(ggplot2)
require(locfdr)

source("csrproject/Code/normalizationFunctions.R")
source("csrproject/Code/fast_lmm.R")

system("mkdir interimData/DEAlmm")

get_fast_lmm=function(myy, ex_sc){
	ee=fast_lmm(my_y=myy,my_X=t(exp_sc))
	ee_tbl=data.table(ee)
	ee_tbl[,pval:=1-pchisq(chi_sq,1)]
	ee_tbl[,chirank:=rank(-abs(chi_sq))]
	rownames(ee_tbl)=rownames(exp_sc)
	return(ee_tbl)
}

get_lm=function(myy, ex_sc){
	i=1
	EE=lapply(c(1:nrow(exp_sc)),function(i){
		ff=summary(lm(myy~exp_sc[i,]))
		coefs=ff$coef[2,]
		tbl=data.table(betas=coefs[1],chi_sq=coefs[1]^2,deltas=NA,all_sigmasqe=NA,pval=coefs[4])
	})
	full_tbl=do.call("rbind",EE)	
	full_tbl[,chirank:=rank(-abs(chi_sq))]
	rownames(full_tbl)=rownames(exp_sc)
	return(full_tbl)
}



DATASET="prot"
CONDITIONAL=FALSE
USECSS=FALSE
USELMM=TRUE
#usecss_vec=c(TRUE,FALSE)
#usecss_vec=c(FALSE)
usecss_vec=c(TRUE,FALSE)
dataset_vec=c("prot","rna","average")
conditional_vec=c(FALSE,TRUE)

#uselmm_vec=c(TRUE,FALSE)
uselmm_vec=c(FALSE,TRUE)

for(DATASET in dataset_vec){
	for(USECSS in usecss_vec){	
		for(CONDITIONAL in conditional_vec){
			for(USELMM in uselmm_vec){
				if(DATASET=="prot"){
					ww=load("interimData/prot_mat.RData")
					exp_mat=prot_mat
					type_used_str="prot"
					locfdrBre=50
					}else if(DATASET=="rna"){
					ww=load("interimData/rnafiltered_mat.RData")
					exp_mat=rnafiltered_mat
					type_used_str="rnaseq"
					locfdrBre=120
				}else{
					ww=load("interimData/rnafiltered_mat.RData")
					ww=load("interimData/prot_mat.RData")
					myintersMMA=intersect(colnames(prot_mat),colnames(rnafiltered_mat))
					protMatSorted=prot_mat[,match(myintersMMA,colnames(prot_mat))]
					rnaMatSorted=rnafiltered_mat[,match(myintersMMA,colnames(rnafiltered_mat))]
					myinters=intersect(rownames(prot_mat),rownames(rnafiltered_mat))
					protMatSorted=protMatSorted[match(myinters,rownames(prot_mat)),]
					rnaMatSorted=rnaMatSorted[match(myinters,rownames(rnafiltered_mat)),]
					protMatSorted=scalerRound(protMatSorted,7)
					rnaMatSorted=scalerRound(rnaMatSorted,7)
					exp_mat=protMatSorted+rnaMatSorted
					type_used_str="avg"
					locfdrBre=50
				}

				pheno0=fread("interimData/phenotypic_tables/meta_info_wcss.txt")

				if(CONDITIONAL){
					pheno0=pheno0[is.element(as.numeric(sub("MMA","",mma_id)),1:150),]
					cond_str="_conditional"
				}else{
					cond_str="_all"
				}
				#pheno0=pheno0[case==1,]
				mmaids=intersect(pheno0[,mma_id],colnames(exp_mat))
				exp_mat_used=exp_mat[,match(mmaids,colnames(exp_mat))]
				pheno=pheno0[match(mmaids,mma_id),]

		
				if(USECSS){
					pheno_str="_wcss"
					pheno[,normedpw:=scale(cssModif)]
				}else{
					pheno0[is.na(OHCblPlus),OHCblPlus:=10000]
					pheno0[,AdoCblPlus:=as.numeric(AdoCblPlus)]
					pheno0[is.na(AdoCblPlus) & case==0,AdoCblPlus:=3000]
					pheno0[is.na(AdoCblPlus) & case==1,AdoCblPlus:=50]
					pheno[,normedpw:=qnorm(rank(OHCblPlus)/(1+nrow(pheno)))]
					pheno_str="_pathwayact"
				}

				mapping_tbl=fread("interimData/protein_coding.bed", header=TRUE, sep="\t")
				intersected_ensg=intersect(rownames(exp_mat_used),mapping_tbl[,ensembl])
				exp_mat_final=exp_mat_used[match(intersected_ensg,rownames(exp_mat_used)),]
				mapping_final=mapping_tbl[match(intersected_ensg,ensembl),]

				myy=scale(pheno[,normedpw])
				exp_sc=scalerRound(exp_mat_final,7)
				if(USELMM){
					method_str=""
					out0_tbl=get_fast_lmm(myy,exp_sc)
				}else{
					method_str="nolmm_"
					out0_tbl=get_lm(myy,exp_sc)
				}
				val=sum(rownames(out0_tbl)!=mapping_final[,ensembl])
				if(val==0){
					out_tbl=cbind(mapping_final,out0_tbl)
				}

				out_tbl[,zsc:=-qnorm(pval/2)*sign(betas)]
				tt=locfdr(out_tbl[,zsc], bre=locfdrBre, df=7, pct=0, pct0=1/4, nulltype=1, type=0, plot=0)
				estM=tt$fp0[3,1]
				estS=tt$fp0[3,2]
				out_tbl[,fdr:=tt$fdr]
				out_tbl[,zscnorm:=(zsc-estM)/estS]
				out_tbl[,pvalnorm:=pnorm(-abs(zscnorm))*2]
				out_tbl[,chiranknorm:=rank(pvalnorm)]
				out_tbl[,ptheornorm:=chiranknorm/nrow(out_tbl)]
				out_str=paste("interimData/DEAlmm/diff_exp_data_",method_str,type_used_str,pheno_str,cond_str,".txt",sep="")
				print(out_str)
				print(out_tbl[fdr<0.2,])
				write.table(out_tbl,file=out_str,quote=F,row.names=F,col.names=TRUE,sep="\t")
			}
		}
	}
}
