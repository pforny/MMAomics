require(data.table)
require(ggplot2)

PROT=FALSE
CONDITIONAL=TRUE
USECSS=TRUE
USELMM=TRUE

system("mkdir Results")
system("mkdir Results/DEA")
system("mkdir Figs")
system("mkdir Figs/DEA")
system("mkdir Figs/DEA/lmm")
system("mkdir interimData/DEA")

usecss_vec=c(TRUE,FALSE)

conditional_vec=c(FALSE,TRUE)

dataset_vec=c("prot","rna","average")

#uselmm_vec=c(TRUE,FALSE)
uselmm_vec=c(FALSE,TRUE)

for(USECSS in usecss_vec){
	for(DATASET in dataset_vec){
		for(CONDITIONAL in conditional_vec){
			for(USELMM in uselmm_vec){
				if(USECSS){
					pheno_str="_wcss"
				}else{
					pheno_str="_pathwayact"
				}

				if(CONDITIONAL){
					cond_str="_conditional"
				}else{
					cond_str="_all"
				}
				if(DATASET=="prot"){
					type_used_str="prot"
				}else if(DATASET=="rna"){
					type_used_str="rnaseq"
					}else{
					type_used_str="avg"
				}
				if(USELMM){
					method_str=""
					}else{
					method_str="nolmm_"
				}

				infile_str=paste("interimData/DEAlmm/diff_exp_data_",method_str,type_used_str,pheno_str,cond_str,".txt",sep="")
				print(infile_str)
				
				outfile_str=paste("Figs/DEA/lmm/diff_exp_",method_str,type_used_str,pheno_str,cond_str,".pdf",sep="")
				outsubsetfile_str=paste("Figs/DEA/lmm/diff_exp_",method_str,type_used_str,pheno_str,cond_str,"_subset.pdf",sep="")
				outtcasubsetfile_str=paste("Figs/DEA/lmm/diff_exp_",method_str,type_used_str,pheno_str,cond_str,"_tcasubset.pdf",sep="")
				outsubsettbl_str=paste("Results/DEA/diff_exp_",method_str,type_used_str,pheno_str,cond_str,"_subset_topout.txt",sep="")
				outtcasubsettbl_str=paste("Results/DEA/diff_exp_",method_str,type_used_str,pheno_str,cond_str,"_tcasubset_topout.txt",sep="")
				outfull_str=paste("interimData/DEA/diff_exp_",method_str,type_used_str,pheno_str,cond_str,"_annotout.txt",sep="")
				plot_title=paste("diff_exp_",method_str,type_used_str,pheno_str, cond_str,sep="")
				plot_title_subset=paste("diff_exp_",method_str,type_used_str,pheno_str,cond_str,"_subset",sep="")
				plot_title_tcasubset=paste("diff_exp_",method_str,type_used_str,pheno_str,cond_str,"_tcasubset",sep="")




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

				res_tbl = fread(infile_str,header=TRUE)

				enzymes = fread("interimData/gene_sets/all_enzyme_ensg.tbl",header=FALSE)
				res_tbl[,isenzyme:=is.element(ensembl,enzymes[,V4])]


				pulld = fread("interimData/gene_sets/pulldown_mmm.txt",header=FALSE)
				myg = pulld[V2<0.05,V4]
				res_tbl[,ispulldown:=is.element(ensembl,myg)]

				# for patrick: You can add new subsetting stratgies here. please make sure that the original plots are kept. (i.e. make additional plots if you want.)

				mitog = fread("interimData/gene_sets/mitoglist.txt",header=FALSE)[,V1]
				
				tcaprot = fread("interimData/gene_sets/tca_proteins.txt", header=FALSE)[,V1]

				res_tbl[,ismitog:=is.element(ensembl,mitog)]
				res_tbl[,istca:=is.element(ensembl,tcaprot)]
				res_tbl[,inpulldown:="not in pulldown"]
				res_tbl[ispulldown==TRUE,inpulldown:="in pulldown"]
				res_tbl[gname=="MMUT",inpulldown:="is MMUT"]
				res_tbl[,p_theor:=rank(pval)/(dim(res_tbl)[1]+1)]
				res_tbl[,pnorm_theor:=rank(pvalnorm)/(dim(res_tbl)[1]+1)]

				res_subset_tbl=res_tbl[ismitog==TRUE & isenzyme==TRUE,]
				res_subset_tbl[,p_theor:=rank(pval)/(dim(res_subset_tbl)[1]+1)]
				
				res_tcasubset_tbl=res_tbl[istca==TRUE]
				res_tcasubset_tbl[,p_theor:=rank(pval)/(dim(res_tcasubset_tbl)[1]+1)]

				ggplot(data=res_tbl,aes(x=-log10(p_theor),y=-log10(pval),colour=inpulldown))+geom_point(size=2)+geom_abline(intercept=0,slope=1)+xlab("-log10(p-value) (theoretical)")+ylab("-log10(p-value) (empirical)")+geom_abline(intercept=1,slope=1,linetype=2)+scale_colour_manual(values=c(myblue,myorange,mygrey))+theme_bw()+ggtitle(plot_title)
				ggsave(outfile_str,width=4,height=3)

				ggplot(data=res_subset_tbl,aes(x=-log10(p_theor),y=-log10(pval),colour=inpulldown))+geom_point(size=2)+geom_abline(intercept=0,slope=1)+xlab("-log10(p-value) (theoretical)")+ylab("-log10(p-value) (empirical)")+geom_abline(intercept=1,slope=1,linetype=2)+scale_colour_manual(values=c(myblue,myorange,mygrey))+theme_bw()+ggtitle(plot_title_subset)
				ggsave(outsubsetfile_str,width=4,height=3)

				ggplot(data=res_tcasubset_tbl,aes(x=-log10(p_theor),y=-log10(pval),colour=inpulldown))+geom_point(size=2)+geom_abline(intercept=0,slope=1)+xlab("-log10(p-value) (theoretical)")+ylab("-log10(p-value) (empirical)")+geom_abline(intercept=1,slope=1,linetype=2)+scale_colour_manual(values=c(myblue,myorange,mygrey))+theme_bw()+ggtitle(plot_title_tcasubset)
				ggsave(outtcasubsetfile_str,width=4,height=3)


				write.table(res_tbl, file=outfull_str,sep="\t",quote=F,row.names=FALSE,col.names=TRUE)

				outtbl=head(res_subset_tbl[order(pval),list(gname, chirank, pval, betas)],5)
				write.table(outtbl, file=outsubsettbl_str,sep="\t",quote=F,row.names=FALSE,col.names=TRUE)

				outtbl=res_tcasubset_tbl[order(pval),list(gname, chirank, pval, betas)]
				write.table(outtbl, file=outtcasubsettbl_str,sep="\t",quote=F,row.names=FALSE,col.names=TRUE)
			}
		}
	}
}


dataset_vec=c("prot","rna","average")
myl=list()
count=0
for(USECSS in usecss_vec){
	for(DATASET in dataset_vec){
		for(CONDITIONAL in conditional_vec){
			for(USELMM in uselmm_vec){
				count=count+1
				if(USECSS){
					pheno_str="_wcss"
				}else{
					pheno_str="_pathwayact"
				}
				if(CONDITIONAL){
					cond_str="_conditional"
				}else{
					cond_str="_all"
				}
				if(DATASET=="prot"){
					type_used_str="prot"
				}else if(DATASET=="rna"){
					type_used_str="rnaseq"
					}else{
					type_used_str="avg"
				}
				if(USELMM){
					method_str=""
					}else{
					method_str="nolmm_"
				}
				outfull_str=paste("interimData/DEA/diff_exp_",method_str,type_used_str,pheno_str,cond_str,"_annotout.txt",sep="")
				ww=fread(outfull_str)
				ww[,uscss:=USECSS]
				ww[,cond:=CONDITIONAL]
				ww[,dataset:=DATASET]
				ww[,lmm:=USELMM]
				ww[,filepath:=outfull_str]
				myl[[count]]=ww
			}
		}
	}
}

fullTbl=do.call("rbind",myl)

write.table(fullTbl,file="interimData/DEA/FullTbl.txt",sep="\t",quote=F,row.names=FALSE,col.names=TRUE)
