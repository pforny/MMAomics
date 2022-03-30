# more /Users/davidlamparter/Documents/workspace/Fcbg/mapMetabolites/testdata/metaboliteDBs/chebi.obo | perl -nle '/(^id:|is_conjugate|is_tautomer)/ && print' | perl -pe 's/\n/\t/mg' > /Users/davidlamparter/Documents/workspace/Fcbg/mapMetabolites/results/chebi_tmp_obo.txt
#more /Users/davidlamparter/Documents/workspace/Fcbg/mapMetabolites/results/chebi_tmp_obo.txt | perl -pe 's/\tid: CHEBI/\nid: CHEBI/mg' | grep "is_conjugate" | perl -nle '/(CHEBI:\d+).+(base|acid).+(CHEBI:\d+)/ && print "$1\t$2\t$3"' > /Users/davidlamparter/Documents/workspace/Fcbg/mapMetabolites/results/chebi_obo.txt 


#more /Users/davidlamparter/Documents/workspace/Fcbg/mapMetabolites/results/chebi_tmp_obo.txt | perl -pe 's/\tid: CHEBI/\nid: CHEBI/mg' | grep "is_conjugate" |  grep "relationship" | perl -ple 's/\t	relationship: is_(conjugate_acid|conjugate_base|tautomer)_of /;/mg' > /Users/davidlamparter/Documents/workspace/Fcbg/mapMetabolites/results/chebi_obo.txt 


#perl -nle '/(CHEBI:\d+).+(base|acid).+(CHEBI:\d+)/ && print "$1\t$2\t$3"' > /Users/davidlamparter/Documents/workspace/Fcbg/mapMetabolites/results/chebi_obo.txt 



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


USE_ATP=FALSE
require(data.table)
require(ggplot2)
require(reshape2)
require(readxl)
### load some helper functions.
source("csrproject/Code/normalizationFunctions.R")
source("csrproject/Code/fast_lmm.R")
############################################################
system("mkdir Figs/metabo")
fig_path <- c("Figs/metabo/")

id_tbl=fread("interimData/metabo_annot/metaboId.tbl")


######################################################################
######################################################################
## get ion_intensity_tbl

## input: intracell_intensities1.csv
##
## output: ion_intensity_tbl
######################################################################
######################################################################
metaboVals=fread("Data/metaboDat/intra/intracell_intensities1.csv")
expid=colnames(metaboVals)[5:ncol(metaboVals)]
sampleid=sub("-\\d","",sub(".+E\\d ","",expid))
repid=sub(".+-","",expid)
batchid=sub("(E\\d+).+","\\1",sub(".+E","E",expid))
expcl=data.table(sampleid,repid,batchid,expid)
mymatcher=match(expcl[,sampleid],id_tbl[,sampleid])
mmaids=id_tbl[mymatcher,mmaid]
expcl[,mmaid:=mmaids]
expcl[,biorep:=sub("\\d+ / ","",expid)]




metaboM=as.matrix(metaboVals[,5:ncol(metaboVals),with=FALSE])

rownames(metaboM)=metaboVals[,ionIdx]

require(reshape2)
ww=data.table(melt(metaboM))
ww[,expid:=Var2]
setkey(ww,expid)
setkey(expcl,expid)
ww2=merge(ww,expcl)

write.table(ww2,file="interimData/intracell_metab_annotated.tbl",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# ggplot(ww2, aes(x=expid,y=log(value),group= expid)) + geom_boxplot()

exp_sc=QQnormalizeMatCol(metaboM)
#exp_sc=log(metaboM)

#####################################
## evaluate replicability
#####################################
## evaluating correlation biological and technical replicate level,
## (technical replicates non collapsed)

ee2=scalerRound(exp_sc,7)
colnames(ee2)=colnames(exp_sc)
rownames(ee2)=rownames(exp_sc)

matcher=match(colnames(ee2),expcl[,expid])
expcl_ordered=expcl[matcher,]

mycor=t(ee2)%*%ee2/nrow(ee2)

expcl_newordered=expcl_ordered[1:nrow(expcl_ordered),]
setkey(expcl_newordered,sampleid)

newmatcher = match(expcl_newordered[,expid],colnames(mycor))
mycor2=mycor[newmatcher,newmatcher]
pdf("Figs/metabo/intracell_cor.pdf")
image(mycor2)
dev.off()


#####################################
#####################################
#####################################
## evaluating correlation on biological replicate level,
## (technical replicates collapsed)

biorepnames = unique(expcl[,biorep])
tmp = sub("\\d+ / ","",colnames(exp_sc))
biorepMat=do.call("cbind",lapply(biorepnames,function(x){rowMeans(exp_sc[,tmp==x])}))
biorepMatSc = scalerRound(biorepMat,7)

colnames(biorepMatSc)=biorepnames
rownames(biorepMatSc)=rownames(exp_sc)

biorep_expcl=unique(expcl[,list(sampleid,biorep)])
setkey(biorep_expcl,sampleid)

matcher=match(biorep_expcl[,biorep],colnames(biorepMatSc))

biorepMatScSorted=biorepMatSc[,matcher]

mycor3=t(biorepMatScSorted)%*%biorepMatScSorted/nrow(biorepMatScSorted)
image(mycor3)

mycorMelted=data.table(melt(mycor3))
mycorMelted=mycorMelted[Var1!=Var2,]
mycorMelted[,Sample1:=sub("-\\d+$","",Var1)]
mycorMelted[,Sample2:=sub("-\\d+$","",Var2)]
#mycorMelted[Sample1==Sample2,]
mycorMelted[,isRep:=Sample1==Sample2]

ggplot(mycorMelted,aes(x=value,fill=isRep))+geom_density(alpha=0.5)
ggsave("Figs/metabo/intracell_hist_BioRep.pdf", width=3.5,height=3.5)


mysamples = unique(mycorMelted[,Sample1])
mysamplesE1 = mysamples[grepl("E1",mysamples)]
mysamplesE2 = mysamples[grepl("E2",mysamples)]


medianval = median(mycorMelted[Sample2==Sample1,value])
randmedians=sapply(c(1:2000), function(dummy){
	matcherTblE1 = data.table(x1=mysamplesE1,x2=sample(mysamplesE1))	
	matcherTblE2 = data.table(x1=mysamplesE2,x2=sample(mysamplesE2))
	matcherTbl = rbind(matcherTblE1, matcherTblE2)
	## matcherTbl =  data.table(x1=mysamples,x2=sample(mysamples))
	var1=mycorMelted[,Sample1]
	matcher=match(var1,matcherTbl[,x2])
	sample1Shuf=matcherTbl[matcher,x1]
	mycorMelted[,sample1Shuffled:=sample1Shuf]
	median(mycorMelted[Sample2==sample1Shuffled,value])
})

ggplot(data.table(val=randmedians),aes(x=val))+geom_histogram()+geom_vline(xintercept=medianval)
ggsave(paste("Figs//metabo/histog_intra_BioRep.pdf",sep=""))







mycorMelted=data.table(melt(mycor2))
mycorMelted=mycorMelted[Var1!=Var2,]
mycorMelted[,tech1:=sub("^\\d+","",Var1)]
mycorMelted[,tech2:=sub("^\\d+","",Var2)]
mycorMelted[,Sample1:=sub("^\\d+","",sub("-\\d+$","",Var1))]
mycorMelted[,Sample2:=sub("^\\d+","",sub("-\\d+$","",Var2))]
#mycorMelted[Sample1==Sample2,]
mycorMelted[,isRep:=(Sample1==Sample2) + (tech1==tech2)]
mycorMelted[,isRepF:=factor(isRep)]
ggplot(mycorMelted,aes(x=value,fill=isRepF))+geom_density(alpha=0.5)
ggsave(paste("Figs//metabo/intracell_hist_TechBioRep.pdf",sep=""))



medianval = median(mycorMelted[Sample2==Sample1,value])


mysamples = unique(mycorMelted[,Sample1])
mysamplesE1 = mysamples[grepl("E1",mysamples)]
mysamplesE2 = mysamples[grepl("E2",mysamples)]


randmedians=sapply(c(1:5000), function(dummy){
	matcherTblE1 = data.table(x1=mysamplesE1,x2=sample(mysamplesE1))	
	matcherTblE2 = data.table(x1=mysamplesE2,x2=sample(mysamplesE2))
	matcherTbl = rbind(matcherTblE1, matcherTblE2)
	## matcherTbl =  data.table(x1=mysamples,x2=sample(mysamples))
	var1=mycorMelted[,Sample1]
	matcher=match(var1,matcherTbl[,x2])
	sample1Shuf=matcherTbl[matcher,x1]
	mycorMelted[,sample1Shuffled:=sample1Shuf]
	median(mycorMelted[Sample2==sample1Shuffled,value])
})

ggplot(data.table(val=randmedians),aes(x=val))+geom_histogram()+geom_vline(xintercept=medianval)
ggsave(paste("Figs//metabo/histog_intra_TechBioRep.pdf",sep=""))


#####################################
#####################################

#####################################
#####################################
FF=lapply(id_tbl[,sampleid],function(x){rowMeans(exp_sc[,sampleid==x])})
ion_intensity_unnormed_mat=do.call("cbind",FF)
colnames(ion_intensity_unnormed_mat)=id_tbl[,mmaid]
rownames(ion_intensity_unnormed_mat)=metaboVals[,ionIdx]

ion_intensity_expand_mat=scalerRound(ion_intensity_unnormed_mat,7)
colnames(ion_intensity_expand_mat)=colnames(ion_intensity_unnormed_mat)
rownames(ion_intensity_expand_mat)=metaboVals[,ionIdx]
#####################################
#####################################
## evaluate  correlation among true biological replicates
colnames(ion_intensity_expand_mat)

expclSample = unique(expcl[,list(sampleid,batchid,mmaid,sid=paste(batchid,sampleid),isc=sub("\\d+","",sampleid))])
setkey(expclSample,isc)
matcher = match(expclSample[,mmaid],colnames(ion_intensity_expand_mat))

ion_intensity_expand_matSorted=ion_intensity_expand_mat[,matcher]
colnames(ion_intensity_expand_matSorted)=expclSample[,sampleid]
mycor4=t(ion_intensity_expand_matSorted)%*%ion_intensity_expand_matSorted/ncol(ion_intensity_expand_matSorted)


mycorMelted=data.table(melt(mycor4))
mycorMelted=mycorMelted[Var1!=Var2,]
mycorMelted[,rep1:=sub("\\d+","",Var1)]
mycorMelted[,rep2:=sub("\\d+","",Var2)]
mycorMelted[,isRep:=(rep1==rep2)]
ggplot(mycorMelted,aes(x=value,fill=isRep))+geom_density(alpha=0.5)
#####################################
#####################################

iscontrol=unique(expcl[grepl("C",sampleid),mmaid])
isc=is.element(colnames(ion_intensity_expand_mat),iscontrol)+0


resdf=data.table(fast_lmm(my_y=scale(isc),my_X=t(ion_intensity_expand_mat)))
res_tbl=data.table(resdf)
res_tbl[,ionIdx:=as.numeric(rownames(ion_intensity_expand_mat))]

setkey(res_tbl,ionIdx)
compoundname_tbl=unique(ion_annot_mat_raw[,list(ionIdx,compoundname=`label (bona fide)`)])
setkey(compoundname_tbl,ionIdx)

out_tbl=merge(res_tbl,compoundname_tbl)


out_tbl[,pval:=1-pchisq(chi_sq,1)]
out_tbl[,p_theor:=rank(pval)/(dim(out_tbl)[1]+1)]
ggplot(data=out_tbl,aes(x=sort(-log10(p_theor)),y=sort(-log10(pval))))+geom_point(size=2)+geom_abline(intercept=0,slope=1)+xlab("-log10(p-value) (theoretical)")+ylab("-log10(p-value) (empirical)")+geom_abline(intercept=1,slope=1,linetype=2)+theme_bw()


######################################################################
######################################################################
## END: get ion_intensity_tbl
######################################################################
######################################################################


