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
### load some helper functions.
source("csrproject/Code/normalizationFunctions.R")
source("csrproject/Code/fast_lmm.R")
############################################################
############################################################
id_tbl=fread("interimData/metabo_annot/metaboId.tbl")


######################################################################
######################################################################
## get ion_intensity_tbl

## input: extracell_intensities1.csv
##
## output: ion_intensity_tbl
######################################################################
######################################################################
metaboVals=fread("Data/metaboDat/extra/extracell_intensities1.csv")
expid=colnames(metaboVals)[4:ncol(metaboVals)]
sampleid=sub("-\\d.+","",sub(".+E\\d ","",expid))
timeid=as.numeric(sub(".+-\\d+ T","",expid))
batchid=sub("(E\\d+).+","\\1",sub(".+E","E",expid))
expcl=data.table(sampleid,timeid,batchid,expid)
mymatcher=match(expcl[,sampleid],id_tbl[,sampleid])
mmaids=id_tbl[mymatcher,mmaid]
expcl[,mmaid:=mmaids]



metaboM=as.matrix(metaboVals[,4:ncol(metaboVals),with=FALSE])
rownames(metaboM)=metaboVals[,ionIdx]
meantbl=data.table(medianval=apply(log(metaboM),2,median),expid=colnames(metaboM))
setkey(meantbl)
setkey(expcl)
tt=merge(meantbl,expcl)
tt[,expnum:=as.numeric(sub(" /.+","",expid))]

#ggplot(tt, aes(x=medianval))+geom_density()+facet_grid(timeid~.)
#ggplot(tt, aes(x=medianval))+geom_density()+facet_grid(batchid~.)
#ggplot(tt, aes(y=medianval,x=sampleid))+geom_boxplot()


require(reshape2)
ww=data.table(melt(metaboM))
ww[,expid:=Var2]
setkey(ww,expid)
setkey(expcl,expid)
ww2=merge(ww,expcl)

#ggplot(ww2, aes(x=expid,y=log(value),group= expid)) + geom_boxplot()

exp_sc=QQnormalizeMatCol(metaboM)

metaboMLog=log(metaboM)
#exp_sc=log(metaboM)

#####################################
## evaluate replicability
#####################################
ee2=scalerRound(metaboMLog,7)

colnames(ee2)=colnames(metaboMLog)
rownames(ee2)=rownames(metaboMLog)

matcher=match(colnames(ee2),expcl[,expid])
expcl_ordered=expcl[matcher,]

mycor=t(ee2)%*%ee2/nrow(ee2)

expcl_newordered=expcl_ordered[1:nrow(expcl_ordered),]
setkey(expcl_newordered,sampleid)


newmatcher = match(expcl_newordered[,expid],colnames(mycor))
mycor2=mycor[newmatcher,newmatcher]
image(mycor2)

selectr=expcl_newordered[,timeid==36]
mycor3 = mycor2[selectr,selectr]
image(mycor3)
#######################
eigv=eigen(mycor3)$vectors[,1]
batchnr=grepl("E2",colnames(mycor3))
plot(eigv,batchnr)
summary(lm(eigv~batchnr))
#######################

pdf("Figs/extraccell.pdf")
image(mycor3)
dev.off()


mycorMelted=data.table(melt(mycor3))

mycorMelted[,Sample1:=sub("-\\d+ T36","",sub("^\\d+","",Var1))]
mycorMelted[,Sample2:=sub("-\\d+ T36","",sub("^\\d+","",Var2))]

mycorMelted2=mycorMelted[Var1!=Var2,]


mysamples = unique(mycorMelted2[,Sample1])
mysamplesE1 = mysamples[grepl("E1",mysamples)]
mysamplesE2 = mysamples[grepl("E2",mysamples)]

set.seed(11)
medianval = median(mycorMelted2[Sample2==Sample1,value])
##ggplot(mycorMelted2,aes(x=value, fill = Sample2==Sample1) + geom_density()
randmedians=sapply(c(1:1000), function(dummy){
	# matcherTblE1 = data.table(x1=mysamplesE1,x2=sample(mysamplesE1))	
	# matcherTblE2 = data.table(x1=mysamplesE2,x2=sample(mysamplesE2))
	# matcherTbl = rbind(matcherTblE1, matcherTblE2)
	 matcherTbl = rbind(x1=mysamples,x2=sample(mysamples))
	var1=mycorMelted2[,Sample1]
	matcher=match(var1,matcherTbl[,x2])
	sample1Shuf=matcherTbl[matcher,x1]
	mycorMelted2[,sample1Shuffled:=sample1Shuf]
	median(mycorMelted2[Sample2==sample1Shuffled,value])
})


#####################################
#####################################


FF=lapply(id_tbl[,sampleid],function(x){rowMeans(ee2[,sampleid==x & grepl("T0",colnames(ee2))])})
ion_intensity_unnormed_mat=do.call("cbind",FF)
colnames(ion_intensity_unnormed_mat)=id_tbl[,mmaid]
rownames(ion_intensity_unnormed_mat)=metaboVals[,ionIdx]

ion_intensity_unnormed_mat=ion_intensity_unnormed_mat[,!is.na(ion_intensity_unnormed_mat[1,])]




ion_intensity_expand_mat=scalerRound(ion_intensity_unnormed_mat,7)
colnames(ion_intensity_expand_mat)=colnames(ion_intensity_unnormed_mat)
rownames(ion_intensity_expand_mat)=metaboVals[,ionIdx]





iscontrol=unique(expcl[grepl("C",sampleid),mmaid])
isc=is.element(colnames(ion_intensity_expand_mat),iscontrol)+0


resdf=data.table(fast_lmm(my_y=scale(isc),my_X=t(ion_intensity_expand_mat)))
res_tbl=data.table(resdf)
res_tbl[,ionIdx:=as.numeric(rownames(ion_intensity_expand_mat))]

setkey(res_tbl,ionIdx)


ion_annot_mat_raw=fread("Data/metaboDat/intra/intracell_annotation.csv")

compoundname_tbl=unique(ion_annot_mat_raw[,list(ionIdx,compoundname=`label (bona fide)`)])
setkey(compoundname_tbl,ionIdx)

out_tbl=merge(res_tbl,compoundname_tbl)


out_tbl[,pval:=1-pchisq(chi_sq,1)]
out_tbl[,p_theor:=rank(pval)/(dim(out_tbl)[1]+1)]
ggplot(data=out_tbl,aes(x=-log10(p_theor),y=-log10(pval)))+geom_point(size=2)+geom_abline(intercept=0,slope=1)+xlab("-log10(p-value) (theoretical)")+ylab("-log10(p-value) (empirical)")+geom_abline(intercept=1,slope=1,linetype=2)+theme_bw()

######################################################################
######################################################################
## END: get ion_intensity_tbl
######################################################################
######################################################################


