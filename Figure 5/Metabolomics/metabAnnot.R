
##ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo

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
require(readxl)

### load some helper functions.
source("csrproject/Code/normalizationFunctions.R")
source("csrproject/Code/fast_lmm.R")



######################################################################
######################################################################
## generate single csv files for each sheet of the raw intracell metabolomics file

## input: froese1 DATA CURATED NORM.xlsx
##
## output: intracell_intensities1.csv
## output: intracell_injections.csv
## output: intracell_ions.csv
## output: intracell_annotation.csv
######################################################################
######################################################################

int_sheet <- read_excel("Data/metaboDat/intra/froese1 DATA CURATED NORM.xlsx", sheet = "intensities1")
write.csv(int_sheet, "Data/metaboDat/intra/intracell_intensities1.csv")

int_sheet <- read_excel("Data/metaboDat/intra/froese1 DATA CURATED NORM.xlsx", sheet = "injections")
write.csv(int_sheet, "Data/metaboDat/intra/intracell_injections.csv")

int_sheet <- read_excel("Data/metaboDat/intra/froese1 DATA CURATED NORM.xlsx", sheet = "ions")
write.csv(int_sheet, "Data/metaboDat/intra/intracell_ions.csv")

int_annot <- read_excel("Data/metaboDat/intra/froese1 DATA CURATED NORM.xlsx", sheet = "annotation")
write.csv(int_annot, "Data/metaboDat/intra/intracell_annotation.csv")



######################################################################
######################################################################
## generate single csv files for each sheet of the raw extracell metabolomics file

## input: froese2 DATA CURATED.xlsx
##
## output: intracell_intensities1.csv
## output: intracell_injections.csv
## output: intracell_ions.csv
## output: intracell_annotation.csv
######################################################################
######################################################################

int_sheet <- read_excel("Data/metaboDat/extra/froese2 DATA CURATED.xlsx", sheet = "intensities1")
write.csv(int_sheet, "Data/metaboDat/extra/extracell_intensities1.csv")

int_sheet <- read_excel("Data/metaboDat/extra/froese2 DATA CURATED.xlsx", sheet = "injections")
write.csv(int_sheet, "Data/metaboDat/extra/extracell_injections.csv")

int_sheet <- read_excel("Data/metaboDat/extra/froese2 DATA CURATED.xlsx", sheet = "ions")
write.csv(int_sheet, "Data/metaboDat/extra/extracell_ions.csv")

int_sheet <- read_excel("Data/metaboDat/extra/froese2 DATA CURATED.xlsx", sheet = "annotation")
write.csv(int_sheet, "Data/metaboDat/extra/extracell_annotation.csv")



############################################################
############################################################
## get gene2compound_tbl : gene names of enzymes and maps to 
## metabolites. 

## input: enzymeTbl.txt
## input: HUMAN_9606_idmapping.dat
## input: metabores.txt
## input: protid_gname.txt
## output:  gene2compound_tbl
############################################################
############################################################

############################################################
## load mapping tbls
############################################################
enz2prot_tbl=fread("Data/metaboAnnot/enzymeTbl.txt",header=FALSE)
protid_map_tbl=fread('cat Data/HUMAN_9606_idmapping.dat | grep "_HUMAN"',header=FALSE)
enz2compound_tbl0=fread("Data/metaboAnnot/metabores.txt",header=FALSE)

conj=fread("Data/metaboAnnot/chebi_obo.txt",header=FALSE)[,V2]

############################################################
## END: load mapping tbls
############################################################
############################################################
## prepare acids and base mappings
############################################################
TT=lapply(c(1:length(conj)),function(i){
	mysplit=strsplit(conj[i],split=";")[[1]]
	ff=kronecker(rep(1,length(mysplit)),c(1:length(mysplit)))
	ff2=kronecker(c(1:length(mysplit)),rep(1,length(mysplit)))
	out=unique(data.table(rh=mysplit[ff],lh=mysplit[ff2]))
	return(out)
})
TT2=do.call("rbind",TT)
TT2=unique(TT2[rh!=lh,])
TT2[,V2:=lh]
###
##TT3=TT2[1:dim(TT2)[1],list(rh,lh)]
getRR4=function(TT3){
	RR1=TT3[,list(v1=rh,v2=lh)]
	RR2=TT3[,list(v2=rh,v1=lh)]
	setkey(RR1,v1)
	setkey(RR2,v1)
	RR3=merge(RR1,RR2,allow.cartesian=TRUE)[,list(v1=v2.x,v2=v2.y),]
	RR4=unique(rbind(RR1,RR3))[,list(rh=v1,lh=v2)]
	return(RR4)
}
TT3=getRR4(getRR4(getRR4(TT2)))
###
TT3[,V2:=lh]
setkey(TT3,V2,physical=FALSE)
setkey(enz2compound_tbl0,V2,physical=FALSE)
RR=merge(TT3,enz2compound_tbl0)

#enz2compound_tbl0[,lis]
RR2=RR[,list(V1,V2=rh,V3,V4)]
enz2compound_tbl=unique(rbind(enz2compound_tbl0,RR2))
############################################################
## END: prepare acids and base mappings
############################################################
gene2protid_raw_tbl=fread("Data/protid_gname.txt")

############################################################
## END: load mapping tbls
############################################################
############################################################
## preprocess mapping  tbls
############################################################
gene2protid_raw_tbl[,protid:=`UniProtKB Gene Name ID`]
gene2protid_raw_tbl[,`UniProtKB Gene Name ID`:=NULL]
gene2protid_tbl=gene2protid_raw_tbl[protid!="",]

protid_map_tbl[,protname:=V3]
protid_map_tbl[,V3:=NULL]
protid_map_tbl[,protid:=V1]
protid_map_tbl[,V1:=NULL]
protid_map_tbl[,V2:=NULL]
enz2prot_tbl[,protname:=V1]
enz2prot_tbl[,V1:=NULL]
enz2prot_tbl[,enzyme:=V2]
enz2prot_tbl[,V2:=NULL]
enz2prot_tbl[,V3:=NULL]

############################################################
## END: preprocess mapping  tbls
############################################################

############################################################
## join tables
############################################################
setnames(enz2compound_tbl,c("enzyme","compound","eq","upst"))

setkey(protid_map_tbl,protname)
setkey(enz2prot_tbl,protname)
enz2protid_tbl=merge(enz2prot_tbl,protid_map_tbl)

setkey(enz2compound_tbl,enzyme)
setkey(enz2protid_tbl,enzyme)
protid2compound_tbl=merge(enz2compound_tbl,enz2protid_tbl)[,list(protid,compound,upst,protname,enzyme)]

setkey(protid2compound_tbl,protid)
setkey(gene2protid_tbl,protid)
gene2compound_tbl=merge(protid2compound_tbl,gene2protid_tbl)[,list(gname=`Gene name`,upst,compound,enzyme)]
write.table(gene2compound_tbl,file="interimData/metabo_annot/gene2compound.tbl",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
############################################################
## END: join tables
############################################################

############################################################
############################################################
## END: get gene2compound_tbl
############################################################
############################################################

############################################################
############################################################
## get ionid2gname_tbl : gene names of enzymes and maps to 
## ionIds 

## input: intracell_annotation.csv
## input: extracell_annotation.csv
## input: gene2compound_tbl

## output: ionid2gname_tbl
############################################################
############################################################


getTbl = function(ion_annot_mat_raw){
	ion_annot_tbl=ion_annot_mat_raw[!duplicated(ionIdx),]
	mylens=c(1:nrow(ion_annot_tbl))
	ion2compounds_list=lapply(mylens,function(i){
		ionIdx=ion_annot_tbl[i,ionIdx]	
		chebs=ion_annot_tbl[i,`other ids`]
		ee=strsplit(chebs,split=";\\s*",perl=T)[[1]]
		if(length(ee)==0){
			return(NULL)
		}else{
			return(data.table(ionIdx,compound=ee))
		}
	})
	ionid2compound_tbl=do.call("rbind",ion2compounds_list)
	ionid2compound_tbl=ionid2compound_tbl[grepl("CHEBI",compound)==TRUE,]

	setkey(ionid2compound_tbl,compound)
	setkey(gene2compound_tbl,compound)
	ionid2gname_raw_tbl=merge(ionid2compound_tbl,gene2compound_tbl)

	ionid2gname_tbl=unique(ionid2gname_raw_tbl[,list(ionIdx,upst,gname,compound)])
	setkey(ionid2gname_tbl,ionIdx)
	labelTbl=ion_annot_mat_raw[,list(ionIdx,compoundname=`label (bona fide)`)]
	setkey(labelTbl,ionIdx)
	ionid2gname_tbl=merge(ionid2gname_tbl,labelTbl)
	return(ionid2gname_tbl)
}

ion_annot_mat_raw=fread("Data/metaboDat/intra/intracell_annotation.csv")
ionid2gname_tbl=getTbl(ion_annot_mat_raw)
write.table(ionid2gname_tbl,file="interimData/metabo_annot/intraIonId2gname.tbl",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

ion_annot_mat_raw=fread("Data/metaboDat/extra/extracell_annotation.csv")
ionid2gname_tbl=getTbl(ion_annot_mat_raw)
write.table(ionid2gname_tbl,file="interimData/metabo_annot/extraIonId2gname.tbl",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
######################################################################
######################################################################
## END: get ionid2gname_tbl
######################################################################
######################################################################

######################################################################
######################################################################
## get id_tbl: hand-coded mapping between individual ids in mma ids vs 
## ids in metabolomics experiment.
######################################################################
######################################################################
aa=c("C4","C5","C6","M4","M7","M10","C8","C9","C10","M1","M3","M6","M8","M9","C1","C2","C3","C7")
bb=c("222","213","230","067","013","036","228","225","215","014","042","104","030","138","219","221","227","226")
matabo_id_tbl=data.table(sampleid=aa,mmaid=paste("MMA",bb,sep=""))

write.table(matabo_id_tbl,file="interimData/metabo_annot/metaboId.tbl",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
######################################################################
######################################################################
## END: get id_tbl
######################################################################
######################################################################

