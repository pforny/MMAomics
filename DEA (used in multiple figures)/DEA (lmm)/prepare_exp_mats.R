PROTEIN_FILE_PATH1=paste("Data/20190305_PHRT-5_MMA_Normalized-Protein-Report_edit_edit_edit_filtered.csv", sep="/")
#PROTEIN_FILE_PATH2=paste("Data/20191011_PHRT-5_MMA_Normalized-Protein-Report_edit_edit.csv",sep="/")

RNASEQ_FILE_PATH=paste("Data/MMA_PHRT_R_gene_rpkm.bed",sep="")
PROTEIN_MAPPING_FILE="Data/HUMAN_9606_idmapping_selected.tab.gz"


##################################################
## prepare prot_mat
##################################################
require(data.table)
prot_tbl0=fread(PROTEIN_FILE_PATH1)
#prot_tbl0=fread(PROTEIN_FILE_PATH2,dec=",")
	#prot_tbl0=prot_tbl0[,4:ncol(prot_tbl0)]
	#prot_tbl0=prot_tbl0[, lapply(.SD, as.numeric), by="PG.ProteinNames"]

ww2=fread(cmd=paste("gunzip -cf ",PROTEIN_MAPPING_FILE,sep=""))
tt=ww2[,list(V2,V19)]
tt2=tt[V19!="" & !grepl(";",V19),]

prot_tbl0_subset=prot_tbl0[is.element(prot_tbl0[,PG.ProteinNames],tt2[,V2]),]
ordered_idmap=tt2[match(prot_tbl0_subset[,PG.ProteinNames],tt2[,V2]),]
sum(ordered_idmap[,V2]!=prot_tbl0_subset[,PG.ProteinNames])

prot_mat0=as.matrix(prot_tbl0_subset[,c(5:ncol(prot_tbl0_subset)),with=FALSE])
#prot_mat0=as.matrix(prot_tbl0_subset[,c(2:ncol(prot_tbl0_subset)),with=FALSE])

rownames(prot_mat0)=ordered_idmap[,V19]
colnames(prot_mat0)=sub(".*MMA_(\\d+).*","MMA\\1",colnames(prot_mat0))
prot_mat=prot_mat0[,!duplicated(colnames(prot_mat0))]##remove technical replicates
save(prot_mat,file="interimData/prot_mat.RData")

####################################
####################################
rnaseq_tbl=fread(cmd=paste("less ",RNASEQ_FILE_PATH, " | perl -ple 's/\\t$//'"))
rnaseq_tbl[,pid:=sub("\\.\\d+","",pid)]
rnaseq_mat=as.matrix(rnaseq_tbl[,c(7:ncol(rnaseq_tbl)),with=FALSE])
rownames(rnaseq_mat)=rnaseq_tbl[,pid]

trimmed_names=sub(".*(MMA\\d+).*","\\1",colnames(rnaseq_mat))
colnames(rnaseq_mat)=trimmed_names
save(rnaseq_mat,file="interimData/rnaseq_mat.RData")

####################################
# filter proteincoding; less than half 0
####################################
require(data.table)
ww2=fread(cmd='gunzip -cf Data/gencode.v31lift37.basic.annotation.gff3.gz  | grep "\tgene\t" | grep "protein_coding"')
myids=ww2[,sub("ID=","",sub(";.+","",V9))]
genenames=ww2[,sub(".*gene_name=(.+?);.+","\\1",V9)]
myids2=sub("\\..+","",myids)
rnafiltered_mat=rnaseq_mat[is.element(rownames(rnaseq_mat),myids2),]
rnafiltered_mat=rnafiltered_mat[rowSums(rnafiltered_mat==0)<dim(rnafiltered_mat)[2]/2,]
save(rnafiltered_mat,file="interimData/rnafiltered_mat.RData")
write.csv(rnafiltered_mat, file="interimData/rnafiltered_mat.csv", row.names=TRUE, quote=FALSE)
####################################
# write gene mapping file
####################################
ww3=ww2[,list(chr=V1,start=V4,end=V5,strand=V7,origin=V2,ensembl=myids2,gname=genenames)]
write.table(ww3,file="interimData/protein_coding.bed",row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

