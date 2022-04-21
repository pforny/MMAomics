##https://ftp.expasy.org/databases/enzyme/enzyme.dat
#more Data/gene_sets/enzyme.dat.txt | perl -ple 's/_HUMAN/_HUMAN\n/mg' | perl -nle '/_HUMAN/ && print' | perl -nle '/(\w+),\s(\w+_HUMAN)/ && print "$1\t$2"' > Data/genesets/all_enzyme_ids.tbl



#more Data/gene_sets/enzyme.dat.txt | perl -ple 's/_HUMAN/_HUMAN\n/mg' | perl -nle '/(ID\s+\d+\.\d+\.\d+\.\d+|_HUMAN)/ && print' | perl -ple 's/^(ID|DR|;)\s+//' | perl -ple 's/\s*;\s+/\n/mg' | perl -nle '/(\d+\.\d+\.\d+\.\d+|_HUMAN)/ && print' > Data/gene_sets/all_enzyme_ids2.tbl
require(data.table)
WW=fread("Data/gene_sets/all_enzyme_ids2.tbl",sep="\t",header=FALSE)[,V1]
system("mkdir interimData/gene_sets")
ee=grepl("\\d+\\.\\d+\\.\\d+\\.+\\d+",WW)
ee2=!grepl("\\d+\\.\\d+\\.\\d+\\.+\\d+",WW)

eee2=c(FALSE,ee2)
eee3=c(ee2,FALSE)
eee4=(eee2!=eee3) & eee3
eee5=(eee2!=eee3) & eee2
tt=cbind(which(eee4),which(eee5)-1)
ee3=intersect(which(ee),which(ee2)-1)
tt2=WW[ee3]
EE=lapply(c(1:length(ee3)),function(i){
	cbind(tt2[i],WW[tt[i,1]:tt[i,2]])
})
EE2=do.call("rbind",EE)
write.table(EE2,file="interimData/gene_sets/enzyme_sec.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

# perl -nle '/(\w+),\s(\w+_HUMAN)/ && print "$1\t$2"' > Data/genesets/all_enzyme_ids2.tbl

#mitochondrian direct
# UniProt search Subcellular localistions: mitochondrion

#http://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity,bioentity_name,qualifier,annotation_class,annotation_extension_json,assigned_by,taxon,evidence_type,evidence_with,panther_family,type,bioentity_isoform,reference,date&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=document_category:%22annotation%22&fq=regulates_closure:%22GO:0005739%22&fq=taxon_subset_closure_label:%22Homo%20sapiens%22&facet.field=aspect&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=evidence_subset_closure_label&facet.field=regulates_closure_label&facet.field=annotation_class_label&facet.field=qualifier&facet.field=annotation_extension_class_closure_label&facet.field=assigned_by&facet.field=panther_family_label&q=*:*
#######



ENZYMES_TBL = "Data/gene_sets/all_enzyme_ids.tbl"
PROTEIN_MAPPING_FILE="Data/HUMAN_9606_idmapping_selected.tab.gz"
OUT_TBL = "interimData/gene_sets/all_enzyme_ensg.tbl"
require(data.table)
#######
#######
#######


ww2=fread(cmd=paste("gunzip -cf ",PROTEIN_MAPPING_FILE,sep=""))
tt=ww2[,list(V1,V2,V19)]
tt2=tt[V19!="" & !grepl(";",V19),]
enzymes_tbl=fread(ENZYMES_TBL,header=FALSE)

setkey(enzymes_tbl,V2)
setkey(tt2,V2)
m3=merge(enzymes_tbl,tt2)
write.table(m3,file=OUT_TBL,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
#######
#######
#######
MITO_TBL = "Data/gene_sets/mitochondrion_uniprot.txt"
mitoprot = sub("UniProtKB:","",fread(MITO_TBL,header=FALSE)[,V1])

tbl3=tt2[is.element(V1,mitoprot),list(V19)]
write.table(tbl3,file = "interimData/gene_sets/mitoglist.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


#######
#######
#######


PULLDOWN_TBL = "Data/gene_sets/Pull_down_MUT_MCEE_MMAB_pooled_edited.csv"
tbl1 = fread(PULLDOWN_TBL,header=TRUE)

tbl2 = tbl1[,list(V2=sub(" .+","",`Accession Number`),pval=as.numeric(sub(".*?(\\d)","\\1",`ANOVA Test (p-value): (p < 0.05)`)))]
setkey(tbl2,V2)
setkey(tt2,V2)
tbl3=merge(tbl2,tt2)
write.table(tbl3, file = "interimData/gene_sets/pulldown_mmm.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


## for Patrick: you can create new subset lists. To make it consistent, use ensembl as identifier 


#######
#######
#######

# import TCA cycle file, pulled from the Human Protein Atlas with search term "protein_class:Citric acid cycle related proteins", resulting in https://www.proteinatlas.org/search/protein_class:Citric+acid+cycle+related+proteins 

TCA_PROTEINS_TBL = "Data/gene_sets/protein_class_Citric.tsv"
tib1 = fread(TCA_PROTEINS_TBL, header = TRUE)[,Ensembl]

# manual TCA gene set

TCA_PROTEINS_TBL_MAN = "Data/gene_sets/tca_geneset_manual.csv"
tib2 = fread(TCA_PROTEINS_TBL_MAN, header = FALSE)[,V1]


write.table(tib2, file = "interimData/gene_sets/tca_proteins.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")