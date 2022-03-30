#'\code{PAEM}
#'
#' Pathway Analysis using Exact Mass (PAEM) - this function annotate using exact mass and  performs a pathway analysis (overrepresentation analysis) using ranked metabolites
#'@param file diffreport file.
#'@param mode 'positive' or 'negative'
#'@param accuracy Mass accuracy: 0.01, 0,001   
#'@param metabolitesDatabase IDs and mass of metabolites taken from database. Default: mref_kegg.xlsx
#'@param pathwayDatabase Metabolites present in each pathway
#'@param exportAnnotation TRUE/FALSE. Export the annotation. File 
#'@param sortList How to sort the input. 'pvalue' 'fold' 'tstat-positive' 'tstat-negative' 'tstat-absolute'
#'@return pathways enrichment. 
#'
PAEM=function(file='diffreport_MMA_liver_koki_vs_wt-pos_original.csv',mode='positive',accuracy=0.01,metabolitesDatabase='mref_kegg.xlsx',pathwayDatabase='hsa_HumanSpecificMetab.csv',exportAnnotation=FALSE,exportPath='Results/',sortList='pvalue'){
  require(ggplot2)
  require(openxlsx)
  require(stringr)
  require(data.table)

  ## Annotation
  # Load kegg metabolites and mass
  mrefKEGG=read.xlsx(metabolitesDatabase,  colNames = FALSE)
  
  ## Mass H
  HMass=1.007276 
  
  # add H+ (Positive mode) H- (Negative mode)
  if(mode=='positive'){mrefKEGG$MassProtonated=mrefKEGG$X2+HMass
  }else if(mode=='negative'){mrefKEGG$MassProtonated=mrefKEGG$X2-HMass
  }else{stop("What is the ionisation mode? ")}
  
  # Load peak list 
  liverMouse_pos=read.csv(file)
  
  # Sort ions
  if(sortList=='pvalue'){liverMouse_pos=liverMouse_pos[order(liverMouse_pos$pvalue),] # Smallest to biggest
  }else if(sortList=='fold'){liverMouse_pos=liverMouse_pos[order(liverMouse_pos$fold,decreasing = T),] # Biggest to smallest
  }else if(sortList=='tstat-positive'){liverMouse_pos=liverMouse_pos[order(liverMouse_pos$tstat,decreasing = T),] # Biggest to smallest
  }else if(sortList=='tstat-negative'){liverMouse_pos=liverMouse_pos[order(liverMouse_pos$tstat),] # Smallest to biggest
  }else if(sortList=='tstat-absolute'){liverMouse_pos=liverMouse_pos[order(abs(liverMouse_pos$tstat),decreasing = T),] # abs(tstat) - biggest to smallest
  }else{stop("How do you want to sort your ions?")}
  
  # Matching mass to KEGG DB
  massIons=as.numeric(liverMouse_pos$mzmed)
  liverMouse_pos$KEGGID=c()
  liverMouse_pos$KEGGName=c()
  
  for(i in 1:length(massIons)){
    # find position de protonated 
    posMref=intersect(which(mrefKEGG$MassProtonated>(massIons[i]-accuracy)), which(mrefKEGG$MassProtonated<(massIons[i]+accuracy)))
    liverMouse_pos$KEGGID[i]=paste0(mrefKEGG$X4[posMref],collapse = ':')
    liverMouse_pos$KEGGName[i]=paste0(mrefKEGG$X3[posMref],collapse = ':')
  }
  
  # Remove all ions that were not annotated
  posWithKEGG=!liverMouse_pos$KEGGID==''
  liverMouse_pos_withKEGG=liverMouse_pos[posWithKEGG,]
  
  # Export Annotation
  if(exportAnnotation){
    fileAnnotationName=paste0(exportPath, gsub(".csv","", gsub(".*diffreport_MMA_", "", file)),'_annotation_accuracy_',accuracy,'.csv')
    write.csv(liverMouse_pos,fileAnnotationName)
  }
  
  ## Pathway Analysis
  # Loading KEGG Pathways
  KEGG_Pathways=read.csv(pathwayDatabase,header = F)
  KEGG_Pathways$V1=gsub('cpd:','',KEGG_Pathways$V1)
  posPathways=which(regexpr('>path:.*',KEGG_Pathways$V1)>0) # get pathways position
  KEGG_PathwaysName=KEGG_Pathways$V1[posPathways]
  KEGG_PathwaysName=gsub('>.*# ','',KEGG_PathwaysName) # clean pathway name
  posPathways=c(posPathways,length(KEGG_Pathways$V1)) # add last pathway
  
  # Convert to data.table
  liverMouse_pos_withKEGG=data.table(liverMouse_pos_withKEGG)
  KEGG_Pathways=data.table(KEGG_Pathways)
  
  # Data to store
  pathwayCovered=c()
  numberOfmetab=length(liverMouse_pos_withKEGG$KEGGID)
  pvalPathways=c()
  positionMetabAllPathways=list()
  
  # For each pathway
  NumPath=length(KEGG_PathwaysName)
  i=1
  while(i<=NumPath){
    # Get the metabolite of that pathway
    metabInPathway=KEGG_Pathways$V1[(posPathways[i]+1):(posPathways[i+1]-1)]
    
    # Check how many metabolites of that pathway you measured
    posMetabMeasured_tmp=sapply(metabInPathway,regexpr,liverMouse_pos_withKEGG$KEGGID)
    posMetabMeasured=apply(posMetabMeasured_tmp!=-1, 1, any)
    numberMetabMeasured=sum(posMetabMeasured)
  
    pvalPathSpec=c()
    if(numberMetabMeasured>0){
      # Find position of measured metabolite 
      posMeasured=which(posMetabMeasured_tmp!=-1,arr.ind = T)
      positionMetab=data.frame(Position=posMeasured[,1],Metab=colnames(posMetabMeasured_tmp)[posMeasured[,2]])
      positionMetabAllPathways[[length(positionMetabAllPathways)+1]]=positionMetab
      
      # run over list of rank metabolite
      j=2
      while(j<=numberOfmetab){
        # Check if metabolite is in the top of list
        posMetab =posMetabMeasured[1:j]
        numberSignMetab=sum(posMetab)
        
        # Running a hypergeometric test
        pval=phyper(numberSignMetab-1,numberMetabMeasured,numberOfmetab-numberMetabMeasured,j,lower.tail=FALSE)
        pvalPathSpec=c(pvalPathSpec,pval)
        j=j+1
      }
      pathwayCovered=c(pathwayCovered,KEGG_PathwaysName[i])
    }
    pvalPathways=rbind(pvalPathways,pvalPathSpec) 
    i=i+1
  }
  
  # New - to export pvalue and leading metabolites
  pathways_pval_min=data.frame()
  for(i in 1:length(pathwayCovered)){
    posMin=which.min(pvalPathways[i,])
    pvalue=pvalPathways[i,posMin]
    positionMetabThisPath=positionMetabAllPathways[[i]]
    metabLeading=positionMetabThisPath[positionMetabThisPath[,1]<=posMin+1,2] # have to add 1 because we only started computing hyergeometric test after 2 metabolites
    metabLeadingName=paste0(mrefKEGG$X3[match(metabLeading,mrefKEGG$X4)],collapse=',')
    metabLeadingKEGGID=paste0(metabLeading,collapse=',')
    pathways_pval_min=rbind(pathways_pval_min,data.frame(Pathway=pathwayCovered[i],pvalue=pvalue,metabLeadingName=metabLeadingName,metabLeadingKEGGID=metabLeadingKEGGID))
  }
  rownames(pathways_pval_min)=c()
  # Old - to export only pvalue
  #pathways_pval_min=apply(pvalPathways, 1, min,na.rm=T)
  #names(pathways_pval_min)=pathwayCovered
  
  return(pathways_pval_min)
}

