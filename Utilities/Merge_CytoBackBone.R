Merge_CytoBackBone <- function(patient,Qnorm=TRUE) {
  
  ########### LOAD PACKAGES #############
  
  library(flowCore) 
  library(data.table)
  library(plyr)
  library(aroma.light)
  library(FlowSOM)
  source("Utilities/get.nann.R")
  
  ########## IMPORT  ##################
  
  filenames=list.files(patient,pattern = ".fcs")
  
  files=list()
  ncells=integer()
  markers=list()
  
  Geom=c("FSC.A","FSC.H","SSC.A","SSC.H")
  
  for (i in 1:length(filenames)){
    
    # Read FCS
    data=read.FCS(paste0(patient,"/",filenames[i]),truncate_max_range = FALSE,column.pattern = "Time",invert.pattern = TRUE) # Read
    files[[i]]=as.data.frame(data@exprs)
  }
  
  ncells=as.integer(lapply(files, nrow))
  markers=lapply(files,colnames)
  
  ########## CREATE MERGED FILE ############
  
  # Get common part
  common=Reduce(intersect,markers)
  common=common[which(common%in%Geom==FALSE)] 
  
  Merged_file=rbindlist(files, fill = TRUE,idcol = TRUE)
  all=colnames(Merged_file)
  left=list()
  for (i in 1:length(files)){
    left[[i]]=all[-c(1,which(all%in%markers[[i]]))]
  }

  
  ########### QUANTILE NORMALIZATION ##########
  
  if (Qnorm==TRUE){
    
    # Group common markers
    files_common=lapply(files,"[",common)
    common_joint=rbindlist(files_common, fill = TRUE,idcol = TRUE)
    common_joint=flowFrame(as.matrix(common_joint))
    
    # Cluster
    fSOM = FlowSOM(common_joint,xdim = 20, ydim = 20, colsToUse = c(2:ncol(common_joint)), nClus = 5,transform = FALSE,compensate = FALSE,scale = FALSE) 
    cluster = GetMetaclusters(fSOM)
    cj_clustered=as.data.frame(cbind(Merged_file,cluster))
    
    # Quantile normalize
    for (i in as.numeric(levels(cluster))){
      keep=which(cj_clustered$cluster==i)
      cl=cj_clustered[keep,]
      for (to_norm in common){
        clm=cl[,to_norm] 
        files_split=split(clm,cl$.id)
        
        short=unname(which(lapply(files_split,length)<2))
        if (length(short)!=0){short_val=files_split[short]}
        if (length(short)!=0){files_split=files_split[-short]}
        
        if (length(files_split)>1){
          marknorm=normalizeQuantile(files_split)
        }else{marknorm=files_split[[1]]}
        
        if (length(short)!=0){for (j in 1:length(short)){marknorm=append(marknorm,short_val[j],short[j])}}
        
        cl[,to_norm]=unlist(marknorm)
      }
      cj_clustered[keep,]=cl
    }
    
    Merged_file=cj_clustered
    Merged_file$cluster=NULL
    Merged_file=as.data.frame(Merged_file)
    Merged_updated=Merged_file
    
  }

  ######### INPUTATION #############
  
  if (length(left)>1){
  
  for (i in sample(1:length(files))){ # for each file
    
    # Initialize list of reference tubes
    reference=list()
    reference[[1]]=0
    
    for (j in 1:length(left[[i]])){ # for each missing marker
      
      # Select marker j in tube i
      marker=left[[i]][j]
      
      # Find which tubes have this marker
      reference[[j+1]]=which(lapply(markers,function(x) marker%in%x)==TRUE)
      
      if (paste(reference[[j+1]],collapse = "")==paste(reference[[j]],collapse = "")){ # Copy from previous iteration
        
        Merged_updated[cells_unmeasured[nn[[1]]],marker]=Merged_file[cells_measured[nn[[2]]],marker]
        
      }else{# compute neighbors
        
        cells_unmeasured=which(Merged_file$.id==i)
        cells_measured=which(Merged_file$.id%in%reference[[j+1]])
        
        df_unmeasured=Merged_file[cells_unmeasured,common]
        df_measured=Merged_file[cells_measured,common]
        
        # Compute neighbors
        nn=get.nann(df_unmeasured,df_measured,1) 
        
        # Impute
        Merged_updated[cells_unmeasured[nn[[1]]],marker]=Merged_file[cells_measured[nn[[2]]],marker]
        
      }
    }
  }
    
  }
  
  
  Merged_updated=na.omit(subset(Merged_updated, select=-.id))
  Merged_fcs=flowFrame(as.matrix(Merged_updated))

  return(Merged_fcs)
  
  
}