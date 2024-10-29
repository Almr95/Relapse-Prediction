EMD_fcs <- function(dir_files,dir_merged) {
  
  ########### LOAD PACKAGES #############
  
  library(flowCore) 
  library(EMDomics)
  source("/Users/alvaromartinezrubio/Desktop/RELAPSE IN B-ALL/Relapse Prediction (gH)/Utilities/extreme_quant_norm.R")
  
  ########## IMPORT  ##################
  
  # Read aliquots
  filenames=list.files(dir_files,pattern = ".fcs",full.names = TRUE)
  files=list()
  ncells=integer()
  markers=list()
  Geom=c("FSC.A","FSC.H","SSC.A","SSC.H")
  
  for (i in 1:length(filenames)){
    
    # Read FCS
    data=read.FCS(filenames[i],truncate_max_range = FALSE,column.pattern = "Time",invert.pattern = TRUE) # Read
    files[[i]]=as.data.frame(data@exprs)
    #files[[i]]=as.data.frame(lapply(files[[i]],extreme_quant_norm))
  }
  
  ncells=as.integer(lapply(files, nrow))
  markers=lapply(files,colnames)
  
  # Get common part
  common=Reduce(intersect,markers)
  common=common[which(common%in%Geom==FALSE)] 
  
  # Read merged file
  Merged_fcs=read.FCS(dir_merged)
  Merged_file=as.data.frame(Merged_fcs@exprs)
 # Merged_file=as.data.frame(lapply(Merged_file,extreme_quant_norm))
  
  
  ################ MEASURE QUALITY #####################
    
    EMDs=integer()
    index=1
    # Common markers
    
    for (mkr in common){
      
      
      EMDs_common=integer()
      for (i in 1:length(files)){
        
        original=files[[i]][,mkr]
        merged=Merged_file[,mkr]
        
        if (length(original)<length(merged)){
          merged=merged[sample(1:length(merged),length(original))]
        }else{
          original=original[sample(1:length(original),length(merged))]
        }
        
        # adapt data for EMDomics package
        values=c(original,merged)
        labels=c(rep("A",length(original)),rep("B",length(merged)))
        names(values) = paste("cell", 1:length(values))
        names(labels)=names(values)
        
        # compute EMD
        EMDs_common[i]=calculate_emd_gene(values, labels, names(values))
        
      }
      
      # Average 
      EMDs[index]=mean(EMDs_common)
      index=index+1
    }
    
    # Imputed markers
    
    imputed=markernames(Merged_fcs)[!markernames(Merged_fcs)%in%common]
    
    for (mkr in imputed){

      # which files contain mkr
      mkr_idx=which(sapply(markers, FUN=function(X) mkr %in% X))
      
      original=numeric()
      for (idx in mkr_idx){
      original=c(original,files[[idx]][,mkr])
      }
      
      merged=Merged_file[,mkr]
      
      if (length(original)<length(merged)){
        merged=merged[sample(1:length(merged),length(original))]
      }else{
        original=original[sample(1:length(original),length(merged))]
      }
      
      # adapt data for EMDomics package
      values=c(original,merged)
      labels=c(rep("Original",length(original)),rep("Merged",length(merged)))
      names(values) = paste("cell", 1:length(values))
      names(labels)=names(values)
      
      # compute EMD
      EMDs[index]=calculate_emd_gene(values, labels, names(values))
      index=index+1
    }
    
    names(EMDs)=c(common,imputed)
    
    quality=as.data.frame(t(c(EMDs,Cells_pre=sum(ncells),Cells_post=nrow(Merged_file))))
    
    return(quality)
    
  
}