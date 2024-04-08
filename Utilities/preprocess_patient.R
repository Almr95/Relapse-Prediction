preprocess_patient <- function(patient) {
  
  ########### LOAD LIBRARIES AND AUXILIARY PROGRAMMES #############
  
  library(flowCore) 
  library(data.table)
  library(plyr)
  library(aroma.light)
  
  # code
  source("Utilities/rename.R")
  source("Utilities/extreme_quant_norm.R")
  
  ########## IMPORT AND PREPROCESS ##################
  
  filenames=as.list(list.files(patient,pattern = ".fcs"))
  
  # initialize structures
  files=list()
  remove=integer()
  ncells=integer()
  markers=list()
  
  # auxiliary variables 
  DEL=c("DEL1","DEL2","DEL3","DEL4","DEL5","DEL6","DEL7","DEL8","DEL9","DEL10") # markers to delete
  Geom=c("FSC.H","FSC.A","SSC.A","SSC.H") # name forward and side scatter parameters
  BB=c("CD19","CD45") # backbone markers
  
  # Logicle transform
  lgcl=logicleTransform(w=0.75,t=262144,m=4.5,a=0) # parameters
  
  for (i in 1:length(filenames)){
    
    # Read FCS
    data=read.FCS(paste0(patient,"/",filenames[[i]]),truncate_max_range = FALSE,column.pattern = "Time",invert.pattern = TRUE) # Read
    
    # Rename
    data=renamed(data) 
    names=unname(data@parameters@data$name)
    
    # Remove margins
    maxRange=data@parameters@data$maxRange[1]
    data=as.data.frame(data@exprs)
    data[data>=0.98*maxRange]=NA
    data=na.omit(data)
    data=flowFrame(as.matrix(data))
    
    # Transform
    trans=names[which(names%in%Geom==FALSE)]
    translgcl=transformList(trans, lgcl)
    data=transform(data, translgcl)
    
    # Store info
    data=as.data.frame(data@exprs)
    files[[i]]=data[,!(names(data) %in% DEL)]
    markers[[i]]=colnames(files[[i]])
    ncells[i]=nrow(files[[i]])
    
    # Normalize (modified min-max)
    files[[i]]=as.data.frame(lapply(files[[i]],extreme_quant_norm))
    
    # Check if it contains backbone markers
    if (all(BB%in%markers[[i]])==FALSE){
      remove=c(remove,i)
    }
    
  }
  
  ########## REMOVE TUBES LACKING BACKBONE MARKERS #########
  
  if (length(remove)!=0){
    filenames[remove]=NULL
    files[remove]=NULL
    ncells=ncells[-remove]
    markers[remove]=NULL
  }
  
  ########## SAMPLE CELLS #################
  
  # Remove tubes with too little cells
  mincells=min(ncells)
  while (mincells*length(files)<10000 && sum(ncells[-which.min(ncells)])>10000 ){
    files=files[-which.min(ncells)]
    filenames=filenames[-which.min(ncells)]
    ncells=as.integer(lapply(files, nrow))
    markers=lapply(files,colnames)
    mincells=min(ncells)
  }
  
  # Sample
  if (mincells>10000){mincells=10000}
  files <- lapply(files, function(x){x[sample(nrow(x), mincells),]})
  
  
  ######### EXPORT ###################
  
  files_fcs = lapply(files, function(x) flowFrame(as.matrix(x)))
  names(files_fcs)=filenames
  return(files_fcs)
  
}