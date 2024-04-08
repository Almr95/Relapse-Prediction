library(FNN)

# Obtains non-ambiguous nearest-neighbors. Adapted from CytoBackBone. 

get.nann <- function(df1, df2,th) {
  
  # Create auxiliary dataframe
  
    df1_iter=df1
    df2_iter=df2
    
  # Get original indices
    
    indices_1=1:nrow(df1)
    indices_2=1:nrow(df2)
    
  # Obtain distances and subset
    
    dist = FNN::knnx.dist(df2_iter,df1_iter,k=1,algorithm="kd_tree")
    df1=df1[dist<th,]
    indices_1=indices_1[dist<th]
    
    dist = FNN::knnx.dist(df1_iter,df2_iter,k=1,algorithm="kd_tree")
    df2=df2[dist<th,]
    indices_2=indices_2[dist<th]
    
    max <- min(nrow(df1),nrow(df2))
  
  
  # Prepare data for first iteration
    
    idxA=rep(NA,max)
    idxB=rep(NA,max)
    
    df1_iter=df1
    df2_iter=df2
    
  # Get respective neighbours
    
    knnx                          <- FNN::get.knnx(df2_iter,df1_iter,k=1,algorithm="kd_tree")
    idx_1                         <- knnx$nn.index
    dist_1                          <- knnx$nn.dist
    
    knnx                          <- FNN::get.knnx(df1_iter,df2_iter,k=1,algorithm="kd_tree")
    idx_2                           <- knnx$nn.index
    dist_2                          <- knnx$nn.dist
  
    idx                           <- paste0(1:nrow(idx_1),"-",idx_1) %in% paste0(idx_2,"-",1:nrow(idx_2))
    
  # Obtain indices
    
    idx_2=idx_1[idx]
    idx_1=(1:nrow(idx_1))[idx]
    
  # Obtain original indices
    
    idxA[1:length(idx_1)]=indices_1[idx_1]
    idxB[1:length(idx_2)]=indices_2[idx_2]
    sig=length(idx_1)
  

    
    if(max-sum(!is.na(idxA))!=0){	

      while(TRUE){
        
        # Adapt data
        
        indices_1=indices_1[-idx_1]
        indices_2=indices_2[-idx_2]
        
        df1=df1[-idx_1,]
        df2=df2[-idx_2,]
        
        df1_iter=df1
        df2_iter=df2
        
        # Run neighbour search
        
        knnx                          <- FNN::get.knnx(df2_iter,df1_iter,k=1,algorithm="kd_tree")
        idx_1                         <- knnx$nn.index
        dist_1                          <- knnx$nn.dist
        idx_1[dist_1>th]=-1
        
        knnx                          <- FNN::get.knnx(df1_iter,df2_iter,k=1,algorithm="kd_tree")
        idx_2                           <- knnx$nn.index
        dist_2                          <- knnx$nn.dist
        idx_2[dist_2>th]=-2
        
        idx                           <- paste0(1:nrow(idx_1),"-",idx_1) %in% paste0(idx_2,"-",1:nrow(idx_2))
        
        idx_2=idx_1[idx]
        idx_1=(1:nrow(idx_1))[idx]
        
        if(sum(idx)==0){break}
        if(sum(idx)==1){
          idxA[1:length(idx_1)+sig]=indices_1[idx_1]
          idxB[1:length(idx_2)+sig]=indices_2[idx_2]
          sig=sig+length(idx_1)
        }else{
          idxA[1:length(idx_1)+sig]=indices_1[idx_1]
          idxB[1:length(idx_2)+sig]=indices_2[idx_2]
          sig=sig+length(idx_1)
        }
        
        
        if(max-sum(!is.na(idxA))==0){break}
          
      }
      
    }
    
    return(list(na.omit(idxA),na.omit(idxB)))
  
}