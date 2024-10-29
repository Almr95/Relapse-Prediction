extract_features <- function(fs,metadata,dir){
  
  ########### LOAD PACKAGES #############
  
  library(moments)
  library(reshape)
  library(dplyr)
  library(scattermore)
  library(pheatmap)
  library(data.table)
  library(corrplot)
  library(Hmisc)
  library(ggplot2)
  
  
  ############## Order 1-4 moments #############
  # Mean, Standard Deviation, Skewness, Kurtosis
  
  Feat_list=list()
  KSp_list=list()
  
  for (i in 1:4){
    
    dir.create(paste0(dir,"/Order_",i), showWarnings = FALSE)
    Mom=as.data.frame(fsApply(fs,each_col,moment,order=i))
    
    Feat_list[[i]]=Mom
    colnames(Feat_list[[i]])=paste0(colnames(Feat_list[[i]])," Order ",i)

    # Compute KStest pvals and print barplot
    Kstest_p=apply(Mom, 2, function(x) ks.test(x[metadata$Status=="R"],x[metadata$Status=="NR"])$p)
    KSp_list[[i]]=Kstest_p
    
    pdf(paste0(dir,"/Order_",i,"/KStest_pvals.pdf"),width=10,height=5)
    print(ggplot(data.frame(Pval=Kstest_p,Marker=colnames(Mom)),aes(Marker,Pval))+geom_col(color="black")+geom_hline(yintercept=0.05, linetype="dashed", color = "red")+theme_bw()+theme(text = element_text(size = 20)))
    dev.off()
    
    # Print boxplot
    Mom$Status=metadata$Status
    meltData=melt(data.table(Mom),id="Status")
    
    pdf(paste0(dir,"/Order_",i,"/Boxplot.pdf"),width=10,height=5)
    print(ggplot(meltData, aes(x=factor(variable), y=value,fill=Status))+geom_boxplot()+ylab(paste0("Order ",i," moment"))+scale_fill_manual(values=c("#00BFC4","#F8766D"))+
          scale_x_discrete(guide = guide_axis(angle = 90)))
    dev.off() 
    
    # Compute correlations and print 
    Mom$Status=as.numeric(as.factor(Mom$Status))-1
    correlations=rcorr(as.matrix(Mom))
    correlations$P[is.na(correlations$P)]=0 # Replace NA with 0
    
    pdf(paste0(dir,"/Order_",i,"/Corr.pdf"),width=7,height=7)
    corrplot(correlations$r, p.mat = correlations$P,type="lower",insig="label_sig")
    dev.off()
    
    # Save values
    Mom$Status=metadata$Status
    write_xlsx(Mom,paste0(dir,"/Order_",i,"/Values.xlsx"))
    
    
    

  }

  ############## Median (MFI) #############
  
  
  i=i+1
  
  
  dir.create(paste0(dir,"/Median"), showWarnings = FALSE)
  Mom=as.data.frame(fsApply(fs,each_col,median))
  
  Feat_list[[i]]=Mom
  colnames(Feat_list[[i]])=paste0(colnames(Feat_list[[i]])," Median")
  
  # Compute KStest pvals and print barplot
  Kstest_p=apply(Mom, 2, function(x) ks.test(x[metadata$Status=="R"],x[metadata$Status=="NR"])$p)
  KSp_list[[i]]=Kstest_p
  
  pdf(paste0(dir,"/Median","/KStest_pvals.pdf"),width=10,height=5)
  print(ggplot(data.frame(Pval=Kstest_p,Marker=colnames(Mom)),aes(Marker,Pval))+geom_col(color="black")+geom_hline(yintercept=0.05, linetype="dashed", color = "red")+theme_bw()+theme(text = element_text(size = 20)))
  dev.off()
  
  # Print boxplot
  Mom$Status=metadata$Status
  meltData=melt(data.table(Mom),id="Status")
  
  pdf(paste0(dir,"/Median","/Boxplot.pdf"),width=10,height=5)
  print(ggplot(meltData, aes(x=factor(variable), y=value,fill=Status))+geom_boxplot()+ylab(paste0("Median"))+scale_fill_manual(values=c("#00BFC4","#F8766D"))+
          scale_x_discrete(guide = guide_axis(angle = 90)))
  dev.off() 
  
  # Compute correlations and print 
  Mom$Status=as.numeric(as.factor(Mom$Status))-1
  correlations=rcorr(as.matrix(Mom))
  correlations$P[is.na(correlations$P)]=0 # Replace NA with 0
  
  pdf(paste0(dir,"/Median","/Corr.pdf"),width=7,height=7)
  corrplot(correlations$r, p.mat = correlations$P,type="lower",insig="label_sig")
  dev.off()
  
  # Save values
  Mom$Status=metadata$Status
  write_xlsx(Mom,paste0(dir,"/Median","/Values.xlsx"))
  
  ############## All metrics together #############
  
  dir.create(paste0(dir,"/All"), showWarnings = FALSE)
  Mom = bind_cols(Feat_list)
  
  # Compute KStest pvals and print barplot
  Kstest_p=apply(Mom, 2, function(x) ks.test(x[metadata$Status=="R"],x[metadata$Status=="NR"])$p)
  
  pdf(paste0(dir,"/All","/KStest_pvals.pdf"),width=20,height=5)
  print(ggplot(data.frame(Pval=Kstest_p,Marker=colnames(Mom)),aes(Marker,Pval))+geom_col(color="black")+geom_hline(yintercept=0.05, linetype="dashed", color = "red")+theme_bw()+theme(text = element_text(size = 20)))
  dev.off()
  
  # Print boxplot
  Mom$Status=metadata$Status
  meltData=melt(data.table(Mom),id="Status")
  
  pdf(paste0(dir,"/All","/Boxplot.pdf"),width=20,height=5)
  print(ggplot(meltData, aes(x=factor(variable), y=value,fill=Status))+geom_boxplot()+scale_fill_manual(values=c("#00BFC4","#F8766D")))
  dev.off() 
  
  # Compute correlations and print 
  Mom$Status=as.numeric(as.factor(Mom$Status))-1
  correlations=rcorr(as.matrix(Mom))
  correlations$P[is.na(correlations$P)]=0 # Replace NA with 0
  
  pdf(paste0(dir,"/All","/Corr.pdf"),width=15,height=15)
  corrplot(correlations$r, p.mat = correlations$P,type="lower",insig="label_sig")
  dev.off()
  
  # Save values
  Mom$Status=metadata$Status
  write_xlsx(Mom,paste0(dir,"/All","/Values.xlsx"))
  
  ############## Patient IDs #############
  
  # Save patient IDs to check consistency between different clusterings
  
  write_xlsx(metadata,paste0(dir,"/IDs.xlsx"))
  

}