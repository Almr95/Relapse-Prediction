classification <- function(dir) {
  
  # Packages ---------------------
  
  library(readxl)
  library(dplyr)
  library(writexl)
  library(Hmisc)
  library(caret)
  library(rsample)
  library(doParallel)
  library(foreach)
  library(entropy)
  library(ggpubr)
  library(precrec)
  
  source("summarize_tune_results.R")
  source("analyse.R")
  source("prSummary_mod.R")
  
  # ------------------------------
  
  # Read data
  file=list.files(dir,pattern="Values")
  features=read_excel(paste0(dir,"/",file))
  features$Status=factor(features$Status,levels=c("R","NR"))
  
  if (nrow(features)>40){
    
  # Create outer loop
  outer <- vfold_cv(features,v=5,repeats=10,strata = "Status")

  # Parallel
  numCores = detectCores()-1
  cl=makeCluster(numCores)
  registerDoParallel(cl)
  
  # Outer fold routine
  finalres=data.frame()
  finalres <- foreach(i=1:nrow(outer), .combine=rbind) %dopar% {
  
    # Run inner loop
    source("summarize_tune_results.R")
    res=suppressWarnings(summarize_tune_results(outer$splits[[i]]))
    res
    
  }
  
  # Label models
  finalres$Model=factor(finalres$Model,levels = c("knn","naive_bayes","rf","svmLinearWeights"))
  levels(finalres$Model)=c("KNN","NB","RF","SVM")
  
  # Compute stability
  Stability=max(table(finalres$Model))/sum(table(finalres$Model))
  Stability_min=max(diff(round(seq(0, sum(table(finalres$Model)), length.out=1+length(table(finalres$Model))))))/sum(table(finalres$Model))
  Stability_norm=(Stability-Stability_min)/(1-Stability_min)
  
  # Summarize outer loop results (stability check) and compare with inner loop (overfitting check)
  png(paste0(dir,"/Outer.png"),width=30,height=30,unit="cm",res=500)
  a=data.frame(Summary=summary(finalres$Model),Model=c("KNN","NB","RF","SVM"))
  p1=ggplot(a,aes(Model,Summary))+geom_col(color="black")+theme_bw()+theme(text = element_text(size = 20))+ylab("")+xlab("Model")+
             ggtitle(paste0("Stability index = ",as.character(Stability_norm,digits=3)))
  p2=ggplot(finalres,aes(AUC))+geom_density(linewidth=1)+geom_density(data=finalres,aes(innerAUC_mean),linetype="dashed",linewidth=1)+theme_bw()+theme(text = element_text(size = 20))+ylab("")+xlab("AUC")+
    ggtitle(paste0("── Outer (Mean ",as.character(format(mean(finalres$AUC),digits=3)),")","\n----- Inner (Mean ",as.character(format(mean(finalres$innerAUC_mean),digits=3)),")"))
  p3=ggplot(finalres, aes(x=1:nrow(finalres), y=AUC)) +
    geom_segment( aes(x=1:nrow(finalres), xend=1:nrow(finalres), y=0, yend=AUC)) +
    geom_point( size=2, color="red", fill=alpha("red", 1), alpha=1, shape=21, stroke=2) +
    geom_segment(data=finalres, aes(x=1:nrow(finalres), xend=1:nrow(finalres), y=0, yend=innerAUC_mean),linetype="dashed") +
    geom_point(data=finalres,aes(x=1:nrow(finalres), y=innerAUC_mean), size=2, color="red", fill=alpha("red", 1), alpha=1, shape=1, stroke=2)+
    xlab("Outer CV folds")+ylab("AUC")+theme_bw()+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),text = element_text(size = 20))+ 
    annotate("text", x = 1:nrow(finalres), y = 1.15, label = finalres$Model, size = 3,angle=90) +
    coord_cartesian(ylim = c(0, 1.2),clip = "off") + geom_hline(yintercept=mean(finalres$AUC), color = "red")+ geom_hline(yintercept=mean(finalres$innerAUC_mean),linetype="dashed", color = "red")
  print(ggarrange(ggarrange(p1,p2,ncol=2),p3,nrow=2))
  dev.off()
  
  # Plot Precision-Recall curve
  scores=lapply(finalres$preds, `[[`, 1)
  labels=lapply(finalres$preds, `[[`, 2)
  labels=lapply(labels, function(x) ifelse(x=="R",1,-1))
  smmdat <- mmdata(scores, labels, dsids = seq(1, length(scores)))
  smcurves <- evalmod(smmdat, raw_curves = TRUE)
  smcurves.df <- as.data.frame(smcurves)
  smcurves.df=smcurves.df[smcurves.df$type=="PRC",]
  p=autoplot(smcurves,"PRC")
  smcurves.df_avg=p$data
    
  png(paste0(dir,"/PRcurve.png"),width=15,height=15,unit="cm",res=300)
  print(ggplot(smcurves.df_avg,aes(x,y,group=modname))+annotate("segment", x = 0, xend = 1, y = 0.202, yend = 0.202,size=2,linetype="dashed")+
    geom_ribbon(aes(ymin=ymin,ymax=ymax),fill = "gray",alpha=0.8)+geom_line(linewidth=3)+
    scale_x_continuous(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1))+
    scale_y_continuous(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1))+
    theme_bw()+xlab("Recall")+ylab("Precision")+
    theme(axis.title = element_text(size=30),
          axis.text = element_text(size=25),
          legend.position = c(0.8, 0.8)))
  dev.off()
  
  
  
  
  # Refit
  models=analyse(features)
  
  # Compare resamples
  resamps <- resamples(list(KNN = models[[1]],
                            NB = models[[2]],
                            RF = models[[3]],
                            SVM = models[[4]]))
  
  png(paste0(dir,"/Resamples.png"),width=25,height=20,unit="cm",res=500)
  print(bwplot(resamps, layout = c(3, 1)))
  dev.off()
  
  # Get best model
  maxAUCs=lapply(models,function(x) x$results$AUC[oneSE(x$results,metric="AUC",num=nrow(x$results),maximize = TRUE)])
  best=which.max(maxAUCs)
  best_model=models[[best]]
  
  png(paste0(dir,"/Best.png"),width=25,height=20,unit="cm",res=500)
  print(plot(best_model)) 
  dev.off()
  
  write_xlsx(best_model$results,path=paste0(dir,"/Best.xlsx"))
  
  # Stop parallel
  stopCluster(cl)
  
  # Output
  classification=list(finalres,models)
  classification
  
  }else{
    
    classification=list()
    classification
    
  }
}