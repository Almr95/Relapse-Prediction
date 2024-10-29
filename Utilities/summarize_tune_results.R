summarize_tune_results <- function(object) {
  # INPUT -------
  # object -> rsplit object from rsample package that contains a split from the original dataset
  # OUTPUT ------
  # result -> list with three elements:
  #   - ROC is the AUC of the ROC curve of the prediction of the best model
  #   - Best is the name of the best model selected by the inner loop
  #   - Params is a data frame with the final parameters selected for the best model
  # -----------------------------
  
  library(caret)
  library(rsample)
  library(pROC)
  library(MLmetrics)
  library(precrec)
  source("Utilities/analyse.R")
  source("Utilities/prSummary_mod.R")
  
  # Get analysis and assessment splits of this particular fold (outer loop)
  
  trn=analysis(object)
  tst=assessment(object)
  
  # ANALYSIS ------------------------------------------------------------------
  
  models=analyse(trn)
  
  # Compare resamples
  
  resamps = resamples(list(KNN = models[[1]],
                            NB = models[[2]],
                            RF = models[[3]],
                            SVM=models[[4]]))
  
  # Get best model

  maxAUCs=lapply(models,function(x) x$results$AUC[oneSE(x$results,metric="AUC",num=nrow(x$results),maximize = TRUE)])
  best=which.max(maxAUCs)
  best_model=models[[best]]
  
  # ASSESSMENT ------------------------------------------------------------------
  
  if (best==1 || best==4){
    # KNN and SVM preprocess data. Test data has to be normalized with train parameters!!
    tst[,1:(ncol(tst)-1)]=scale(tst[,1:(ncol(tst)-1)],center=apply(trn[,1:(ncol(trn)-1)],2,mean),scale=apply(trn[,1:(ncol(trn)-1)],2,sd))
  }
  prd=predict(best_model,tst,type="prob")
  
  # precrec
  sscurves=evalmod(scores=prd$R,labels=ifelse(tst$Status == "R", 1, -1))
  prauc=precrec::auc(sscurves)
  auc=subset(prauc,curvetypes=="PRC")$aucs
  
  # Extract results
  
  result=data.frame(AUC=auc,Model=best_model$method)
  result$Params=list(best_model$bestTune)
  result$innerAUC_mean=best_model$results[as.numeric(rownames(best_model$bestTune)),]$AUC
  result$preds=list(data.frame(scores=prd$R,labels=tst$Status))
  result

}