prSummary_mod <- function (data, lev = NULL, model = NULL) 
{

  # Modified from base caret to employ precrec package to compute AUCPR
  
  require("precrec")
  if (length(levels(data$obs)) > 2) 
    stop(paste("Your outcome has", length(levels(data$obs)), 
               "levels. `prSummary2`` function isn't appropriate.", 
               call. = FALSE))
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
    stop("Levels of observed and predicted data do not match.", 
         call. = FALSE)
  if (!lev[1] %in% colnames(data)) 
    stop(paste("Class probabilities are needed to score models using the", 
               "area under the PR curve. Set `classProbs = TRUE`", 
               "in the trainControl() function."), call. = FALSE)
  
  sscurves=precrec::evalmod(scores=data[,lev[1]],labels=ifelse(data$obs == lev[1], 1, -1))
  prauc=precrec::auc(sscurves)

  out=c(subset(prauc,curvetypes=="PRC")$aucs)
  names(out) <- c("AUC")
  return(out)
                                                                                                                                                                             
}