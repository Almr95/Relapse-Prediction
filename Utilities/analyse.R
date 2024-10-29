analyse = function(train) {
  
  source("Utilities/prSummary_mod.R")
  
  models=list()
  cvIndex = createMultiFolds(factor(train$Status), k=9,times=20)
  
  # Specify inner loop
  train_control = trainControl(method="repeatedcv",number=9, repeats = 20, returnResamp = "final",selectionFunction = "oneSE",summaryFunction =prSummary_mod, savePredictions = "all",classProbs = TRUE,allowParallel = FALSE,index = cvIndex)
  
  # Fit models
  
  seed=round(runif(1)*1000)
  
  # K-nearest-neighbors
  set.seed(seed)
  models[[1]] = train(Status~., data=train, trControl=train_control, method="knn",metric = "AUC",preProcess = c("center", "scale"), tuneGrid = expand.grid(k = round(seq(round(0.2*length(which(train$Status=="R"))), round(0.8*length(which(train$Status=="R"))), length.out=6))))
 
  # Naive Bayes
  set.seed(seed)
  models[[2]] = train(Status~., data=train, trControl=train_control, method="naive_bayes",metric = "AUC", tuneGrid = expand.grid(laplace = 1,usekernel=c(TRUE,FALSE),adjust=c(0.01,0.1,1)))
  
  # Random Forest
  set.seed(seed)
  models[[3]] = train(Status~., data=train, trControl=train_control, method="rf",metric = "AUC",tuneGrid=expand.grid(mtry=round(seq(round(0.2*ncol(train)),round(0.8*ncol(train)),length.out=6))),ntree=200)
  
  # Support Vector Machine
  set.seed(seed)
  models[[4]] = train(Status~.,data=train,trControl=train_control,method="svmLinearWeights",metric="AUC",preProcess = c("center", "scale"),tuneGrid = expand.grid(cost=c(1, 10),weight=c(0.01,0.1,1)))
  
  models
  
}