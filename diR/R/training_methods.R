# trains a model with some preprocessing, stores all the data in a list (model), uses 1/6 of data as holdout set
trainModel <- function(x, 
                       y, 
                       method, 
                       tune.grid, 
                       cores = 1, 
                       fold.index = 1,
                       partition.seed = 1,
                       training.seed = 2,
                       freqCut = 30/1,
                       correlationCut = 0.95,
                       ... ) {
  ### data partitioning
  set.seed(partition.seed)
  folds <- createFolds(y, k=6)
  in.test <- folds[[fold.index]]
  in.train <- (1:length(y))[-in.test]
  
  x.train <- x[in.train,] 
  x.test <- x[-in.train,]
  
  y.train <- y[in.train]
  y.test <- y[-in.train]
  
  # remove near zero variance descriptors
  nzv.columns <- nearZeroVar(x.train, freqCut = freqCut)
  if(length(nzv.columns) != 0) {
    nzv.names <- names(x.train)[nzv.columns]
    x.train <- x.train[, -nzv.columns]
    x.test <- x.test[, -nzv.columns]
  }
  
  # remove highly correlated descriptors
  is.highly.correlated <- findCorrelation(cor(x.train), correlationCut)
  if(length(is.highly.correlated) != 0) {
    x.train <- x.train[,-is.highly.correlated]
    x.test <- x.test[,-is.highly.correlated]
  }
  
  # apply a centering and scaling transformation to the data
  transformation <- preProcess(x.train, method = c("center", "scale"))
  x.train <- predict(transformation, x.train)
  x.test <- predict(transformation, x.test)

  # build a model and return it
  registerDoMC(cores)
  set.seed(training.seed)
  trControl <- trainControl(method = "cv", number = 5, returnData = FALSE, 
                            savePredictions=TRUE, verboseIter=TRUE, 
                            allowParallel=TRUE,
                            index=createMultiFolds(y.train, k=5, times=1))
  print(method)
  print(cores)
  print(tune.grid)
  fit <- train(x.train, y.train, method = method, trControl = trControl, tuneGrid = tune.grid)

  model <- list()
  model$fit <- fit
  model$transformation <- transformation
  model$x.train <- x.train
  model$y.train <- y.train
  model$x.test <- x.test
  model$y.test <- y.test
  model$in.test <- in.test
  model$in.train <- in.train
  model
}

# trains a model with some preprocessing, stores all the data in a list (model), uses all the data for training
trainModelOnAll <- function(x, 
                          y, 
                          method, 
                          tune.grid, 
                          cores = 1, 
                          training.seed = 2,
                          freqCut = 30/1,
                          correlationCut = 0.95,
                          ... ) {
  ### data partitioning 
  x.train <- x
  y.train <- y
  
  # remove near zero variance descriptors
  nzv.columns <- nearZeroVar(x.train, freqCut = freqCut)
  if(length(nzv.columns) != 0) {
    nzv.names <- names(x.train)[nzv.columns]
    x.train <- x.train[, -nzv.columns]
  }
  
  # remove highly correlated descriptors
  is.highly.correlated <- findCorrelation(cor(x.train), correlationCut)
  if(length(is.highly.correlated) != 0) {
    x.train <- x.train[,-is.highly.correlated]
  }
  
  # apply a centering and scaling transformation to the data
  transformation <- preProcess(x.train, method = c("center", "scale"))
  x.train <- predict(transformation, x.train)
  
  # build a model and return it
  registerDoMC(cores)
  trControl <- trainControl(method = "cv", number = 2) 
  set.seed(training.seed)
  print(method)
  print(cores)
  print(tune.grid)
  fit <- train(x.train, y.train, method = method, trControl = trControl, tuneGrid = tune.grid, ...)

  model <- list()
  model$fit <- fit
  model$transformation <- transformation
  model$x.train <- x.train
  model$y.train <- y.train
  model
}

