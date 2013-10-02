# the following functions assume that you already have features calculated and that the features are seperate from the target values
ViewTargets <- function(y, bw) {
  plot(density(y, bw = bw), main="Targets")
}

# TBD write function to view a feature by name in the same way as the targets above

ImputeFeatures <- function(x) {
  library(impute)
  as.data.frame(impute.knn(as.matrix(x), k = 10)$data)
}

SplitSet <- function(x, y, percentage = 20, seed = 1) {
  holdout.size <- round(nrow(x) * (percentage/100))
  set.seed(seed)
  holdout.indexes <- sample(1:nrow(x), holdout.size, replace=FALSE)
  train.indexes <- (1:length(y))[-holdout.indexes]
  x.train <- x[train.indexes, ]
  x.holdout <- x[holdout.indexes, ]
  y.train <- y[train.indexes]
  y.holdout <- y[holdout.indexes]
  l <- list()
  l$holdout.indexes <- holdout.indexes
  l$train.indexes <- train.indexes
  l$x.train <- x.train
  l$x.holdout <- x.holdout
  l$y.train <- y.train
  l$y.holdout <- y.holdout
  l
}

RemoveNearZeroVarianceFeatures <- function(ss, frequencyCutoff = 30/1) {
  nzv.columns <- nearZeroVar(ss$x.train, freqCut = frequencyCutoff)
  if (length(nzv.columns) != 0) {
    print(paste(length(nzv.columns), "features removed with variance below cutoff"))
    ss$x.train <- ss$x.train[, -nzv.columns]
    ss$x.holdout <- ss$x.holdout[, -nzv.columns]
  }
  else {
    print("no features removed")
  }
  ss
}

RemoveHighlyCorrelatedFeatures <- function(ss, correlationCutoff = 0.95) {
  hc.columns <- findCorrelation(cor(ss$x.train), correlationCutoff)
  if (length(hc.columns) != 0) {
    print(paste(length(hc.columns), "features removed with correlation above cutoff"))
    ss$x.train <- ss$x.train[, -hc.columns]
    ss$x.holdout <- ss$x.holdout[, -hc.columns]
  }
  else {
    print("no features removed")
  }
  ss
}

PreProcess <- function(ss, steps = c("center", "scale")) {
  transformation <- preProcess(ss$x.train, method = steps)
  ss$x.train <- predict(transformation, ss$x.train)
  ss$x.holdout <- predict(transformation, ss$x.holdout)
  ss$transformation <- transformation
  ss
}

GetCVTrainControl <- function(ss, seed = 1, folds = 5, repeats = 1) {
  set.seed(seed)
  ss$trControl <- trainControl(method='cv', number=folds, repeats=repeats, returnResamp='none',
                               returnData=FALSE, savePredictions=TRUE,
                               verboseIter=TRUE, allowParallel=TRUE,
                               index=createMultiFolds(ss$y.train, k=folds, times=repeats))
  ss
}