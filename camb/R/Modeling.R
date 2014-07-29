#################################################################################
## Modeling 
#################################################################################
asNumeric <- function(x) as.numeric(as.character(x))
factorsNumeric <- function(d) modifyList(d, lapply(d[, sapply(d, is.factor)],   
                                                   asNumeric))

#' Remove infinite values from the descriptor data.frame
#' 
#' Any infinites found in the descriptor data.frame are replaced with NA.
#' 
#' @param d A data.frame
#' @export
#' @return A data.frame with infinite values replaced by NA.
#' @author Daniel Murrell <dsmurrell@@gmail.com> and Isidro Cortes <isidrolauscher@@gmail.com>
ReplaceInfinitesWithNA <- function(d) {
  do.call(data.frame,lapply(d, function(x) replace(x, is.infinite(x), NA)))
}
                         
#' Impute missing descriptor values using knn.impute
#' 
#' A nearest neighbour based 
#' 
#' @param d A data.frame
#' @export
#' @return A data.frame with infinite values replaced by NA.
#' @author Daniel Murrell <dsmurrell@@gmail.com> and Isidro Cortes <isidrolauscher@@gmail.com>
ImputeFeatures <- function(d, k=10,...) {
  suppressWarnings(require(impute)) || stop("Package impute is required. Install from CRAN or Bioconductor -depending on the R version you are using-.")
  as.data.frame(impute.knn(as.matrix(factorsNumeric(d)), k = k, ...)$data)
} 

#' Split into training and hold out data.
#' 
#' Creates a training/holdout split for the data. The training data is used to train models and the holdout data is used to check model 
#' performance on data that is pulled from the same distribution as the training data.
#' 
#' @param ids The names of the molecules
#' @param x The descriptors
#' @param y The target values
#' @param percentage The percentage of data to put into the holdout set
#' @param seed The seed for randomization so that the split can be does in a repeatable way if required
#' @export
#' @return a list designed to be passed between functions in the camb workflow containing the variables required for training: 
#' ids, holdout.indexes, train.indexes, x.train, x.holdout, y.train, y.holdout
#' @author Daniel Murrell <dsmurrell@@gmail.com> and Isidro Cortes <isidrolauscher@@gmail.com>
SplitSet <- function(ids, x, y, percentage = 20, seed = 1) {
  holdout.size <- round(nrow(x) * (percentage/100))
  set.seed(seed)
  holdout.indexes <- sample(1:nrow(x), holdout.size, replace=FALSE)
  train.indexes <- (1:length(y))[-holdout.indexes]
  x.train <- x[train.indexes, ]
  x.holdout <- x[holdout.indexes, ]
  y.train <- y[train.indexes]
  y.holdout <- y[holdout.indexes]
  l <- list()
  l$ids <- ids
  l$holdout.indexes <- holdout.indexes
  l$train.indexes <- train.indexes
  l$x.train <- x.train
  l$x.holdout <- x.holdout
  l$y.train <- y.train
  l$y.holdout <- y.holdout
  l
}

#' Calls the \code{\link{nearZeroVar}} function of the caret package and removes descriptos with near zero variance in the
#' appropriate columns from the training and holdout sets.
#' 
#' \code{\link{nearZeroVar}} diagnoses predictors that have one unique value (i.e. are zero variance predictors) or predictors 
#' that are have both of the following characteristics: they have very few unique values relative to the number 
#' of samples and the ratio of the frequency of the most common value to the frequency of the second most common 
#' value is large. checkConditionalX looks at the distribution of the columns of x conditioned on the levels of 
#' y and identifies columns of x that are sparse within groups of y.
#' 
#' @param dataset The training dataset returned by \code{\link{SplitSet}}
#' @param frequencyCutoff The cutoff for the ratio of the most common value to the second most common value
#' @export
#' @return A dataset with the appropriate columns cut out of the x.train and x.holdout dataframes
#' @author Daniel Murrell <dsmurrell@@gmail.com> and Isidro Cortes <isidrolauscher@@gmail.com>
RemoveNearZeroVarianceFeatures <- function(dataset, frequencyCutoff = 30/1,...) {
  nzv.columns <- nearZeroVar(dataset$x.train, freqCut = frequencyCutoff,...)
  if (length(nzv.columns) != 0) {
    message(paste(length(nzv.columns), "features removed with variance below cutoff"))
    dataset$x.train <- dataset$x.train[, -nzv.columns]
    dataset$x.holdout <- dataset$x.holdout[, -nzv.columns]
  }
  else {
    message("no features removed")
  }
  dataset
}

#' Calls the \code{\link{findCorrelation}} function of the caret package which finds the highly correlated descriptors
#' and removes them from the training and holdout sets.
#' 
#' The absolute values of pair-wise correlations are considered. If two variables have a high correlation, 
#' the function looks at the mean absolute correlation of each variable and removes the variable with the 
#' largest mean absolute correlation.
#' 
#' @param dataset The training dataset returned by \code{\link{SplitSet}}
#' @param correlationCutoff A numeric value for the pair-wise absolute correlation cutoff
#' @export
#' @return A dataset with the appropriate columns cut out of the x.train and x.holdout dataframes
#' @author Daniel Murrell <dsmurrell@@gmail.com> and Isidro Cortes <isidrolauscher@@gmail.com>
RemoveHighlyCorrelatedFeatures <- function(dataset, correlationCutoff = 0.95,...) {
  hc.columns <- findCorrelation(cor(dataset$x.train), correlationCutoff,...)
  if (length(hc.columns) != 0) {
    message(paste(length(hc.columns), "features removed with correlation above cutoff"))
    dataset$x.train <- dataset$x.train[, -hc.columns]
    dataset$x.holdout <- dataset$x.holdout[, -hc.columns]
  }
  else {
    message("no features removed")
  }
  dataset
}

#' Cals the \code{\link{preProcess}} function of the caret package which handles data transformation before training.
#' 
#' Pre-processing transformation (centering, scaling etc.) can be estimated from the training data and applied to 
#' any data set with the same variables. The simplest form of this is to center and scale all the variables so that
#' each of their means are 0 and each of their standard devitations is 1. See \code{\link{preProcess}} arguments for 
#' more control.
#' 
#' @param dataset The training dataset returned by \code{\link{SplitSet}}
#' @param steps a character vector specifying the type of processing. Possible values are "BoxCox", "YeoJohnson", 
#' "expoTrans", "center", "scale", "range", "knnImpute", "bagImpute", "medianImpute", pca", "ica" and "spatialSign" 
#' for details see the \code{\link{preProcess}} function for more details.
#' 
#' @export
#' @return A dataset with the preprocessed training and holdout set that also stores the transformation
#' @author Daniel Murrell <dsmurrell@@gmail.com> and Isidro Cortes <isidrolauscher@@gmail.com>
PreProcess <- function(dataset, steps = c("center", "scale"),...) {
  transformation <- preProcess(dataset$x.train, method = steps,...)
  dataset$x.train <- predict(transformation, dataset$x.train)
  dataset$x.holdout <- predict(transformation, dataset$x.holdout)
  dataset$transformation <- transformation
  dataset
}

#' Sets up the control parameters of caret's \code{\link{train}} function.
#' 
#' Calls caret's \code{\link{trainControl}} function to set up the parameters of the \code{\link{train}} function. This control variable is saved
#' into the \code{dataset} list as \code{dataset$trControl}.
#' 
#' @param dataset The training dataset returned by \code{\link{SplitSet}}
#' @param seed The seed for randomization so that the fold selection can be done in a repeatable way if desired
#' @param folds The number of folds to use during cross-validation
#' @param repeats For repeated k-fold cross-validation only: the number of complete sets of folds to compute
#' @param returnResamp A character string indicating how much of the resampled summary metrics should be saved. 
#' Values can be “final”, “all” or “none”
#' @param returnData A logical for saving the data
#' @param savePredictions A logical to save the hold-out predictions for each resample
#' @param verboseIter A logical for printing a training log.
#' @param allowParallel If a parallel backend is loaded and available, should the function use it?
#' 
#' @export
#' @return a dataset with the traincontrol saved within for future training
#' @author Daniel Murrell <dsmurrell@@gmail.com> and Isidro Cortes <isidrolauscher@@gmail.com>
GetCVTrainControl <- function(dataset, seed = 1, folds = 5, repeats = 1, method='cv', returnResamp='none', returnData=TRUE, savePredictions=TRUE, verboseIter=TRUE, allowParallel=TRUE,...) {
  set.seed(seed)
  dataset$trControl <- trainControl(method='cv', number=folds, repeats=repeats, returnResamp='none',
                                    returnData=TRUE, savePredictions=TRUE,
                                    verboseIter=TRUE, allowParallel=TRUE,
                                    index=createMultiFolds(dataset$y.train, k=folds, times=repeats),...)
  dataset
}

#' Creates an exponential grid
#' 
#' @param power.from The first power to use
#'
#' @export
#' @author Daniel Murrell <dsmurrell@@gmail.com> and Isidro Cortes <isidrolauscher@@gmail.com>
expGrid <- function(power.from, power.to, power.by, base){
  grid <- c()
  for (i in seq(power.from, power.to, power.by)){
    grid <- append(grid,base^i)
  }
  return(grid)
}

#' @export
#' @author Daniel Murrell <dsmurrell@@gmail.com> and Isidro Cortes <isidrolauscher@@gmail.com>
YScrambling <- function(y,percent){
if (percent < 0 || percent > 1){stop("The percent value needs to be between 0 and 1")}
inds_to_resamp <- sample.int(length(y), length(y)*percent)
resamp_vals <- y[inds_to_resamp]
resamp_vals2 <- sample(resamp_vals,length(resamp_vals),replace = FALSE, prob = NULL)
y[inds_to_resamp] <- resamp_vals2
return(y)
}


