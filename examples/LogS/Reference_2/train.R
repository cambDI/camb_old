library(caret)
library(kernlab)
library(doMC)

expGrid <- function(power.from, power.to, power.by, base){
  grid <- c()
  for (i in seq(power.from, power.to, power.by)){
    grid <- append(grid,base^i)
  }
  return(grid)
}

registerDoMC(cores=10)
dataset <- readRDS("dataset.rds")
method <- "svmRadial"
tune.grid <- expand.grid(.sigma = expGrid(-8, 4, 2, 2), .C = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100))
model <- train(dataset$x.train, dataset$y.train, method, tuneGrid=tune.grid, trControl=dataset$trControl)
saveRDS(model, file=paste(method,".rds",sep=""))