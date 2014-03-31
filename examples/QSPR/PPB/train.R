ss <- readRDS("ss.rda")
library(caret)
library(doMC)
library(gbm)
registerDoMC(cores=12)

method <- "gbm"
#tune.grid <- expand.grid(.n.trees=500, .interaction.depth=c(25), .shrinkage = c(0.0025, 0.005, 0.01, 0.02))

#method <- "rf"
#tune.grid <- expand.grid(.mtry = seq(5,100,5))

#method <- "svmRadial"
#tune.grid <- expand.grid(.sigma = c(0.001, 0.002, 0.004, 0.008, 0.012), .C = c(0.5, 1, 2, 3, 4))

#method <- "gaussprRadial"
#tune.grid <- expand.grid(.sigma = c(0.001))

#method <- "rvmRadial"
#tune.grid <- expand.grid(.sigma = c(0.001, 0.0015, 0.002))

model <- train(ss$x.train, ss$y.train, method, tuneGrid=tune.grid, trControl=ss$trControl)
saveRDS(model, file=paste(method,".rds",sep=""))

# show some plots in a pdf file
pdf(paste(method,"_results.pdf", sep=""))
CVRMSE <- signif(min(as.vector(na.omit(model$results$RMSE))), digits=3)
holdout.predictions <- as.vector(predict(model, newdata = ss$x.holdout))
holdoutRMSE <- signif(RMSE(holdout.predictions, ss$y.holdout), digits=3)
plot(model, metric = "RMSE", main=paste("CV RMSE:", CVRMSE, "    Holdout RMSE", holdoutRMSE))
plot(ss$y.holdout, holdout.predictions)
plot(density(holdout.predictions, bw = 0.5))
dev.off()

# save out the ids, smiles, predictions and observed values so that interactive software may be used
origin <- read.csv("from_chris.csv")
# do the holdout set results
is <- ss$holdout.indexes
ids <- ss$ids[is]
smiles <- origin$original_format[match(ids, origin$ChemSpiderID)]
obs <- ss$y.holdout
preds <- as.vector(predict(model, newdata = ss$x.holdout)) 
odf <- data.frame(ID = ids, SMILES = smiles, OBSERVED = obs, PREDICTED = preds)
write.csv(odf, file=paste(method,"_holdout.csv", sep=""), row.names = F)

# do the CV results
bt <- model$bestTune
RowMatch <- function(x, y) {
  identical(as.numeric(x), as.numeric(y))
}
btindexes <- which(apply(as.matrix(model$pred[,colnames(bt)]), 1, RowMatch, as.numeric(bt)))
cv <- model$pred[btindexes, ]
ids <- ss$ids[ss$train.indexes[cv$rowIndex]]
smiles <- origin$original_format[match(ids, origin$ChemSpiderID)]
odf <- data.frame(ID = ids, SMILES = smiles, OBSERVED = cv$obs, PREDICTED = cv$pred, fold=cv$Resample)
write.csv(odf, file=paste(method,"_cv.csv", sep=""), row.names = F)