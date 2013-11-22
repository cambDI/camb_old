''' Cambridge Workshop '''
''' November 2013 '''
''' Daniel Murrell and Isidro Cortes '''
''' QSPR example with camb '''

library(camb)
library(ggplot2)
library(doMC)

#########################################
# Standardise molecules
#########################################

StandardiseMolecules(structures.file="solubility_2007_ref2.sdf", 
                     standardised.file="standardised.sdf", 
                     removed.file="removed.sdf",
                     target.field.name = "EXPT",
                     remove.inorganic=TRUE, 
                     fluorine.limit=3, 
                     chlorine.limit=3, 
                     bromine.limit=3, 
                     iodine.limit=3, 
                     min.mass.limit=20, 
                     max.mass.limit=900)

#########################################
# Calculate descriptors for the molecules
#########################################
#, ExtendedFingerprinter

descriptors <- GeneratePadelDescriptors(standardised.file = "standardised.sdf", types=c("2D"), threads = 1)
dim(descriptors)
head(descriptors)

descriptors <- RemoveStandardisedPrefix(descriptors)
saveRDS(descriptors, file="descriptors.rds")
descriptors <- readRDS("descriptors.rds")

#########################################
# Target Visualization
#########################################

# Having a look at the response variable
targets <- read.csv("targets.csv", header=FALSE)
names(targets) <- c("Name", "target")
p <- DensityResponse(targets$Target)
p + labs(title="LogS target value distribution")

#########################################
# Merge the target values together with the descriptors
#########################################

all <- merge(x=targets, y=descriptors, by="Name")
# check the number of rows are the same
dim(all)
dim(targets)
dim(descriptors)
ids <- all$Name
x <- all[3:ncol(all)]
y <- all$target

#########################################
# Replace infinite values with NA and impute the NAs
#########################################
x.finite <- ReplaceInfinitesWithNA(x)
x.imputed <- ImputeFeatures(x.finite)

#########################################
# Split the dataset into a training and holdout set
#########################################
ss <- SplitSet(ids, x.imputed, y, percentage=20)

#########################################
# Remove the descriptors that are highly correlated or have low variance 
#########################################
ss <- RemoveNearZeroVarianceFeatures(ss)
ss <- RemoveHighlyCorrelatedFeatures(ss)

#########################################
# Preprocess the data (center and scale are the defaults used here)
#########################################
ss <- PreProcess(ss)

#########################################
# Generate a 5 folds for cross validation setup
# and save the data that is prepared for model training
#########################################
ss <- GetCVTrainControl(ss)
saveRDS(ss, file="ss.rds")

#########################################
# Training an Random Forest
#########################################

# Set the number of cores for parallelization of the training
registerDoMC(cores=1)

method <- "rf"
tune.grid <- expand.grid(.mtry = seq(5,100,5))

model <- train(ss$x.train, ss$y.train, method, tuneGrid=tune.grid, trControl=ss$trControl)
saveRDS(model, file=paste(method,".rds",sep=""))

#########################################
# Training an SVM
#########################################

method <- "svmRadial"
tune.grid <- expand.grid(.sigma = c(0.001, 0.002, 0.004, 0.008, 0.012), .C = c(0.5, 1, 2, 3, 4))

model <- train(ss$x.train, ss$y.train, method, tuneGrid=tune.grid, trControl=ss$trControl)
saveRDS(model, file=paste(method,".rds",sep=""))

library(errorestimatoR)
estimator = BuildCaretErrorEstimator(model, Nmax=20, cores=2)

indexes <- which(apply(model$pred[,names(model$bestTune)], 1, function(x,y) {identical(as.numeric(x), as.numeric(y))}, model$bestTune))
best <- model$pred[indexes,]
order <- order(best$rowIndex, decreasing=FALSE)
best <- best[order,]
obs <- best$obs
preds <- best$pred
fold.strings <- unique(best$Resample)
folds <- list()
for(fold in 1:length(fold.strings)) {
  folds[[fold]] <- which(best$Resample == fold.strings[fold])
}
estimator <- BuildEstimator(ss$x.train, folds, obs, preds, Nmax = 20, cores = 2)
cores <- 2
Nmax <- 20
errors <- preds - obs
dms <- CreateDistanceMatricesMC(x, folds, cores)
m <- CreateEstimatorMatrix(Nmax, folds, dms, errors, obs, preds)
sl <- GetSigmaSplineList(m, errors)
weights <- GetWeights(m, errors, sl, cores, optFunc)
list(x = x, sl = sl, weights = weights, obs = obs, preds = preds, errors = errors, Nmax = Nmax)

dim(ss$x.train)
    


#########################################
# Assesing Model Performance
#########################################

# Cross Validation Metrics.
# We assume the metric used for the choice of the best combination of hyperparameters is 'RMSE'.
# This can e chacke by: _my_model_$metric
RMSE_CV = signif(min(as.vector(na.omit(modelCoxRF$results$RMSE))), digits=3)
Rsquared_CV = modelCoxRF$results$Rsquared[which( modelCoxRF$results$RMSE %in% min(modelCoxRF$results$RMSE, na.rm=TRUE))]

# Predict the values of the hold-out (external) set
holdout.predictions <- as.vector(predict(modelCoxRF, newdata = dataset$x.holdout))

# Statistics for Model Validation
MetricsRf <- Validation(holdout.predictions,dataset$y.holdout)

# Correlation between observed and predicted
ObsPred(pred=holdout.predictions,obs=dataset$y.holdout,PointSize=3,ColMargin='green',
        margin=1,PointColor="black",PointShape=16,MarginWidth=2)