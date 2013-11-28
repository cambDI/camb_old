#########################################
# Error estimation
#########################################
model <- readRDS("svmRadial.rds")
library(errorestimatoR)

estimator = BuildCaretErrorEstimator(dataset$x.train, model, Nmax=20, cores=1)

sigpred <- PredictSigmas(x=ss$x.holdout, estimator)

preds <- predict(model, newdata = ss$x.holdout)
errors <- preds - ss$y.holdout

model

k <- rbfdot(sigma=0.00391)

k(x=0, y=10)

k(3)
k(1,2)

plot(errors, sigpred$sigmas)
plot(errors, sigpred$sigma.matrix[,2])
x <- ss$x.test
CEC(sigpred$sigmas, abs(errors))

CEC(sigpred$sigma.matrix[,2], abs(errors))

trainsig <- PredictSigmas(x=ss$x.train, estimator)
trainpreds <- predict(model, newdata = ss$x.train)
trainerrors <- trainpreds - ss$y.train

plot(trainerrors, trainsig$sigmas)

CEC(abs(trainerrors), trainsig$sigmas)

CEC(abs(trainerrors), trainsig$sigma.matrix[,2])

length(abs(trainerrors))
length(trainsig$sigma.matrix[,2])

dim(estimator$x)
dim(x)

dm <- CreateNewDistanceMatrix(estimator$x, x)
nm <- CreateNewEstimatorMatrix(estimator$Nmax, dm, estimator$errors, 
                               estimator$obs, estimator$preds)
sigma.matrix <- GetNewSigmaMatrix(nm, estimator$sl)
sigmas <- GetNewSigmas(sigma.matrix, estimator$weights)
prediction <- list()
prediction$dm <- dm
prediction$nm <- nm
prediction$sigma.matrix <- sigma.matrix
prediction$sigmas <- sigmas
prediction