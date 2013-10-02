PredictLogPFromDescriptors <- function(descriptors, error.variance=FALSE) {
  print("Predicting LogP")
  load(file = system.file("extdata", "data.rda", package="camb"))
  load(file = system.file("extdata", "svm.rda", package="camb"))
  load(file = system.file("extdata", "gbm.rda", package="camb"))
  names <- descriptors[, 1]
  used.descriptors <- names(x.train)
  descriptors <- descriptors[, used.descriptors]
  suppressWarnings(descriptors <- apply(descriptors, 2, as.numeric))
  nrows <- nrow(descriptors)
  to.impute <- rbind(descriptors, x.train)
  set.seed(777)
  imputed <- as.data.frame(impute.knn(as.matrix(to.impute), k = 10)$data)
  to.predict <- imputed[1:nrows, ]
  x <- predict(transformation, to.predict)
  svm_pred <- predict(svm$finalModel, newdata = x)
  gbm_pred <- predict(gbm$finalModel, newdata = x, n.trees = 500)
  greedy_pred <- svm_pred*0.712 + gbm_pred*0.288
  r <- data.frame(ID = names, smLogP = greedy_pred)
  if(error.variance) {
    n <- 5
    p <- 1.5
    load(file = system.file("extdata", "spline.rda", package="camb"))
    load(file = system.file("extdata", "as.rda", package="camb"))
    outset <- to.predict
    for(i in 1:length(as)) {
      outset[,i] <- outset[,i] * (as[i]^p)
    }
    inset <- x.train
    for(i in 1:length(as)) {
      inset[,i] <- inset[,i] * (as[i]^p)
    }
    distances <- findMinDistances(outset, inset, n)
    sf <- 15/mean(distances)
    distances <- distances*sf
    variances <- predict(spline, distances)$y
    r$variance <- variances
  }
  r
}

PredictLogPFromDescriptorsFile <- function(descriptors.file, error.variance=FALSE) {
  descriptors <- read.csv(file = descriptors.file)
  PredictLogPFromDescriptors(descriptors, error.variance)
}

Train <- function(structures.file, error.variance=FALSE, threads = -1) {
  standardised.file <- tempfile("standardised", fileext=".sdf")
  name.file <- tempfile("name", fileext=".txt") # used to recover the original ordering of molecules (multicore descriptor generation messes this up)
  descriptors.file <- tempfile("descriptors", fileext=".csv")
  StandardiseMolecules(structures.file, standardised.file, name.file = name.file, limit = -1)
  GenerateDescriptors.internal(standardised.file, descriptors.file, name.file, threads)
  PredictLogPFromDescriptorsFile(descriptors.file, error.variance)
}

connectToDatabase <- function() {
  con <- dbConnect(MySQL(), user="root", password="leaves", dbname= "property", host = "localhost", port = 3306)
  con
}

getMeasurements <- function(type) {
  # retrieve the SMILES, compound id and measured value for compounds measurement sets of 'type'
  getMeasurementsQ = paste("select c.id as compoundID, s.id as measurement_setID, s.original_smiles as original_SMILES, c.SMILES as SMILES, m.value as Value, source.name as Source ",
                           "from source, compound as c, measurement as m, measurement_set as s, measurement_type as t ", 
                           "where s.source_id = source.id and m.measurement_set_id = s.id and s.compound_id = c.id and s.measurement_type_id = t.id and t.type = '", type, "'", sep = "")
  measurements = dbGetQuery(con, getMeasurementsQ)
  measurements
}

getReliableMeasurements <- function(measurements) {
  count <- count(measurements[,1]) 
  N <- nrow(count)
  reliable <- data.frame(compoundID = rep(0,N), smiles = rep("",N), value = rep(0,N), stringsAsFactors=FALSE)
  c <- 0
  for (i in 1:nrow(count)) {
    same <- measurements[which(measurements[,1] == count[i,1]), ]
    sd <- sd(same$Value)
    # reliable rule, if more than one measurment or sd is <= 0.5 then accept the mean value
    if(nrow(same) == 1 || sd <= 0.5) {
      c <- c+1
      reliable[c,] <- c(as.numeric(same$compoundID[1]), same$SMILES[1], mean(as.numeric(same$Value)))
    }
  }
  reliable <- reliable[1:c,]
  reliable <- transform(reliable, compoundID = as.numeric(compoundID), value = as.numeric(value))
  reliable
}

