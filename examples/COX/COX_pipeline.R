''' Cambridge Workshop '''
''' November 2013 '''
''' Isidro Cortes and Daniel Murrell '''
''' PCM Example with camb '''

library(camb)
setwd('/Users/icortes/Desktop/camb_final/camb/examples/COX')

#########################################
# Read, preprocess and calculate descriptors for the compounds
#########################################
smiles <- read.table("smiles_COX.smi", header=FALSE)
#StandardiseMolecules(structures.file="smiles_COX.smi", standardised.file="smiles_COX_processed.sdf", is.training=TRUE)
descriptors_COX <- GeneratePadelDescriptors(standardised.file = "smiles_COX.smi", threads = 1)

descriptors <- RemoveStandardisedPrefix(descriptors)
saveRDS(descriptors, file="descriptors.rds")
descriptors <- readRDS("descriptors.rds")

#########################################
# Calculate Circular Morgan Fingerprints
#########################################
Sys.setenv(RDBASE="/usr/local/share/RDKit")
Sys.setenv(PYTHONPATH="/usr/local/lib/python2.7/site-packages")
#fps_COX_512 <- MorganFPs(bits=512,radius=2,type='smi',mols='smiles_COX.smi',output='COX',keep='hashed_counts')
#saveRDS(fps_COX_512,file="fps_COX_512.rds")
fps_COX_512 <- readRDS("fps_COX_512.rds")

#########################################
# Read and calculate target descriptors
#########################################

amino_accompound_compound_IDs <- read.table("AAs_COX.csv",sep=",",header=TRUE,colClasses=c("character"),row.names=1)
amino_accompound_IDs <- amino_accompound_IDs[,2:ncol(amino_accompound_IDs)]
amino_accompound_IDs_zscales <- AA_descs(Data=amino_accompound_IDs,type="Z3")
saveRDS(amino_accompound_IDs_zscales,file="Z3_COX.rds")

#########################################
# Read the file with the info about the dataset: {target names, bioctivities, etc..}
#########################################

setwd('/Users/icortes/Desktop/camb_final/camb/examples/COX')
dataset <- readRDS("COX_dataset_info.rds")
bioactivity <- dataset$standard_value

# The bioactivity is in nM, we convert it to pIC50
bioactivity <- bioactivity * 10^-9
bioactivity <- - log(bioactivity,base=10)

#########################################
# Dataset Visualization
#########################################

# Having a look at the response variable
DensityResponse(bioactivity)

# Let's have a look at the PCA of the target descriptors
target_PCA <- PCAProt(amino_accompound_IDs_zscales,SeqsName=dataset$accession)
PCAProtPlot(target_PCA,main="PCA COX dataset",PointSize=8,TitleSize=25)

# Pairwise Compound Similarity
pw_dist_comp_fps <- PairwiseDist(fps_COX_512,method="jaccard")
#saveRDS(pw_dist_comp_fps,file="pairiwse_dist_COX.rds")
readRDS("pairiwse_dist_COX.rds")
PairwiseDistPlot(pw_dist_comp_fps)

# Which is the maximum achievable performance?
MaxPerf(meanNoise=0,sdNoise=0.68,meanResp=mean(bioactivity),sdResp=sd(bioactivity),lenPred=800)

#########################################
# Removing repeated bioactivity datapoints (more than one annotation for the same compound-target combination)
#########################################

# We run the file "remove_duplicates.R"
source("remove_duplicates.R")

#########################################
# Preprocessing of the Dataset
#########################################

dataset <- readRDS("Whole_dataset_NO_REP.rds")
# Removing columns not containing descriptors
killset <- expression(c(tid,pref_name,accession,organism,chembl_id,standard_value,standard_units,
                        standard_type,chembl_id.1,Name,Name.1,Name.2,rows))
# We get the processed bioactivity (it was already converted to pIC50 units when running 'remove_duplicates.R')
bioactivity <- dataset$standard_value
compound_IDs <- dataset$chembl_id.1
dataset <- subset(dataset,select=-eval(killset))

# split the dataset into a training and holdout set
dataset <- SplitSet(compound_IDs, dataset, bioactivity, percentage=30)

# remove the descriptors that are highly correlated or have low variance
dataset <- RemoveNearZeroVarianceFeatures(dataset,frequencyCutoff=30/1)
dataset <- RemoveHighlyCorrelatedFeatures(dataset)

# center and scale all the descriptors (conversion to z-scores)
dataset <- PreProcess(dataset)

# generate 5 folds for cross validation
dataset <- GetCVTrainControl(dataset)

# save out the split set
saveRDS(dataset, file="dataset_COX_preprocessed.rda")

#########################################
# Training a SVM
#########################################
dataset <- readRDS("dataset_COX_preprocessed.rda")

# Set the number of cores for parallelization of the training
registerDoMC(cores=4)

method <- "rvmRadial"
# We define an exponential grid to optimize the hyperparameters
exp_grid <- expGrid(ini=-8,end=-6,stride=2,base=2)
tune.grid <- expand.grid(.sigma = exp_grid)

modelCoxSVMrad <- train(dataset$x.train, dataset$y.train, method, tuneGrid=tune.grid, trControl=dataset$trControl)
saveRDS(modelCoxSVMrad, file=paste("COX",method,".rds",sep=""))

#########################################
# Training a Random Foresy
#########################################

method <- "rf"
tune.grid <- expand.grid(.mtry = seq(5,100,5))

modelCoxRF<- train(dataset$x.train, dataset$y.train, method, tuneGrid=tune.grid, trControl=dataset$trControl)
saveRDS(modelCoxRF, file=paste("COX",method,".rds",sep=""))

dataset <- readRDS("dataset_COX_preprocessed.rda")


#########################################
# Assesing Model Performance
#########################################
CVRMSE <- signif(min(as.vector(na.omit(modelCoxSVMrad$results$RMSE))), digits=3)
holdout.predictions <- as.vector(predict(modelCoxSVMrad, newdata = dataset$x.holdout))


# Correlation between observed and predicted
ObsPred()



