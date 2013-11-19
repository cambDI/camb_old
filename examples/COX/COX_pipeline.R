#### Cambridge Workshop ####
#### November 2013 ####

# PCM Example with 'camb' #

library(camb)
setwd('/Users/icortes/Desktop/camb_final/camb/examples/COX')

#########################################
# Read, preprocess and calculate descriptors for the compounds
#########################################
smiles <- read.table("smiles_COX.smi", header=FALSE)
#StandardiseMolecules(structures.file="smiles_COX.smi", standardised.file="smiles_COX_processed.sdf", is.training=TRUE)
descriptors_COX <- GeneratePadelDescriptors(standardised.file = "smiles_COX.smi", threads = 1)

# imputation etc..

descriptors <- RemoveStandardisedPrefix(descriptors)
saveRDS(descriptors, file="descriptors.rds")
descriptors <- readRDS("descriptors.rds")


#########################################
# Calculate Circular Morgan Fingerprints
#########################################
Sys.setenv(RDBASE="/usr/local/share/RDKit")
Sys.setenv(PYTHONPATH="/usr/local/lib/python2.7/site-packages")

#########################################
# Read and calculate target descriptors
#########################################

amino_acids <- read.table("AAs_COX.csv",sep=",",header=TRUE,colClasses=c("character"),row.names=1)
amino_acids <- amino_acids[,2:ncol(amino_acids)]
amino_acids_zscales <- AA_descs(Data=amino_acids,type="Z3")
saveRDS(amino_acids_zscales,file="Z3_COX.rds")

#########################################
# Read the file with the info about the dataset: {target names, bioctivities, etc..}
#########################################

setwd('/Users/icortes/Desktop/camb_final/camb/examples/COX')
dataset <- read.table("COX_dataset_info.csv",sep=",",header=TRUE)
bioactivity <- dataset$standard_value

# The bioactivity is in nM, we convert it to pIC50
bioactivity <- bioactivity * 10^-9
bioactivity <- - log(bioactivity,base=10)

# Having a look at the response variable
DensityResponse(bioactivity)

# Let's have a look at the PCA of the target descriptors
target_PCA <- PCAProt(amino_acids_zscales,SeqsName=dataset$accession)
PCAProtPlot(target_PCA,main="PCA COX dataset",PointSize=8,TitleSize=25)


#########################################
# Removing repeated bioactivity datapoints (more than one annotation for the same compound-target combination)
#########################################

# We run the file "remove_duplicates.R"
source("remove_duplicates.R")
