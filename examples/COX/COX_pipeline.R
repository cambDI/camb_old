#### Cambridge Workshop ####
#### November 2013 ####

# PCM Example with 'camb' #

library(camb)
setwd('/Users/icortes/Desktop/camb_final/camb/examples/COX')

#########################################
# Read and preprocess compounds
#########################################
smiles <- read.table("smiles_COX.smi", header=FALSE)
StandardiseMolecules(structures.file="smiles_COX.smi", standardised.file="smiles_COX_processed.sdf", is.training=TRUE)
descriptors <- GeneratePadelDescriptors(standardised.file = "smiles_COX_processed.sdfdim(des.sdf", threads = 1)

# Generate descriptors for the standardised list
descriptors <- GeneratePadelDescriptors(standardised.file = "standardised.sdf", threads = 1)
descriptors <- RemoveStandardisedPrefix(descriptors)
saveRDS(descriptors, file="descriptors.rds")
descriptors <- readRDS("descriptors.rds")


#########################################
# Read and preprocess target descriptors
#########################################

amino_acids <- read.table("AAs_COX.csv",sep=",",header=TRUE)



# Merge the target values together with the descriptors
targets <- data.frame(Name = origin$ChemSpiderID, target = origin$Value)
all <- merge(x=targets, y=descriptors, by="Name")
ids <- all$Name
x <- all[3:ncol(all)]
y <- all$target

