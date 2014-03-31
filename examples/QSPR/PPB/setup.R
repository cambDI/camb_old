# read in the original data
origin <- read.csv("from_chris.csv")

# write out a SMILES file
df <- data.frame(SMILES = origin$original_format, Name = origin$ChemSpiderID)
write.table(df, file="drugbank.smi", sep=" ", col.names = F, row.names = F, quote=F)

library(camb)
# Standardise the original molecule set
StandardiseMolecules(structures.file="drugbank.smi", standardised.file="standardised.sdf", is.training=TRUE)

# Generate descriptors for the standardised list
descriptors <- GeneratePadelDescriptors(standardised.file = "standardised.sdf", threads = 1)
descriptors <- RemoveStandardisedPrefix(descriptors)
saveRDS(descriptors, file="descriptors.rds")
descriptors <- readRDS("descriptors.rds")

# Merge the target values together with the descriptors
targets <- data.frame(Name = origin$ChemSpiderID, target = origin$Value)
all <- merge(x=targets, y=descriptors, by="Name")
ids <- all$Name
x <- all[3:ncol(all)]
y <- all$target

# look at the distribution of the target values
ViewTargets(y, bw=00.1)

# replace infinite values with NA and impute NAs
x.finite <- ReplaceInfinitesWithNA(x)
x.imputed <- ImputeFeatures(x.finite)

# split the dataset into a training and holdout set
ss <- SplitSet(ids, x.imputed, y, percentage=20)

# remove the descriptors that are highly correlated or have low variance 
ss <- RemoveNearZeroVarianceFeatures(ss)
ss <- RemoveHighlyCorrelatedFeatures(ss)

# center and scale all the descriptors
ss <- PreProcess(ss)
# generate a 5 folds for cross validation
ss <- GetCVTrainControl(ss)

# save out the split set
saveRDS(ss, file="ss.rda")


