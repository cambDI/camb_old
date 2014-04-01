library(camb)

descriptors <- GeneratePadelDescriptors(standardised.file="../testing_standardistaion/standardised.sdf", types=c("2D", "ExtendedFingerprinter"))
dim(descriptors)
names(descriptors)