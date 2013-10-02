# format data
load("Whole_dataset.RData")
x <- final_data[,20:ncol(final_data)]
y <- final_data$standard_value
yis4 <- which(y==4)
yis5 <- which(y==5)
y <- y[-c(yis4,yis5)]
x <- x[-c(yis4,yis5), ]
dim(x)

library(caret)
library(camb)
# training workflow
ViewTargets(y, bw=0.01)
x.imputed <- ImputeFeatures(x)
ss <- SplitSet(x.imputed, y, percentage=20)
ss <- RemoveNearZeroVarianceFeatures(ss)
ss <- RemoveHighlyCorrelatedFeatures(ss)
ss <- PreProcess(ss)
ss <- GetCVTrainControl(ss)
saveRDS(ss, file="ss.rda")





