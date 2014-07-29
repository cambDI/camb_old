library(devtools)
setwd("~/Dropbox/projects/camb/roxygen")
document('../camb')
setwd("~/Dropbox/projects/camb/examples/QSPR/LogS/Reference_2")

library(camb)
ls("package:camb")

?StandardiseMolecules
?DensityResponse
?GeneratePadelDescriptors
?SplitSet
?RemoveNearZeroVarianceFeatures
?RemoveHighlyCorrelatedFeatures
?PreProcess
?GetCVTrainControl
?CorrelationPlot
?caretEnsemble
?caretStack
?PredictExternal