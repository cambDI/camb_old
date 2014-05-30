library(devtools)

setwd("~/Dropbox/projects/camb/roxygen")
document('../camb')
setwd("~/Dropbox/projects/camb/examples/QSPR/LogS/Reference_2")

#?StandardiseMolecules
#.C("R_drawMoleculeInSDF", system.file("test_structures", "structures_10.sdf", package = "camb"), 1, "output.png", FALSE)