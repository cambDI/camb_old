library(devtools)
document('../camb')

.C("R_drawMoleculeInSDF", system.file("test_structures", "structures_10.sdf", package = "camb"), 1, "output.png", FALSE)