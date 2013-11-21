library(camb)

StandardiseMolecules(structures.file="structures.sdf", 
                     standardised.file="standardised.sdf", 
                     removed.file="removed.sdf", 
                     remove.inorganic=TRUE, 
                     fluorine.limit=3, 
                     chlorine.limit=3, 
                     bromine.limit=3, 
                     iodine.limit=3, 
                     min.mass.limit=20, 
                     max.mass.limit=900)
                     