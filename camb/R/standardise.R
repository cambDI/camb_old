StandardiseMolecules <- function(structures.file, standardised.file, is.training = FALSE, limit = -1) {
  # handle non-existant file
  if (!file.exists(structures.file)) {
    print("File does not exist")
  }
  
  # deal with sdf or smi difference
  split <- strsplit(structures.file, "\\.")[[1]]
  filetype <- split[length(split)]
  if(tolower(filetype) == "sdf") {
    print("Standardising Structures: Reading SDF (R)")
    sink(file="standardisation.log", append=FALSE, split=FALSE)
    .C("R_standardiseMolecules", structures.file, standardised.file, as.integer(1), as.integer(is.training), as.integer(limit)) 
    sink()
  }
  else if(tolower(filetype) == "smi") {
    print("Standardising Structures: Reading SMILES (R)")
    sink(file="standardisation.log", append=FALSE, split=FALSE)
    .C("R_standardiseMolecules", structures.file, standardised.file, as.integer(0), as.integer(is.training), as.integer(limit))
    sink()
  }
  else {
    print("Unrecognised file type")
  }
}

