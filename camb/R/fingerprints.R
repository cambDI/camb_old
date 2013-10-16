Fingerprints <- function(structures.file, structure.number, output.file, use.name.as.title) {
  
  # handle non-existant file
  if (!file.exists(structures.file)) {
    print("File does not exist")
  }
  
  split <- strsplit(structures.file, "\\.")[[1]]
  filetype <- split[length(split)]
  if(tolower(filetype) == "sdf") {
    .C("R_Fingerprints", structures.file, as.integer(structure.number), output.file, as.integer(use.name.as.title)) 
  }
  else if(tolower(filetype) == "smi") {
    print("Smiles doesn't work at the moment")
  }
  else {
    print("Unrecognised file type")
  }
}