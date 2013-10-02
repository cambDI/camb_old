GeneratePadelDescriptors <- function(structures.file, threads = -1, limit = -1) {
  standardised.file <- tempfile("standardised", fileext=".sdf")
  name.file <- tempfile("name", fileext=".txt")
  descriptors.file <- tempfile("descriptors", fileext=".csv")
  StandardiseMolecules(structures.file, standardised.file, name.file = name.file, limit = limit)
  GenerateDescriptors.internal(standardised.file, descriptors.file, name.file, threads)
  read.csv(descriptors.file)
}

GeneratePadelDescriptorsFile <- function(structures.file, descriptors.file, threads = -1, limit = -1) {
  standardised.file <- tempfile("standardised", fileext=".sdf")
  name.file <- tempfile("name", fileext=".txt")
  StandardiseMolecules(structures.file, standardised.file, name.file = name.file, limit = limit)
  GenerateDescriptors.internal(standardised.file, descriptors.file, name.file, threads)
}

GeneratePadelDescriptors.internal <- function(structures.file, descriptors.file, name.file, threads = -1) {
  print("Generating Descriptors")
  .jinit()
  .jcall("java/lang/System","S","setProperty","java.awt.headless","true")
  # add all JARs
  .jaddClassPath(Sys.glob("lib/*.jar"))
  # call the main() method on the main class
  padel_config_file <- system.file("extdata", "config.txt", package="camb")
  readCon  <- file(padel_config_file, open = "r")
  writefile <- tempfile("config", fileext=".txt")
  writeCon <- file(writefile)
  lines <- c()
  while (length(line <- readLines(readCon, n = 1, warn = FALSE)) > 0) {
    vec <- unlist(strsplit(line, "="))
    if(vec[1]=="MaxThreads") {
      lines <- c(lines, paste("MaxThreads=", threads, sep=""))
    }
    else if(vec[1]=="DescriptorFile") {
      lines <- c(lines, paste("DescriptorFile=", descriptors.file, sep=""))
    } 
    else if(vec[1]=="Directory") {
      lines <- c(lines, paste("Directory=", structures.file, sep=""))
    }
    else {
      lines <- c(lines, line)
    }
  }
  #print(lines)
  writeLines(lines, writeCon)
  close(readCon)
  close(writeCon)
      
  .jcall("padeldescriptor.PaDELDescriptorApp", , "main", c("-config", writefile))
  
  # reads the name.file and puts the molecules back in their original order using the names in the SMILES/SDF    
  descs <- read.csv(descriptors.file)
  df <- data.frame(number = 1:nrow(descs))
  readCon2  <- file(name.file, open = "r")
  names <- c()
  while (length(line <- readLines(readCon2, n = 1, warn = FALSE)) > 0) {
    names <- c(names,line)
  }
  close(readCon2)
  df$Name <- names
  m <- merge(df, descs, by.x="number", by.y="Name")
  write.csv(m[, 2:ncol(m)], descriptors.file, row.names=FALSE)    
}