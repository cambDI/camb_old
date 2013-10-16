RemoveStandardisedPrefix <- function(descriptors) {
  descriptors$Name <- sapply(descriptors$Name, function(x) {strsplit(as.character(x), "Standardised_")[[1]][2]})
  descriptors
}

GeneratePadelDescriptors <- function(standardised.file, threads = -1, limit = -1) {
  descriptors.file <- tempfile("descriptors", fileext=".csv")
  GeneratePadelDescriptors.internal(standardised.file, descriptors.file, threads)
  read.csv(descriptors.file)
}

GeneratePadelDescriptorsFile <- function(standardised.file, descriptors.file, threads = -1, limit = -1) {
  GeneratePadelDescriptors.internal(standardised.file, descriptors.file, threads)
}

GeneratePadelDescriptors.internal <- function(structures.file, descriptors.file, threads = -1) {
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
  writeLines(lines, writeCon)
  close(readCon)
  close(writeCon)
      
  .jcall("padeldescriptor.PaDELDescriptorApp", , "main", c("-config", writefile))
}