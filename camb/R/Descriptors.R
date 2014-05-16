#################################################################################
## Standardization and Descriptor Calculation of Molecules and Peptides / Sequences
#################################################################################

#' Convert molecules to a standard representation
#' 
#' Details about the standardisation procedure
#' 
#' @param structures.file A character, vector, matrix or data.frame containing the amino acids in either one-letter or three-letter format. 
#' Amino acids symbols are valid in capitals or in lower-case.
#' @export
#' @return A data.frame with the properties of the molecules read from the original structure file.
#' @references \url{http://www.ggasoftware.com/opensource/indigo}
#' @examples
#' test_mols <- system.file("test_structures", "structures_10.sdf", package = "camb")
#' StandardiseMolecules(structures.file=test_mols, standardised.file="st.sdf", removed.file="removed.sdf", properties.file="properties.csv")
#' @author Daniel Murrell <dsmurrell@@gmail.com> and Isidro Cortes <isidrolauscher@@gmail.com>
StandardiseMolecules <- function(structures.file, 
                                 standardised.file,
                                 removed.file = "",
                                 properties.file="standardisation_info.csv",
                                 remove.inorganic = FALSE, 
                                 fluorine.limit = -1,
                                 chlorine.limit = -1,
                                 bromine.limit = -1,
                                 iodine.limit = -1,
                                 min.mass.limit = -1,
                                 max.mass.limit = -1,                      
                                 number.processed = -1) {
  # handle non-existant file
  if (!file.exists(structures.file)) {
    stop("File does not exist")
  }
  if (file.info(structures.file)$size  == 0) {stop("Input file is empty")}
  
  # deal with sdf or smi difference
  split <- strsplit(structures.file, "\\.")[[1]]
  filetype <- split[length(split)]
  if(tolower(filetype) == "sdf") {
    print("Standardising Structures: Reading SDF (R)")
    sink(file="standardisation.log", append=FALSE, split=FALSE)
    .C("R_standardiseMolecules", 
       structures.file, 
       standardised.file, 
       removed.file,
       properties.file,
       as.integer(1), # process SDF
       as.integer(remove.inorganic), 
       as.integer(fluorine.limit),
       as.integer(chlorine.limit),
       as.integer(bromine.limit),
       as.integer(iodine.limit),
       as.integer(min.mass.limit),
       as.integer(max.mass.limit),
       as.integer(number.processed)) 
    sink()
  }
  else if(tolower(filetype) == "smi") {
    print("Standardising Structures: Reading SMILES (R)")
    sink(file="standardisation.log", append=FALSE, split=FALSE)
    .C("R_standardiseMolecules", 
       structures.file, 
       standardised.file, 
       removed.file,
       properties.file,
       as.integer(0), # process SMILES
       as.integer(remove.inorganic), 
       as.integer(fluorine.limit),
       as.integer(chlorine.limit),
       as.integer(bromine.limit),
       as.integer(iodine.limit),
       as.integer(min.mass.limit),
       as.integer(max.mass.limit),
       as.integer(number.processed)) 
    sink()
  }
  else {
    print("Unrecognised file type")
  }
}

GetPropertiesSDF <- function(structures.file,number_processed=-1, type=1){ ## 1 refers to not smiles
##print("Get properties from SDF")
if (!file.exists(structures.file)) {stop("File does not exist")}
if (file.info(structures.file)$size  == 0) {stop("Input file is empty")}
sink(file="getPropertiesSDF.log", append=FALSE, split=FALSE)
properties.file <- tempfile("properties",fileext=".csv")
.C("R_GetPropertiesSDF",structures.file,as.integer(number_processed),as.integer(type),output)
print("Reading..")
properties <- read.table(output,sep="\t",header=TRUE)
properties <- properties[, !apply(is.na(properties), 2, all)]
sink()
return(properties)
}

ShowPropertiesSDF <- function(structures.file,type=1){ ## 1 refers to not smiles
if (!file.exists(structures.file)) {stop("File does not exist")}
if (file.info(structures.file)$size  == 0) {stop("Input file is empty")}
output <- tempfile("props_temp",fileext=".csv")
.C("R_ShowPropertiesSDF",structures.file,output,as.integer(type))
if (file.info(output)$size  == 0) {stop("The molecules in the file provided do not contain any property")}
props <- as.vector(read.table(output,header=FALSE)$V1)
return(props)
}

GetPropertySDF <- function(structures.file, property="", number_processed=-1, type=1 ){ ## 1 refers to not smiles
##print("Get properties from SDF")
if (!file.exists(structures.file)) {stop("File does not exist")}
if (file.info(structures.file)$size  == 0) {stop("Input file is empty")}
sink(file="getPropertySDF.log", append=FALSE, split=FALSE)
output <- tempfile("prop_temp",fileext=".csv")
.C("R_GetPropertySDF",structures.file,property,as.integer(number_processed),as.integer(type),output)
sink()
if (file.info(output)$size  == 0) {stop("The molecules in the file provided do not contain the specified property")}
prop <- read.table(output)
names(prop) <- property
return(prop)
}
##############
## Compound Descriptors
RemoveStandardisedPrefix <- function(descriptors) {
  descriptors$Name <- sapply(descriptors$Name, function(x) {strsplit(as.character(x), "Standardised_")[[1]][2]})
  descriptors
}

GeneratePadelDescriptors <- function(standardised.file, types = c("2D"), threads = -1, limit = -1) {
  if (file.info(standardised.file)$size  == 0) {stop("Input file is empty")}
  descriptors.file <- tempfile("descriptors", fileext=".csv")
  GeneratePadelDescriptors.internal(standardised.file, descriptors.file, types, threads)
  read.csv(descriptors.file)
}

GeneratePadelDescriptorsFile <- function(standardised.file, descriptors.file, types = c("2D"), threads = -1, limit = -1) {
  if (file.info(standardised.file)$size  == 0) {stop("Input file is empty")}
  GeneratePadelDescriptors.internal(standardised.file, descriptors.file, types, threads)
}

# options for types "2D", "Fingerprinter", "ExtendedFingerprinter", "EStateFingerprinter", "GraphOnlyFingerprinter",
# "MACCSFingerprinter", "PubchemFingerprinter", "SubstructureFingerprinter", "SubstructureFingerprintCount", 
# "KlekotaRothFingerprinter", "KlekotaRothFingerprintCount"
GeneratePadelDescriptors.internal <- function(structures.file, descriptors.file, types, threads = -1) {
  print("Generating Descriptors")
  .jinit()
  .jcall("java/lang/System","S","setProperty","java.awt.headless","true")
  # add all JARs
  .jaddClassPath(Sys.glob("lib/*.jar"))
  # call the main() method on the main class
  
  # replace the options in the padel_config.txt file
  padel_config_file <- system.file("extdata", "padel_config.txt", package="camb")
  readCon  <- file(padel_config_file, open = "r")
  conffile <- tempfile("padel_config", fileext=".txt")
  confCon <- file(conffile)
  lines <- c()
  while (length(line <- readLines(readCon, n = 1, warn = FALSE)) > 0) {
    vec <- unlist(strsplit(line, "="))
    if(vec[1]=="Compute2D") {
      if("2D" %in% types) {
        lines <- c(lines, "Compute2D=true")
      }
    }
    else if(vec[1]=="MaxThreads") {
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
  writeLines(lines, confCon)
  close(readCon)
  close(confCon)
  
  # replace the options in the padel_types.xml file
  padel_types_file <- system.file("extdata", "padel_types.xml", package="camb")
  readCon  <- file(padel_types_file, open = "r")
  descfile <- tempfile("padel_types", fileext=".xml")
  descCon <- file(descfile)
  lines <- c()
  while (length(line <- readLines(readCon, n = 1, warn = FALSE)) > 0) {
    vec <- unlist(strsplit(line, "\""))
    if(vec[2] %in% types) {
      vec[4] <- "true"
      lines <- c(lines, paste(vec, collapse="\""))
    }
    else {
      lines <- c(lines, line)
    }
  }
  writeLines(lines, descCon)
  close(readCon)
  close(descCon)
  
  .jcall("padeldescriptor.PaDELDescriptorApp", , "main", c("-config", conffile, "-descriptortypes", descfile))
}

##############
# Check if an AA is natural
checkAA <- function(x) {
  AADict = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", 
             "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
             "ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS",
             "LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL",
             "TRP","TYR")
  return(x %in% AADict)
}

##############
## Three letter to one letter AA code
convert31 <- function(AA) {  
  threeL <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
              "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
              "PRO", "SER", "THR", "TRP", "TYR", "VAL")  
  names(threeL) <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
                     "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")  
  res <- names(threeL[which(threeL==toupper(AA))])
  return(res)
}

#' Amino Acid Descriptor Calculation
#' 
#' The function calculates amino acid descriptors for natural amino acids.
#' Currently available descriptors are: 3 and 5 z-scales ("Z3" and "Z5"),
#' T-scales ("TScales"), ST-scales ("STScales"),
#' principal components score Vectors of Hydrophobic, Steric, and Electronic properties ("VHSE"),
#' BLOSUM ("BLOSUM"), FASGAI ("FASGAI"), MSWHIM ("MSWHIM"), and ProtFP ("ProtFP8").
#' See references for further information on these descriptors.
#' 
#' @param data   A character, vector, matrix or data.frame containing the amino acids in either one-letter or three-letter format. 
#' Amino acids symbols are valid in capitals or in lower-case.
#' @param type   Type of descriptors to be calculated. Default value is 5 z-scales. 
#' Any combination of descriptors is valid. 
#' A vector containing the abbreviation of the desired descriptors is taken as argument.
#' @export
#' @return   A data.frame which columns are indexed by the descriptors, and
#' rows by the rows amino acid sequence.
#' If the input data is a matrix or data.frame, the number of rows in the original data.frame or matrix,
#' and the number of rows of the ouput data.frame are equal.
#' If several descriptor types are chosen, descriptors are concatenated for the ease of further modeling.
#' Column names indicate the amino acid position in the original input data, and the type of descriptor. 
#' @references \url{http://www.jcheminf.com/content/5/1/41}
#' @references \url{http://www.jcheminf.com/content/5/1/42}
#' @examples
#' AA_descs(c("A","A"))
#' @author Isidro Cortes <isidrolauscher@@gmail.com> and Daniel Murrell <dsmurrell@@gmail.com>
AA_descs <- function(data, type="Z5", ..) {
  if (!is.vector(data) && !is.character(data) && !is.data.frame(data) && !is.matrix(data)){
    stop("Input must be a character, vector, data frame or matrix")
  } else {    
    match_AA1 <- function(AA,sel.=sel){
      if (checkAA(toupper(AA))){
        AA <- toupper(AA)
        if (nchar(AA) == 1){
          d <- descs[descs$AA1 == AA,sel.+2]
        } else {
          AA <- convert31(AA)
          d <- descs[descs$AA1 == AA,sel.+2]
        }
        return(d)
      } else {
        stop("Non natural Amino acid provided")
      }
    }
    
    match_AA1_vec <-function(v,sel.=sel){
      if (is.vector(v)){
        des <- sapply(v,match_AA1,sel.=sel)
        namesDes <- row.names(des)
        des <- unlist(des)
        namesDes <- c(t(sapply(namesDes,paste0,"_",v,seq(1,length(v)))))
        names(des) <- namesDes
        return(des)
      } else {
        stop("Input is not a vector or character")
      }
    }
    
    match_AA1_df <- function(df,colNames,sel.){
      if (is.data.frame(df) || is.matrix(df)){
        des <- t(apply(df,1,match_AA1_vec))
        row.names(des) <- seq(1,nrow(df))
        tableDescs <- table(colNames)
        colNamesPos=c()
        for(i in 1:length(tableDescs)){
          now <- which(colNames %in% names(tableDescs)[i])
          colNamesPos <- append(colNamesPos,paste0(colNames[now],"_",seq(1,length(now))))
        }
        colnames(des) <- c(t(sapply(colNamesPos,paste0,"_Aa",seq(1,ncol(df)))))
        return(des)
      } else {
        stop("Input is neither a data.frame nor a matrix")
      }
    }
    
    if  (is.vector(data)){
	  descs_path <- system.file("extdata", "aa_descs.rds", package="camb")
	  descs <- readRDS(descs_path)
      types <- c("ProtFP8","TScales","VHSE","STScales","BLOSUM","FASGAI","MSWHIM","Z5","Z3")
      type <- match.arg(type,types,several.ok=TRUE)
      root <- strsplit(names(descs)[3:ncol(descs)],"_")
      root <- unlist(root)[seq(1,length(unlist(root)),2)]
      sel <- which(root %in% type)
      return(match_AA1_vec(data,sel))
    } else {
	  descs_path <- system.file("extdata", "aa_descs.rds", package="camb")
	  descs <- readRDS(descs_path)
	  types <- c("ProtFP8","TScales","Tscales","VHSE","STScales","BLOSUM","FASGAI","MSWHIM","Z5","Z3")
      type <- match.arg(type,types,several.ok=TRUE)
      root <- strsplit(names(descs)[3:ncol(descs)],"_")
      root <- unlist(root)[seq(1,length(unlist(root)),2)]
      typeExt <- root[which(root %in% type)]
      sel <- which(root %in% type)
      return(match_AA1_df(data,colNames=typeExt,sel=sel))
    }
  }
}





#' Whole Protein Sequence Descriptor Calculation
#' 
#' Calculation of the following 12 whole sequence protein descriptors:
#'   - Amino Acid Composition ('AAC')
#'   - Dipeptide Composition ('DC')
#'   - Tripeptide Composition ('TC')
#'   - Normalized Moreau-Broto Autocorrelation ('MoreauBroto')
#'   - Moran Autocorrelation ('Moran')
#'   - Geary Autocorrelation ('Geary')
#'   - CTD (Composition/Transition/Distribution) ('CTD')
#'   - Conjoint Traid ('CTriad')
#'   - Sequence Order Coupling Number ('SOCN')
#'   - Quasi-sequence Order Descriptors ('QSO')
#'   - Pseudo Amino Acid Composition ('PACC')
#'   - Amphiphilic Pseudo Amino Acid Composition ('APAAC')
#' 
#' @param data   One or more protein sequences, or one or several UniProt IDs.
#' @param UniProtID If TRUE the argument calculates the descriptors for the proteins 
#' which UniProt IDs have been indicated in the argument 'data'.
#' @param type     The type of protein descriptors to be calculated (see above).
#' Any combination thereof is valid. A vector containing the abbreviation of the desired descriptors
#' is taken as argument. Default value 'AAC'.
#' @export
#' @return     A numeric matrix which rows are indexed by proteins and the columns by descriptors.
#' If multiple descriptors are chosen, the function returns a matrix where descriptors are concatenated per row
#' for the ease of modeling. 
#' @references R package \code{\link{protr}}.
#' @author Isidro Cortes <isidrolauscher@@gmail.com> and Daniel Murrell <dsmurrell@@gmail.com>
SeqDescs <- function(data,UniProtID=TRUE,type="AAC",..){
  if (UniProtID){
    if (!is.vector(data)) stop("The input UNIPROT identifiers must be in a vector")
    d <- data.frame(unlist(lapply(data,getUniProt)))
    types <- c("AAC","DC","TC","MoreauBroto","Moran","Geary",
               "CTDC","CTDT","CTDD","CTriad","SOCN","QSO",
               "PACC","APAAC")
    type <- match.arg(type,types,several.ok=TRUE)
    type <- paste0("extract",type,"(x)")
    type <- paste("c(",paste(type,collapse=","),")",sep="")
    t <- paste("t(apply(d,1,FUN=function(x)",type,"))",sep="") 
    des=eval(parse(text=t))
    row.names(des) <- data
    return(des) 
  } else {
    if (!is.data.frame(data) && !is.matrix(data)) stop("The input sequences must be a dataframe or a matrix")
    types <- c("AAC","DC","TC","MoreauBroto","Moran","Geary",
               "CTDC","CTDT","CTDD","CTriad","SOCN","QSO",
               "PACC","APAAC")
    type <- match.arg(type,types,several.ok=TRUE)
    type <- paste0("extract",type,"(x)")
    type <- paste("c(",paste(type,collapse=","),")",sep="")
    t <- paste("t(apply(data,1,FUN=function(x)",type,"))",sep="") 
    des=eval(parse(text=t))
    row.names(des) <- rownames(data)
    return(des) 
  }
}

#' Circular Morgan Fingerprints as specified in RDkit
#' 
#' The function calculates circular Morgan fingerprints for chemical compounds
#' using the RDkit python library (Greg Landrum).
#' Hashed fingeprints are calculated in binary format or with counts.
#' In addition, it also calcualates unhashed fingerprints, both binary and with counts.
#' 
#' @param bits Number of bits of the hashed fingerprints.
#' @param radius Radius of the hashed fingerprints. A radius of 2 is equivalent to ECFP-4 fingerprints.
#' @param type File format containing the input molecules.
#' @param mols File containing the input molecules.
#' @param output Labels that will be appended to all ouput files.
#' @param keep The fingeprints that will be kept after the calculation.
#' Apart from calculating different types of fingerprints,
#' the function returns a data.frame with the type of fingerprints specified here.
#' Possible types are: hashed_binary, hashed_counts, unhashed_binary, unhashed_counts, 
#' hashed_binaryEXT, hashed_countsEXT, unhashed_binaryEXT and unhashed_countsEXT.
#' @param images   If TRUE individual .pdf files containing the image of each substructure present in the input file,
#' and for each molecule, are created. Be aware that the number of fingerprints can be large depending on the number 
#' and diversity of the molecules present in the input file. Thus, allow for sufficient memory.
#' @param If TRUE, unhashed fingeprints -both in binary format and with counts- are calculated.
#' @param verbose If TRUE, information about the progression of the calculation is printed.
#' @param RDkitPath The path to the folder containing the RDkit library in your computer.
#' On mac, the is equal to the environment variable $RDBASE.
#' @param PythonPath Path to python ($PYTHONPATH).
#' @param extFileExtension File extension for the file containing the molecules for which unhashed fingeprints are to be calculated
#' with respect to the pool of substructures in the molecules present in the file specified in 'mols'.
#' @param extMols   File containing the molecules for which unhashed fingeprints are to be calculated
#' with respect to the pool of substructures in the molecules present in the file specified in 'mols'.
#' @param unhashedExt If TRUE, unhashed fingerprints are calcualted for the molecules in 'extMols'.
#' @param logFile File where the log messages will be dropped.
#' @export
#' @return In the working directory, .csv files will be created containing the different fingeprint types
#' specified with the function arguments. By default, hashed fingerprint in binary format and with counts 
#' will be created. In addition, the function returns a data.frame with the fingerprint types defined in 
#' the argument 'keep'. In the data.frame, rows are indexed by the molecules in the file containing the 
#' molecules, and columns by the bits in the fingerprint vector. 
#' @references FingeprintCalculator.py. Isidro Cortes. 2013/2014. \url{http://github.com/isidroc/FingerprintCalculator}.
#' @examples   
#' test_mols <- system.file("test_structures", "structures_10.sdf", package = "camb")
#' MorganFPs(bits=28,radius=4,type="sdf",mols=test_mols,output="test_mols")
#' @author Isidro Cortes <isidrolauscher@@gmail.com> and Daniel Murrell <dsmurrell@@gmail.com>
MorganFPs <- function (bits=512,radius=2,type="smi",mols,output,keep="hashed_binary",images=FALSE,unhashed=FALSE,
                       verbose=FALSE,RDkitPath="/usr/local/share/RDKit",PythonPath="/usr/local/lib/python2.7/site-packages",
                       extFileExtension=FALSE,extMols=FALSE,unhashedExt=FALSE,logFile=FALSE) {
  types <- c("smi","smiles","sdf")
  type <- match.arg(type,types,several.ok=FALSE)
  if (is.na(type)){
    stop("Input formats currently supported are .smi, .smiles and .sdf")
  }
  if (is.null(RDkitPath)) RDkitPath <- "/usr/local/share/RDKit"
  if (is.null(PythonPath)) RDkitPath <- "/usr/local/lib/python2.7/site-packages"
  
  pyfpath <- system.file("extdata", "FingerprintCalculator.py", package="camb")
  t <- paste(pyfpath," --bits",bits,"--rad",radius,"--f",type,
             "--mols",mols,"--output",output,"--RDkitPath",RDkitPath,sep=" ")
  if (images) t <- paste(t,"--image",sep=" ")
  if (verbose) t <- paste(t,"--v",sep=" ")
  if (unhashed) t <- paste(t,"--unhased",sep=" ")
  if (extFileExtension && extMols) t <- paste(t,"--extF", extFileExtension, "--molsEXT",extMols,sep=" ")
  if (unhashedExt) t <- paste(t, "--unhashedEXT",unhashedExt,sep=" ")
  if (logFile) t <- paste(t, " > ",output,".log",sep="")
  system(t) 
  types_keep <- c("hashed_binary","hashed_counts","unhashed_binary","unhashed_counts",
                  "hashed_binaryEXT","hashed_countsEXT","unhashed_binaryEXT","unhashed_countsEXT")
  keep <- match.arg(keep,types_keep,several.ok=FALSE)
  if (is.na(keep)){
    stop("Cannot load the fingerprints file. Possible types are:
         hashed_binary, hashed_counts, unhashed_binary, unhashed_counts, hashed_binaryEXT, hashed_countsEXT, unhashed_binaryEXT and unhashed_countsEXT")
  }  
  loadFile=paste(output,"_",keep,".csv",sep="")
  p=read.table(loadFile,sep=",")
  return(p)
}

##############
# TBD: not sure if this function is called externally... I don't think it's called from somewhere else in this package... should we remove it?
loadMorganFPs <- function (type="hashed_binary",output){
  types_keep <- c("hashed_binary","hashed_counts","unhashed_binary","unhashed_counts",
                  "hashed_binaryEXT","hashed_countsEXT","unhashed_binaryEXT","unhashed_countsEXT")
  type <- match.arg(type,types_keep,several.ok=FALSE)
  if (is.na(type)){
    stop("Cannot load the fingerprints file. Possible types are:
         hashed_binary, hashed_counts, unhashed_binary, unhashed_counts, hashed_binaryEXT, hashed_countsEXT, unhashed_binaryEXT and unhashed_countsEXT")
  }  
  loadFile=paste(output,type,".csv",sep="_")
  p=read.table(loadFile,sep=",")
  return(p)
}

#' Merge descriptors blocks.
#' 
#' The function merges blocks of descriptors by columns.
#' 
#' @param a1..a6 Descriptors blocks with the same number of rows to be merged.
#' @export
#' @return The merged block of descriptors.
#' @author Isidro Cortes <isidrolauscher@@gmail.com> and Daniel Murrell <dsmurrell@@gmail.com>
mergeData <- function(a1,a2,a3,a4=NULL,a5=NULL,a6=NULL){
  return(cbind(a1,a2,a3,a4,a5,a6,deparse.level=0))
}

