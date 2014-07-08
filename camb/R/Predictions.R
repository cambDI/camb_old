#################################################################################
## Predictions on new molecules
#################################################################################

#' Make predictions on new molcules using a saved standardisation procedure and a saved model
#' 
#' Molecules are converted to a standard representation in the same way as during model training.
#' A saved model is used to make predictions on new molecules.
#' 
#' @param structures.file The name of the file containing the chemical structures. SMILES and SDF are currently supported formats.
#' @param standardisation.options The options saved during the standardisation procedure. These options are returned from the StandardiseMolecules function.
#' @param descriptor.types A named list of the types of descriptors used in model training.
#' @param dataset The dataset used in model training. This is used for the preprocessing applied to the training data as well as the descriptors used in training.
#' @param model The trained model.
#' @export
#' @return A data.frame containing the original ids of the molecules as well as their predicted values.
#' @examples
#' test_structures_file <- system.file("test_structures", "structures_10.sdf", package = "camb")
#' predictions <- PredictExternal(test_structures_file, standardisation.options, descriptor.types, dataset, readRDS("rf.rds"))
#' @author Daniel Murrell <dsmurrell@@gmail.com> and Isidro Cortes <isidrolauscher@@gmail.com>
PredictExternal <- function(structures.file, standardisation.options, descriptor.types, dataset, model) {
  model = model # this is done to make sure model gets evaluated successfully first
  
  standardised.file <- tempfile("standardised", fileext=".sdf")
  removed.file <- tempfile("removed", fileext=".sdf")
  properties.file <- tempfile("properties", fileext=".csv")
  
  StandardiseMolecules(structures.file = structures.file, 
                       standardised.file = standardised.file, 
                       removed.file = removed.file,
                       properties.file = properties.file,
                       standardisation.options)
  
  descriptors <- GeneratePadelDescriptors(standardised.file = standardised.file, types=c("2D"), threads = 1)
  descriptors <- RemoveStandardisedPrefix(descriptors)
  
  ids <- descriptors[,1]
  x <- descriptors[,2:ncol(descriptors)]
  
  x <- x[, names(dataset$x.train)]
  x.finite <- ReplaceInfinitesWithNA(x)
  x.imputed <- ImputeFeatures(x.finite)
  
  x.preprocessed <- predict(dataset$transformation, x.imputed)
  predictions <- predict(model, newdata=x.preprocessed)
  return(data.frame(id=ids, prediction=predictions))
}