# Author: Isidro Cortes Ciriano
# Functions to analyse of the results of learning algorithms
# calculates the RMSE between two vectors
RMSE <- function(v1, v2) {
  i1 <- which(!is.na(v1))
  i2 <- which(!is.na(v2))
  is <- intersect(i1, i2)
  v1 <- v1[is]
  v2 <- v2[is]
  residuals <- abs(v1-v2)
  return(sqrt( (residuals%*%residuals)/length(v1) ))
}

# calculates the MAE between two vectors
MAE <- function (v1, v2) {
  i1 <- which(!is.na(v1))
  i2 <- which(!is.na(v2))
  is <- intersect(i1, i2)
  v1 <- v1[is]
  v2 <- v2[is]
  residuals <- abs(v1 - v2)
  return(sum(residuals)/length(v1))
}

# Tropsha's Statistics for Model Assessment
# Tropsha, A.; Golbraikh, A. Predictive Quantitative Structure–Activity Relationships Modeling: 
#Development and Validation of QSAR Models. In: Handbook of Chemoinformatics Algorithms 
#(Faulon, J.-L.; Bender, A., Eds.), Chapter 7, pp. 213-233, Chapman & Hall / CRC, London, UK, 2010.

# Calculates the slope between two vector (k')
slope <- function(v1,v2){ # v1=z.test v2=y.test
  return(sum(v2*v1)/sum(v1*v1))
}

# Calculates the regression coefficient through the origin
Rsquared0 <- function(v1,v2) { #v1=z.test (y), v2=y.test (x)
  if (is.vector(v1) && is.vector(v2) && length(v1)==length(v2)){
  y_obs_mean <- mean(v2)
  yr0 = v1 * slope(v1,v2)
  first_term = (v2 - yr0)*(v2 - yr0)
  second_term= (v2-y_obs_mean)*(v2-y_obs_mean)
  return(1-(sum(first_term)/sum(second_term)))
  }
  else {print("Wrong input: input arguments are not vector or have unequal length")}
}


# Calculates the regression coefficient 
Rsquared <- function(v1,v2) { # v1=z.test (y), v2=y.test (x)
  if (is.vector(v1) && is.vector(v2) && length(v1)==length(v2)){
  y_obs_mean <- mean(v2)
  y_pred_mean <- mean(v1)
  first_term <- sum((v2-y_obs_mean) * (v1 - y_pred_mean))
  second_term <- sqrt(sum((v2-y_obs_mean)*(v2-y_obs_mean)) * sum((v1 - y_pred_mean)*(v1 - y_pred_mean)))
  division <- first_term / second_term
  return(division * division)
  }
  else {print("Wrong input: input arguments are not vector or have unequal length")}
  
}

# Calculates the Q squared 
#Qsquared (z.test,y.test) (predicted vs observed)
Qsquared <- function(v1, v2) {
  if (is.vector(v1) && is.vector(v2) && length(v1)==length(v2)){
  y_obs_mean <- mean(v2)
  first_term <- abs(v1-v2)*abs(v1-v2)
  second_term <- abs(v2-y_obs_mean)*abs(v2-y_obs_mean)
  return(1-(sum(first_term)/sum(second_term)))
  }
  else {print("Wrong input: input arguments are not vector or have unequal length")}
}

################################


