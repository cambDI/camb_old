# Author: Isidro Cortes Ciriano
# Plot erros bars
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


# Euclidean Distance of two vectors
euc.dist <- function(x1,x2){
  if (is.vector(v1) && is.vector(v2) && length(v1)==length(v2)){
  return(sqrt(sum((x1 - x2) ^ 2)))
  }
  else {print("Wrong input: input arguments are not vector or have unequal length")}
}

# Tanimoto Distance between two vectors
tani <- function(a,b){  
  if (length(a) != length(b)){print("vectors of unqual length");break}
  differentes = sum((a == b)*1) 
  comunes = length(a)
  tani =  differentes / comunes
  return(tani)
}

