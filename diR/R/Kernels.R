# Author: Isidro Cortes Ciriano
###############################################################
############## Kernels ########################################
###############################################################

require(kernlab)
# Pearson VII PUK Kernel
puk<- function(sigma=1, omega =1) 
{
  rval <- function(x,y=NULL)   
  {    
    if(!is(x,"vector")) stop("x must be a vector")    
    if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")   
    if (is(x,"vector") && is.null(y)){      
      return(1)     
    }    
    if (is(x,"vector") && is(y,"vector")){      
      if (!length(x)==length(y))        
        stop("number of dimension must be the same on both data points")            
      return(1/((1 + ((2*(sqrt(((2*crossprod(x,y) - crossprod(x) - crossprod(y))^2 * sqrt(2^(1/omega) -1)))))/sigma)^2))^ omega)      
    }   
  }  
  return(new("puk",.Data=rval,kpar=list(sigma=sigma,omega=omega)))  
}
setClass("puk",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))
############################################################
