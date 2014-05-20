\name{Validation}
\alias{Validation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Statistical Metrics to Assess Model Peformance
%%  ~~function to do ... ~~
}
\description{
	The function calculates correlation and error metrics between two 
	numeric vectors.
	These metrics are used to evaluate the performance of predictive models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Validation(pred, obs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred}{
	  A numeric vector with the predicted values.
%%     ~~Describe \code{pred} here~~
}
  \item{obs}{
	  A numeric vector with the observed values.
%%     ~~Describe \code{obs} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A list containing the values for the followin metrics:
	(i) R2, (ii) R02, (iii) three formualtions of Q2 (see the PCM camb tutorial),
	(iv) RMSE, (v) Slope, and (vi) MAE.
	Further details can be found in the references below.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
	\cite{beware}
	\cite{earnest}
	\cite{cortesGP}
	\cite{consonni}
%% ~put references to the literature/web site here ~
}
\author{
	Isidro Cortes and Daniel Murrell.
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

Validation(pred=seq(1,10),obs=seq(1,10))

}
\keyword{ ~Validation }
