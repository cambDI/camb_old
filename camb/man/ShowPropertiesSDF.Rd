\name{ShowPropertiesSDF}
\alias{ShowPropertiesSDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{List Properties in a .sdf file.
%%  ~~function to do ... ~~
}
\description{
	The function returns a vector with the different fields found in a .sdf file.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ShowPropertiesSDF(structures.file, type = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structures.file}{
File containing the molecules in .sdf format.
%%     ~~Describe \code{structures.file} here~~
}

  \item{type}{
	  Format type. Default value is one, which corresponds to .sdf format.
      The only format available currently is .sdf.
%%     ~~Describe \code{type} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A vector with all the fields found in a .sdf file.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
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
	}
\keyword{ ~ShowPropertiesSDF }
