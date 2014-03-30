\name{GetPropertiesSDF}
\alias{GetPropertiesSDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Retrieval of Properties from .sdf file
%%  ~~function to do ... ~~
}
\description{
	The function retrieves the values of all fields in a .sdf file for a user-defined
	number of molecules.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetPropertiesSDF(structures.file, number_processed = -1, type = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structures.file}{
	  File in .sdf format with the compounds.
%%     ~~Describe \code{structures.file} here~~
}
  \item{number_processed}{
	  Number of molecules from which the properties are to be extracted.
	  The function starts at the beginning of the .sdf file till the molecule 
	  which ordinal position is equal to the argument.
	  The default value is -1. In that case, the properties are extracted for all molecules in the .sdf file.
%%     ~~Describe \code{number_processed} here~~
}
  \item{type}{
	  1 corresponds to .sdf format, which is the only format currentyl supported.
%%     ~~Describe \code{type} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A data.frame with all properties which rows are indexed by the compounds in the .sdf file 
	and columns by the properties.
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
	GetPropertySDF
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
\keyword{ ~GetPropertiesSDF }
