\name{GetPropertySDF}
\alias{GetPropertySDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Retrieval of a given Property from a .sdf File.
%%  ~~function to do ... ~~
}
\description{
	The function extracts the values of one .sdf field for a user-defined number of molecules.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetPropertySDF(structures.file, property = "", number_processed = -1, type = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structures.file}{
	  The .sdf file containing the molecules.
%%     ~~Describe \code{structures.file} here~~
}
  \item{property}{
	  The property which values are to be extracted.
%%     ~~Describe \code{property} here~~
}
  \item{number_processed}{
	  Number of molecules for which the value are to be extracted.
If this argument is set to -1 (default) the property will be extracted for all compounds.
%%     ~~Describe \code{number_processed} here~~
}
  \item{type}{
	  -1 corresponds to .sdf file format, which is the only format currently supported.
%%     ~~Describe \code{type} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A one-field data.frame with the values of the property.
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
	GetPropertiesSDF
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
\keyword{ ~GetPropertySDF}
