\name{PlotMolecules}
\alias{PlotMolecules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Compounds from a .sdf File.
%%  ~~function to do ... ~~
}
\description{
	The function plots the chemical structures provided in a .sdf file.
	The plots can also be sent to a 2x2 grid in a .pdf file.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PlotMolecules(sdf.file, IDs, pdf.file = NULL, PDFMain = NULL, useNameAsTitle = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdf.file}{
	  The .sdf file with the molecules.
%%     ~~Describe \code{sdf.file} here~~
}
  \item{IDs}{
	  The IDs of the molecules to be depicted (the ordinal position of the molecules in the .sdf file). 
	  Currently, a maximum of four IDs is supported.
%%     ~~Describe \code{IDs} here~~
}
  \item{pdf.file}{
	  If not NULL, the .pdf where the molecules will be depicted.
%%     ~~Describe \code{pdf.file} here~~
}
  \item{PDFMain}{
If not NULL, the title of the molecule depiction in the .pdf file.
%%     ~~Describe \code{PDFMain} here~~
}
  \item{useNameAsTitle}{
	  If TRUE, the names of the molecules as especified in the .sdf file are used as molecules names in the depiction.
%%     ~~Describe \code{useNameAsTitle} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A list with the plots of the molecules.
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
\keyword{ ~PlotMolecules }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
