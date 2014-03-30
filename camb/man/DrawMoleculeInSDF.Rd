\name{DrawMoleculeInSDF}
\alias{DrawMoleculeInSDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Draw a molecule from 
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DrawMoleculeInSDF(structures.file, structure.number, file.name, useNameAsTitle)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structures.file}{
%%     ~~Describe \code{structures.file} here~~
}
  \item{structure.number}{
%%     ~~Describe \code{structure.number} here~~
}
  \item{file.name}{
%%     ~~Describe \code{file.name} here~~
}
  \item{useNameAsTitle}{
%%     ~~Describe \code{useNameAsTitle} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (structures.file, structure.number, file.name, useNameAsTitle) 
{
    if (file.info(structures.file)$size == 0) {
        stop("Input file is empty")
    }
    print(structure.number)
    .C("R_drawMoleculeInSDF", structures.file, as.integer(structure.number), 
        file.name, useNameAsTitle)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
