\name{YScrambling}
\alias{YScrambling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Randomization of the Response Variable
%%  ~~function to do ... ~~
}
\description{
	The function randomizes a user defined percentage of the response variable.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
YScrambling(y, percent)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	  A numeric vector with the response varible.
%%     ~~Describe \code{y} here~~
}
  \item{percent}{
	  The fraction of the response variable to be randomized.
	  0 equals 0\%, whereas 1 equals 100\%.
%%     ~~Describe \code{percent} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A numeric vector with the response variable randimized as indicated by the user.
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
# 0% randomized
YScrambling(seq(1,10),0)
# 50% randomized
YScrambling(seq(1,10),0.5)
# 100% randomized
YScrambling(seq(1,10),1)
}
\keyword{ ~YScrambling }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
