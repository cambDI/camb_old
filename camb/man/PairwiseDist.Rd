\name{PairwiseDist}
\alias{PairwiseDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pairwise Distance (Similarity) Matrix
%%  ~~function to do ... ~~
}
\description{
	The function is based on the vegdist function from the vegan package.
	It calculated the pairwise distance similarity matrix for all vectors
	input in a matrix or data.frame.
	The functions operates on a row basis.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PairwiseDist(Data, method = "jaccard", ..)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
	  A numeric data.frame or matrix containing compound, protein or amino acid descriptors (or any combination thereof).
%%     ~~Describe \code{Data} here~~
}
  \item{method}{
Available distance metrics are: "manhattan", "euclidean", "canberra", "bray", 
	            "kulczynski", "jaccard", "gower", "altGower", "morisita", 
				            "horn", "mountford", "raup", "binomial", "chao", 
							            "cao".
See the documentation of the R package vegan for details.
%%     ~~Describe \code{method} here~~
}
  \item{..}{
%%     ~~Describe \code{..} here~~
}
}
\details{
	For further detials see the documentation in the R package vegan.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A data.frame with the all pairwise distances.
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
	PairwiseDistPlot
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	m = matrix(abs(rnorm(20)),4,4) 
	mDist = PairwiseDist(m)
	head(mDist)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~PairwiseDist }
