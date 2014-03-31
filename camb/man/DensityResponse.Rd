\name{DensityResponse}
\alias{DensityResponse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Distribution of the Response Variable
%%  ~~function to do ... ~~
}
\description{
This function creates plots the distribution of the response variable.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DensityResponse(Data, xlab = "", ylab = "", main = "", alpha = 0.2, binwidth = NULL, histFill = "white", histCol = "black", densityFill = "#FF6666", TitleSize = 15, TextSize = 15, XAxisSize = 15, YAxisSize = 15, AngleLab = 30, LegendPosition = "right", TitleAxesSize = 15, tmar = 1, bmar = 1, rmar = 1, lmar = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
	 A numeric vector.
%%     ~~Describe \code{Data} here~~
}
  \item{xlab}{
	  Title of the x axis.
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
	  Title of the y axis.
%%     ~~Describe \code{ylab} here~~
}
  \item{main}{
	  Title of the plot.
%%     ~~Describe \code{main} here~~
}
  \item{alpha}{
	  Alpha for the fill color of the distribution. Default value 0.2.
%%     ~~Describe \code{alpha} here~~
}
  \item{binwidth}{
	  Width of the histogram bins. Default value NULL.
%%     ~~Describe \code{binwidth} here~~
}
  \item{histFill}{
	  Fill color of the histogram bars. Default value 'white'.
%%     ~~Describe \code{histFill} here~~
}
  \item{histCol}{
	  Color of the histogram lines. Default value 'black'.
%%     ~~Describe \code{histCol} here~~
}
  \item{densityFill}{
	  Fill color of the distribution. Default value "#FF6666".
%%     ~~Describe \code{densityFill} here~~
}
  \item{TitleSize}{
	  Title font size. Default value 15.
%%     ~~Describe \code{TitleSize} here~~
}
  \item{TextSize}{
	  Text font size. Default value 15.
%%     ~~Describe \code{TextSize} here~~
}
  \item{XAxisSize}{
	  Size of the text on the X axis. Default value 15.
%%     ~~Describe \code{XAxisSize} here~~
}
  \item{YAxisSize}{
	  Size of the text on the Y axis. Default value 15.
%%     ~~Describe \code{YAxisSize} here~~
}
  \item{AngleLab}{
	  Angle of the labels in the X axis. Default value 30. 
%%     ~~Describe \code{AngleLab} here~~
}

  \item{LegendPosition}{
	  Position of the legend. Default value 'right'.
%%     ~~Describe \code{LegendPosition} here~~
}
  \item{TitleAxesSize}{
	  Font size of the axes lables. Default value 15.
%%     ~~Describe \code{TitleAxesSize} here~~
}

  \item{tmar}{
	  Top margin size. Default values is 1.
%%     ~~Describe \code{tmar} here~~
}
  \item{bmar}{
	  Bottom margin size. Default values is 1.
%%     ~~Describe \code{bmar} here~~
}
  \item{rmar}{
	  Right margin size. Default values is 1.
%%     ~~Describe \code{rmar} here~~
}
  \item{lmar}{
	  Left margin size. Default values is 1. 
%%     ~~Describe \code{lmar} here~~
}
}
\details{
	Additional ggplot2 layers can be added with "+".
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	Returns a ggplot object.
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

## Example:

DensityResponse(rnorm(100),xlab = "Random values", ylab = "Frequency", main = "Exampl DensityResponse function", alpha = 0.2, 
    binwidth = NULL, histFill = "white", histCol = "black", densityFill = "#FF6666", 
	TitleSize = 15, TextSize = 15, XAxisSize = 15, YAxisSize = 15, 
	AngleLab = 30, LegendPosition = "right", TitleAxesSize = 15, 
	tmar = 1, bmar = 1, rmar = 1, lmar = 1
}
\keyword{ ~DensityResponse }
