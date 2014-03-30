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
	  Alpha for the fill color of the distribution.
%%     ~~Describe \code{alpha} here~~
}
  \item{binwidth}{
	  Width of the histogram bins.
%%     ~~Describe \code{binwidth} here~~
}
  \item{histFill}{
	  Fill color of the histogram bars.
%%     ~~Describe \code{histFill} here~~
}
  \item{histCol}{
	  Color of the histogram lines.
%%     ~~Describe \code{histCol} here~~
}
  \item{densityFill}{
	  Fill color of the distribution.
%%     ~~Describe \code{densityFill} here~~
}
  \item{TitleSize}{
	  Title font size.
%%     ~~Describe \code{TitleSize} here~~
}
  \item{TextSize}{
	  Text size.
%%     ~~Describe \code{TextSize} here~~
}
  \item{XAxisSize}{
	  Size of the text on the X axis.
%%     ~~Describe \code{XAxisSize} here~~
}
  \item{YAxisSize}{
	  Size of the text on the Y axis.
%%     ~~Describe \code{YAxisSize} here~~
}
  \item{AngleLab}{
	  Angle of the labels in the X axis.
%%     ~~Describe \code{AngleLab} here~~
}

  \item{LegendPosition}{
	  Position of the legend.
%%     ~~Describe \code{LegendPosition} here~~
}
  \item{TitleAxesSize}{
	  Font size of the axes lables.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Data, xlab = "", ylab = "", main = "", alpha = 0.2, 
    binwidth = NULL, histFill = "white", histCol = "black", densityFill = "#FF6666", 
    TitleSize = 15, TextSize = 15, XAxisSize = 15, YAxisSize = 15, 
    AngleLab = 30, LegendPosition = "right", TitleAxesSize = 15, 
    tmar = 1, bmar = 1, rmar = 1, lmar = 1) 
{
    if (!is.vector(Data)) 
        stop("Input data must be a numeric matrix or data.frame")
    Data <- data.frame(Values = Data)
    if (is.null(binwidth)) 
        binwidth <- abs(range(Data)[1] - range(Data)[2])/20
    p <- ggplot(Data, aes(x = Values)) + theme_bw() + geom_histogram(aes(y = ..density..), 
        binwidth = binwidth, colour = histCol, fill = histFill) + 
        geom_density(alpha = alpha, fill = densityFill) + ylab(ylab) + 
        xlab(xlab) + ggtitle(main) + theme(text = element_text(size = TextSize), 
        axis.text.x = element_text(size = XAxisSize, angle = AngleLab, 
            hjust = 1), axis.title.x = element_text(size = TitleAxesSize), 
        axis.title.y = element_text(size = TitleAxesSize), axis.text.y = element_text(size = YAxisSize), 
        legend.position = LegendPosition, plot.title = element_text(size = TitleSize), 
        legend.key = element_blank(), plot.margin = unit(c(tmar, 
            rmar, bmar, lmar), "cm"))
    return(p)
  }
}

## Example:

DensityResponse(rnorm(100),xlab = "Random values", ylab = "Frequency", main = "Exampl DensityResponse function", alpha = 0.2, 
    binwidth = NULL, histFill = "white", histCol = "black", densityFill = "#FF6666", 
	TitleSize = 15, TextSize = 15, XAxisSize = 15, YAxisSize = 15, 
	AngleLab = 30, LegendPosition = "right", TitleAxesSize = 15, 
	tmar = 1, bmar = 1, rmar = 1, lmar = 1
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
