\name{CorrelationPlot}
\alias{CorrelationPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scatterplot of the Observed against the Predicted Values.
%%  ~~function to do ... ~~
}
\description{
The function depicted a scatterplot of the observed against the predicted values with a machine learning model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CorrelationPlot(pred, obs, margin = NULL, main = "", ylab = "Predicted", xlab = "Observed", PointSize = 4, ColMargin = "blue", TextSize = 15, TitleSize = 15, XAxisSize = 15, YAxisSize = 15, TitleAxesSize = 15, tmar = 1, bmar = 1, rmar = 1, lmar = 1, AngleLab = 30, LegendPosition = "right", PointColor = "black", PointAlpha = 1, PointShape = 16, MarginWidth = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred}{
  Predicted values.
%%     ~~Describe \code{pred} here~~
}
  \item{obs}{
  Observed values.
%%     ~~Describe \code{obs} here~~
}
  \item{margin}{
  Bioactivity margin centered in the diagonal of the correlation plot. 
Default value NULL.
%%     ~~Describe \code{margin} here~~
}
  \item{main}{
  Plot title.
%%     ~~Describe \code{main} here~~
}
  \item{ylab}{
  Title of the Y axis.
%%     ~~Describe \code{ylab} here~~
}
  \item{xlab}{
  Title of the X axis.
%%     ~~Describe \code{xlab} here~~
}
  \item{PointSize}{
  Size of the points.
%%     ~~Describe \code{PointSize} here~~
}
  \item{ColMargin}{
  Color of the bioactivity margin ('margin').
%%     ~~Describe \code{ColMargin} here~~
}
  \item{TextSize}{
  Text font size. Default value 15.
%%     ~~Describe \code{TextSize} here~~
}
  \item{TitleSize}{
  Title font size. Default value 15.
%%     ~~Describe \code{TitleSize} here~~
}
  \item{XAxisSize}{
  Size of the text on the X axis. Default value 15.
%%     ~~Describe \code{XAxisSize} here~~
}
  \item{YAxisSize}{
  Size of the text on the Y axis. Default value 15.
%%     ~~Describe \code{YAxisSize} here~~
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
  \item{AngleLab}{
  Angle of the labels in the X axis. Default value 30.
%%     ~~Describe \code{AngleLab} here~~
}
  \item{LegendPosition}{
  Position of the legend. Default value 'right'.
%%     ~~Describe \code{LegendPosition} here~~
}
  \item{PointColor}{
  Color of the points in the scatterplot. Default value 'black'.
%%     ~~Describe \code{PointColor} here~~
}
  \item{PointAlpha}{
  Color alpha of the points in the scatterplot. Default value 1.
%%     ~~Describe \code{PointAlpha} here~~
}
  \item{PointShape}{
  Shape of the points in the scatterplot. Default value 16.
%%     ~~Describe \code{PointShape} here~~
}
  \item{MarginWidth}{
%%     ~~Describe \code{MarginWidth} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list (ggplot2 plot) with the scatterplot.
Further layers can be added to the plot with the ggplot function "+".
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
CorrelationPlot(pred=seq(1,10)+rnorm(10), obs = seq(1,10),
margin = NULL, main = "", ylab = "Predicted", 
    xlab = "Observed", PointSize = 4, ColMargin = "blue", TextSize = 15, 
    TitleSize = 15, XAxisSize = 15, YAxisSize = 15, TitleAxesSize = 15, 
    tmar = 1, bmar = 1, rmar = 1, lmar = 1, AngleLab = 30, LegendPosition = "right", 
    PointColor = "black", PointAlpha = 1, PointShape = 16, MarginWidth = 1) 
}
\keyword{ ~CorrelationPlot }
