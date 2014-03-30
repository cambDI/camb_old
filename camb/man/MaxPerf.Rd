\name{MaxPerf}
\alias{MaxPerf}
\title{Distribution of Maximum Theoretical Values Achievable given the Dataset and its Uncertainty
}
\description{
	Calculates the ditribution of model validation metrics that are achievable given the size of the dataset,
	the uncertainty in the response variable, and the distribution of the responsible variable quantified by its mean and standard deviation.
	Therefore, these distributions help to assess models overfitting; e.g. a model trained on a dataset with high uncertainty exhibiting
	high correlation values might be overoptimistic.
}
\usage{
MaxPerf(meanNoise = 0, sdNoise, meanResp, sdResp, lenPred, iters = 1000, filename = NULL, pdfW = 10, pdfH = 10, TextSize = 15, TitleSize = 15, XAxisSize = 15, YAxisSize = 15, TitleAxesSize = 15, tmar = 1, bmar = 1, rmar = 1, lmar = 1, AngleLab = 30, LegendPosition = "right")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meanNoise}{
Mean value of the noise in the data. Default value 0.
}
  \item{sdNoise}{
Standard deviation of the noise in the data. See the work by Kramer et al. about uncertainty in public bioactivity databases.
}
  \item{meanResp}{
Mean of the response variable.
}
  \item{sdResp}{
Standard deviation of the response variable.
}
  \item{lenPred}{
Number of datapoints of the external (hold-out) set.
%%     ~~Describe \code{lenPred} here~~
}
  \item{iters}{
	  Number of iterations. Default value 1000.
%%     ~~Describe \code{iters} here~~
}
  \item{filename}{
If not NULL, file where the plot will be saved. Default value NULL.
%%     ~~Describe \code{filename} here~~
}
  \item{pdfW}{
Width of the .pdf file, in centimeters, where the plot will be saved. Default value 10.
%%     ~~Describe \code{pdfW} here~~
}
  \item{pdfH}{
Height of the .pdf file, in centimeters, where the plot will be saved. Default value 10.
%%     ~~Describe \code{pdfH} here~~
}
  \item{TextSize}{
Fontsize of the text in the plot. Default value 15.
%%     ~~Describe \code{TextSize} here~~
}
  \item{TitleSize}{
Fontsize of the title. Default value 15.
%%     ~~Describe \code{TitleSize} here~~
}
  \item{XAxisSize}{
%%     ~~Describe \code{XAxisSize} here~~
Fontsize of the X axis. Default value 15.
}
  \item{YAxisSize}{
%%     ~~Describe \code{YAxisSize} here~~
Fontsize of the Y axis. Default value 15.
}
  \item{TitleAxesSize}{
%%     ~~Describe \code{TitleAxesSize} here~~
Fontsize of the axes titles. Default value 15.
}
  \item{tmar}{
	  Top margin size. Default value 1.
%%     ~~Describe \code{tmar} here~~
}
  \item{bmar}{
	  Bottom margin size. Default value 1.
%%     ~~Describe \code{bmar} here~~
}
  \item{rmar}{
	  Right margin size. Default value 1.
%%     ~~Describe \code{rmar} here~~
}
  \item{lmar}{
	  Left margin size. Default value 1.
%%     ~~Describe \code{lmar} here~~
}
  \item{AngleLab}{
	  Angle of the labels of the X axis. Default value 30.
%%     ~~Describe \code{AngleLab} here~~
}
  \item{LegendPosition}{
	  Position of the legend. Default value 'right'.
%%     ~~Describe \code{LegendPosition} here~~
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
Cortes-Ciriano et al. 'Proteochemometric Modeling in a Bayesian Framework'. 2014
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
\keyword{ ~MaxPerf }
\keyword{ ~MaximumPerformance }
