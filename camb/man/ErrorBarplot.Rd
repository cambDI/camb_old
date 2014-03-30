\name{ErrorBarplot2}
\alias{ErrorBarplot2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Barplot with error bars.
%%  ~~function to do ... ~~
}
\description{
	The function creates easily customizable barplots with error bars.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ErrorBarplot2(X, Y, err, fill = as.factor(rep(4, length(X))), main = "", ylab = "", xlab = "", minn = NULL, maxx = NULL, TextSize = 15, TitleSize = 15, XAxisSize = 15, YAxisSize = 15, TitleAxesSize = 15, AngleLab = 35, barcol = "red", barSize = 1, barWidth = 0.3, LegendName = "Legend", ColLegend = 1, RowLegend = NULL, LegendPosition = "right", tmar = 1, bmar = 1, rmar = 1, lmar = 1, stat = "identity")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	  A vector containing the aesthetics corresponding to the X axis.

%%     ~~Describe \code{X} here~~
}
  \item{Y}{
	  The values for the ordenate axis.
%%     ~~Describe \code{Y} here~~
}
  \item{err}{
%%     ~~Describe \code{err} here~~
}
  \item{fill}{
%%     ~~Describe \code{fill} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{minn}{
%%     ~~Describe \code{minn} here~~
}
  \item{maxx}{
%%     ~~Describe \code{maxx} here~~
}
  \item{TextSize}{
%%     ~~Describe \code{TextSize} here~~
}
  \item{TitleSize}{
%%     ~~Describe \code{TitleSize} here~~
}
  \item{XAxisSize}{
%%     ~~Describe \code{XAxisSize} here~~
}
  \item{YAxisSize}{
%%     ~~Describe \code{YAxisSize} here~~
}
  \item{TitleAxesSize}{
%%     ~~Describe \code{TitleAxesSize} here~~
}
  \item{AngleLab}{
%%     ~~Describe \code{AngleLab} here~~
}
  \item{barcol}{
%%     ~~Describe \code{barcol} here~~
}
  \item{barSize}{
%%     ~~Describe \code{barSize} here~~
}
  \item{barWidth}{
%%     ~~Describe \code{barWidth} here~~
}
  \item{LegendName}{
%%     ~~Describe \code{LegendName} here~~
}
  \item{ColLegend}{
%%     ~~Describe \code{ColLegend} here~~
}
  \item{RowLegend}{
%%     ~~Describe \code{RowLegend} here~~
}
  \item{LegendPosition}{
%%     ~~Describe \code{LegendPosition} here~~
}
  \item{tmar}{
%%     ~~Describe \code{tmar} here~~
}
  \item{bmar}{
%%     ~~Describe \code{bmar} here~~
}
  \item{rmar}{
%%     ~~Describe \code{rmar} here~~
}
  \item{lmar}{
%%     ~~Describe \code{lmar} here~~
}
  \item{stat}{
%%     ~~Describe \code{stat} here~~
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
function (X, Y, err, fill = as.factor(rep(4, length(X))), main = "", 
    ylab = "", xlab = "", minn = NULL, maxx = NULL, TextSize = 15, 
    TitleSize = 15, XAxisSize = 15, YAxisSize = 15, TitleAxesSize = 15, 
    AngleLab = 35, barcol = "red", barSize = 1, barWidth = 0.3, 
    LegendName = "Legend", ColLegend = 1, RowLegend = NULL, LegendPosition = "right", 
    tmar = 1, bmar = 1, rmar = 1, lmar = 1, stat = "identity") 
{
    data <- data.frame(X = X, Y = Y, err = err, Group = fill)
    errbar <- aes(ymin = Y - err, ymax = Y + err)
    p <- ggplot(data, aes(x = X, y = Y, fill = Group)) + theme_bw() + 
        geom_bar(position = "dodge", stat = stat) + geom_errorbar(mapping = errbar, 
        position = position_dodge(0.9), width = barWidth, color = barcol, 
        size = barSize) + theme(text = element_text(size = TextSize), 
        axis.text.x = element_text(size = XAxisSize, angle = AngleLab, 
            hjust = 1), axis.title.x = element_text(size = TitleAxesSize), 
        axis.title.y = element_text(size = TitleAxesSize), axis.text.y = element_text(size = YAxisSize), 
        legend.position = LegendPosition, plot.title = element_text(size = TitleSize), 
        legend.key = element_blank(), plot.margin = unit(c(tmar, 
            rmar, bmar, lmar), "cm"))
    return(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
