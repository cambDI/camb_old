#################################################################################
## Visualization
#################################################################################

#' Plot Distribution of the Response Variable
#' 
#' This function creates plots the distribution of the response variable.
#' 
#' @param Data A numeric vector. TBD: explain more
#' @param xlab Title of the x axis.
#' @param ylab Title of the y axis.
#' @param main Title of the plot.
#' @param alpha Alpha for the fill color of the distribution. Default value 0.2.
#' @param binwidth Width of the histogram bins. Default value NULL.
#' @param histFill Fill color of the histogram bars. Default value 'white'.
#' @param histCol Color of the histogram lines. Default value 'black'.
#' @param densityFill Fill color of the distribution. Default value "#FF6666".
#' @param TitleSize Title font size. Default value 15.
#' @param TextSize Text font size. Default value 15.
#' @param XAxisSize Size of the text on the X axis. Default value 15.
#' @param YAxisSize Size of the text on the Y axis. Default value 15.
#' @param AngleLab Angle of the labels in the X axis. Default value 30. 
#' @param LegendPosition Position of the legend. Default value 'right'.
#' @param TitleAxisSize Font size of the axes lables. Default value 15.
#' @param tmar Top margin size. Default values is 1.
#' @param bmar ottom margin size. Default values is 1.
#' @param rmar Right margin size. Default values is 1.
#' @param lmar Left margin size. Default values is 1. 
#' @details Additional ggplot2 layers can be added with "+".
#' @export
#' @return Returns a ggplot object.
#' @author Isidro Cortes <isidrolauscher@@gmail.com> and Daniel Murrell <dsmurrell@@gmail.com>
DensityResponse <- function(Data,xlab="",ylab="",main="",alpha=0.2,
                            binwidth=NULL,histFill="white",histCol="black",
                            densityFill="#FF6666",TitleSize=15,TextSize=15,
                            XAxisSize=15,YAxisSize=15,AngleLab=30,LegendPosition="right",
                            TitleAxesSize=15,tmar=1,bmar=1,rmar=1,lmar=1){
  if (!is.vector(Data)) stop("Input data must be a numeric vector")
  Data <- data.frame(Values=Data)
  if (is.null(binwidth)) binwidth <- abs(range(Data)[1] - range(Data)[2]) / 20
  p <- ggplot(Data, aes(x=Values)) + theme_bw() + 
    geom_histogram(aes(y=..density..),binwidth=binwidth,colour=histCol, fill=histFill) + 
    geom_density(alpha=alpha, fill=densityFill)+ ylab(ylab) + xlab(xlab) + ggtitle(main) +
    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),axis.text.y=element_text(size=YAxisSize),
          legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar), "cm"))
return(p)
}

#################################################################################
## PCA analysis 

isnot.null <- function(x) ! is.null(x)
isnot.vector <- function(x) ! is.vector(x)

##############
PCA <- function (Data, RowNames = NULL,cor=TRUE, scale = TRUE, center = TRUE,...) {
  if (is.matrix(Data) || is.data.frame(Data)) {
      ana <- prcomp(t(Data), cor = cor, scale = scale, center = center,...)
      PC1 <- as.vector(ana$rotation[, 1])
      PC2 <- as.vector(ana$rotation[, 2])
      if (is.null(RowNames)) {
        Data <- data.frame(PC1, PC2)
      out <- list(Data=Data,PCs_all=ana$rot,Std=ana$sdev,Info=summary(ana))
      }
      else {
        Data <- data.frame(PC1, PC2, RowNames)
      out <- list(Data=Data,PCs_all=ana$rot,Std=ana$sdev,Info=summary(ana))
      }
      return(out)
    }
  else {
      stop("Input data must be a numeric matrix or data.frame")
    }
}

##############
# Plot the towo first PC of the sequence descriptors
PCAPlot <- function (Data,main="",ylab="PC2",xlab="PC1",labels=NULL,PointSize=4,
                         LegendPosition="right",LegendName="",ColLegend=1,RowLegend=NULL,
                         TitleSize=15,TextSize=15,XAxisSize=15,YAxisSize=15,AngleLab=30,
                         TitleAxesSize=15,LegendTitleSize=15,LegendTextSize=15,tmar=1,bmar=1,rmar=1,lmar=1) 
{
	isnot.null <- function(x) ! is.null(x)
	isnot.vector <- function(x) ! is.vector(x)
  if (length(names(Data)) < 2 || length(names(Data)) > 3){
    stop("Two PCs required. The Data.frame provided has less than two columns (PCs) or more than 3")
  } else if (names(Data)[1] != "PC1" || names(Data)[2] != "PC2"){
    stop("Column names have to be {PC1, PC2}")
  } else if (length(names(Data)) == 2 && is.null(labels)){
    print("No names provided")
    p <- ggplot(Data, aes(x=PC1, y=PC2)) +
      geom_point(size=PointSize) + theme_bw() + ggtitle(main) + ylab(ylab) + xlab(xlab) +
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,
            plot.title=element_text(size=TitleSize),legend.key=element_blank(),
            legend.text = element_text(size=LegendTextSize),legend.title= element_text(size=LegendTitleSize),
            plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
      guides(colour = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend),
             shape = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend))
  } else if (length(names(Data)) == 2 && isnot.null(labels)){
    if (length(labels) != nrow(Data) || isnot.vector(labels)) {
      stop("Either the names are not in a vector, or its length is not equal to the number of datapoints (rows of the input data)")
    } else {
      print("Names provided")
      Data <- data.frame(Data,labels=labels)
      a <- length(unique(labels))
      p <- ggplot(Data, aes(x=PC1, y=PC2, color=labels,shape=labels)) +
        geom_point(size=PointSize)  + scale_shape_manual(values=1:a) + theme_bw() + ggtitle(main) + ylab(ylab) + xlab(xlab) +
        theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
              axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
              axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,
              legend.text = element_text(size=LegendTextSize),legend.title= element_text(size=LegendTitleSize),
              plot.title=element_text(size=TitleSize),legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
        guides(colour = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend),
               shape = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend))
    } 
  } else {
    print("Names provided in the third column of the data.frame")
    labels <- unlist(Data[names(Data)[3]])
    a <- length(unique(labels))
    namee <- names(Data)[3]
    p <- ggplot(Data, aes_string(x="PC1", y="PC2",colour=namee,shape=namee)) +
      geom_point(size=PointSize) +scale_shape_manual(values=1:a) + theme_bw() + ggtitle(main) + ylab(ylab) + xlab(xlab) +
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
            legend.text = element_text(size=LegendTextSize),legend.title= element_text(size=LegendTitleSize),
            legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
      guides(colour = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend), shape = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend))
  }
  return(p) 
}

##############
PairwiseDist <- function(Data,method="jaccard",...){
  if (is.matrix(Data) || is.data.frame(Data)){
    Data <- unique(Data)
    methods <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard",
                 "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial",
                 "chao", "cao")
    method <- match.arg(method,methods)
    pwdist <- vegdist(Data, method = method,...)
    pwdist <- data.frame(as.vector(pwdist))
    names(pwdist) <- "Distance"
    return(pwdist)
  } else {
    stop("Input data must be a numeric matrix or data.frame")
  }
}

##############
PairwiseDistPlot <- function(Data,xlab="",ylab="",main="",TextSize=15,TitleSize=15,XAxisSize=15,YAxisSize=15,
                             TitleAxesSize=15,tmar=1,bmar=1,rmar=1,lmar=1,AngleLab=30,
                             binwidth=NULL,fillCol="white",Colour="black",DensityFill="#FF6666",DensityAlpha=.2){
  if (is.matrix(Data) || is.data.frame(Data)){
    if (is.null(binwidth)) binwidth <- abs(range(Data)[1] - range(Data)[2]) / 20
    p <- ggplot(Data, aes(x=Distance)) + theme_bw() +
      geom_histogram(aes(y=..density..),binwidth=binwidth,colour=Colour, fill=fillCol) + 
      geom_density(alpha=DensityAlpha, fill=DensityFill)+ ylab(ylab) + xlab(xlab) + ggtitle(main) +
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),plot.title=element_text(size=TitleSize),
            legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) 
    return(p)
  } else {
    stop("Input data must be a numeric matrix or data.frame")
  }
}

##############
##############
## Maximum Model Performance
#
#MaxPerf <- function(meanNoise=0,sdNoise,meanResp,sdResp,lenPred,iters=1000,
#                    filename=NULL,pdfW=10,pdfH=10,TextSize=15,TitleSize=15,
#                    XAxisSize=15,YAxisSize=15,TitleAxesSize=15,tmar=1,bmar=1,
#                    rmar=1,lmar=1,AngleLab=30,LegendPosition="right"){
#  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#  R2 <- c()
#  R02 <- c()
#  Q2 <- c()
#  rmsep <- c()
#  for (i in 1:iters){
#    x <- rnorm(lenPred,mean=meanResp,sd=sdResp)
#    noise <- rnorm(length(x),mean=meanNoise,sd=sdNoise)
#    y <- x+noise
#    R2[i] <- Rsquared(y,x)
#    Q2[i] <- Qsquared2(y,x)
#    R02[i] <- Rsquared0(y,x)
#    rmsep[i] <- RMSE(y,x)
#  }
#  
#  pp <- data.frame(R2=R2)
#  title <- expression(paste("R"^{2}))
#  p1 <- ggplot(pp, aes(x=R2)) + 
#    theme_bw() + 
#    geom_density(alpha=.2, fill="#FF6666") + ggtitle(title)+ aes(y = ..count..)+ 
#    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
#          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
#          axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
#          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
#    geom_vline(aes(xintercept=mean(R2, na.rm=T)), color="blue", linetype="dashed", size=1)
#  
#  pp <- data.frame(rmsep=rmsep)
#  p2 <- ggplot(pp, aes(x=rmsep)) + 
#    theme_bw() + 
#    geom_density(alpha=.2, fill="#FF6666") + ggtitle('RMSEP') +  aes(y = ..count..)+ 
#    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
#          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
#          axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
#          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
#    geom_vline(aes(xintercept=mean(rmsep, na.rm=T)), color="blue", linetype="dashed", size=1)
#  
#  pp <- data.frame(R02=R02)
#  title <- expression(paste("R"[0]^{2}))
#  
#  p3 <- ggplot(pp, aes(x=R02)) + 
#    theme_bw() + 
#    geom_density(alpha=.2, fill="#FF6666") + ggtitle(title) + aes(y = ..count..)+ 
#    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
#          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
#          axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
#          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
#    geom_vline(aes(xintercept=mean(R02, na.rm=T)), color="blue", linetype="dashed", size=1)
#  
#  title <- expression(paste("Q"^{2}))
#  pp <- data.frame(Q2=Q2)
#  p4 <- ggplot(pp, aes(x=Q2)) + 
#    theme_bw() + 
#    geom_density(alpha=.2, fill="#FF6666") + ggtitle(title) + aes(y = ..count..)+ 
#    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
#          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
#          axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
#          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
#    geom_vline(aes(xintercept=mean(Q2, na.rm=T)), color="blue", linetype="dashed", size=1)
#  
#  if (isnot.null(filename)){
#    pdfname=paste(filename,".pdf",sep="")
#    pdf(pdfname,width=pdfW,height=pdfH)
#    grid.newpage()
#    pushViewport(viewport(layout = grid.layout(2, 2)))
#    print(p1, vp = vplayout(1,1))
#    print(p2, vp = vplayout(1,2))
#    print(p3, vp = vplayout(2,1))
#    print(p4, vp = vplayout(2,2))
#    dev.off()
#  } #else {
#    #grid.newpage()
#    #pushViewport(viewport(layout = grid.layout(2, 2)))
#    #print(p1, vp = vplayout(1,1))
#    #print(p2, vp = vplayout(1,2))
#    #print(p3, vp = vplayout(2,1))
#    #print(p4, vp = vplayout(2,2))
#  #}
#  p <- list()
#  p$p1 <- p1
#  p$p2 <- p2
#  p$p3 <- p3
#  p$p4 <- p4
#  return(p)
#}
#
##########
# Get Legend
GetLegend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
