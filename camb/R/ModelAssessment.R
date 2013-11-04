#################################################################################
## Model Assessment and Results Visualization
#################################################################################

##############
ErrorBarplot <- function(Data,X,Y,std,colour=NULL,shape=NULL,fill=NULL,main="",ylab="",xlab="",
                         minn=NULL,maxx=NULL,TextSize=15,TitleSize=15,XAxisSize=15,YAxisSize=15,
                         TitleAxesSize=15,AngleLab=35,barcol="red",
                         barWidth=0.3, LegendName="Legend",ColLegend=1,
                         RowLegend=NULL,LegendPosition="right",
                         tmar=1,bmar=1,rmar=1,lmar=1,stat="identity"){
  low <- Data[,Y] - Data[,std]
  high <- Data[,Y] + Data[,std]
  err <- data.frame(low,high)
  err_rang <- range()
  if (is.null(minn)) { minn <- min(err) - 0.20*diff((range(err)))}
  if (is.null(maxx)) { maxx <- max(err) + 0.20*diff((range(err)))}
  
  p <- ggplot(Data, aes_string(x=X, y=Y,fill=fill,colour=colour,shape=shape)) + theme_bw() + 
    ggtitle(main)+ ylab(ylab) + xlab(xlab) +   
    coord_cartesian(ylim=c(minn,maxx))+
    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
          axis.text.y=element_text(size=YAxisSize),
          legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar), "cm")) +
    geom_bar(position="dodge",stat=stat) + 
    geom_errorbar(data=err, aes(ymin=low, ymax=high),
                  position=position_dodge(0.9), width=barWidth,color=barcol) +   
    guides(colour = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend),
           shape=guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend),
           fill=guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend)) 
  return(p)
}

##############
plotGrid <- function(plots,MyLegend,RowLegend,heightPlots,heightLegend,filename=NULL,height=10,width=10){
  if (length(plots) < RowLegend) {warning("WARNING: More rows than plots! An optimal number of rows should be equal or smaller to the number of plots")}
  if (is.null(filename)){
    t <- c("grid.arrange(arrangeGrob(")
    for (i in 1:length(plots)){
      t <- paste(t,plots[i]," + theme(legend.position='none'),",sep="")
    }
    t <- paste(t,"nrow=",RowLegend,"),arrangeGrob(MyLegend,nrow=1),heights=c(",heightPlots,",",heightLegend,"))" )
    return(eval(parse(text=t)))
  } else {
    t <- c("grid.arrange(arrangeGrob(")
    for (i in 1:length(plots)){
      t <- paste(t,plots[i]," + theme(legend.position='none'),",sep="")
    }
    t <- paste(t,"nrow=",RowLegend,"),arrangeGrob(MyLegend,nrow=1),heights=c(",heightPlots,",",heightLegend,"))" )
    pdfname <- paste(filename,".pdf",sep="")
    pdf(file=pdfname,width=width,height=height)
    eval(parse(text=t))
    dev.off()
    return(1) 
  }
}

##############
## Validation of models
Validation <- function(pred,obs){
  if (is.vector(pred) && is.vector(obs) && length(pred)==length(obs)){
    metrics <- list(R2 = Rsquared(pred,obs), R02 = Rsquared0(pred,obs), Q2 = Qsquared(pred,obs), RMSEP = RMSE(pred,obs), Slope=slope(pred,obs))
  } else {
    stop("Wrong input: input arguments are not vector or have unequal length")
  }
  return(metrics)
}

##############
## Observed vs predicted
ObsPred <- function (pred,obs,margin=NULL,main="",ylab="Observed",xlab="Predicted",
                   PointSize=4,ColMargin="blue",TextSize=15,TitleSize=15,
                   XAxisSize=15,YAxisSize=15,TitleAxesSize=15,tmar=1,bmar=1,
                   rmar=1,lmar=1) 
{ 
  if (isnot.vector(obs) || isnot.vector(pred)){
    stop("The input data must be two vectors")
  } else if ( length(obs) != length(pred) ){
    stop("Both vectors have to be of equal length")
  } else if(isnot.null(margin)) {
    Data <- data.frame(Observed=obs,Predicted=pred)
    p <- ggplot(Data, aes(x=Observed, y=Predicted)) +
      geom_point(size=PointSize) +  geom_abline(slope=1,intercept=margin/2,colour=ColMargin) + 
      geom_abline(slope=1,intercept=-(margin/2),colour=ColMargin) + theme_bw() + 
      ggtitle(main) + ylab(ylab) + xlab(xlab)+
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
            legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) 
  } else {
    Data <- data.frame(Observed=obs,Predicted=pred)
    p <- ggplot(Data, aes(x=Observed, y=Predicted)) + geom_point(size=PointSize) + theme_bw() + 
      ggtitle(main) + ylab(ylab) + xlab(xlab) +
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
            legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm"))   
  }
  return(p) 
}
