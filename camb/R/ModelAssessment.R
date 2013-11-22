#################################################################################
## Model Assessment and Results Visualization
#################################################################################

##############
ErrorBarplot <- function(Data,X,Y,err,colour=NULL,shape=NULL,fill=NULL,main="",ylab="",xlab="",
                         minn=NULL,maxx=NULL,TextSize=15,TitleSize=15,XAxisSize=15,YAxisSize=15,
                         TitleAxesSize=15,AngleLab=35,barcol="red",barSize=1,
                         barWidth=0.3, LegendName="Legend",ColLegend=1,
                         RowLegend=NULL,LegendPosition="right",
                         tmar=1,bmar=1,rmar=1,lmar=1,stat="identity"){
  
  yerr_names <- names(Data)[c(Y,err)]
  yerrbar <- aes_string(ymin = paste(yerr_names, collapse = '-'), 
                        ymax = paste(yerr_names,collapse='+'))
  
  p <- ggplot(Data, aes_string(x=names(Data)[X], y=names(Data)[Y],fill=fill,colour=colour,shape=shape)) + theme_bw() + 
    geom_bar(position="dodge",stat=stat)+
    geom_errorbar(mapping=yerrbar,
                  position=position_dodge(0.9), width=barWidth,color=barcol,size=barSize) 
  
  return(p)
}

##############
plotGrid <- function(plots,NRows,NCols,HeightBlocks,MyLegend=NULL,LegendRight=NULL,filename=NULL,PDFheight=10,PDFwidth=10){ 
  if(is.null(MyLegend) && length(HeightBlocks) != NRows){stop("The length of each column is given in HeightBlocks. Thus, the length of HeightBlocks should be equal to the number of columns")}
  if(isnot.null(MyLegend) && length(HeightBlocks) != 2 && is.null(LegendRight)){stop("HeightBlocks defines the height of the plots and the legend. Therefore, its length has to be equal to 2")}
  if(isnot.null(MyLegend) && isnot.null(LegendRight) && length(HeightBlocks) != NRows){stop("The length of each column is given in HeightBlocks. Thus, the length of HeightBlocks should be equal to the number of columns")}
  
  if (is.null(filename)){
    t <- c("grid.arrange(arrangeGrob(")
    for (i in 1:length(plots)){
      t <- paste(t,plots[i]," + theme(legend.position='none'),",sep="")
    }
    if (isnot.null(MyLegend)){
      if (length(HeightBlocks)>1){
        jj <-  c(paste0(HeightBlocks[1:length(HeightBlocks)-1],","),HeightBlocks[length(HeightBlocks)])
      } else {
        jj <- HeightBlocks[1]
      }
      jj <- paste(jj,collapse=" ")
      t <- paste(t,"nrow= ",NRows ,", ncol= ",NCols,"),arrangeGrob(MyLegend,nrow=1),heights=c(",jj,"))" )
    } else {
      if (length(HeightBlocks)>1){
        jj <-  c(paste0(HeightBlocks[1:length(HeightBlocks)-1],","),HeightBlocks[length(HeightBlocks)])
      } else {
        jj <- HeightBlocks[1]
      }
      jj <- paste(jj,collapse=" ")
      t <- paste(t,"nrow= ",NRows ,", ncol= ",NCols ,", heights=c(",jj,")))" )
    }
    if(isnot.null(LegendRight)){
      substr(t,nchar(t),nchar(t)) = ""
      t <- paste(t,",ncol=2)")
    }
    return(eval(parse(text=t)))
  } else {
    t <- c("grid.arrange(arrangeGrob(")
    for (i in 1:length(plots)){
      t <- paste(t,plots[i]," + theme(legend.position='none'),",sep="")
    }
    if (isnot.null(MyLegend)){
      if (length(HeightBlocks)>1){
        jj <-  c(paste0(HeightBlocks[1:length(HeightBlocks)-1],","),HeightBlocks[length(HeightBlocks)])
      } else {
        jj <- HeightBlocks[1]
      }
      jj <- paste(jj,collapse=" ")
      t <- paste(t,"nrow= ",NRows ,", ncol= ",NCols,"),arrangeGrob(MyLegend,nrow=1),heights=c(",jj,"))" )
    } else {
      if (length(HeightBlocks)>1){
        jj <-  c(paste0(HeightBlocks[1:length(HeightBlocks)-1],","),HeightBlocks[length(HeightBlocks)])
      } else {
        jj <- HeightBlocks[1]
      }
      jj <- paste(jj,collapse=" ")
      t <- paste(t,"nrow= ",NRows ,", ncol= ",NCols ,", heights=c(",jj,")))" )
    }
    if(isnot.null(LegendRight)){
      substr(t,nchar(t),nchar(t)) = ""
      t <- paste(t,",ncol=2)")
    }
    pdfname <- paste(filename,".pdf",sep="")
    pdf(file=pdfname,width=PDFwidth,height=PDFheight)
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
                   rmar=1,lmar=1,AngleLab=30,LegendPosition="right",PointColor="black",
                   PointShape=16,MarginWidth=1) 
{ 
  if (isnot.vector(obs) || isnot.vector(pred)){
    stop("The input data must be two vectors")
  } else if ( length(obs) != length(pred) ){
    stop("Both vectors have to be of equal length")
  } else if(isnot.null(margin)) {
    Data <- data.frame(Observed=obs,Predicted=pred)
    p <- ggplot(Data, aes(x=Observed, y=Predicted)) +
      geom_point(size=PointSize,colour=PointColor,shape=PointShape) +  
      geom_abline(slope=1,intercept=margin/2,colour=ColMargin,size=MarginWidth) + 
      geom_abline(slope=1,intercept=-(margin/2),colour=ColMargin,size=MarginWidth) + theme_bw() + 
      ggtitle(main) + ylab(ylab) + xlab(xlab)+
      ylim(c(min(c(obs,pred)), max(c(obs,pred)))) + xlim(c(min(c(obs,pred)), max(c(obs,pred)))) +
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
            legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) 
  } else {
    Data <- data.frame(Observed=obs,Predicted=pred)
    p <- ggplot(Data, aes(x=Observed, y=Predicted)) + geom_point(size=PointSize,colour=PointColor,shape=PointShape) + theme_bw() + 
      ggtitle(main) + ylab(ylab) + xlab(xlab) +
      ylim(c(min(c(obs,pred)), max(c(obs,pred)))) + xlim(c(min(c(obs,pred)), max(c(obs,pred)))) +
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
            legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm"))   
  }
  return(p) 
}

# Functions to evaluate models performance
# Tropsha, A.; Golbraikh, A. Predictive Quantitative Structure–Activity Relationships Modeling: 
#Development and Validation of QSAR Models. In: Handbook of Chemoinformatics Algorithms 
#(Faulon, J.-L.; Bender, A., Eds.), Chapter 7, pp. 213-233, Chapman & Hall / CRC, London, UK, 2010.

# calculates the RMSE between two vectors
RMSE <- function(v1, v2) {
  i1 <- which(!is.na(v1))
  i2 <- which(!is.na(v2))
  is <- intersect(i1, i2)
  v1 <- v1[is]
  v2 <- v2[is]
  residuals <- abs(v1-v2)
  return(as.numeric(sqrt( (residuals%*%residuals)/length(v1) )))
}

# calculates the MAE between two vectors
MAE <- function (v1, v2) {
  i1 <- which(!is.na(v1))
  i2 <- which(!is.na(v2))
  is <- intersect(i1, i2)
  v1 <- v1[is]
  v2 <- v2[is]
  residuals <- abs(v1 - v2)
  return(sum(residuals)/length(v1))
}

# Calculates the slope between two vector (k')
slope <- function(v1,v2){ # v1=z.test v2=y.test
  return(sum(v2*v1)/sum(v1*v1))
}

# Calculates the regression coefficient through the origin
Rsquared0 <- function(v1,v2) { #v1=z.test (y), v2=y.test (x)
  if (is.vector(v1) && is.vector(v2) && length(v1)==length(v2)){
    y_obs_mean <- mean(v2)
    yr0 = v1 * slope(v1,v2)
    first_term = (v2 - yr0)*(v2 - yr0)
    second_term= (v2-y_obs_mean)*(v2-y_obs_mean)
    return(1-(sum(first_term)/sum(second_term)))
  }
  else {print("Wrong input: input arguments are not vector or have unequal length")}
}


# Calculates the regression coefficient 
Rsquared <- function(v1,v2) { # v1=z.test (y), v2=y.test (x)
  if (is.vector(v1) && is.vector(v2) && length(v1)==length(v2)){
    y_obs_mean <- mean(v2)
    y_pred_mean <- mean(v1)
    first_term <- sum((v2-y_obs_mean) * (v1 - y_pred_mean))
    second_term <- sqrt(sum((v2-y_obs_mean)*(v2-y_obs_mean)) * sum((v1 - y_pred_mean)*(v1 - y_pred_mean)))
    division <- first_term / second_term
    return(division * division)
  }
  else {print("Wrong input: input arguments are not vector or have unequal length")}
  
}

# Calculates the Q squared 
#Qsquared (z.test,y.test) (predicted vs observed)
Qsquared <- function(v1, v2) {
  if (is.vector(v1) && is.vector(v2) && length(v1)==length(v2)){
    y_obs_mean <- mean(v2)
    first_term <- abs(v1-v2)*abs(v1-v2)
    second_term <- abs(v2-y_obs_mean)*abs(v2-y_obs_mean)
    return(1-(sum(first_term)/sum(second_term)))
  }
  else {print("Wrong input: input arguments are not vector or have unequal length")}
}

################################


