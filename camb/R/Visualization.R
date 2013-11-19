#################################################################################
## Visualization
#################################################################################

##############
## Density of the response variable
DensityResponse <- function(Data,xlab="",ylab="",main="",alpha=0.2,
                            binwidth=NULL,histFill="white",histCol="black",
                            densityFill="#FF6666",TitleSize=15,TextSize=15,
                            XAxisSize=15,YAxisSize=15,AngleLab=30,LegendPosition="right",
                            TitleAxesSize=15,tmar=1,bmar=1,rmar=1,lmar=1){
  if (!is.vector(Data)) stop("Input data must be a numeric matrix or data.frame")
  Data <- data.frame(Values=Data)
  if (is.null(binwidth)) binwidth <- abs(range(Data)[1] - range(Data)[2]) / 10
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
## PCA analysis of the Proteins 

isnot.null <- function(x) ! is.null(x)
isnot.vector <- function(x) ! is.vector(x)

##############
## Get the two first PC of the sequence descriptors
PCAProt <- function(Data,SeqsName=NULL){ 
  if (is.matrix(Data) || is.data.frame(Data)){
    ana <- prcomp(t(Data), cor=TRUE,scale=TRUE,center=TRUE)  
    PC1 <- as.vector(ana$rotation[,1])   
    PC2 <- as.vector(ana$rotation[,2])
    if (is.null(SeqsName)){
      Data <- data.frame(PC1,PC2)
    } else {
      Data <- data.frame(PC1,PC2,SeqsName)
    }
    return(Data)
  }else{
    stop("Input data must be a numeric matrix or data.frame")
  }
}

##############
# Plot the towo first PC of the sequence descriptors
PCAProtPlot <- function (Data,main="",ylab="PC2",xlab="PC1",Seqs=NULL,PointSize=4,
                         LegendPosition="right",LegendName="Sequences",ColLegend=1,RowLegend=NULL,
                         TitleSize=15,TextSize=15,XAxisSize=15,YAxisSize=15,
                         TitleAxesSize=15,tmar=1,bmar=1,rmar=1,lmar=1) 
{
  if (length(names(Data)) < 2 || length(names(Data)) > 3){
    stop("Two PCA required. The Data.frame provided has less than two columns (PCA) or more than 3")
  } else if (names(Data)[1] != "PC1" || names(Data)[2] != "PC2"){
    stop("Column names have to be {PC1, PC2}")
  } else if (length(names(Data)) == 2 && is.null(Seqs)){
    print("No sequence names provided")
    p <- ggplot(Data, aes(x=PC1, y=PC2)) +
      geom_point(size=PointSize) + theme_bw() + ggtitle(main) + ylab(ylab) + xlab(xlab) +
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,
            plot.title=element_text(size=TitleSize),legend.key=element_blank(),
            plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
      guides(colour = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend),
             shape = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend))
  } else if (length(names(Data)) == 2 && isnot.null(Seqs)){
    if (length(Seqs) != nrow(Data) || isnot.vector(Seqs)) {
      stop("Either the sequences are not in a vector, or its length is not equal to the number of datapoints (rows of the input data)")
    } else {
      print("Sequence names provided")
      Data <- data.frame(Data,Seqs=Seqs)
      a <- length(unique(Seqs))
      p <- ggplot(Data, aes(x=PC1, y=PC2, color=Seqs,shape=Seqs)) +
        geom_point(size=PointSize)  + scale_shape_manual(values=1:a) + theme_bw() + ggtitle(main) + ylab(ylab) + xlab(xlab) +
        theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
              axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
              axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
              legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
        guides(colour = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend),
               shape = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend))
    } 
  } else {
    print("Sequence names provided in the third column of the data.frame")
    Seqs <- unlist(Data[names(Data)[3]])
    a <- length(unique(Seqs))
    namee <- names(Data)[3]
    p <- ggplot(Data, aes_string(x="PC1", y="PC2",colour=namee,shape=namee)) +
      geom_point(size=PointSize) +scale_shape_manual(values=1:a) + theme_bw() + ggtitle(main) + ylab(ylab) + xlab(xlab) +
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
            legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
      guides(colour = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend), shape = guide_legend(LegendName,ncol=ColLegend,nrow=RowLegend))
  }
  return(p) 
}

##############
PairwiseDist <- function(Data,method="jaccard",..){
  if (is.matrix(Data) || is.data.frame(Data)){
    Data <- unique(Data)
    methods <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard",
                 "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial",
                 "chao", "cao")
    method <- match.arg(method,methods)
    pwdist <- vegdist(Data, method = method)
    pwdist <- data.frame(as.vector(pwdist))
    names(pwdist) <- "Distance"
    return(pwdist)
  } else {
    stop("Input data must be a numeric matrix or data.frame")
  }
}

##############
PairwiseDistPlot <- function(Data,xlab="",ylab="",main="",TextSize=15,TitleSize=15,XAxisSize=15,YAxisSize=15,
                             TitleAxesSize=15,tmar=1,bmar=1,rmar=1,lmar=1,LegendPosition="right"){
  if (is.matrix(Data) || is.data.frame(Data)){
    p <- ggplot(Data, aes(x=Distance)) + 
      geom_histogram(aes(y=..density..),binwidth=.05,colour="black", fill="white") + 
      geom_density(alpha=.2, fill="#FF6666")+ ylab(ylab) + xlab(xlab) + ggtitle(main) +
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
            legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) 
    return(p)
  } else {
    stop("Input data must be a numeric matrix or data.frame")
  }
}

##############
# Diversity selection of compounds
