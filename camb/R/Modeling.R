#################################################################################
## Modeling 
#################################################################################

##############
## The following functions assume that you already have features calculated and that
## the features are seperate from the target values
ViewTargets <- function(y, bw) {
  plot(density(y, bw = bw), main="Targets")
}

# TBD write function to view a feature by name in the same way as the targets above
##############
ReplaceInfinitesWithNA <- function(d) {
  do.call(data.frame,lapply(d, function(x) replace(x, is.infinite(x),NA)))
}

##############
ImputeFeatures <- function(d) {
  library(impute)
  as.data.frame(impute.knn(as.matrix(d), k = 10)$data)
}

##############
SplitSet <- function(ids, x, y, percentage = 20, seed = 1) {
  holdout.size <- round(nrow(x) * (percentage/100))
  set.seed(seed)
  holdout.indexes <- sample(1:nrow(x), holdout.size, replace=FALSE)
  train.indexes <- (1:length(y))[-holdout.indexes]
  x.train <- x[train.indexes, ]
  x.holdout <- x[holdout.indexes, ]
  y.train <- y[train.indexes]
  y.holdout <- y[holdout.indexes]
  l <- list()
  l$ids <- ids
  l$holdout.indexes <- holdout.indexes
  l$train.indexes <- train.indexes
  l$x.train <- x.train
  l$x.holdout <- x.holdout
  l$y.train <- y.train
  l$y.holdout <- y.holdout
  l
}

##############
RemoveNearZeroVarianceFeatures <- function(ss, frequencyCutoff = 30/1) {
  nzv.columns <- nearZeroVar(ss$x.train, freqCut = frequencyCutoff)
  if (length(nzv.columns) != 0) {
    print(paste(length(nzv.columns), "features removed with variance below cutoff"))
    ss$x.train <- ss$x.train[, -nzv.columns]
    ss$x.holdout <- ss$x.holdout[, -nzv.columns]
  }
  else {
    print("no features removed")
  }
  ss
}

##############
RemoveHighlyCorrelatedFeatures <- function(ss, correlationCutoff = 0.95) {
  hc.columns <- findCorrelation(cor(ss$x.train), correlationCutoff)
  if (length(hc.columns) != 0) {
    print(paste(length(hc.columns), "features removed with correlation above cutoff"))
    ss$x.train <- ss$x.train[, -hc.columns]
    ss$x.holdout <- ss$x.holdout[, -hc.columns]
  }
  else {
    print("no features removed")
  }
  ss
}

##############
PreProcess <- function(ss, steps = c("center", "scale")) {
  transformation <- preProcess(ss$x.train, method = steps)
  ss$x.train <- predict(transformation, ss$x.train)
  ss$x.holdout <- predict(transformation, ss$x.holdout)
  ss$transformation <- transformation
  ss
}

##############
GetCVTrainControl <- function(ss, seed = 1, folds = 5, repeats = 1) {
  set.seed(seed)
  ss$trControl <- trainControl(method='cv', number=folds, repeats=repeats, returnResamp='none',
                               returnData=FALSE, savePredictions=TRUE,
                               verboseIter=TRUE, allowParallel=TRUE,
                               index=createMultiFolds(ss$y.train, k=folds, times=repeats))
  ss
}

##############
expGrid <- function(ini,end,stride,base){
  grid <- c()
  for (i in seq(ini,end,stride)){
    grid <- append(grid,base^i)
  }
  return(grid)
}

##############
## Maximum Model Performance
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

MaxPerf <- function(meanNoise=0,sdNoise,meanResp,sdResp,lenPred,iters=1000,
                    filename=NULL,pdfW=10,pdfH=10,TextSize=15,TitleSize=15,
                    XAxisSize=15,YAxisSize=15,TitleAxesSize=15,tmar=1,bmar=1,
                    rmar=1,lmar=1){
  R2 <- c()
  R02 <- c()
  Q2 <- c()
  rmsep <- c()
  for (i in 1:iters){
    x <- rnorm(lenPred,mean=meanResp,sd=sdResp)
    noise <- rnorm(length(x),mean=meanNoise,sd=sdNoise)
    y <- x+noise
    R2[i] <- Rsquared(y,x)
    Q2[i] <- Qsquared(y,x)
    R02[i] <- Rsquared0(y,x)
    rmsep[i] <- RMSE(y,x)
  }
  
  pp <- data.frame(R2=R2)
  title <- expression(paste("R"^{2}))
  p1 <- ggplot(pp, aes(x=R2)) + 
    theme_bw() + 
    geom_density(alpha=.2, fill="#FF6666") + ggtitle(title)+ aes(y = ..count..)+ 
    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
          axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
    geom_vline(aes(xintercept=mean(R2, na.rm=T)), color="blue", linetype="dashed", size=1)
  
  pp <- data.frame(rmsep=rmsep)
  p2 <- ggplot(pp, aes(x=rmsep)) + 
    theme_bw() + 
    geom_density(alpha=.2, fill="#FF6666") + ggtitle('RMSEP') +  aes(y = ..count..)+ 
    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
          axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
    geom_vline(aes(xintercept=mean(rmsep, na.rm=T)), color="blue", linetype="dashed", size=1)
  
  pp <- data.frame(R02=R02)
  title <- expression(paste("R"[0]^{2}))
  
  p3 <- ggplot(pp, aes(x=R02)) + 
    theme_bw() + 
    geom_density(alpha=.2, fill="#FF6666") + ggtitle(title) + aes(y = ..count..)+ 
    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
          axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
    geom_vline(aes(xintercept=mean(R02, na.rm=T)), color="blue", linetype="dashed", size=1)
  
  title <- expression(paste("Q"^{2}))
  pp <- data.frame(Q2=Q2)
  p4 <- ggplot(pp, aes(x=Q2)) + 
    theme_bw() + 
    geom_density(alpha=.2, fill="#FF6666") + ggtitle(title) + aes(y = ..count..)+ 
    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
          axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
    geom_vline(aes(xintercept=mean(Q2, na.rm=T)), color="blue", linetype="dashed", size=1)
  
  if (isnot.null(filename)){
    pdfname=paste(filename,".pdf",sep="")
    pdf(pdfname,width=pdfW,height=pdfH)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 2)))
    print(p1, vp = vplayout(1,1))
    print(p2, vp = vplayout(1,2))
    print(p3, vp = vplayout(2,1))
    print(p4, vp = vplayout(2,2))
    dev.off()
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 2)))
    print(p1, vp = vplayout(1,1))
    print(p2, vp = vplayout(1,2))
    print(p3, vp = vplayout(2,1))
    print(p4, vp = vplayout(2,2))
  }
  return(list(p1,p2,p3,p4))
}
