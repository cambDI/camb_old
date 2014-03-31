library(grid)
library(camb)
setwd("/Users/icortes/Desktop/camb_final/camb/examples/testing/testing_drawing_molecules")

base <- theme(axis.line = element_blank(), 
              axis.text.x = element_blank(), 
              axis.text.y = element_blank(),
              axis.ticks = element_blank(), 
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.margin=unit(c(0,0,0,0), "cm"))

img <- readPNG("temp.png")
g <- rasterGrob(img, interpolate=TRUE)

p <- qplot(1, 1, geom="blank") + theme_bw() 
tt = p +annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + base
tr = p +annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + base
grid.arrange(tt,tr,tt,tr,nrow=2,main="Test Drawing")

pdf("test.pdf")
grid.arrange(tt,tr,tt,tr,nrow=2,main="Test Drawing")
dev.off()


DrawMoleculeInSDFbyID <- function(structures.file, structure.ID, file.name, useNameAsTitle) {
  .C("R_drawMoleculeInSDFbyID", structures.file, as.character(structure.ID), file.name, useNameAsTitle)
}

DrawMoleculeInSDF <- function(structures.file, structure.number, file.name, useNameAsTitle) {
  print(structure.number)
  .C("R_drawMoleculeInSDF", structures.file, as.integer(structure.number), file.name, useNameAsTitle)
}

PlotMolecules <- function(sdf.file, IDs,pdf.file="aaa.pdf",PDFMain="") {
  temp.png <- tempfile("temp", fileext=".png")
  base <- theme(axis.line = element_blank(), 
                axis.text.x = element_blank(), 
                axis.text.y = element_blank(),
                axis.ticks = element_blank(), 
                axis.title.x = element_blank(), 
                axis.title.y = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                plot.margin=unit(c(0,0,0,0), "cm"))
  
  for(id in 1:IDs) {
    DrawMoleculeInSDF(structures.file=sdf.file, structure.number=id, temp.png, TRUE)
    img <- readPNG(temp.png)
    g <- rasterGrob(img, interpolate=TRUE)
    p <- qplot(1, 1, geom="blank") + theme_bw() 
    assign(paste("p",id,sep=""), p +annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + base)
  }
  pdf(pdf.file)
  grid.arrange(p1,p2,p3,p4,nrow=2,main=PDFMain)
  dev.off()
  
  return(list(p1,p2,p3,p4))
}


