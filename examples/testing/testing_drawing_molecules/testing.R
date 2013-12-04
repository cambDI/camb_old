library(camb)
library(png)

DrawMoleculeInSDFbyID <- function(structures.file, structure.ID, file.name, useNameAsTitle) {
  .C("R_drawMoleculeInSDFbyID", structures.file, as.character(structure.ID), file.name, useNameAsTitle)
}

DrawMoleculeInSDF <- function(structures.file, structure.number, file.name, useNameAsTitle) {
  print(structure.number)
  .C("R_drawMoleculeInSDF", structures.file, as.integer(structure.number), file.name, useNameAsTitle)
}

PlotMolecules <- function(sdf.file, IDs) {
  ggplot() + annotation + ylim(0,1) + xlim(0,1)
  temp.png <- tempfile("temp", fileext=".png")
  for(id in IDs) {
    DrawMoleculeInSDFbyID(structures.file=sdf.file, structure.ID=id, temp.png, TRUE)
    img <- readPNG(temp.png)
    plot(1:2, type='n')
    rasterImage(img, xleft=1, xright=2, ybottom=1, ytop=2)
  }
  dev.off()
}

molX + theme(legend.position="none"),
library(ggplot2)
df <- data.frame()
g <- ggplot(df) + xlim(0, 1) + ylim(0, 1)
g + annotation

g <- ggplot(df) + geom_blank() + xlim(0, 10) + ylim(0, 10) + theme_bw()
g

img <- readPNG("temp.png")
g <- rasterGrob(img, interpolate=TRUE)
annotation <- annotation_custom(g, xmin=0, xmax=1, ymin=0, ymax=1)
base <- qplot(0:1, 0:1, geom = "blank")  +
base <- theme(axis.line = element_blank(), 
                       axis.text.x = element_blank(), 
                       axis.text.y = element_blank(),
                       axis.ticks = element_blank(), 
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       plot.margin=unit(c(0,0,0,0), "cm"))
base + annotation

library(gridExtra)
grid.arrange(arrangeGrob(base + annotation, base + annotation, base + annotation, base + annotation), nrow=1)
  
arrangeGrob(molX + theme(legend.position="none"),molX2 + theme(legend.position="none"),molX3 + theme(legend.position="none"),molX4 + theme(legend.position="none"),nrow=2)
  
img <- readPNG("temp.png")
g <- rasterGrob(img, interpolate=TRUE)

ttt= list(g,g)

qplot(1, 1, geom="blank") +annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + theme_bw()

p <- qplot(1, 1, geom="blank") + theme_bw() 
tt = p +annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + base
tr = p +annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + base
grid.arrange(tt,tr,tt,tr,nrow=2,main="kkk")

library(png)
library(grid)
img <- readPNG("temp.png")
g <- rasterGrob(img, interpolate=TRUE)
annotation <- annotation_custom(g, xmin=0, xmax=1, ymin=0, ymax=1)

annotation_custom(

MoleculesToPDF <- function(pdf.file, sdf.file, IDs) {
  temp.png <- tempfile("temp", fileext=".png")
  pdf(pdf.file)
  
  for(id in IDs) {
    DrawMoleculeInSDFbyID(structures.file=sdf.file, structure.ID=id, temp.png, TRUE)
    img <- readPNG(temp.png)
    g <- rasterGrob(img, interpolate=TRUE)
    
    annotation <- annotation_custom(g, xmin=position[1]-size/2, xmax=position[1]+size/2, 
                                    ymin=position[2]-size/2, ymax=position[2]+size/2)
    #plot(1:2, type='n')
    rasterImage(img, xleft=1, xright=2, ybottom=1, ytop=2)
  }
  dev.off()
}

MoleculesToPDF(pdf.file="test.pdf", sdf.file="structures.sdf", IDs=c("Tolazamide", "Oseltamivir", "Granisetron"))

library(ggplot2)

MoleculesToPDF("test.pdf", "structures.sdf", IDs=1:10)
temp.png <- tempfile("temp", fileext=".png")
pdf(pdf.file)
for(id in IDs) {
  DrawMoleculeInSDFbyID(structures.file=sdf.file, structure.ID=id, temp.png, TRUE)
  img <- readPNG(temp.png)
  img <- readPNG(file.name)
  g <- rasterGrob(img, interpolate=TRUE)
  
  annotation <- annotation_custom(g, xmin=position[1]-size/2, xmax=position[1]+size/2, 
                                  ymin=position[2]-size/2, ymax=position[2]+size/2)
  ggplot(
  plot(1:2, type='n')
  rasterImage(img, xleft=1, xright=2, ybottom=1, ytop=2)
}
dev.off()

DrawMoleculeInSDFbyID(


DrawMoleculeInSDF(structures.file=sdf.file, structure.number=1, file.name="temp.png", useNameAsTitle=FALSE)