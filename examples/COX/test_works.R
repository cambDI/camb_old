library(grid)
library(camb)

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
