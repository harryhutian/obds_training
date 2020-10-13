library(tidyverse)
library(gridExtra)
library(cowplot)
library(patchwork)
library(biomaRt)
library(pheatmap)
library(scales)

setwd("C:/Users/surface/github/obds_training/project_tidyverse/data")
coding_gene_region <- read.table("coding_gene_region.bed")
colnames(coding_gene_region) <- c("chrom","startpos","endpos","ID","score","strand")
View(coding_gene_region)
coding_gene_region$genomic_interval <- coding_gene_region$endpos - coding_gene_region$startpos
View(coding_gene_region)
# Plot a histogram of the lengths using ggplot2:
# Add a plot title
# Change the x and y axis titles and sizes
# Change the size and angle of the x tick labels
# Change the x axis upper limit to 500,000
# Change the number of bins or the bin width
# Change the fill and border colour of the bars

ggplot(coding_gene_region,aes(x = genomic_interval))+
        geom_histogram(bins = 100, colour = "blue", fill = "yellow")+
        labs(title = "Histograms of genomic interval",
             x = "Genomic Interval Length",
             y = "Count")+
        theme(axis.title.y = element_text(size = 10, face = "bold"),
              axis.title.x = element_text(size = 10, face = "bold"),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(angle = 45, hjust = 1, vjust = 0.5),
              plot.title = element_text(hjust = 0.5, size = 10)) + 
        xlim(0, 50000) + 
        scale_x_continuous(labels = comma)

        
        







data("mtcars")
str(mtcars)
par(mar=c(6,6,5,5), lwd=5, mfrow = c(1,2))   
# mfrow sets the plotting area into a 1*2 array

table(mtcars$gear)
barplot(table(mtcars$gear),
       xlab = "Numbers of gears",
       ylab = "Numbers of cars",
       main = "A main title\nspread across two lines",
       col = "grey",
       border = "red")
abline(h = 6, lwd= 1)
# Making vector of colours, named by gears

colours <- c("red","green","blue")
names(colours) <- c(3,4,5)
gears_column <- as.character(mtcars$gear)
colours[gears_column]


colors
par(mar=c(6,6,5,5), lwd=1)
plot(mtcars$mpg, mtcars$hp,
     col = colours[gears_column],
     pch = 9,
     cex = 2,
     xlab = "Miles per gallon",
     ylab = "Horse power",
     cex.lab = 1.5)

legend(legend = sort(unique(mtcars$gear)),
       x = "right",
       fill = c("red","green","blue"))


