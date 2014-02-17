#!/usr/bin/Rscript

library(ggplot2)
library(grid)
library(MASS)

usage <-function() {
cat("
USAGE: 1.73-plot_tracebacks.R <data_divid.tab> <native_tracebacks.tab> [<secondSieve.tab>]

        data_divid.tab		- full dataset to plot (output of 1.3-merge_clustal)
	native_tracebacks.tab	- subset to highlight (output of 1.72-merge_tracebacks) (in blue)
	secondSieve.tab		- second subset to highlight (red)

Created by Chaim A Schramm 2014-01-07
Copyright (c) 2014 Columbia University Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

")
}


####################################
# CORE PLOTTING FUNCTION
####################################

TBplotter <- function (data, trace, native, plot_title="") {

	g <- kde2d(data$germ_div,data[[native]],n=100)
	densf = data.frame(expand.grid(x=g$x, y=g$y), z=as.vector(g$z))
	r=10^seq(-3.5,log10(max(g$z))+1,.5)

	my_y<-paste("% ID to", native)
	my_x<-"% Divergence from germline V"
	
	p<-ggplot(data,aes_string(x="germ_div",y=native)) +
		stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
		scale_fill_gradientn(colours=rev(gray.colors(15)), trans="log", limits=c(1e-4,12*max(g$z)),na.value="white", guide="none") +
		geom_contour(data=densf, aes(x=x,y=y,z=z), colour='black', lwd=0.5, breaks=r,guide="none") + 
		theme_bw(18) + scale_x_continuous(expand=c(0,0),limits=c(0,40)) +
		scale_y_continuous(expand=c(0,0),limits=c(60,101)) +
		#opts(axis.text=element_blank()) +
		labs(x=my_x, y=my_y, title=plot_title)

	p <- p + geom_point(data=trace,colour='blue', shape=16,guide="none")

}



##########################################################################################
#         MAIN PROGRAM
##########################################################################################

#plot_tracebacks <- function(file1, file2) {

input = commandArgs(T)
#print(input)

#check for parameters
if(length(input)<2){
  stop(usage())
}

#now check to make sure the files are there
if (! all(file.exists(input))) {
  stop("Please check input parameters and make sure all files exist!\n\n")
}

#load full dataset
bigdata <- read.table(input[1],h=T)

#load tracebacks
sieve <- read.table(input[2],h=T)
natName <- colnames(sieve)[3]

#plot
myplot <- TBplotter(bigdata, sieve, natName, paste("CDR3 tracebacks of",natName))


#second sieve?
if (length(input)>2) {
   second <- read.table(input[3],h=T)
   myplot <- myplot + geom_point(data=second,colour='red', shape=16,guide="none")
}

#and print to file
outname=paste0("analysis/figures/",natName,"_tracebacks.png")
#outname=paste0(natName,"_tracebacks.png")
png(outname)
print(myplot)
dev.off()

