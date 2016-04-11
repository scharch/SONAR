#!/usr/bin/env Rscript


########################
#
# ADD CUSTOMIZATION OPTIONS:
#     - figure size
#     - dpi
#     - page title
#     - colors and contours
#     - show/hide color guide
#
# GENERAL FEATURES STILL NEEDED:
#     - implement showNames
#     - correct default handling of axis labels
#     - row and column titles
#     - what should default guide behavior be?
#
########################


"Usage: 4.3-plot_identity_divergence.R <idDivFile>... [ --plot <ids.list>... ] [ --sieve1 <ids.list>... --sieve2 <ids.list> ] [ --mab <abID>... ] [ --xaxis <label> --outdir </path> --output <output.pdf> --labels <dates.list> --reference <idDivFile> --showNames <0|1|2> --transpose]

Options:
  -h --help		   Show this documentation.
  <idDivFile>	           Output from 2.1-calulate_id-div.pl with data for plotting.
			       Multiple files will be treated as longitudinal time
		    	       points and displayed in the order entered.
  --plot <ids.list>    	   Specify a subset of sequences to be plotted (eg those from
   	       	      	       a particular V gene). One id per line. If used, must
			       match number and order of data files.
  --sieve1 <ids.list>	   A subset of reads to overlay on the I-D plot as points
   	       		       (eg those with CDR3 homology to the antibody of interest).
			       If used, must match number and order of data files.
  --sieve2 <ids.list>	   A second subset of reads to overlay on the I-D plot as points.
			       If used, must match number and order of data files.
  --mab <abID>		   Identity referent to use on y-axis. Default is the first mAb
  			       in the first id-div.tab file. If multiple referents are
			       specified, each will be shown in its own row of plots.
  --xaxis <label>	   X-axis label is '% divergence from <label>'.
  	  		       [default: germline V]
  --outdir </path>	   Directory for saving output. [default: output/plots]
  --output <output.pdf>	   File name of output. [default: id-div.pdf]
  --labels <dates.list>    List of labels for the longitudinal time points, one label
  	   		       per line. Must match number and order of data files.
  --reference <idDivFile>  Data points to display on ALL plots, eg other members of the
  	      		       antibody lineage of interest.
  --showNames <0|1|2>	   Show text labels for native antibodies supplied with --reference.
  	      		       0 = none; 1 = first plot in each row; 2 = all
			       [default: 0]
  -t --transpose	   Show time in rows and different antibody referents in columns.
  			       [default: FALSE]
			       

Created by Chaim A Schramm 2013-07-11.
Edited and commented by CAS 2015-03-09.
Overhauled for generality CAS 2016-04-11.

Copyright (c) 2013-2016 Columbia University and Vaccine Research Center, National
                               Institutes of Health, USA. All rights reserved.

" -> usage



####################################
# PLOTTING FUNCTION
####################################

plot_all <- function (data, native,  heavy, plot_title=NULL, color=TRUE, guide=TRUE, xlabel=TRUE, ylabel=TRUE, contour=TRUE, conCol = 'black') {

	my_x<-paste("% divergence from", heavy)
	my_y<-paste("% ID to", native)

	if (color) {
	   my_colors=rev(rainbow(15,end=4/6))
	} else {
	   my_colors=rev(gray.colors(5))
	}

	g     <- kde2d(data$germ_div,data[[native]],n=100,h=1,lim=c(0,60,40,100))
	densf <- data.frame(expand.grid(x=g$x, y=g$y), z=as.vector(g$z))
	b     <- (sum(g$z) / length(data$germ_div))/2
	t     <- b * 10^4
	if ( max(g$z) > t ) {
	      t <- b * 10^ceiling( log10( max(g$z)/b ) )
	}
	r     <- 10^seq(log10(b), log10(t), 1)

	p<-ggplot(densf,aes(x,y,z=z)) +
	 		geom_tile(aes(fill = z)) +
			scale_fill_gradientn(colours=rev(rainbow(15,end=4/6)), trans="log10", limits=c(b,t),
				na.value="white", breaks=r, labels=signif(r/b,1),
			        guide = guide_colorbar( title="number of\nsequences", title.theme=element_text(size=4,angle=0),
				        barheight = unit(.5,"in"), barwidth = unit(.1,"in"), label.theme=element_text(size=3,angle=0,) ) )+
			theme_bw() + scale_x_continuous(expand=c(0,0),limits=c(0,50)) +
			scale_y_continuous(expand=c(0,1),limits=c(50,100)) +
	     		theme(plot.background = element_blank(),panel.grid.major = element_blank(),
				axis.ticks.length = unit(.02,"in"), axis.ticks = element_line(size = .5),
		          	panel.grid.minor = element_blank(), axis.text = element_text(size = 6),
				axis.title = element_text(size = 8), plot.margin = unit(c(.1,.1,.1,.1),"in"),
				plot.title = element_text(size = 8) )

	if ( contour ) 		     { p <- p + stat_contour(colour=conCol, size=.25, breaks=10*r) }
	if ( ! is.null(plot_title) ) { p <- p + labs( title=plot_title ) }
	if ( xlabel ) 		     { p <- p + labs( x=my_x ) }
	if ( ylabel ) 		     { p <- p + labs( y=my_y ) }

	p
}




####################################
# DISPLAY FUNCTION
####################################

# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




##########################################################################################
#         MAIN PROGRAM
##########################################################################################

lds <- function(dataFileList, subsetFileList, firstSieveList, secondSieveList, natAbList, xLabel, outDir, outFile, timepointList, refPointsData, showNames, transpose) {

#initiate
myplots=list()

if ( ! is.null(refPointsData) ) {
   refData <- read.table( refPointsData, h=T )
}

#process each time point
for ( time in seq(length(dataFileList)) ) {

    	#main plot
	bigdata   <- read.table( dataFileList[time], h=T )
	smalldata <- bigdata

	if ( is.null(natAbList) )  natAbList <- c( names(bigdata)[3] )

	for (mab in natAbList) {

	    if ( ! is.null(subsetFileList) ) {
	       subset    <- scan(subsetFileList[time], what='character', quiet=T)
	       smalldata <- bigdata[which(bigdata$ID %in% subset),]
	    }

	    ####
	    # ADD OPTIONS (plot color, contours, yes/no to axis labels, scale bar)
	    ####
	    p <- plot_all(smalldata, mab, xLabel)

	    if ( ! is.null(firstSieveList) ) {
	       subset <- scan(firstSieveList[time], what='character', quiet=T)
	       s1data <- bigdata[which(bigdata$ID %in% subset),]
	       p <- p + geom_point(data=s1data,aes_string(x="germ_div",y=mab,z=NA),colour='goldenrod', shape=16, size=.5)
	    }

	    if ( ! is.null(secondSieveList) ) {
	       subset <- scan(secondSieveList[time], what='character', quiet=T)
	       s2data <- bigdata[which(bigdata$ID %in% subset),]
	       p <- p + geom_point(data=s2data,aes_string(x="germ_div",y=mab,z=NA),colour='magenta', shape=16, size=.5)
	    }

	    if ( ! is.null(refPointsData) ) {
	       #####
	       # NEED TO HANDLE SHOW_NAMES OPTION HERE
	       #####
	       p <- p + geom_point(data=refData,aes_string(x="germ_div",y=mab,z=NA),colour='black', shape=4, size=1)
	    }

	    myplots[[length(myplots)+1]] <- p

	} #next mab in list

} #next time point


numCols <- length(dataFileList)
numRows <- length(natAbList)

if ( transpose ) {
   numCols <- length(natAbList)
   numRows <- length(dataFileList)
}

height  <- 2*numRows
width   <- 2.5*numCols

ggsave(
	paste(outDir,outFile,sep="/"),
	multiplot( plotlist=myplots, cols=numCols ),
	dpi = 600,
	height = height,
	width = width
)


}



##########################################################################################
#         ENTRY POINT FROM COMMAND LINE
##########################################################################################

library(docopt)
library(ggplot2)
library(grid)
library(MASS)


opts<-docopt(usage, strip=T, help=T)


# check to make sure the files are there
if (! all(file.exists(opts$idDivFile))) {
  stop("Cannot find 1 or more input files...\n\n")
}

if (! is.null(opts$plot)) {
   if ( length(opts$plot) != length(opts$idDivFile) ) { stop("Please make sure that the number of plot files matches the number of input files and that they are in the same order\n\n") }
   if (! all(file.exists(opts$plot))) {
      stop("Cannot find 1 or more plot files...\n\n")
   }
}

if (! is.null(opts$sieve1)) {
   if ( length(opts$sieve1) != length(opts$idDivFile) ) { stop("Please make sure that the number of sieve1 files matches the number of input files and that they are in the same order\n\n") }
   if (! all(file.exists(opts$sieve1))) {
      stop("Cannot find 1 or more sieve1 files...\n\n")
   }
}

if (! is.null(opts$sieve2)) {
   if ( length(opts$sieve2) != length(opts$idDivFile) ) { stop("Please make sure that the number of sieve2 files matches the number of input files and that they are in the same order\n\n") }
   if (! all(file.exists(opts$sieve2))) {
      stop("Cannot find 1 or more sieve2 files...\n\n")
   }
}

if (! is.null(opts$reference) && ! file.exists(opts$reference)) { stop("Cannot find reference data file...\n\n") }
 

lds(opts$idDivFile, opts$plot, opts$sieve1, opts$sieve2, opts$mab, opts$xaxis, opts$outdir, opts$output, opts$labels, opts$reference, opts$showNames, opts$transpose)
