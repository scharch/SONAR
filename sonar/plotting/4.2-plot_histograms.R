#!/usr/bin/env Rscript


usage <-function() {
cat("
4.2-plot_histogram.R

This is a generic script for plotting nice histograms. Expected usage
      is to be called directly by 4.1-setup_plots.pl, which parses out
      specific desired data from annotation of the repertoire and 
      reformats it to be compatible with this script.

Usage: 4.2-plot_histogram.R data.txt outfile.png
       			    [ title=title bars=dodge percent=T xlab=category 
			      ylab=percent xlim=c(0,100) ylim=c(0,100)
       			      logx=F logy=F mids=F magnify=1 showlegend=T
			      legendpos=right h=height w=width d=dpi ]

    Invoke with -h or --help to print this documentation.

    data.txt      Tab-delimited file with data to be plotted. Expects 3
    		        unlabled columns: group, x, y
    outfile.png	  Where to save the image. Filetype will be parsed
    		  	automatically from the extension.
    title	  Display title of plot. Default = none.
    bars	  Display bars for multiple groups side-by-side ('dodge')
    		  	or stacked ('stack'). Default = dodge.
    percent	  Boolean indicating whether the y axis should be plotted
    		  	as percent of total for each group or raw counts.
			Default = T.
    xlab	  Display label for x axis. Default = 'category'.
    ylab	  Display label for y axis. Default = 'percent' or 'number'
    		  	depending on the value of the percent flag.
    xlim          Display limits for x axis, specified as an R list.
    		  	Default = chosen automatically.
    ylim          Display limits for y axis, specified as an R list.
    		  	Default = chosen automatically.
    logx          Boolean indicating if the x axis should be log-scaled. If
    		  	true, plot type is switched from bar to step.
			Default = F.
    logy          Boolean indicating if the y axis should be log-scaled.
			Default = F.
    mids          Boolean indicating if position of bars should be shifted
    		  	so that they are centered on the midpoint of the
			bins instead of the lower bound. Default = F.
    magnify       Magnification factor for labels and titles. Default = 1.
    showlegend    Boolean indicating whether legend should be displayed.
    		  	Default = T.
    legendpos     Position for placement of legend. Default = 'right'.
    h             Desired height of final plot, in inches. Default = 4.
    w		  Desired width of final plot, in inches. Default = 6.
    dpi 	  Desired resolution of output image. Default = 150.

Created by Chaim A Schramm 2015-11-19.

Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National
                               Institutes of Health, USA. All rights reserved.

")
}



##########################################################################################
#         MAIN PROGRAM
##########################################################################################

main <- function(infile, outfile, title, bars, percent, xlab, ylab, xlim, ylim, logx, logy, mids, magnify, showlegend, legendpos, h, w, d) {

    data<-read.table( infile, h=F )

    #make sure R plots things in the correct order:
    data$V1 <- ordered(data$V1,levels=unique(as.character(data$V1)))
    if (! is.numeric(data$V2)) { data$V2 <- ordered(data$V2,levels=unique(data$V2)) }

    
    #adjust X to midpoint of bin, if desired
    if (mids) {
       	      if (is.numeric(data$V2)) {
	      	 	movedData <- data
			data <- data.frame( V1=c(), V2=c(), V3=c() )
			for ( g in levels(movedData$V1) ) {
			      	subData <- movedData[ movedData$V1==g, ]
			      	subData$oldX <- subData$V2
       	      		      	for (i in seq(length(subData[,1]))) { 
				       subData$V2[i] <- subData$oldX[i] + ( (subData$oldX[i+1]-subData$oldX[i]) / 2 )
				}
				data <- rbind( data, subData[1:length(subData[,1])-1, c(1,2,3)] )
			}
	      } else {
			message("Option 'mids' not compatible with category plots; ignoring...")      
	      }
    }


    #fill in missing points to control bar width
    if (bars == "dodge") {
          if (is.numeric(data$V2)) {
       		data <- rbind( data, expand.grid(V1=levels(data$V1), V2=as.numeric(levels(as.factor(data$V2))), V3=NA) )
		if (xlab=="category") { xlab <- "number" } #check default while we're here
          } else {
 		data <- rbind( data, expand.grid(V1=levels(data$V1), V2=levels(data$V2), V3=NA) )
          }
    }


    #add floor bin for step plot if x-axis is log-scaled
    # shouldn't really matter, but add to top in reverse order so levels stay in same order when done
    # changing data must be done before setting up plot, hence 2 if(logx) statements
    if (logx) {
       	      for ( g in rev(levels(data$V1)) ) {
       	      	    newMin = 0.999 * data$V2[data$V1 == g][1]
	      	    data <- rbind( data.frame(V1=g, V2=newMin, V3=0), data )
	      }
    }

    # change zeroes to NAs for ylog plots
    if (logy) {
       	      data$V3[which(data$V3==0)] <- NA
    }


    if (percent) {
       	 for (g in levels(data$V1) ) {
		total = sum( data$V3[data$V1 == g], na.rm=T )
		data$V3[ data$V1==g ] <- data$V3[data$V1==g] * 100 / total
	 }
    } else if (ylab == "percent") {
      	      	 ylab <- "number"
    }


    #which column is the x axis and which is the grouping depends on use of 'stack'
    myGroupColumn <- "V1"
    myXposColumn  <- "V2"
    if (bars=="stack") {
            myGroupColumn <- "V2"
    	    myXposColumn  <- "V1"
    }


    #colors generated stochastically via 
    #http://tools.medialab.sciences-po.fr/iwanthue/
    #and sorted for increasing hue
    fullList = c("#BE4229", "#E74721", "#8C431C", "#CB6A27", "#E98C25", "#946E13",
    	       	 "#D0A620", "#8B8A22", "#A9B71D", "#555C10", "#839A21", "#82B531",
		 "#4D8426", "#55C42A", "#34992A", "#2B6B2E", "#4FC456", "#33B66D",
		 "#4296CB", "#5A8CE2", "#3E5988", "#656CE2", "#524EA0", "#8F83CC",
		 "#A57CE4", "#8E46AD", "#C056EB", "#CA6BE4", "#7B4D87", "#D186D7")
    subset = fullList[floor( seq(1,length(levels(data[[myGroupColumn]]))) * 30 / length(levels(data[[myGroupColumn]])) )]
    

    #start plot
    p <- ggplot(data, aes_string(x=myXposColumn,y="V3"))

    #logx is proxy for plot type, so check first
    if (logx) {
       	   p <- p + geom_step(aes_string(colour=myGroupColumn), show.legend=showlegend) + 
       	      	    scale_colour_manual( values=subset ) + scale_x_log10( )
    } else {
	   p <- p + geom_bar(aes_string(fill=myGroupColumn), stat='identity', position=bars, show.legend=showlegend) + 
       	      	    scale_fill_manual( values=subset )
    }
        

    #check y axis
    if (logy) {
       	      p <- p + scale_y_log10( )
    }

    myXrotation <- 0
    myHjust     <- 0.5
    if (nlevels(data[[myXposColumn]]) > 6) {
            myXrotation <- 90
	    myHjust     <- 1
    }
    #do all other prettifying here
    p <- p + labs(x=xlab, y=ylab, title=title) + theme_bw() +
      	     coord_cartesian( xlim=eval(parse(text=xlim)), ylim=eval(parse(text=ylim)) ) +
             theme( axis.text=element_text(size=10*magnify), axis.text.x=element_text(angle=myXrotation, hjust=myHjust, vjust=0.5), 
	     	    	axis.title=element_text(size=12*magnify),
			legend.title=element_blank(), legend.position=legendpos, 
			legend.text=element_text(size=10*magnify), legend.key.size=unit(0.2*magnify,"in"),
                        legend.background=element_blank() )


    #and output
    ggsave(outfile, p, height=h, width=w, units="in", dpi=dpi)
	
}


##########################################################################################
#         ENTRY POINT FROM COMMAND LINE
##########################################################################################

library(ggplot2)
library(grid)

#args are:
#  data.txt outfile.png 
#  [ title=title bars=dodge percent=T xlab=category ylab=percent 
#    xlim=c(0,100) ylim=c(0,100) logx=F logy=F mids=F magnify=1
#    showlegend=T legendpos="top" h=height w=width d=dpi ]
input = commandArgs(T)

#check for parameters
if(length(input)<2){
  stop(usage())
}

#now check to make sure the data file is there
if (! file.exists(input[1])) {
  stop("Please check and make sure data file exists!\n\n")
}


#set defaults and process optional parameters
infile     <- input[1]
outfile    <- input[2]
title      <- ""
bars       <- "dodge"
percent    <- T
xlab       <- "category"
ylab       <- "percent"
xlim       <- NULL
ylim       <- NULL
logx       <- F
logy       <- F
mids       <- F
magnify    <- 1
showlegend <- T
legendpos  <- "right"
h          <- 4
w          <- 6
dpi        <- 150

for (i in seq(3,length(input))) { assign(sub("=.*$", "", input[i]), sub("^[^=]*=", "", input[i])) }

#check input type for numerics and booleans
magnify    <- as.numeric(magnify)
h          <- as.numeric(h)
w          <- as.numeric(w)
dpi        <- as.numeric(dpi)
percent    <- as.logical(percent)
logx       <- as.logical(logx)
logy       <- as.logical(logy)
mids       <- as.logical(mids)
showlegend <- as.logical(showlegend)

main(infile, outfile, title, bars, percent, xlab, ylab, xlim, ylim, logx, logy, mids, magnify, showlegend, legendpos, h, w, d)

