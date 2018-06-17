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
#     - what should default guide behavior be?
#     - change default sizing based on guides and presence/absence of labels/titles
#     - fix title size (adjust based on num plots)
#     - fix y-axis label in transpose mode
#     - "catch" and supress this dumb: Error in UseMethod("grid.draw") : no applicable method for 'grid.draw' applied to an object of class "NULL"
#
########################


"Usage: 4.3-plot_identity_divergence.R <idDivFile>... [ --plot <ids.list>... ] [ --sieve1 <ids.list>... --sieve2 <ids.list>... ] [ --mab <abID>... ] [ --xaxis <vgene> --outdir </path> --output <output.pdf> --title <title> --labels <timecode>... --reference <idDivFile> --showNames <0|1|2> --transpose]

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
  --xaxis <vgene>	   X-axis label is '% divergence from <label>'.
  	  		       [default: germline V]
  --outdir </path>	   Directory for saving output. [default: output/plots]
  --output <output.pdf>	   File name of output. [default: id-div.pdf]
  --title <title>          Title to display at the top of page. (MUST escape any spaces
                               eg \"'Like this'\" or \"This\ is\ my\ title\"!!)
  --labels <timecode>      Labels for the longitudinal time points. Must match number and
  	   		       order of data files.
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
Moved several functions to plottingFunctions.R so they can be imported by
      	      multiple scripts CAS 2016-08-30.

Copyright (c) 2013-2016 Columbia University and Vaccine Research Center, National
                               Institutes of Health, USA. All rights reserved.

" -> usage


##########################################################################################
#         MAIN PROGRAM
##########################################################################################

lds <- function(dataFileList, subsetFileList, firstSieveList, secondSieveList, natAbList, xLabel, outDir, outFile, pageTitle, timepointList, refPointsData, showNames, transpose) {

#initiate
myplots=list()

#page title
if ( ! is.null(pageTitle) ) {
    p <- ggplot(data.frame(x=1,y=1,text=pageTitle)) + geom_text(aes(x,y,label=text),size=15)+theme_void()
    myplots[[length(myplots)+1]] <- p
}

#timepoint labels
if ( ! is.null(timepointList) ) {

   #generate a data frame to plot these from
   x<-seq(0.5, length(timepointList)-.5)
   y<-rep(1, length(timepointList))
   d<-data.frame(x=x,y=y,label=timepointList)

   p<-ggplot(d)+geom_text(aes(x,y,label=label),size=6)+theme_void()+scale_x_continuous(lim=c(0,length(timepointList)),expand=c(0,0))
   if (transpose) { p<-ggplot(d)+geom_text(aes(x=y,y=x,label=label),size=6,angle=90)+theme_void()+scale_y_continuous(lim=c(0,length(timepointList)),expand=c(0,0)) }

   myplots[[length(myplots)+1]] <- p
}


if ( ! is.null(refPointsData) ) {
   refData <- read.table( refPointsData, h=T )
}

#process each time point
for ( time in seq_along(dataFileList) ) {

    	#main plot
        #specify the "ID" column as "character" so that we don't have to worry about mismatches between strings and ints
	bigdata   <- read.table( dataFileList[time], h=T, colClasses=c("ID"="character") )
	smalldata <- bigdata
        if ( ! is.null(subsetFileList) ) {
            subset    <- scan(subsetFileList[time], what='character', quiet=T)
            smalldata <- bigdata[which(bigdata$ID %in% subset),]
        }


        #get a default reference antibody if one was not provided
	if ( is.null(natAbList) )  natAbList <- c( names(bigdata)[3] )


	#only show axis labels on edges
	if (transpose) {
	   if ( time == length(dataFileList) ) {
	      showX = T
	   } else {
	      showX = F
	   }
	} else {
	   if ( time == 1 ) {
	      showY = T
	   } else {
	      showY = F
	   }
	}


	for ( mab_num in seq_along(natAbList) ) {

	    mab   <- natAbList[mab_num]
            mab.R <- make.names(mab)

            if (! mab.R %in% names(bigdata)) {
                stop(paste0("Antibody ",mab," was not found in file ", dataFileList[time]))
            }
            
            #only show axis labels on edges
	    if (transpose) {
	       if ( mab_num == 1 ) {
	       	  showY = T
	       } else {
	       	  showY = F
	       }   
	    } else {
	       if ( mab_num == length(natAbList) ) {
	       	  showX = T
	       } else {
	       	  showX = F
	       }
	    }

	    ####
	    # ADD OPTIONS (plot color, contours, yes/no to axis labels, scale bar)
	    ####
	    p <- plot_all(smalldata, mab.R, mab, xLabel, xlabel=showX, ylabel=showY)

	    if ( ! is.null(firstSieveList) ) {
	       subset <- scan(firstSieveList[time], what='character', quiet=T)
	       s1data <- bigdata[which(bigdata$ID %in% subset),]
	       p <- p + geom_point(data=s1data,aes_string(x="germ_div",y=mab.R,z=NA),colour='goldenrod', shape=16, size=.5)
	    }

	    if ( ! is.null(secondSieveList) ) {
	       subset <- scan(secondSieveList[time], what='character', quiet=T)
	       s2data <- bigdata[which(bigdata$ID %in% subset),]
	       p <- p + geom_point(data=s2data,aes_string(x="germ_div",y=mab.R,z=NA),colour='magenta', shape=16, size=.5)
	    }

	    if ( ! is.null(refPointsData) ) {
	       p <- p + geom_point(data=refData,aes_string(x="germ_div",y=mab.R,z=NA),colour='black', shape=4, size=1)
	       if( showNames > 1 || (showNames == 1 &&  time == 1) ) {
	       	   p <- p + geom_text(data=refData, aes_string(x="germ_div",y=mab.R,z=NA,label="ID"), size=2, colour='black', hjust=2, vjust=1)
	       }
	    }

	    myplots[[length(myplots)+1]] <- p

	} #next mab in list

} #next time point


#figure out how to arrange everything on the page
myLayout <- layoutGrid( !is.null(pageTitle), !is.null(timepointList), length(dataFileList), length(natAbList), transpose )

#figure size (use the boolean transpose as 0 or 1 to get the correct dimensions)
height  <- 2 * (transpose*length(dataFileList) + (!transpose)*length(natAbList))
width   <- 2.5 * ((!transpose)*length(dataFileList) + transpose*length(natAbList))

ggsave(
	paste(outDir,outFile,sep="/"),
	multiplot( plotlist=myplots, layout=myLayout ),
	dpi = 200,
	height = height/2,
	width = width/2,
	limitsize = FALSE
)


}



##########################################################################################
#         ENTRY POINT FROM COMMAND LINE
##########################################################################################

library(docopt)

cmd.line <- commandArgs()

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", cmd.line[grep(file.arg.name, cmd.line)])
script.basename <- dirname(script.name)

source(file.path(script.basename, "..", "plottingFunctions.R"))
source(file.path(script.basename, "..", "logging.R"))

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

if (! is.null(opts$labels)) {
   if ( length(opts$labels) != length(opts$idDivFile) ) { stop("Please make sure that the number of timepoint labels matches the number of input files and that they are in the same order\n\n") }
}

if (! is.null(opts$reference) && ! file.exists(opts$reference)) { stop("Cannot find reference data file...\n\n") }


if (dir.exists("output/logs")) { 
    saveCommandLine( cmd.line )
    options("error"=my.error.fun)
} else { print("SONAR log directory not found; command line and output will not be saved") }

lds(opts$idDivFile, opts$plot, opts$sieve1, opts$sieve2, opts$mab, opts$xaxis, opts$outdir, opts$output, opts$title, opts$labels, opts$reference, opts$showNames, opts$transpose)

if (dir.exists("output/logs")) { logSuccess() }
