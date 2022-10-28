#!/usr/bin/env Rscript


"Usage: 2.2-get_island_interactive.R <idDivFile> [ --plot <ids.list> ] [ --mab <abID>... ] [ --outdir <output/tables> --output <islandSeqs.txt> --reference <idDivFile> --plotmethod <plotmethod>]

Options:
  -h --help  	             Show this documentation.
  <idDivFile>	             Output from 2.1-calulate_id-div.pl with data for plotting.
  --plot <ids.list>          Specify a subset of sequences to be plotted (eg those from
                                 a particular V gene). One id per line.
  --mab <abID>		     Identity referent to use on y-axis. Default is the first mAb
                                 in the first id-div.tab file. If multiple referents are
                                 specified, each will be processed in order. Use '--mab =a'
				 to look at all antibodies in the idDivFile.
  --outdir <output/tables>   Directory for saving output. [default: output/tables]
  --output <islandSeqs>      File name of output in outdir. Automatically appends '.txt' for
  	   		     	 list of sequences and '.png' for image. 
				 [default: islandSeqs]
  --reference <idDivFile>    Other data points to display on plots, eg other members of
                                 the antibody lineage of interest.
  --plotmethod <plotmethod>  Plotting method to use. 'original' uses kernel density
                                 estimation to show smoothed distribution of sequences.
                                 'binned' directly plots counts of sequences
                                 for tiled ID/DIV regions [default: original]


Created by Chaim A Schramm 2016-08-30.

Copyright (c) 2016 Columbia University and Vaccine Research Center, National
                              Institutes of Health, USA. All rights reserved.

" -> usage


#### WRAPPER FOR gglocator() TO LIMIT TO POINTS WITHIN PLOT ####
getLocation <- function(xmin=-1,xmax=50,ymin=50,ymax=101) {
  
  myPoint <- data.frame( x=c(xmin-1), y=c(ymax+1) )
  while( myPoint$x<xmin || myPoint$x>xmax || myPoint$y<ymin || myPoint$y>ymax ) {
    myPoint <- gglocator( n=1,xexpand=c(0,0),yexpand=c(0,1) )
    if (is.na(myPoint$x)) { myPoint <- data.frame( x=c(xmin-1), y=c(ymax+1) ) }
  }

  return( myPoint )
  
}


#manually adding @mvkorpel's proposed fixes to ggmap's gglocator function
# https://github.com/dkahle/ggmap/pull/119/commits/9fd4c29d26f0a08d99608a502f90c472f44d3af2
gglocator <- function (n = 1, message = FALSE, xexpand = c(0.05, 0), yexpand = c(0.05, 0)) 
{
    if (n > 1) {
        df <- NULL
        for (k in 1:n) {
            df <- rbind(df, gglocator(message = message, xexpand = xexpand, 
                yexpand = yexpand))
        }
        return(df)
    }

    object <- last_plot()
    if(is.null(object)){
        stop("no plots available")
    }
    
    x <- unlist(current.vpTree())
    x <- unname(x[grep("\\.name$", names(x))])
    x <- grep("panel", x, fixed = TRUE, value = TRUE)
    n_panels <- length(x)
    if(n_panels == 0){
        stop("ggmap plot not detected in current device")
    }
    if(n_panels > 1){
        x <- x[1]
        warning(gettextf("multiple plots detected, choosing one (\"%s\")",
                          x), domain = NA)
    }
    previous_viewport <- current.vpPath()
    seekViewport(x, recording = FALSE)
    
    # when exiting function, return to previous position in viewport tree
    on.exit(upViewport(0, recording = FALSE))
    if(!is.null(previous_viewport)){
        on.exit(downViewport(previous_viewport, strict = TRUE, recording = FALSE),
                 add = TRUE)
    }
    
    # get the position relative to that viewport
    loc <- as.numeric(grid.locator("npc"))

    # get the x.range and y.range from ggplot
    plot_info <- ggplot_build(object)

    if("layout" %in% names(plot_info)){
        ranges <- plot_info$layout$panel_ranges[[1]]
    } else{
        ranges <- plot_info$panel$ranges[[1]]
    }

    xrng <- ranges$x.range
    yrng <- ranges$y.range

    if(is.null(ranges)) {
    	#ggplot has changed this yet again
	# https://gist.github.com/tomhopper/9076152 seems like a good resource to keep track of
	xrng <- plot_info$layout$panel_params[[1]]$x.range
	yrng <- plot_info$layout$panel_params[[1]]$y.range
    }
    
    xrng <- expand_range(range = xrng, mul = xexpand[1], add = xexpand[2])
    yrng <- expand_range(range = yrng, mul = yexpand[1], add = yexpand[2])
    point <- data.frame(x = xrng[1] + loc[1] * diff(xrng), y = yrng[1] + loc[2] * diff(yrng))

    #this is more general, but it's now returning '~x' instead of x, so just fixing it manually
    #names(point) <- with(object, c(deparse(mapping$x), deparse(mapping$y)))

    point
}


#### SURVEYORS FORMULA ####
edgeArea <- function( coords ) {
  return( (coords['x2'] - coords['x1'])*(coords['y2'] + coords['y1']) )
}


#### DETERMINE CW vs CCW ####
checkPoly <- function( point_list ) {
  numPoints    <- dim(point_list)[1]
  edges        <- cbind(point_list[c(numPoints,1:numPoints-1),], point_list)
  names(edges) <- c('x1','y1','x2','y2')
  return(sum( apply(edges,1,edgeArea) ))
}


#### CHECK FOR LEFT CLICK VS RIGHT CLICK ####
getClick <- function ( buttons, x, y ) { 
  if ( buttons[1] == 1 ) {
        NULL
  } else {
    	myClick <<- buttons[1]
        return( invisible(1) )
  }
}


#### GET DISTANCE BETWEEN POINTS (DID USER CLOSE SELECTION POLYGON?) ####
distance <- function (point1, point2) {
  if ( is.null(point1$x) || is.null(point2$x) ) {
      return (100)
  } else {
      return ( sqrt( (point1$x-point2$x)^2 + (point1$y-point2$y)^2 ) )
  }
}


#### MAIN FUNCTION ####
getIsland <- function (dataFile, subsetFile, natAbList, outDir, outFile, refPointsData ) {

  if (.Platform$OS.type == "windows") { 
     win.graph()
  } else {
  misc<-X11Font("-misc-fixed-%s-%s-*-*-%d-*-*-*-*-*-*-*")
  X11Fonts(misc=misc)
  X11(type="Xlib",family="misc")
  }
  
  setGraphicsEventHandlers(prompt = "Left click to accept, right click to start over",
                           onMouseDown = getClick,
                           onMouseUp = NULL,
                           onKeybd = NULL) 
  eventEnv <- getGraphicsEventEnv()
  

  if ( ! is.null(refPointsData) ) {
    refData <- read.table( refPointsData, h=T )
  }

  #specify the ID column as "character" so that we don't have to worry about mismatches between strings and ints
  bigdata   <- read.table( dataFile, h=T, colClasses=c("sequence_id"="character") )
  smalldata <- bigdata
  if ( ! is.null(subsetFile) ) {
    subset    <- scan(subsetFile, what='character', quiet=T)
    smalldata <- bigdata[which(bigdata$sequence_id %in% subset),]
  }
  
  #get a default reference antibody if one was not provided
  if ( is.null(natAbList) )  natAbList <- c( names(bigdata)[4] )
  if ( "=a" %in% natAbList ) natAbList <- c( names(bigdata)[4:length(names(bigdata))] )
  
  lineageReads          <- smalldata[ 0, ] #just get column headers
  lineageReads$referent <- factor(levels=natAbList)
  lineageReads$round    <- factor(levels=c("current","previous"))
  idsOnly               <- c()
  allPlots              <- list() #to avoid recalculating density plots at the end
  
  for ( mab_num in seq_along(natAbList) ) {

      mab   <- natAbList[mab_num]
      mab.R <- make.names(mab)
    
      if ( ! mab.R %in% names(bigdata) ) stop( paste0("Antibody ",mab," was not found in file ", dataFile) )
    
      #generate initial plot; supress color bar and increase size of plot title from default
      #   keep title separate, because we'll want to use different titles at different stages
      pp <- plot_all(smalldata, mab.R, mab, "germline V", plotmethod=opts$plotmethod) + guides(fill=F) + theme( plot.title=element_text(size = 18) )

      #want referents to look different on interactive and final figure (mostly about size)
      #   so save this first, then add
      allPlots[[mab_num]] <- pp

      if ( ! is.null(refPointsData) ) {
      	 p  <- pp + geom_point( data=refData, aes_string(x="germ_div",y=mab.R, z=NA), shape=5, stroke=2, colour="black" )
      } else {
      	 p  <- pp
      }
      t <- labs( title=sprintf("Click to draw a border around an island for %s\nClick original point to complete",mab) )

      myClick <<- 2 #value returned by right click 
    
      while( myClick > 1) {

        print(p + t)

        # get first point
        first_point <- data.frame(x=c(),y=c())
        point_list  <- data.frame(x=c(),y=c())
        last_point  <- getLocation()
        
        # get more points until selected region is closed
        while(distance(first_point,last_point) > 2) { #probably need to test/adjust this tolerance
          point_list <- rbind( point_list, last_point )
          if( length(point_list$x)<2 ) { #special treatment for first round
            first_point <- last_point
	    print( p + t + geom_point(inherit.aes=F,data=point_list,aes(x,y), size=2) )
          } else {
            print( p + t + geom_path(inherit.aes=F,data=point_list,aes(x,y), colour="red", size=2) + 
                           geom_point(inherit.aes=F,data=point_list,aes(x,y), size=2) )
          }
	  temp_point <- getLocation()
	  while(distance(last_point,temp_point) < 2) { #prevent double clicks from blowing things up
	  	temp_point <- getLocation()
	  }
	  last_point <- temp_point
         }

        # ptinpoly closes path automatically, but geom_path doesn't
        box <- rbind( point_list, first_point )
        
        #ptinpoly only works correctly in 2D if points are selected in counterclockwise order...
        if( checkPoly( point_list ) > 0 ) {
          point_list <- point_list[rev(rownames(point_list)),]
        }
	
        included           <- smalldata[ which(pip2d( as.matrix(point_list), as.matrix(cbind(smalldata$germ_div,smalldata[[mab.R]])) )==1), ]
	included$referent  <- rep( mab, length(included$sequence_id) )
	included$round     <- as.factor( rep( "current", length(included$sequence_id) ) )
	lineageReads$round <- rep( "previous", length(lineageReads$sequence_id) )
	tempLineageReads   <- rbind( included, lineageReads )


	if (length(tempLineageReads$germ_div)==0) {
	   #no reads included from the first mAb
	   pp <- p + geom_path(inherit.aes=F, data=box, aes(x,y), col="red", size=2 ) + 
                     geom_point(inherit.aes=F, data=box, aes(x,y) )
	} else {
           pp <- p + geom_point( inherit.aes=F, data=tempLineageReads, aes_string(x="germ_div",y=mab.R,colour="round") ) +
                     theme( legend.position="bottom", legend.key=element_blank() ) + 
              	     scale_colour_manual( values=c("previous"="darkgoldenrod", "current"='magenta'),
                                       	  breaks=c("previous","current"), labels=c("previous","new"), name="" ) + 
                     geom_path(inherit.aes=F, data=box, aes(x,y), col="red", size=2 ) + 
                     geom_point(inherit.aes=F, data=box, aes(x,y) )
	}


	#pick appropriate title
	if ( length(included$sequence_id)==0 ) {
	    print( pp + labs( title="WARNING: NO POINTS FOUND IN CURRENT REGION\nLeft click to accept anyway, right click to reselect region" ) + 
                        theme( plot.title=element_text(colour="red") ) )
	} else if( length(idsOnly)==0 ) {
	    print( pp + labs( title=sprintf("Found %d transcripts in region\nLeft click to accept, right click to reselect region for %s",length(included$sequence_id), mab) ) +
	      	        guides(colour=F) )
	} else {
	    print( pp + labs( title=sprintf("%d transcripts previously identified\nFound %d transcripts (%d new) in current region\nLeft click to accept, right click to re-select region for %s", length(idsOnly), length(included$sequence_id), sum( ! included$sequence_id %in% idsOnly ), mab) ) )
	}
	
	#getGraphicsEvent doesn't work because it checks interactive() which is always FALSE when using Rscript
	#   However, we can just call the core function directly
 	.External2(grDevices:::C_getGraphicsEvent, "Left click to accept, right click to re-select region")

    }

    # Add current selection to total and remove duplicates (for multiple Abs)
    lineageReads<-tempLineageReads
    idsOnly <- unique( c(idsOnly, included$sequence_id) )

  }

  #generate a final set of plots for inspection
      dev.off()
      pl <- list()
      titlePlot <- ggplot( data.frame(x=1,y=1,text=sprintf( "Final Selections (total %d transcripts)",length(idsOnly) )) ) +
      		   	   geom_text(aes(x,y,label=text), size=4) + scale_x_continuous(limits=c(0,2)) + scale_y_continuous(limits=c(0,2)) +
			   geom_point( data=lineageReads, aes_string(x="germ_div",y=make.names(natAbList[1]),colour="referent"),size=0 ) +
			   guides( colour=guide_legend(nrow=1, title="", override.aes=list(size=3), order=1) ) + theme_void() + 
			   theme( legend.position="bottom", legend.key=element_blank(), legend.text=element_text(size=8),
			   	  legend.margin=unit(0,"in") )
      if ( ! is.null(refPointsData) ) {
      	   titlePlot <- titlePlot + geom_point( inherit.aes=F, data=data.frame(x=0,y=0,shape="reference antibodies"), aes(x,y,shape=shape), size=0, stroke=0  ) +
	   	     		    scale_shape_manual( values=c(5) ) +
	   	     		    guides( shape=guide_legend(nrow=1, title="", override.aes=list(size=2,stroke=.5), order=2) )
      }
      pl[[1]] <- titlePlot
      for( i in seq_along(allPlots) ) {
      	   thisAb  <- natAbList[i]
	   current <- which( lineageReads$referent == thisAb )	   
	   lineageReads$round <- rep( "previous", length(lineageReads$sequence_id) )
	   lineageReads$round[current] <- rep( "current", length(current) )
	   lineageReads <- lineageReads[ order(lineageReads$round, decreasing=T), ] #put reads from current mAb at the bottom so they plot on top
	   pl[[i+1]] <- allPlots[[i]] + geom_point( inherit.aes=F, data=lineageReads, aes_string(x="germ_div",y=make.names(thisAb),colour="referent"), size=0.5 ) +
		  		        guides(colour=F) + labs( title=thisAb ) + theme( plot.title=element_text(size = 10) )
	   if ( ! is.null(refPointsData) ) {
      	      	  pl[[i+1]] <- pl[[i+1]] + geom_point( inherit.aes=F, data=refData, aes_string(x="germ_div",y=make.names(thisAb)), shape=5, stroke=.5, size=.75, colour="black" )
      	   }
      }
      myLayout <- layoutGrid( TRUE, FALSE, length(natAbList), 1, FALSE )

      if (.Platform$OS.type == "windows") { 
      	 win.graph()
      } else {
         X11()
      }
      
      png(sprintf("%s/%s.png",outDir,outFile),h=3,w=2*length(natAbList),unit="in",res=300)
      multiplot( plotlist=pl, layout=myLayout )
      #ggsave(sprintf("%s/%s.png",outDir,outFile),multiplot( plotlist=pl, layout=myLayout ),h=3,w=2*length(natAbList),dpi=300)
      #Sys.sleep(10)
      dev.off()
  
  #output
  write.table( idsOnly, file=sprintf("%s/%s.txt",outDir,outFile), quote=F, row.names=F, col.names=F )
  print( sprintf( "Saved %d transcript IDs to %s", length(idsOnly), file.path(outDir,outFile) ) )
  
}



##### ENTRY FROM COMMAND LINE #####

library(docopt)
library(ptinpoly)
#library(ggmap)
library(scales)
library(ggplot2)

cmd.line <- commandArgs()

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", cmd.line[grep(file.arg.name, cmd.line)])
script.basename <- dirname(script.name)

source(file.path(script.basename, "..", "plottingFunctions.R"))
source(file.path(script.basename, "..", "logging.R"))

opts<-docopt(usage, strip=T, help=T)

# check to make sure the files are there
if ( ! file.exists(opts$idDivFile) )                              stop("Cannot find input file...\n\n")
if ( ! is.null(opts$plot)      && ! file.exists(opts$plot) )      stop("Cannot find plot file...\n\n")
if ( ! is.null(opts$reference) && ! file.exists(opts$reference) ) stop("Cannot find reference data file...\n\n")


if (dir.exists("output/logs")) { 
  saveCommandLine( "output/logs/command_history.log", cmd.line )
} else {
  saveCommandLine( "SONAR_command_history.log", cmd.line )
}
options("error"=my.error.fun)

getIsland(opts$idDivFile, opts$plot, opts$mab, opts$outdir, opts$output, opts$reference)

if (dir.exists("output/logs")) { logSuccess() }
