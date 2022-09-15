			       
#split from 4.3 by Chaim A Schramm 2016-08-30.

library(ggplot2)
library(grid)
library(MASS)


####################################
# PLOTTING FUNCTION
####################################

plot_all <- function (data, native, pretty, heavy, plot_title=NULL, color=TRUE, guide=TRUE, xlabel=TRUE, ylabel=TRUE, contour=TRUE, conCol = 'black') {

	my_x<-paste("% divergence from", heavy)
	my_y<-paste("% ID to", pretty)

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
			theme_bw() + scale_x_continuous(expand=c(0,0),limits=c(-1,50)) +
			scale_y_continuous(expand=c(0,1),limits=c(50,101)) +
	     		theme(plot.background = element_blank(),panel.grid.major = element_blank(),
				axis.ticks.length = unit(.02,"in"), axis.ticks = element_line(size = .5),
		          	panel.grid.minor = element_blank(), axis.text = element_text(size = 6),
				axis.title = element_text(size = 8), plot.margin = unit(c(.1,.1,.1,.1),"in"),
				plot.title = element_text(size = 8) )

	if ( contour ) 		     { p <- p + stat_contour(colour=conCol, size=.25, breaks=10*r) }
	if ( ! is.null(plot_title) ) { p <- p + labs( title=plot_title ) }
	if ( xlabel ) 		     { p <- p + labs( x=my_x ) } else { p <- p + labs( x="" ) }
	if ( ylabel ) 		     { p <- p + labs( y=my_y ) } else { p <- p + labs( y="" ) }

	p
}



####################################
# LAYOUT FUNCTION
####################################

layoutGrid <- function( title, label, columns, rows, transpose ) {

    #do the easy case first
    if (!title && !label) {
        if (transpose) { layout <- matrix(seq(columns*rows), columns, byrow=T) }
        else { layout <- matrix(seq(columns*rows), rows) }

    } else {
        #some housekeeping to figure out relative sizes
        hExp <- 1
        vExp <- 1
        if (title) { vExp <- 3 }
        if (label) {
            if (transpose) {
                hExp <- 3
            } else {
                vExp <- 3
            }
        }

        #figure out offset for counting graphs
        #  since titles are placed in myplots list first
        offset <- title + label #booleans eval'ed to 0 or 1
        layout <- NULL
        for (x in seq(0,columns-1)) {
            current <- NULL
            for (y in seq(rows)) {

                #get the position of this plot in the multiplot list
                num <- x*rows + y + offset

                if (transpose) {
                    current <- cbind( current, matrix(rep(num,hExp*vExp), vExp) )
                } else {
                    current <- rbind( current, matrix(rep(num,hExp*vExp), vExp) )
                }

            }
            
            if (transpose) { layout <- rbind(layout, current) }
            else           { layout <- cbind(layout, current) }
        }

    }

    #add time labels
    if (label) {
        plotNum <- 1 + title
        if (transpose) { layout <- cbind( matrix( rep(plotNum, columns*vExp) ), layout ) }
        else           { layout <- rbind( matrix( rep(plotNum, columns*hExp), nrow=1 ), layout ) }
    }

    if (title) {
        if (transpose) { layout <- rbind( matrix( rep(1, rows*hExp+label), nrow=1 ), layout ) }
        else           { layout <- rbind( matrix( rep(1, columns*hExp), nrow=1 ), layout ) }
    }

    layout
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
    return(plots[[1]])

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


