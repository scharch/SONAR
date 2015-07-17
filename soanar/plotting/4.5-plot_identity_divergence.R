#!/usr/bin/Rscript


usage <-function() {
cat("
USAGE: longitudinal_double_sieve.r <div.list> <sieve1.list> <sieve2.list> [unique=unique.list 
       				   mab=VRC08 gene=VH1-2 output=lds.pdf labels=dates.list 
				   black=sieve1b.list]

The first 3 parameters are lists of other file names.
	<div.list>	List of the divid.tab pipeline outputs for each time point, in 
			     chronological order
	<sieve1.list>	List of the files containing the lists of reads to plot in the 
			     first sieve (chrono order!)
			     **FIRST LINE IS FORMAT OF THE FILES, NOT A FILE NAME**
			     1 = 'naked' list of reads in single column format
			     2 = tabulated list, with column headers and ids in (for now) 
			       	 the 2nd column (1-indexed!)
	<seive2.list>	List of the files containing the lists of reads to plot in the 
			     first sieve (chrono order!)
			     **FIRST LINE IS FORMAT OF THE FILES, NOT A FILE NAME** (as above)
			     **MULTIPLE FILES CAN BE COMBINED FOR EACH TIME POINT separated by
			     	 a comma (no space!)**

Optional parameters **MUST** be specified in the 'var=value' format.
	<unique>	List of files containing lists of unique reads at each timepoint (to 
			     	plot instead of using raw/bulk data).
				Rules are same as sieve1.list above.
	<black>		List of files containing lists of partially sieved reads at each 
			     	timepoint to plot as black dots under the red dots in sieve2.
				Rules are same as sieve2.list above.
	<mab>		Identity referent to use on y-axis. Default = 1st one in the divid.tab
				files
	<gene>		Specify x-axis label as '% divergence from <gene>'.
				Default = '% divergence from germline V'.
	<output>	File name to save plots to. Default = lds.pdf.
	<labels>	List of timepoints for labeling the longitudinal plots.

")
}


####################################
# PLOTTING FUNCTION
####################################

plot_all <- function (data, native, pretty, my_x, subset=c(), plot_title=NULL, xlims=c(0,30), ylims=c(70,101), 
	 	    	     color=TRUE, number=TRUE, xlabel=FALSE, ylabel=FALSE, contour = FALSE, conCol = 'black') {

	#print(c(native, pretty, my_x, plot_title, xlims, ylims, color, number, xlabel, ylabel))

	my_y<-paste("% ID to", pretty)
	if (color) my_colors=rev(rainbow(15,end=4/6))
	else my_colors=rev(gray.colors(5))

	if (length(subset) > 0)	 {
		data<-data[which(data$read_id %in% subset),]
	}
	data.labels <- data.frame(x=45,y=85, label=format(length(data$germ_div),big.mark=","))

	p<-ggplot(data,aes_string(x="germ_div",y=native))

	if (! contour) {
		p <- p + stat_bin2d(binwidth=c(1,1)) +
		scale_fill_gradientn(colours=my_colors, trans="log10", limits=c(1,25000),na.value="white", guide="none")
	} else {
	        library(MASS)
		#add a smidge to default bandwidth to handle occassional datasets with super redundancy
	        g <- kde2d(data$germ_div,data[[native]],n=101,h=1,lim=c(0,100,0,100))
		#print(c(min(g$z[g$z>10^-10]),log10(max(g$z))))
		densf = data.frame(expand.grid(x=g$x,y=g$y), z=as.vector(g$z))
		r<-10^seq(-4,log10(max(g$z))+1,1.5)
		p <- p + stat_contour(data=densf, aes(x=x,y=y,z=z), colour=conCol, size=.25, breaks=r)
	}

	p <- p + theme_bw() + labs(x=NULL, y=NULL, title=plot_title) +
	     theme(plot.background = element_blank(),panel.grid.major = element_blank(),
				axis.ticks.length = unit(.02,"in"), axis.ticks = element_line(size = .5),
		          	panel.grid.minor = element_blank(), axis.text = element_blank(),
				plot.margin = unit(c(.02,.02,.02,.02),"in") )

	if ( xlabel ) p <- p + scale_x_continuous(expand=c(0,0),limits=xlims, my_x)
	else p <- p + scale_x_continuous(expand=c(0,0),limits=xlims)

     	if ( ylabel ) p <- p + scale_y_continuous(expand=c(0,0),limits=ylims, my_y)
	else p <- p + scale_y_continuous(expand=c(0,0),limits=ylims)

	if (number) p <- p + geom_text(data=data.labels, aes(x=x, y=y, label=label), angle=90, size=3, col="black")

	p
}


##########################################################################################
#         MAIN PROGRAM
##########################################################################################

lds <- function(file1, file2, file3, unique=NULL, nat="", gene="germline V", output="lds.pdf", labels=NULL, black=NULL) {

#get the names of the data files
div_files = scan(file1,what='character',quiet=T)
s1 = scan(file2,what='character',quiet=T)
s2 = scan(file3,what='character',quiet=T,blank.lines.skip=F)

#process formats
s1_type  = as.numeric(s1[1])
s1_files = s1[-1]
s2_type  = as.numeric(s2[1])
s2_files = strsplit(s2[-1],",")

if (length(div_files) != length(s1_files) || length(div_files) != length(s2_files)) {
	stop("Error: number of timepoints does not match across sieve levels. Please check file formats and try again")
}

if (! is.null(unique) && file.exists(unique)) {
	u  = scan(unique, what='character',quiet=T)
	ut = as.numeric(u[1])
	uf = u[-1]
	if (length(div_files) != length(uf)) { stop ("Error: number of timepoints in unique file does not match. Please check file format and try again") }
} else { ut = NULL }

if (! is.null(black) && file.exists(black)) {
	b  = scan(black, what='character',quiet=T,blank.lines.skip=F)
	bt = as.numeric(b[1])
	bf = strsplit(b[-1],",")
	if (length(div_files) != length(bf)) { stop ("Error: number of timepoints in black file does not match. Please check file format and try again") }
} else { bt = NULL }


#initiate
myplots=list()


#process each time point
for (time in seq(length(div_files))) {

    	#main plot
	data <- read.table( div_files[time], h=T )
	setU <- c()

	if (! is.null(ut)) setU <- switch( ut, scan(uf[time],quiet=0), read.table(uf[time],h=T)[,2] )

	if (mab == "") mab <- names(data)[3] 

	myplots[[length(myplots)+1]] <- plot_all(data, mab, mab, gene, subset=setU, xlims=c(0,50), ylims=c(50,101))#, ylabel=(time==1))

	#sieve1 plot
	set1 <- switch( s1_type, scan(s1_files[time],quiet=0), read.table(s1_files[time],h=T,fill=T)[,2] )
	ind1 <- which(data$read_id %in% set1)
	data.labels <- data.frame(x=45,y=85, label=format(length(ind1),big.mark=","))
	p <- plot_all(data, mab, mab, gene, subset=setU, color=FALSE, number=FALSE, contour=TRUE, xlims=c(0,50), ylims=c(50,101))#, ylabel=(time==1))
	myplots[[length(myplots)+1]] <- p + geom_point(data=data[ind1,],colour='blue', shape=16, size=.5) + 
				     geom_text(data=data.labels, aes(x=x, y=y, label=label), angle=90, size=3, col="blue")

	#sieve2 plot
	set2<-c()
	for (group in s2_files[[time]]) set2 <- c(set2, switch( s2_type, scan(group,quiet=0), read.table(group,h=T,fill=T)[,2] ))
	ind2 <- which(data$read_id %in% set2)
	data.labels <- data.frame(x=45,y=90, label=format(length(ind2),big.mark=","))
	p <- plot_all(data, mab, mab, gene, subset=set1, color=FALSE, number=FALSE, contour=TRUE, conCol='cornflowerblue', xlims=c(0,50), ylims=c(50,101))#, ylabel=(time==1), xlabel=TRUE)

	if (! is.null(bt)) {
	      setB <- c()
	      for (group in bf[[time]]) setB <- c(setB, switch( bt, scan(group,quiet=0), read.table(group,h=T,fill=T)[,2] ))
	      indB <- which(data$read_id %in% setB)
	      data.Blabels <- data.frame(x=45,y=70, label=format(length(indB),big.mark=","))
	      p <- p + geom_point(data=data[indB,],colour='goldenrod', shape=16, size=1) +
	      	       		geom_text(data=data.Blabels, aes(x=x, y=y, label=label), angle=90, size=3, col="goldenrod")
	}

	myplots[[length(myplots)+1]] <- p + geom_point(data=data[ind2,],aes_string(x="germ_div", y=mab),colour='darkorchid4', shape=16, size=1) + 
				     geom_text(data=data.labels, aes(x=x, y=y, label=label), angle=90, size=3, col="darkorchid4")

}

	row_titles <- c("All", "Sieve1", "Sieve2")
	if (! is.null(labels)) col_titles = scan(labels,quiet=T)

	#and print to file	
	pdf(output, height=3, width=length(div_files))
	multiplot(plotlist=myplots,rows=3,layout=matrix(seq(2, 1+3*length(div_files)),nrow=3, byrow=F))#, row_titles=row_titles)
	dev.off()
#break

}



##########################################################################################
#         ENTRY POINT FROM COMMAND LINE
##########################################################################################

library(ggplot2)
library(grid)
library(MASS)
source('/ifs/home/c2b2/ls_nih/cs3037/scripts/multiplot.R')

# args are: <div.list> <sieve1.list> <sieve2.list> [unique=unique.list mab=VRC08 gene=VH1-2 output=lds.pdf labels=dates.list black=sieve1b.list]
input = commandArgs(T)

#check for parameters
if(length(input)<3){
  stop(usage())
} else if (length(input)>8) {
  warning("I only know what to do with the first 8 parameters; ignoring the rest.....\n\n")
}

#now check to make sure the files are there
if (! all(file.exists(input[1:3]))) {
  stop("Please check input parameters and make sure all files exist!\n\n")
}


#set defaults and process optional parameters
unique<-NULL
mab<-""
gene<-"germline V"
output<-"lds.pdf"
labels<-NULL
black<-NULL

for (i in seq(4,length(input))) { assign(sub("=.*$", "", input[i]), sub("^[^=]*=", "", input[i])) }

lds(input[1], input[2], input[3], unique, mab, gene, output, labels, black)

