#!/usr/bin/Rscript

# USAGE: plotDivHist.r <divergence.txt> <germ_stat.txt>
# 
# now makes Jiang-style usage/divergence plots

input = commandArgs(T)


#start with some housekeeping
b <- basename(input[1])
title <- paste(strsplit(b,"_")[[1]][1],strsplit(b,"_")[[1]][2], sep=", ")
onename = sub("data/(.*)_germ_divergence.txt", "figures/\\1-divergence.png", input[1])
twoname = sub("data/(.*)_germ_divergence.txt", "figures/\\1-usage.png", input[1])


#read in divergence data
cn = scan(input[1],what='character', nlines=1)
dat = read.delim(input[1],header=T,row.names=NULL)
colnames(dat)<-cn

#process each gene
maxdiv=c()
mindiv=c()
avgdiv=c()
alldivs=c()
for (gene in cn[order(cn)] ){
    if (is.na(gene)) next
    maxdiv[gene] = max(as.numeric(dat[[gene]]), na.rm=T)
    mindiv[gene] = min(as.numeric(dat[[gene]]), na.rm=T)
    avgdiv[gene] = mean(as.numeric(dat[[gene]]), na.rm=T)
    alldivs = c(alldivs,as.numeric(dat[[gene]]))
}

#plot divergence histogram
png(onename, height=480, width=640) # For png output
par(mar=c(5,5,4,4)+.1, mgp=c(4,1,0), fig=c(0,1,0,1),cex=1)
hist(alldivs, breaks=seq(0,100,2), xlim=c(0,40),xlab="Germline divergence (%)", ylab="Number of sequences", las=1, col="green", main=title)

#inset to show high divergence bump **ADD FLAG TO MAKE THIS OPTIONAL**
par(mar=c(5,4,4,4)+.1, mgp=c(3,1,0), fig=c(.5,1,.35,.85),new=T,cex=.8)
hist(alldivs, breaks=seq(0,100,2), xlim=c(0,40),ylim=c(0,500),xlab="Germline divergence (%)", ylab="Number of sequences", las=1, col="green", main="")


#now a usage histogram
#start by collecting data
usage=c()
con <- file(input[2], "rt")
line <- readLines(con,1) #dump header
repeat{
	line <- readLines(con, 1)
	if (!length(line)) break
	toks = strsplit(line,"[\t*]")
	usage[toks[[1]][1]] = sum(c(usage[toks[[1]][1]],as.numeric(toks[[1]][3])),na.rm=T)
}
close(con)
genes=names(maxdiv) #because there might be extra genes in the usage data (!?)

#and plot
png(twoname, height=480, width=640) # For png output
par(mar=c(5,5,4,5)+.1,mgp=c(4,1,0),fig=c(0,1,0,1),cex=1)
barplot(usage[genes], axisnames=F, width=.8, space=.25, col="green", ylab="Number of sequences",las=1, ylim=c(0,ceiling((max(usage)+100)/500)*500))
text(seq(0.6,length(genes)-.4,1),par("usr")[3]-2, srt=45, adj=c(1.1,1.1), labels=genes,xpd=T)
box()
par(new=T)
newMaxX=1.25*length(genes) - .25
plot(seq(0.5, newMaxX, 1.25), maxdiv, axes=F, type="o", pch=16, lty=2, xlab="", ylab="", main="", xaxs="r", yaxs="i", xlim=c(0,newMaxX), ylim=c(0,ceiling((max(maxdiv)+1)/5)*5))
lines(seq(.5, newMaxX, 1.25), mindiv, type="o", pch=15, lty=2, col="blue")
lines(seq(.5, newMaxX, 1.25), avgdiv, type="o", pch=17, lty=2, col="red")
axis(side=4, las=1)
mtext("Germline divergence (%)", side=4, line=3)

dev.off() # For png output

