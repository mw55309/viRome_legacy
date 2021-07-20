
size.strand.bias.plot <- function(bp=NULL, minlen=1, maxlen=37, line.col="red", sym.axes=TRUE, title="Strand bias plot", xlab="+ strand counts", ylab="- strand counts", lty=1, lwd=2, cextxt=1, mar=c(5,4,4,1), tpos=1, ...) {

	# limit the data according to user input
	bp <- bp[bp$lgth>=minlen&bp$lgth<=maxlen,]

	# figure out the axes we will need
	# depending on whether we want
	# symmetrical X and Y axes
	if (sym.axes==TRUE) {
		miny <- min(bp[,c("pos","neg")])
		maxy <- max(bp[,c("pos","neg")])
		minx <- miny
		maxx <- maxy
	} else {
		miny <- min(bp[,c("neg")])
		maxy <- max(bp[,c("neg")])
		minx <- min(bp[,c("pos")])
		maxx <- max(bp[,c("pos")])
	}

	# basic scatter plot
	par(mar=mar)
	plot(bp$pos, bp$neg, main=title, xlab=xlab, ylab=ylab, xlim=c(minx,maxx), ylim=c(miny, maxy), ...)

	# add the line of y=x
	abline(0,1,col=line.col, lwd=lwd, lty=lty)

	# label the points
	text(bp$pos, bp$neg, bp$lgth, cex=cextxt, pos=tpos)


}
