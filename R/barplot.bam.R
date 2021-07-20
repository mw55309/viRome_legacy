barplot.bam <- function(vdf=NULL, minlen=1, maxlen=37, poscol="red", negcol="green", main="Sequence length distribution",
				xlab="Map length", ylab="Count", legend=c("+ve strand","-ve strand"), legendx=NULL, legendy=NULL, plot=TRUE, down=FALSE, sym.axes=TRUE, ...) {


	# count frequencies in the negative strand
	neg <- .sum.bam.length(vdf, "-", minlen, maxlen)
	

	# count frequencies on the positive strand
	pos <- .sum.bam.length(vdf, "+", minlen, maxlen)
	
	# join the two on column "length" and give the result
	# sensible column names
	bavdf <- merge(pos, neg, by="length", sort=FALSE)
	colnames(bavdf) <- c("lgth","pos","neg")

	# if the user wants a plot
	if (plot==TRUE) {


		# draw the barplot of pos and neg frequencies
		# and with the user specified options
		if (down==TRUE) {

			# if the user wants symmetrical y-axes
			if (sym.axes==TRUE) {
				miny <- -1 * max(bavdf[,2:3])
				maxy <- max(bavdf[,2:3])
			} else {
				# else let the data decide
				miny <- -1 * max(bavdf[,3])
				maxy <- max(bavdf[,2])
			}

			# plot the poistive counts
			barplot(as.matrix(t(bavdf[,2])), beside=TRUE, names=bavdf[,1], col=c(poscol), main=main, xlab=xlab, ylab=ylab, ylim=c(miny,maxy), ...)

			# add zero minus the negative counts
			barplot(as.matrix(0-t(bavdf[,3])), beside=TRUE, col=c(negcol), add=TRUE, ...)
		} else {
			barplot(as.matrix(t(bavdf[,2:3])), beside=TRUE, names=bavdf[,1], col=c(poscol,negcol), main=main, xlab=xlab, ylab=ylab, ...)
		}


		# set defaults for legendx and legendy
		# if no values are provided
		if (is.null(legendx)) {
			legendx <- 5
		}

		if (is.null(legendy)) {
			legendy <- max(bavdf[,2:3])/2
		}


		# draw the legend
		legend(legendx,legendy,legend=legend,fill=c(poscol,negcol))
	}

	# return the frequencies
	return(bavdf)
	
}


.sum.bam.length <- function(vdf=NULL, str=NULL, minlen=NULL, maxlen=NULL) {

	# first limit the input to the strand we are interested in
	vdf <- vdf[vdf$strand==str,]

	# count the occurrence of each length of read
	# and convert to a data.frame
	tbl <- table(vdf$cliplen, useNA="always")
	tvdf <- as.data.frame(tbl)

	# The above may not contain every number between
	# minlen and maxlen, so we create a dummy
	# data.frame that does and "outer join" it 
	# to the above result, replacing the resulting
	# NA with zero counts
	yvdf <- data.frame(Counts=minlen:maxlen)
	avdf <- merge(yvdf, tvdf, by.x="Counts", by.y="Var1", all.x=TRUE, sort=FALSE)
	avdf[is.na(avdf)] <- 0

	# give it sensible column and row names
	colnames(avdf) <- c("length","freq")
	rownames(avdf) <- avdf$length
	
	# return 
	return(avdf[order(avdf$length),])
}
