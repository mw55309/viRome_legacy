size.position.heatmap <- function(dm=NULL, minlen=1, maxlen=37, start=1, end=1e+07, scale=TRUE, col.fun=colorRampPalette(c("black","red"), space="rgb"), log=FALSE, mar=c(5,4,4,2), main=NULL) {
	
	# filter input columns
	mr <- as.integer(rownames(dm))
	mc <- as.integer(colnames(dm))
	dm <- dm[,mc>=minlen&mc<=maxlen]
	dm <- dm[mr>=start&mr<=end,]

	# reset minlen and maxlen in case
	# they were inside the range anyway
	minlen <- min(colnames(dm))
	maxlen <- min(colnames(dm))


	if (log==TRUE) {
		dm <- log(dm+1)
	}

	# scale it if the user wants
	if (scale==TRUE) {
		dm <- t(scale(t(dm)))
		dm[is.na(dm)] <- min(dm[!is.na(dm)])
	}


	# get rownames as vector of integers
	r <- as.integer(rownames(dm))

	# heatmap of read lengths
	par(mar=mar)
	image(z=as.matrix(dm), col=col.fun(32), xaxt="n", yaxt="n", y=1:ncol(dm), x=min(r):max(r), xlab="", ylab="")
	axis(side=2, at=1:ncol(dm),  labels=colnames(dm))
	axis(side=1)
	if (! is.null(main)) {
		title(main)
	}

}
