stacked.barplot <- function(dm=NULL, minlen=1, maxlen=37, start=1, end=1e+07, internal.margins=c(0,0,0,1), skip.x=2, bicol=NULL, col.fun=rainbow, axis.col="black", main.col="black", main.adj=1, samey=FALSE, ...) {

	# filter input columns
	mr <- as.integer(rownames(dm))
	mc <- as.integer(colnames(dm))
	dm <- dm[,mc>=minlen&mc<=maxlen]
	dm <- dm[mr>=start&mr<=end,]

	# reset minlen and maxlen in case
	# they were inside the range anyway
	minlen <- min(colnames(dm))
	maxlen <- min(colnames(dm))

	# figure out how many graphs we have
	ngraph <- ncol(dm)

	# warn if it is lots!
	if (ngraph>7) {
		print("WARNING: you have chosen to plot a lot of graphs to the current device")
		print("WARNING: this may be too many for the device to cope with")
		print("WARNING: if you get an error message, then consider using a bigger device")
		print("WARNING: with bigger dimensions: see ?png ?jpeg ?pdf")
	}

	# create ngraph rows in the device
	split.screen(c(ngraph,1))

	# set the colours
	if (is.null(bicol)) {
		cols <- col.fun(ncol(dm))
	} else {
		cols <- rep(bicol, ncol(dm))
	}

	# iterate over the number of graphs
	for (i in 1:ngraph) {


		maxy <- max(dm[,i])
		miny <- min(dm[,i])

		if (samey==TRUE) {
			maxy <- max(dm)
			miny <- min(dm)
		}

		# select the relevant screen
		screen(i)
	
		# set the margins inside each graph
		par(mar=internal.margins)

		# plot the graph
		plot(rownames(dm), dm[,i], col=cols[i], type="l", cex.axis=0.7, bty="n", xaxt="n", ylim=c(miny, maxy), ...)

		# only plot x-axis every skip.x graph
		if (i%%skip.x==0) {
			axis(side=1, cex.axis=0.7, col=axis.col, col.axis=axis.col)
		}

		# plot a title for each graph
		title(colnames(dm)[i], line=-1, cex.main=0.7, col.main=main.col, adj=main.adj)
	}

	# close the screens
	close.screen(all=TRUE)

}