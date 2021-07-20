position.barplot <- function(vdf=NULL, minlen=1, maxlen=37, reflen=10000, samp="", plot=TRUE, poscol="red", negcol="green") {


	# limit the input to only contain reads
	# of size between minlen and maxlen
	vdf <- vdf[vdf$cliplen>=minlen&vdf$cliplen<=maxlen,]

	# get the reference name and size range as strings
	refname = vdf$rname[1]
	maps <- paste(minlen, "-", maxlen, sep="")

	
	# get counts of reads at each position for 
	# negative strand
	npos <- .count.reads.by.position(vdf, "-")

	# get counts of reads at each position for 
	# positive strand
	ppos <- .count.reads.by.position(vdf, "+")

	# both of the above may "miss" positions out
	# due to zero counts so we create a dummy
	# data.frame with all positions in it
	# and outer-join each of npos and ppos to it
	# filling in resulting NS's with 0 counts
	lgths <- data.frame(position=1:reflen)

	lpos <- merge(lgths, ppos, all.x=TRUE, sort=TRUE)
        lneg <- merge(lgths, npos, all.x=TRUE, sort=TRUE)
      
        lpos[is.na(lpos)] <- 0
        lneg[is.na(lneg)] <- 0

	# sort each of the new merged data.frames
	# by position
	lpos <- lpos[order(lpos$position),]
	lneg <- lneg[order(lneg$position),]
      
	# if the users has decided to plot
	if (plot == TRUE) {

		# calculate the range to be plotted
       		mx <- max(c(lpos[,2],lneg[,2]))
        	rng <- c(-1*mx,mx)

		# create a title for the plot
		tit   <- paste(samp,"Distribution of",maps,"bp maps along",refname, sep=" ")

		# create the barplot as counts on the positive strand and 0-counts on the negative strand
		barplot(as.matrix(lpos[,2]), beside=TRUE, names=lpos[,1], col=c(poscol), xlab="Nucleotide position", ylab="Count", border=poscol, ylim=rng, xaxt="n", main=tit)
        	barplot(as.matrix(0-lneg[,2]), beside=TRUE, add=TRUE, col=c(negcol),  border=negcol)

		# add an axis
        	axis(side=1)
	}

	# merge the positive and negative count data.frames
	# and give them useful column names and return to the user
	ret <- merge(lpos, lneg, by="position", sort=FALSE)
	colnames(ret) <- c("position","poscount","negcount")
	return(ret)
}


.count.reads.by.position <- function(vdf=NULL, str=NULL) {

	# subset by strand
	svdf <- vdf[vdf$strand==str,]

	if (nrow(svdf) > 0) {
		# if there are any mappings on this strand
		# count the number of reads by position
		bn <- by(svdf$qname, svdf$pos, length)
		
		# create data.frame and return it
		ret <- data.frame(position=as.integer(names(bn)), count=as.vector(bn))
		return(ret)

	} else {
		# otherwise just return a dummy data.frame
		return(data.frame(position=1,poscount=0))
	}

}
