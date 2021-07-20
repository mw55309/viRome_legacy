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
clip.bam <- function(vdf = NULL) {

	# set up extra columns with null values
	vdf$softclip <- rep(0,nrow(vdf))
	vdf$leftclip <- rep(0,nrow(vdf))
	vdf$rightclip <- rep(0,nrow(vdf))
	vdf$clipseq   <- as.vector(vdf$seq)

	# we take a copy of the cigar string
	# and remove hard clipping information
	# for convenience
	vdf$cigar2 <- gsub("\\d+H","",vdf$cigar)

	# Record the column indices for cigar columns
	cidx <- grep("cigar", colnames(vdf))[1]
	cidx2 <- grep("cigar", colnames(vdf))[2]

	# get the subset of rows with soft-clipping
	clipr <- grep("S", vdf$cigar)
	svdf <- vdf[clipr,]

	# calculate the total number of soft clipped bases
	# by summing up numbers that are adjacent to an S
	# in the cigar string
	svdf$softclip  <- apply(svdf,1,function(x){sum(as.numeric(strapply(as.character(x[cidx]),"(\\d+)S", simplify=c)))})

	# calculate the number of bases soft clipped from 
	# the left of the read by extracting the number
	# adjacent to any S at the beginning of the string
	svdf$leftclip  <- apply(svdf,1,function(x){sum(as.numeric(strapply(as.character(x[cidx2]),"^(\\d+)S",simplify=c)))})

	# calculate the number of bases soft clipped from 
	# the right of the read by extracting the number
	# adjacent to any S at the end of the string
	svdf$rightclip <- apply(svdf,1,function(x){sum(as.numeric(strapply(as.character(x[cidx2]),"(\\d+)S$",simplify=c)))})

	# extract the sequence in the middle of leftclip and rightclip
	# this is the sequence that has actually aligned
	svdf$clipseq   <- substr(svdf$seq, as.numeric(svdf$leftclip)+1, nchar(as.vector(svdf$seq))-as.numeric(svdf$rightclip))

	# substitute in the values
	vdf$softclip[clipr] <- svdf$softclip
	vdf$leftclip[clipr] <- svdf$leftclip
	vdf$rightclip[clipr] <- svdf$rightclip
	vdf$clipseq[clipr] <- svdf$clipseq

	# calculate the length of the clipped sequence
	vdf$cliplen <- nchar(vdf$clipseq)

	# remove the second cigar string
	# we calculated above
	vdf <- vdf[,-cidx2]

	# return the data
	return(vdf)

}
dna.complement <- function(x=NULL, reverse=TRUE) {

	# convert the sequence to upper case
	x <- toupper(x)

	# split sequence into a vector
	vec <- strsplit(x,split="")[[1]]

	# prepare a vector to hold the
	# revcom data
	cvec <- vector(mode="character",length(vec))

	# for each base....
	for (i in 1:length(vec)) {

		# complement the relevant bases
		if (vec[i] == "A") {
			cvec[i] <- "T"
		}
		if (vec[i] == "G") {
			cvec[i] <- "C"
		}
		if (vec[i] == "C") {
			cvec[i] <- "G"
		}
		if (vec[i] == "T") {
			cvec[i] <- "A"
		}
		if (vec[i] == "U") {
			cvec[i] <- "A"
		}
		if (vec[i] == "N") {
			cvec[i] <- "N"
		}
	}

	# reverse if requested (default)
	if (reverse == TRUE) {
		cvec <- rev(cvec)
	}

	# return the sequence as a string
	return(paste(cvec, sep="", collapse=""))
}
make.pwm <- function(vdf=NULL, minlen=1, maxlen=37, scaled=TRUE, strand="pos", revcom=FALSE, ttou=FALSE) {


	# set the defaults strand and change it
	# if we want to look at the negative strand
	strnum <- "+"
	if (strand == "neg") {
		strnum <- "-"
		
		# negative strand reads need to be reverse
		# complemented as BAM always reports alignments
		# to the positive strand
		revcom <- TRUE
	} 

	# filter the data.frame according to the input
	myvdf <- vdf[vdf$cliplen>=minlen&vdf$cliplen<=maxlen&vdf$strand==strnum,]


	# create a list to hold the output
	lst <- list()
	lst[["A"]] <- vector(mode="numeric", length=maxlen)
	lst[["C"]] <- vector(mode="numeric", length=maxlen)
	lst[["G"]] <- vector(mode="numeric", length=maxlen)
	lst[["T"]] <- vector(mode="numeric", length=maxlen)

	# iterate through the aligned sequences
	# in the filtered input
	for (s in as.vector(myvdf$clipseq)) {

		# reverse complement the read if needed		
		if (revcom == TRUE) {
			s <- dna.complement(s)
		}

		# convert T to U if the user requests it
		if (ttou == TRUE) {
			s <- gsub("T","U", toupper(s))
		}

		# split the sequence into a vector
		# iterate over each base and record counts
		# in the appropriate position using the idx
		# counter
		vec <- strsplit(s,split="")[[1]]
		idx <- 1
		for (v in vec) {
			lst[[v]][idx]=lst[[v]][idx]+1
			idx <- idx+1
		}
	}

	# create a matrix and convert
	# the counts from the list
	out <- matrix(nrow=4, ncol=maxlen)
	out[1,] <- lst[["A"]] 
	out[2,] <- lst[["C"]] 
	out[3,] <- lst[["G"]] 
	out[4,] <- lst[["T"]] 

	# provide appropriate row and column names
	rownames(out) <- c("A","C","G","T")
	colnames(out) <- 1:maxlen

	# vector to hold indices of columns with
	# zero counts
	skipcol <- vector(mode="numeric")

	# if we want to scale the data
	# into proportions (we generally do)
	if (scaled == TRUE) {

		# iterate over the columns
		for(j in 1:ncol(out)) {

			# if the col sums > 0 divide by the total
			# otherwise add it to the columns to be skipped
			if (sum( out[,j] ) > 0) {
				out[,j] <- out[,j] / sum( out[,j] )
			} else {
				skipcol <- append(skipcol, j)
			}
		}
	}

	# remove zero sum columns
	if (length(skipcol) > 0) {
		out <- out[,-skipcol]
	}
	
	# return the data
	return(out)

}
make.simple.consensus <- function(vdf=NULL, reflen=12000) {

	# create a vector for the consensus
	cons <- vector(mode="character", length=reflen)

	# create a list to count bases
	# at each position
	cl <- vector("list", reflen)

	# initialise the list with zeros
	for (i in 1:reflen) {
		cl[[i]]["A"]<-0
		cl[[i]]["G"]<-0
		cl[[i]]["C"]<-0
		cl[[i]]["T"]<-0
		cl[[i]]["N"]<-0
	}

	# iterate through the input data.frame
	for (i in 1:nrow(vdf)) {

		# store the row in r for
		# convenience
		r <- vdf[i,]

		# record the 5p position for convenience
		cpos <- r$pos

		# iterate over each base in the clipped sequence
		for (base in strsplit(toupper(r$clipseq), split="")[[1]]) {

			# record the base at position cpos
			# and move to the next position
			cl[[cpos]][base] <- cl[[cpos]][base] + 1
			cpos <- cpos + 1
		}
	}

	# iterate over the psoitions in the reference
	for (i in 1:reflen) {

		# if there are no bases, record an N
		if (sum(cl[[i]]) == 0) {
			cons[i] <- "N"
		} else {
			# otherwise sort and record the base that occurs
			# the most
			# NB does not deal with ties
			cons[i] <- names(cl[[i]])[order(cl[[i]])[5]]
		}
	}

	# return the consensus as a string	
	fas <- paste(cons, sep="", collapse="")
	return(fas)
}
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
pwm.heatmap <- function(pwm=NULL, col.fun=colorRampPalette(c("black","red"), space="rgb"), mar=c(3,2,1,1), cex.axis=0.8) {

	# set the margins
	par(mar=mar)

	# draw the image
	image(z=t(pwm), y=1:nrow(pwm), x=1:ncol(pwm), xaxt="n", yaxt="n", col=col.fun(10), xlab="",ylab="")

	# add axes
	axis(side=2, at=1:nrow(pwm),  labels=rownames(pwm), cex.axis=cex.axis)

	# add axes
	axis(side=1, at=1:ncol(pwm), labels=colnames(pwm), cex.axis=cex.axis)
}read.bam <- function (bamfile = NULL, chr = NULL, start = 1, end = 1e+07, 
    what = c("qname", "flag", "rname", "strand", "pos", "qwidth", 
        "mapq", "cigar", "mrnm", "mpos", "isize", "seq"), tag = c("NM"), 
    removeN = TRUE) 
{


    # read the data using scanBam from Rsamtools() package
    if (is.null(chr)) {

	# if no chromosome if given, read everything
	param <- ScanBamParam(what = what, tag = tag)
    	bam <- scanBam(bamfile, param = param)

    } else {
	# otherwise limit the data by chromosome and the
	# given start and end
    	which <- RangesList(chr = IRanges(start, end))
    	names(which) <- chr
    	param <- ScanBamParam(which = which, what = what, tag = tag)
    	bam <- scanBam(bamfile, param = param)
    }

    # the result is in a rather difficult
    # to use S4 class/list and so we go through 
    # the following steps to coerce the object
    # to a native data.frame

    # get the names of the object (the columns/"what")
    # and apply lapply
    elts <- names(bam[[1]])
    lst <- lapply(elts, function(elt) .unlist(lapply(bam, "[[", elt)))
    names(lst) <- elts

    # coerce to DataFrame using the in-built Rsamtools function
    vdf <- do.call("DataFrame", lst)
    colnames(vdf) <- names(lst)

    # remove blanks
    vdf <- vdf[!is.na(vdf$qwidth), ]

    # convert important columns to character
    vdf$seq <- as.character(vdf$seq)
    vdf$cigar <- as.character(vdf$cigar)

    # if the user wants to remove sequences
    # that contain N's, remove them
    if (removeN == TRUE) {
        ns <- grep("N", vdf$seq)
        if (length(ns) > 0) {
            vdf <- vdf[-ns, ]
        }
    }

    # call the native as.data,frame function 
    # to create a native data.frame and return it
    vdf <- as.data.frame(vdf, stringsAsFactors = FALSE)
    return(vdf)
}



.unlist <- function (x) {
	# code provided by Martin Morgan (author of Rsamtools)
	# to assist in coercion of Rsamtools classes to data.frames
	## do.call(c, ...) coerces factor to integer, which is undesired
	x1 <- x[[1L]]
	if (is.factor(x1)) {
		structure(unlist(x), class = "factor", levels = levels(x1))
	} else {
		do.call(c, x)
	}
}read.dist.plot <- function (sr = NULL, minlen = 1, maxlen = 37, method = "add", pad = 30, primary = "pos", plot=TRUE, title="5' read distance plot", xlab="Distance", ylab="Count") {

    # check the input paramters
    # check method argument is valid
    if (! method %in% c("pos", "neg", "mean", "multiply", "add", "min", "max")) {
        print("Method should be one of:")
        print(c("pos", "neg", "mean", "multiply", "add", "min", "max"))
	print(paste("You gave:", method))
        print("Setting method to mean")
        method = "mean"
    }

    # check primary is one of "pos" or "neg"
    if (! primary %in% c("pos", "neg")) {
        print("Primary should be one of:")
        print(c("pos", "neg"))
        print("Setting Primary to pos")
        pri = "pos"
    }

    # subset the input according to the
    # given size range
    psr <- sr[sr$len >= minlen & sr$len <= maxlen, ]

    # calculate the end position of each read
    psr$posend <- psr$pos + (nchar(psr$seq) - 1)

    # calculate the 5 prime and 3 prime ends
    # of each read
    psr$p5 <- psr$pos
    psr$p3 <- psr$pos
    psr$p5[psr$strand == "+"] <- psr$pos[psr$strand == "+"]
    psr$p3[psr$strand == "+"] <- psr$posend[psr$strand == "+"]
    psr$p3[psr$strand == "-"] <- psr$pos[psr$strand == "-"]
    psr$p5[psr$strand == "-"] <- psr$posend[psr$strand == "-"]

    # set the "primary" strand
    # it is this strand that we will iterate over
    # and then we will count overlaps on the opposite strand
    if (primary == "pos") {
        pri <- psr[psr$strand == "+", ]
        alt <- psr[psr$strand == "-", ]
    }   else {
        pri <- psr[psr$strand == "+", ]
        alt <- psr[psr$strand == "-", ]
    }

    # create an empty list for the output
    lst <- list()

    # for each entry in the sequence report
    # for the primary strand
    for (i in 1:nrow(pri)) {

	# r contains the row as a vector
        r <- pri[i, ]

        # set the beginning of the range in which we will
	# count reads on the opposite strand.  Set to 1
	# if it is <= 0
        begin <- r$p5 - pad
        if (begin <= 0) 
            begin <- 1

	# set the ending of the range
        ending <- r$p5 + pad

	# subset the alternative strand
	# according to the defined range
        sn <- alt[alt$p5 >= begin & alt$p5 <= ending, ]

	# proceed only if there are actually reads....
        if (nrow(sn) > 0) {

	    # cum up the counts for reads, by position, within the range
            ssn <- aggregate(sn$count, by = list(pos = sn$p5), sum)

	    # give the result some useful column names
            colnames(ssn) <- c("p5", "count")

	    # for every rown in the result
            for (j in 1:nrow(ssn)) {

		# e contains the 5p end of the data
                e <- ssn[j, ]$p5

		# cnt is the summed counts
                cnt <- ssn[j, ]$count

		# idx is the distance between the query read
		# and the read on the opposite strand
                idx <- as.character(e - r$p5)

		# initialise the output list
		# if no value has been set yet
                if (is.null(lst[[idx]])) 
                  lst[[idx]] <- 0

		# depending on which method has been chosen
		# record the result appropriately
                switch(method, 
		   mean = {
		  # add the average of counts on the primary and alternative strand
                  lst[[idx]] <- lst[[idx]] + (mean(c(cnt, r$count)))
                }, pos = {
		  # add just the count of reads on the primary strand
                  lst[[idx]] <- lst[[idx]] + r$count
                }, neg = {
		  # add just the count of reads on the alternate strand
                  lst[[idx]] <- lst[[idx]] + cnt
                }, multiply = {
		  # add the product of counts on the primary and alternate strands
                  lst[[idx]] <- lst[[idx]] + (cnt * r$count)
                }, add = {
		  # add the sum of counts on the primary and alternate strands
                  lst[[idx]] <- lst[[idx]] + sum(c(cnt,r$count))
                }, min = {
		  # add the min of counts on the primary and alternate strands
                  lst[[idx]] <- lst[[idx]] + min(c(cnt,r$count))
		}, max = {
		  # add the max of counts on the primary and alternate strands
                  lst[[idx]] <- lst[[idx]] + max(c(cnt,r$count))
		})
            }
        }
    }

    # unlist the result and create a data.frame output
    ul <- unlist(lst)
    du <- data.frame(loc = as.numeric(names(ul)), count = as.numeric(ul))

    # order by location
    du <- du[order(du$loc), ]

    # plot if the user wants it
    if (plot == TRUE) {

	plot(du$loc, du$count, type="l", main=title, xlab=xlab, ylab=ylab)

    }

    # return the result
    return(du)
}
sequence.report <- function(df=NULL, minlen=1, maxlen=37) {

	# filter the input according to the provided
	# size range
	mydf <- df[df$cliplen>=minlen&df$cliplen<=maxlen,]

	# use ddply to count the occurence of reads
	# grouped by reference (rname), position (pos), 
	# the actual sequence (clipseq), the sequencelength (cliplen)
	# and strand (strand)
	cts <- ddply(mydf, c("rname","pos","clipseq","cliplen","strand"), nrow)

	# order the result by position, strand, length and count
	cts <- cts[order(cts$pos, cts$strand, cts$cliplen, cts$V1),]

	# reorder the columns and rename them with sensible names
	cts <- cts[,c("rname","pos","strand","clipseq","cliplen","V1")]
	colnames(cts) <- c("ref","pos","strand","seq","len","count")

	# return
	return(cts)
}
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

}summarise.by.length <- function(vdf=NULL, minlen=1, maxlen=37, start=1, end=1e+07, strand=NULL) {

	# sort out end which has been given
	# a silly high default number so as not to miss
	# anything
	if (end > max(vdf$pos)) {
		end <- max(vdf$pos)
	}

	# filter according to input
	vdf <- vdf[vdf$cliplen>=minlen&vdf$cliplen<=maxlen&vdf$pos>=start&vdf$pos<=end,]

	# summarise data for each length of read
	# count the occurrence at each position
	if (is.null(strand)) {
		bypos <- acast(vdf, pos~cliplen, length)
	} else if (strand=="pos") {
		bypos <- acast(vdf[vdf$strand=="+",], pos~cliplen, length)

	} else if (strand=="neg") {
		bypos <- acast(vdf[vdf$strand=="-",], pos~cliplen, length)
	} else {
		print("WARNING: don't understand your strand argument")
		print("WARNING: sending back data for both strands")
		bypos <- acast(vdf, pos~cliplen, length)
	}

	# in case of blanks, create dummy data.frame
	dummy <- data.frame(pos=start:end)

	# and merge it with the results from cast(...)
	dm <- merge(dummy, bypos, by.x="pos", by.y="row.names", all.x=TRUE, sort=FALSE)

	# deal with NAs which are zero counts
	dm[is.na(dm)] <- 0

	# order by position
	dm <- dm[order(dm$pos),]

	# give the result sensible rownames
	rownames(dm) <- dm$pos

	# remove superfluous column
	dm <- dm[,-1]

	# return data
	return(dm)

}