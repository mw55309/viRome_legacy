read.bam <- function (bamfile = NULL, chr = NULL, start = 1, end = 1e+07, 
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
}