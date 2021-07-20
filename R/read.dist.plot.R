read.dist.plot <- function (sr = NULL, minlen = 1, maxlen = 37, method = "add", pad = 30, primary = "pos", plot=TRUE, title="5' read distance plot", xlab="Distance", ylab="Count") {

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
