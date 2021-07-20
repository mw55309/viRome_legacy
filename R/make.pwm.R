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
