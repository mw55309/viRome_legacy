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
