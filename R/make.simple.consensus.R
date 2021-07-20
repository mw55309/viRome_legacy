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
