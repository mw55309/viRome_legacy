summarise.by.length <- function(vdf=NULL, minlen=1, maxlen=37, start=1, end=1e+07, strand=NULL) {

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