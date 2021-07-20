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
