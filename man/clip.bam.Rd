\name{clip.bam}
\alias{clip.bam}
\title{
Deals with soft-clipping encoded in the CIGAR string and calculates statistics on the clipped sequences
}
\description{
The function investigates the CIGAR string and performs  soft-clipping of the sequences within the "seq" column, storing the results (and other various statistics) within additional columns.
}
\usage{
clip.bam(vdf = NULL)
}
\arguments{
  \item{vdf}{
A data frame as read from \code{\link{read.bam}}
}
}
\details{
Aligned sequences are stored in BAM files with the caveat that they may have been soft-clipped - i.e. the clipped sequence was not used in the alignment, but is included in the BAM file.  This function looks at the CIGAR string and clips data in the "seq" column appropriately.

For example, conside the CIGAR string "12S25M" and the sequence read "CACCCGAGAATACCCCAGAACCATTATGCTGTGACTT".  The sequence read is 37bp long; however, the cigar string tells us that only 25 "matched" (25M) i.e. only 25 were used in the alignment.  The CIGAR string also tells us that 12bp were soft-clipped (12S).  We can tell that they were soft-clipped from the start of the read as 12S occurs before 25M.  

This function first removes from the CIGAR string any notation of "hard-clipping" - we do not need this information.  The function then finds all reads that are marked as "soft-clipped", and examines whether they are marked as clipped from the beginning or end of the read.  The total amount of clipped sequence, the amount clipped from the left- and right- sides, the clipped sequence and the length of the clipped sequence are all then stored in additional columns.
}
\value{
The same data.frame you put it with the additional fields:
  \item{softclip }{Total number of bases soft clipped from the read}
  \item{leftclip }{Number of bases clipped from the left of the read}
  \item{rightclip }{Number of bases clipped from the right of the read}
  \item{clipseq }{The sequence left after clipping}
  \item{cliplen }{The length of the clipped sequence}
}
\references{

}
\author{
Mick Watson
}
\note{

}
\seealso{
\code{\link{read.bam}}
}
\examples{
 \dontrun{infile <- system.file("data/SRR389184_vs_SINV_sorted.bam", package="viRome")}
 \dontrun{bam <- read.bam(bamfile=infile, chr="SINV", start=1, end=11703, removeN=TRUE)}
 \dontrun{bam[1:10,]}
 \dontrun{bamc <- clip.bam(bam)}
 \dontrun{bamc[1:10,]}
}
\keyword{ data }
\keyword{ manip }
