\name{make.simple.consensus}
\alias{make.simple.consensus}
\title{
Make Simple Consensus
}
\description{
From aligned data read in by \code{\link{read.bam}} and processed by \code{\link{clip.bam}} create a consensus sequence
}
\usage{
make.simple.consensus(vdf = NULL, reflen = 11703)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vdf}{
A data frame as produced by \code{\link{read.bam}} and processed by \code{\link{clip.bam}}.  The fields the data.frame absolutely must have are "pos" (left-most position of the alignment on the positive strand) and clipseq (sequence of the aligned read on the positive strand
}
  \item{reflen}{
The reference length
}
}
\details{
This is a *very simple* way to calculate the consensus sequence, and in general we recommend other tools, such as samtools, for doing this.

Basically, for each location in the reference, this function counts the most common aligned base.  These are then returned as a string.  By definition the return is the same length as the reference.  Indels not handled.
}
\value{
A string of the consensus
}
\references{
}
\author{
Mick Watson
}
\note{
}
\seealso{
}
\examples{
 \dontrun{infile <- system.file("data/SRR389184_vs_SINV_sorted.bam", package="viRome")}
 \dontrun{bam <- read.bam(bamfile=infile, chr="SINV", start=1, end=12000, removeN=TRUE)}
 \dontrun{bamc <- clip.bam(bam)}
 \dontrun{seq <- make.simple.consensus(bamc)}
}
\keyword{ data }
\keyword{ manip }
