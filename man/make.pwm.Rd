\name{make.pwm}
\alias{make.pwm}
\title{
Construct position-weight matrix for subsets of reads from the output of \code{\link{read.bam}} and \code{\link{clip.bam}}
}
\description{
This function takes a data frame as input, which should be identical to the output of \code{\link{read.bam}} and \code{\link{clip.bam}}. It filters according to the length of the mapped reads, and then countsthe occurrence of each base along the length of the read(s).  The data are then optionally scaled.
}
\usage{
make.pwm(vdf = NULL, minlen = 1, maxlen = 37, scaled = TRUE, strand = "pos", revcom = FALSE, ttou = FALSE)
}
\arguments{
  \item{vdf}{
Data frame, should be the output of the output of \code{\link{clip.bam}}
}
  \item{minlen}{
The minimum length of mapped read to include
}
  \item{maxlen}{
The maximum length of mapped read to include
}
  \item{scaled}{
Whether or not to scale each column of base counts to the total number of bases in that column.  Default: TRUE
}
  \item{strand}{
The strand to calculate the PWM on: either "pos" or "neg"
}
  \item{revcom}{
Whether or not to reverse-complement the sequence.  We recommend you leave this as default.  The default should work except those times when, for example, a negative strand virus has been published as a positive strand, and you have aligned your data to the positive strand.
}
  \item{ttou}{
Whether or not to convert T to U (Uracil).  Currently not recommended.  In all liklihood, you measured cDNA anyway, not RNA, and therefore you should report a T :)
}
}
\details{

}
\value{
A matrix with the c("A","G","C","T") as rows, the position as columns and the (scaled) counts as values
}
\references{

}
\author{
Mick Watson
}
\examples{
\dontrun{infile <- system.file("examples/SRR389184_vs_SINV_sorted.bam", package="viRome")}
 \dontrun{bam <- read.bam(bamfile=infile, chr="SINV", minlen=1, maxlen=11703, removeN=TRUE)}
 \dontrun{bamc <- clip.bam(bam)}
 \dontrun{make.pwm(bamc, minlen=25, maxlen=37)}

}
\keyword{ data }
\keyword{ manip }
