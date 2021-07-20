\name{sequence.report}
\alias{sequence.report}
\title{
Produce a summary of the entire BAM alignment, grouped by reference, position, sequence and strand
}
\description{
For the positive and negative strand, group the data by reference, position, sequence length and strand and count the occurrence, producing a data.frame.  Basically counts the occurrence of each unique read/mapping location.
}
\usage{
sequence.report(df = NULL, minlen = 1, maxlen = 37)
}
\arguments{
  \item{df}{
A data.frame as output by \code{\link{read.bam}} and \code{\link{clip.bam}}
}
  \item{minlen}{
The minimum length of reads/alignments to include
}
  \item{maxlen}{
The maximum length of reads/alignments to include
}
}
\details{
}
\value{
A data.frame with columns:
  \item{ref }{The name of the reference sequence}
  \item{pos }{Position in the reference}
  \item{strand }{The strand}
  \item{seq }{The sequence of the aligned read}
  \item{len }{The length of the alignment}
  \item{count }{The count of the occurrence of this sequence at this location/strand}
}
\references{
}
\author{
Mick Watson
}
\note{
}

\examples{
 \dontrun{infile <- system.file("data/SRR389184_vs_SINV_sorted.bam", package="viRome")}
 \dontrun{bam <- read.bam(bamfile=infile, chr="SINV", minlen=1, maxlen=12000, removeN=TRUE)}
 \dontrun{bamc <- clip.bam(bam)}
 \dontrun{sr <- sequence.report(bamc)}
 \dontrun{sr}
}
\keyword{ data }
\keyword{ manip }
