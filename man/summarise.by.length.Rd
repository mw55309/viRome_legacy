\name{summarise.by.length}
\alias{summarise.by.length}
\title{
Summarise reads mapped by length of the read and position in the genome
}
\description{
This function takes the results of \code{\link{read.bam}} and \code{\link{clip.bam}} and counts the occurrence of the length of reads mapping to each position in the reference genome
}
\usage{
summarise.by.length(vdf = NULL, minlen = 1, maxlen = 37, start = 1, end = 1e+07, strand = NULL)
}
\arguments{
  \item{vdf}{
the results of \code{\link{read.bam}} and \code{\link{clip.bam}} 
}
  \item{minlen}{
The minimum length of aligned read to consider
}
  \item{maxlen}{
The maximum length of aligned read to consider
}
  \item{start}{
The start position from which alignments should be extracted
}
  \item{end}{
The end position from which alignments should be extracted
}
  \item{strand}{
Choose which strand to count reads on.  One of either "pos" or "neg".  For both strands, set to NULL (the default).  
}
}
\details{

}
\value{
The result is a data.frame, the columns being read lengths, the rows being position in the genome, and the values being counts.
}
\author{
Mick Watson <mick.watson@roslin.ed.ac.uk>
}
\note{
}

\examples{
 \dontrun{infile <- system.file("data/SRR389184_vs_SINV_sorted.bam", package="viRome")}
 \dontrun{bam <- read.bam(bamfile=infile, chr="SINV", start=1, end=11703, removeN=TRUE)}
 \dontrun{bamc <- clip.bam(bam)}
\dontrun{dm <- summarise.by.length(bamc)}
\dontrun{dmp <- summarise.by.length(bamc, strand="pos")}
\dontrun{dmn <- summarise.by.length(bamc, strand="neg")}
}
\keyword{ manip }
