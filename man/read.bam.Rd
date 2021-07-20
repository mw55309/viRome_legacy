\name{read.bam}
\alias{read.bam}

\title{
Reads a BAM file using the Rsamtools package and manipulates the data into a more useable form
}
\description{
The read.bam function reads all alignments between two points from a particular chromosome/reference sequence.  A variety of manipulations are performed, including conversion of the read data into a \code{\link{data.frame}}.
}
\usage{
read.bam(bamfile = NULL, chr = NULL, start = 1, end = 1e+07, what = c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar", "mrnm", "mpos", "isize", "seq"), tag = c("NM"), removeN = TRUE)
}
\arguments{
  \item{bamfile}{
Full path to the indexed BAM file
}
  \item{chr}{
The name of the chromosome from which alignments should be extracted
}
  \item{start}{
The start position from which alignments should be extracted
}
  \item{end}{
The end position from which alignments should be extracted
}
  \item{what}{
A vector of column names to be extracted.  See \code{\link{scanBamWhat}}.  At present allowed values are any subset of c("qname","flag","rname","strand","pos","qwidth","mapq","cigar","mrnm","mpos","isize","seq","qual").  We recommend sticking with the default.
}
  \item{tag}{
Tags to be exctracted from the BAM file.  See \code{\link{scanBam}}.  Many different tags are used by the various aligners that exist, and we recommend referring to the manual of your chosen aligner to discover which tags are used, and which you may find useful.
}
  \item{removeN}{
Logical.  Whether or not to remove all sequences that contain Ns
}
}
\details{
A wrapper around \code{\link{scanBam}} to read in a subset of aliignments from a BAM file and return those as a \code{\link{data.frame}}.
}
\value{
A data frame with reads/alignments as rows and columns defined by the "what" and "tag" columns.
}
\references{

}
\author{
Mick Watson
}
\note{

}
\seealso{
\code{\link{scanBam}}, \code{\link{clip.bam}}
}
\examples{

\dontrun{infile <- system.file("data/SRR389184_vs_SINV_sorted.bam", package="viRome")}
\dontrun{bam <- read.bam(bamfile=infile, chr="SINV", start=1, end=11703, removeN=TRUE)}
\dontrun{bam[1:10,]}
}
\keyword{ data }
\keyword{ manip }
\keyword{ IO }

