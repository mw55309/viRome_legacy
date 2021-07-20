\name{read.dist.plot}
\alias{read.dist.plot}
\title{
Plot distances between 5' ends of reads on complementary strands
}
\description{
Plot a graph and produce a data.frame, ready for plotting, that counts the distances between the 5' end of reads on opposite strands
}
\usage{
read.dist.plot(sr = NULL, minlen = 1, maxlen = 37, method = "add", pad = 30, primary = "pos", plot=TRUE, title="5' read distance plot", xlab="Distance", ylab="Count")
}
\arguments{
  \item{sr}{
Sequence report as output by \code{\link{sequence.report}}
}
  \item{start}{
The minimum length of mapped read to consider
}
  \item{end}{
The maximum length of mapped read to consider
}
  \item{method}{
Must be one of c("pos", "neg", "mean", "multiply", "add", "min", "max"); mean (takes the mean of counts on opposite strands), pos (counts reads on the pos strand only), neg (counts reads on the neg strand only), multiply (multiplies count of pos with count of neg), add (adds counts of pos and neg), min (minimum of pos and neg counts) and max (maximum of pos and neg counts)
}
  \item{pad}{
Distance to pad either side of the 5' end to look for overlaps
}
  \item{primary}{
Whether the positive or negative strand is iterated over first.  Acceptes "pos" or "neg"
}
  \item{plot}{
Whether or not to plot to the current device.  Default TRUE
}
  \item{title}{
Title for the plot.  Ignored is plot==FALSE
}
 \item{xlab}{
X axis labels for the plot.  Ignored is plot==FALSE
}
  \item{ylab}{
Y axis labels for the plot.  Ignored is plot==FALSE
}
}
\details{
This function can be used to look for enrichment of particular distances between reads on opposite strands: for example, it has been shown in quite a few publications that the piRNA response (25-29bp) is characterised by reads on complementary strands whose 5' ends overlap by 10bp.

This function is more generic, and just summarises distances between the 5' ends of reads on opposite strands.

The main issue is which summary function one chooses.  For example, let us say that we have 10 reads on the positive strand whose 5' end map at position 50; we also have 20 reads on the negative strand whose 5' ends map at position 59.  The distance between the 5' ends is +10: but how many reads do we count as having that distance?  Is it 10 (pos)?  20 (neg)?   min(10,20)?  max(10,20)?  mean(10,20)?  10 * 20 (multiply)?  You may decide, though we set the default as "add".
}
\value{
Data frame with values:
  \item{loc  }{Distance between 5' ends }
  \item{count }{Cumulative function (see "method") of reads that are "loc" distance apart at the 5' end}
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
 \dontrun{sr <- sequence.report(bamc)}
 \dontrun{sr}
 \dontrun{rd <- read.dist.plot(sr, minlen=25, maxlen=29)}
 }
\keyword{ data }
\keyword{ manip }
