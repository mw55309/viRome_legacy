\name{pwm.heatmap}
\alias{pwm.heatmap}
\title{
Draw a heatmap of a position-weight matrix (from \code{\link{make.pwm}}
}
\description{
Draw a heatmap of a position-weight matrix (from \code{\link{make.pwm}}.  Uses the \code{\link{image}} function
}
\usage{
pwm.heatmap(pwm = NULL, col.fun = colorRampPalette(c("black", "red"), space = "rgb"), mar = c(3, 2, 1, 1), cex.axis = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pwm}{
Output from \code{\link{make.pwm}}.  Rows are bases/nucleotides, columns are position and values are proportions
}
  \item{col.fun}{
A colour function such as \code{\link{rainbow}}.  Defaults to colorRampPalette(c("black", "red"), space = "rgb")
}
  \item{mar}{
The margins for the plot.  See \code{\link{par}}.  A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot.
}
  \item{cex.axis}{
A scaling factor for the axes labels
}
}
\value{
An image is plotted to the current device
}
\author{
Mick Watson <mick.watson@roslin.ed.ac.uk>
}
\examples{
\dontrun{infile <- system.file("examples/SRR389184_vs_SINV_sorted.bam", package="viRome")}
 \dontrun{bam <- read.bam(bamfile=infile, chr="SINV", minlen=1, maxlen=11703, removeN=TRUE)}
 \dontrun{bamc <- clip.bam(bam)}
 \dontrun{pwm1 <- make.pwm(bamc, minlen=25, maxlen=29, strand="neg")}
\dontrun{pwm.heatmap(pwm1)}
\dontrun{pwm.heatmap(pwm1, col.fun=colorRampPalette(c("green","red"), space="rgb"))}
}
\keyword{ aplot }

