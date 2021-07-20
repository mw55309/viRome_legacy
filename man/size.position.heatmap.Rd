\name{size.position.heatmap}
\alias{size.position.heatmap}
\title{
Produce a heatmap of counts from the results of \code{\link{summarise.by.length}}
}
\description{
A heatmap will be drawn, the rows being read lengths, the columns being position in the genome, and the colours related to counts of read lengths mapping to those locations
}
\usage{
size.position.heatmap(dm = NULL, minlen = 1, maxlen = 37, start = 1, end = 1e+07, scale = TRUE, col.fun = colorRampPalette(c("black", "red"), space = "rgb"), log = FALSE, mar = c(5, 4, 4, 2), main = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dm}{
The output of \code{\link{summarise.by.length}}.  A data.frame, the columns will be read lengths, the rows will be genome position and the values counts
}
  \item{minlen}{
The minimum length of aligned read to consider
}
  \item{maxlen}{
The maximum length of aligned read to consider
}
  \item{start}{
The minimum start position to consider
}
  \item{end}{
The maximum start position to consider
}
  \item{scale}{
Whether to scale the data for each genome position
}
  \item{col.fun}{
A function that will provide colours; could be \code{\link{rainbow}} or \code{\link{heat.colors}} etc.  Defaults to \code{\link{colorRampPalette}} scaling from black to red
}
  \item{log}{
Whether to log() the data first (prior to scaling)
}
  \item{mar}{
The margins for the plot.  See \code{\link{par}}.  A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot.
}
  \item{main}{
A title for the plot
}
}
\value{
A heatmap is drawn on the current graphical device
}

\author{
Mick Watson <mick.watson@roslin.ed.ac.uk>
}
\examples{
 \dontrun{infile <- system.file("data/SRR389184_vs_SINV_sorted.bam", package="viRome")}
 \dontrun{bam <- read.bam(bamfile=infile, chr="SINV", start=1, end=11703, removeN=TRUE)}
 \dontrun{bamc <- clip.bam(bam)}
\dontrun{dm <- summarise.by.length(bamc)}
\dontrun{dmp <- summarise.by.length(bamc, strand="pos")}
\dontrun{dmn <- summarise.by.length(bamc, strand="neg")}

\dontrun{size.position.heatmap(dm)}
\dontrun{size.position.heatmap(dm, log=TRUE, scale=FALSE)}

}
\keyword{ aplot }

