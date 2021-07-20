\name{stacked.barplot}
\alias{stacked.barplot}
\title{
Produce a stacked barplot from the output of \code{\link{summarise.by.length}}
}
\description{
The function produces a barplot for each read length showing the number of mapped reads of that read length at each position of the reference genome.  One barplot will be plotted for each read length.  The current graphical device will be used, and the function \code{\link{split.screen}} is used to produce the graphic
}
\usage{
stacked.barplot(dm = NULL, minlen = 1, maxlen = 37, start = 1, end = 1e+07, internal.margins = c(0, 0, 0, 1), skip.x = 2, bicol=NULL, col.fun = rainbow, axis.col = "black", main.col = "black", main.adj = 1, samey=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dm}{
The output of \code{\link{summarise.by.length}}
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
  \item{internal.margins}{
The margins for each plot.  See \code{\link{par}}.  A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot.
}
  \item{skip.x}{
Only plot an x-axis every skip.x plots
}
  \item{bicol}{
Set this equal to a vector of two colours e.g. bicol=c("red","blue") and alterate graphs will take on the alternate colours.  Unless this is NULL (the default), col.fun is ignored.
}
  \item{col.fun}{
A colour function to produce colours for the barplots.  Defaults to \code{\link{rainbow}}
}
  \item{axis.col}{
A colour for the x-axis
}
  \item{main.col}{
A colour for the title of each plot (which will be the read length)
}
  \item{main.adj}{
Whether to plot the title left-justified (main.adj=0), centre-justified (main.adj=0.5) or right-justified (main.adj=1)
}
  \item{samey}{
Whether or not all barplots should have the same y-axis.  If set to TRUE, all barplots will have ylim=c(min(dm),max(dm)).  If set to FALSE, each barplot will have its y-axis set individually
}
  \item{\dots}{
Further graphical paramters to be passed to the plot() function
}
}
\value{
A plot is drawn on the current graphical device
}
\author{
Mick Watson <mick.watson@roslin.ed.ac.uk>
}
\note{
When one draws multiple plots stacked one on top of another, one rapidly runs out of space vertically.  Warnings are issued to this effect.  If you run out of space, or receive errors to that effect, please use a plotting device such as \code{\link{png}}, \code{\link{jpeg}}, or \code{\link{pdf}}, and be sure to provide enough height to the device
}
\examples{
 \dontrun{infile <- system.file("data/SRR389184_vs_SINV_sorted.bam", package="viRome")}
 \dontrun{bam <- read.bam(bamfile=infile, chr="SINV", start=1, end=11703, removeN=TRUE)}
 \dontrun{bamc <- clip.bam(bam)}
 \dontrun{dm <- summarise.by.length(bamc)}

\dontrun{stacked.barplot(dm)}

\dontrun{stacked.barplot(dm, internal.margins=c(0,2,0,1))}

\dontrun{stacked.barplot(dm, internal.margins=c(0,2,0,1), skip.x=4, main.adj=0.5)}
	
}

\keyword{ aplot }

