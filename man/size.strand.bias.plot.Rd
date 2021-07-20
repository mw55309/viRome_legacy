\name{size.strand.bias.plot}
\alias{size.strand.bias.plot}
\title{
Plot counts of read lengths on + and - strands
}
\description{
This simple function takes the output of \code{\link{barplot.bam}} and plots counts of read length on the + and - strand against one another.  It can be used to visually inspect the data to see if there is a strand bias for particular read lengths
}
\usage{
size.strand.bias.plot(bp = NULL, minlen = 1, maxlen = 37, line.col = "red", sym.axes = TRUE, title = "Strand bias plot", xlab = "+ strand counts", ylab = "- strand counts", lty = 1, lwd = 2, cextxt = 1, mar = c(5, 4, 4, 1), tpos = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bp}{
The output from \code{\link{barplot.bam}} 
}
   \item{minlen}{
The minimum length of aligned read to consider
}
  \item{maxlen}{
The maximum length of aligned read to consider
}
  \item{line.col}{
The colour of the y=x line to plot (default: "red")
}
  \item{sym.axes}{
Whether or not the X- and Y- axes should be symmetrical
}
  \item{title}{
A title for the plot
}
  \item{xlab}{
A label for the x-axis
}
  \item{ylab}{
A label for the y-axis.
}
  \item{lty}{
The line type. Line types can either be specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) or as one of the character strings "blank", "solid", "dashed", "dotted", "dotdash", "longdash", or "twodash", where "blank" uses ‘invisible lines’ (i.e., does not draw them). 

}
  \item{lwd}{
The line width, a positive number, defaulting to 1. The interpretation is device-specific, and some devices do not implement line widths less than one.
}
  \item{cextxt}{
A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.  Use this to adjust the size of the labels in the plot
}
  \item{mar}{
A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot.

}
  \item{tpos}{
Control the position of the text label relative to the point. Values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the specified coordinates.
}
  \item{\dots}{
Further graphical paramters to be passed to the plot() function
}
}
\details{
}
\value{
A plot is created on the current device.
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
 \dontrun{bp <- barplot.bam(bamc)}
 \dontrun{size.strand.bias.plot(bp)}

 \dontrun{size.strand.bias.plot(bp, mar=c(4,4,3,1), pch=".", tpos=3, cextxt=0.8)}
}

\keyword{ aplot }
\keyword{ manip }
