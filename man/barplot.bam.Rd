\name{barplot.bam}
\alias{barplot.bam}
\title{
Draw/calculate a barplot showing the length distribution of reads aligned to the reference.
}
\description{
This function counts the occurrence of reads of a given size range mapping to the positive and negative strands of a reference sequence.  

Generally, users of viRome may be interested in how many 21-22bp (viRNA) and 24-29bp (piRNA) reads align to the reference genome, on both the positive and negative strand.  
}
\usage{
barplot.bam(vdf = NULL, minlen = 1, maxlen = 37, poscol="red", negcol="green", main = "Sequence length distribution", xlab = "Map length", ylab = "Count", legend = c("+ve strand", "-ve strand"),
legendx = NULL, legendy = NULL, plot=TRUE, down=FALSE, sym.axes=TRUE,  ...)
}
\arguments{
  \item{vdf}{
A \code{\link{data.frame}} as produced by \code{\link{read.bam}} and \code{\link{clip.bam}}
}
  \item{minlen}{
The minimum length of aligned read to consider
}
  \item{maxlen}{
The maximum length of aligned read to consider
}
  \item{poscol}{
The colour of the "positive" bars.  Defaults to "red"
}
 \item{negcol}{
The colour of the "negative" bars.  Defaults to "green"
}
  \item{main}{
Title for the plot
}
  \item{xlab}{
Labels for the x-axis
}
  \item{ylab}{
Labels for the y-axis
}
  \item{legend}{
Vector of text for the positive and negative strand in the legend
}  
  \item{legendx}{
The x location for the legend
}
  \item{legendy}{
The y location for the legend
}
  \item{plot}{
Whether or not to draw a plot.  Set this to "FALSE" to suppress plotting.
}
  \item{down}{
Whether or not the plot for the negative strand should be plotted downwards.  Defaults to FALSE.
}
  \item{sym.axes}{
If down is TRUE, this determines whether the x-axis is symmetrical around zero.  Defaults to TRUE
}
  \item{...}{
Further parameters to pass to the \code{\link{barplot}} function
}

}
\details{
A barplot is plotted, and the summarised data returned as a data.frame.  Set plot=FALSE to suppress plotting.
}
\value{
A barplot is plotted, and the data returned as a data.frame
  \item{lgth }{The length of the aligned read}
  \item{pos }{Count of alignments on the positive strand}
  \item{neg }{Count of alignments on the negative strand}
}
\references{

}
\author{
Mick Watson
}
\note{

}
\seealso{
\code{\link{read.bam}}, \code{\link{clip.bam}} 
}
\examples{

 \dontrun{infile <- system.file("data/SRR389184_vs_SINV_sorted.bam", package="viRome")}
 \dontrun{bam <- read.bam(bamfile=infile, chr="SINV", start=1, end=11703, removeN=TRUE)}
 \dontrun{bamc <- clip.bam(bam)}

 \dontrun{bi <- barplot.bam(vdf=bamc, minlen=1, maxlen=37)}
 \dontrun{bi}

 \dontrun{bp <- barplot.bam(bamc)}
 \dontrun{bp <- barplot.bam(bamc, down=TRUE)}


 \dontrun{bp <- barplot.bam(bamc, space=c(0,0))}
 \dontrun{bp <- barplot.bam(bamc, down=TRUE, space=c(0,0))}


 \dontrun{bp <- barplot.bam(bamc, down=TRUE, space=c(0,0), minlen=17, maxlen=36, legendx=17)}
}
\keyword{ aplot }
\keyword{ manip }
