\name{dna.complement}
\alias{dna.complement}
\title{
Complement (and reverse) DNA sequence
}
\description{
Complement (and optionally reverse) a DNA or RNA sequence
}
\usage{
dna.complement(x = NULL, reverse = TRUE)
}
\arguments{
  \item{x}{
A single character string consisting of characters from [AGCTagctNn]
}
  \item{reverse}{
Logical - whether to reverse the result too.  Default: TRUE
}
}
\details{
This function will complement DNA and RNA sequence according to the rules of DNA complementarity, and optionally reverse the result
}
\value{
A string is returned
}
\references{
http://en.wikipedia.org/wiki/Complementarity_%28molecular_biology%29
}
\author{
Mick Watson (mick.watson@roslin.ed.ac.uk)
}
\note{
}
\examples{

	seq <- "AGCTAGCTAGCT"
	dna.complement(seq)
	dna.complement(seq, reverse=FALSE)
	
}
\keyword{ data }
\keyword{ manip }
