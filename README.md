# viRome_legacy
Legacy R code for the R package viRome DOI:10.1093/bioinformatics/btt297

Install packages:

```sh
install.packages(c("BiocManager","Rcpp","seqinr","plyr","gsubfn","Rsamtools","reshape2","seqLogo", "motifStack", "S4Vectors"))

BiocManager::install("seqLogo")
BiocManager::install("motifStack")
```

Load the libraries

```sh
library("seqinr")
library("plyr")
library("gsubfn")
library("Rsamtools")
library("reshape2")
library("seqLogo")
library("motifStack")
library("S4Vectors")
library("Rcpp")
```

Load the code
```sh
source("https://raw.githubusercontent.com/mw55309/viRome_legacy/main/R/viRome_functions.R")
```

There is a BAM file and associated .bai file called "SRR389184_vs_SINV_sorted.bam" in the data directory of this repo. Download them as examples

Run viRome

```sh
infile <- "SRR389184_vs_SINV_sorted.bam"

bam <- read.bam(infile, chr="SINV")
# requires only the output of read.bam()
bamc <- clip.bam(bam)
# requires only the output of clip.bam()
bpl <- barplot.bam(bamc)
# requires only the output of barplot.bam()
ssp <- size.strand.bias.plot(bpl)
# requires only the output of clip.bam()
dm <- summarise.by.length(bamc)
sph <- size.position.heatmap(dm)
# requires only the output of summarise.by.length()
sbp <- stacked.barplot(dm)
# requires only the output of clip.bam()
# though one should alter minlen, maxlen
# and reflen
sir <- position.barplot(bamc)
# requires only the output of clip.bam()
sr <- sequence.report(bamc)
# requires only the output of clip.bam()
pwm <- make.pwm(bamc)
# requires only the output of make.pwm()
pmh <- pwm.heatmap(pwm)
# requires only the output of sequence.report()
rdp <- read.dist.plot(sr)
```
