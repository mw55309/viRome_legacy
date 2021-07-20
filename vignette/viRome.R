### R code from vignette source 'viRome.Rnw'

###################################################
### code chunk number 1: viRome.Rnw:32-45
###################################################
# get the necessary Bioconductor packages
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("Rsamtools")

# install optional Bioconductor packages
# biocLite("seqLogo")
# biocLite("motifStack")

# install the necessary R packages
# install.packages(c("plyr","gsubfn","seqinr","reshape2"))

# install optional R packages
# install.packages("ggplot2")


###################################################
### code chunk number 2: viRome.Rnw:75-87
###################################################
?read.bam
?clip.bam
?barplot.bam
?size.strand.bias.plot
?summarise.by.length
?size.position.heatmap
?stacked.barplot
?position.barplot
?sequence.report
?make.pwm
?pwm.heatmap
?read.dist.plot


###################################################
### code chunk number 3: viRome.Rnw:92-137
###################################################
# load the library
# find example data
library(viRome)
infile <- system.file("data/SRR389184_vs_SINV_sorted.bam", package="viRome")

# minimal commands

# requires the full path to a bam file, 
# and the name of the reference the data are aligned to
bam    <- read.bam(infile, chr="SINV") 

# requires only the output of read.bam()
bamc   <- clip.bam(bam)  

# requires only the output of clip.bam()              
bpl    <- barplot.bam(bamc)            

# requires only the output of barplot.bam()
ssp    <- size.strand.bias.plot(bpl)   

# requires only the output of clip.bam()
dm     <- summarise.by.length(bamc)    

# requires only the output of summarise.by.length()
sph    <- size.position.heatmap(dm)    

# requires only the output of summarise.by.length()
sbp    <- stacked.barplot(dm)          

# requires only the output of clip.bam()
# though one should alter minlen, maxlen 
# and reflen
sir    <- position.barplot(bamc)       
              
# requires only the output of clip.bam()                         
sr     <- sequence.report(bamc)        

# requires only the output of clip.bam()
pwm    <- make.pwm(bamc)               

# requires only the output of make.pwm()
pmh    <- pwm.heatmap(pwm)             

# requires only the output of sequence.report()
rdp    <- read.dist.plot(sr)           


###################################################
### code chunk number 4: viRome.Rnw:146-153
###################################################
library(viRome)
infile <- system.file("data/SRR389184_vs_SINV_sorted.bam", package="viRome")
bam    <- read.bam(bamfile = infile, chr = "SINV", start = 1, end = 11703,
                   what = c("qname", "flag", "rname", "strand", 
			    "pos", "qwidth", "mapq", "cigar", "mrnm", 
                            "mpos", "isize", "seq"),
                   tag = c("NM"), removeN = TRUE)


###################################################
### code chunk number 5: viRome.Rnw:158-159
###################################################
?read.bam


###################################################
### code chunk number 6: viRome.Rnw:176-177
###################################################
bamc <- clip.bam(bam)


###################################################
### code chunk number 7: viRome.Rnw:182-183
###################################################
?clip.bam


###################################################
### code chunk number 8: viRome.Rnw:198-204
###################################################
b <- barplot.bam(vdf = bamc, minlen = 1, maxlen = 37, 
                 poscol="red", negcol="green",
                 main = "Sequence length distribution", 
                 xlab = "Map length", ylab = "Count",
                 legend = c("+ve strand", "-ve strand"),
                 legendx = NULL, legendy = NULL)


###################################################
### code chunk number 9: fig1
###################################################
b <- barplot.bam(vdf = bamc, minlen = 1, maxlen = 37, poscol="red", negcol="green",
                 main = "Sequence length distribution", xlab = "Map length", ylab = "Count",
                 legend = c("+ve strand", "-ve strand"),
                 legendx = NULL, legendy = NULL)


###################################################
### code chunk number 10: viRome.Rnw:219-220
###################################################
b


###################################################
### code chunk number 11: viRome.Rnw:228-235
###################################################
b <- barplot.bam(vdf = bamc, minlen = 17, maxlen = 37, 
		 poscol="blue", negcol="yellow",
                 main = "Sequence length distribution", 
		 xlab = "Map length", ylab = "Count",
                 legend = c("+ve strand", "-ve strand"),
                 legendx = 17, legendy = NULL, down=TRUE, 
		 space=c(0,0))


###################################################
### code chunk number 12: fig3
###################################################
b <- barplot.bam(vdf = bamc, minlen = 17, maxlen = 37, 
		 poscol="blue", negcol="yellow",
                 main = "Sequence length distribution", 
		 xlab = "Map length", ylab = "Count",
                 legend = c("+ve strand", "-ve strand"),
                 legendx = 17, legendy = NULL, down=TRUE, 
		 space=c(0,0))


###################################################
### code chunk number 13: viRome.Rnw:264-265
###################################################
size.strand.bias.plot(b, mar=c(4,4,3,1), pch=".", tpos=3, cextxt=0.8)


###################################################
### code chunk number 14: fig4
###################################################
size.strand.bias.plot(b, mar=c(4,4,3,1), pch=".", tpos=3, cextxt=0.8)


###################################################
### code chunk number 15: viRome.Rnw:287-290
###################################################
dm <- summarise.by.length(bamc)
dmp <- summarise.by.length(bamc, strand="pos")
dmn <- summarise.by.length(bamc, strand="neg")


###################################################
### code chunk number 16: viRome.Rnw:297-298
###################################################
dm[1:10,2:19]


###################################################
### code chunk number 17: viRome.Rnw:305-306
###################################################
size.position.heatmap(dm, mar=c(3,3,2,1))


###################################################
### code chunk number 18: fig5
###################################################
size.position.heatmap(dm, mar=c(3,3,2,1))


###################################################
### code chunk number 19: viRome.Rnw:320-321
###################################################
size.position.heatmap(dm, scale=FALSE, log=TRUE, mar=c(3,3,2,1))


###################################################
### code chunk number 20: fig6
###################################################
size.position.heatmap(dm, scale=FALSE, log=TRUE, mar=c(3,3,2,1))


###################################################
### code chunk number 21: viRome.Rnw:338-346
###################################################
split.screen(c(3,1))
screen(1)
size.position.heatmap(dm, log=TRUE, scale=FALSE, mar=c(2,2,2,1), main="Both")
screen(2)
size.position.heatmap(dmp, log=TRUE, scale=FALSE, mar=c(2,2,2,1), main="Pos")
screen(3)
size.position.heatmap(dmn, log=TRUE, scale=FALSE, mar=c(2,2,2,1), main="Neg")
close.screen(all=TRUE)


###################################################
### code chunk number 22: fig7
###################################################
split.screen(c(3,1))
screen(1)
size.position.heatmap(dm, log=TRUE, scale=FALSE, mar=c(2,2,2,1), main="Both")
screen(2)
size.position.heatmap(dmp, log=TRUE, scale=FALSE, mar=c(2,2,2,1), main="Pos")
screen(3)
size.position.heatmap(dmn, log=TRUE, scale=FALSE, mar=c(2,2,2,1), main="Neg")
close.screen(all=TRUE)


###################################################
### code chunk number 23: viRome.Rnw:371-372
###################################################
stacked.barplot(dm, internal.margins=c(0,2,0,1), skip.x=4, main.adj=0.5)


###################################################
### code chunk number 24: fig8
###################################################
stacked.barplot(dm, internal.margins=c(0,2,0,1), skip.x=4, main.adj=0.5)


###################################################
### code chunk number 25: viRome.Rnw:386-390
###################################################
stacked.barplot(dm, internal.margins=c(0,2,0,1), skip.x=4, 
		main.adj=0.5, 
		col.fun=colorRampPalette(c("green", "red"), 
			space = "rgb"))


###################################################
### code chunk number 26: fig9
###################################################
stacked.barplot(dm, internal.margins=c(0,2,0,1), skip.x=4, 
		main.adj=0.5, 
		col.fun=colorRampPalette(c("green", "red"), 
			space = "rgb"))


###################################################
### code chunk number 27: viRome.Rnw:407-410
###################################################
stacked.barplot(dm, internal.margins=c(0,2,0,1), 
		skip.x=4, main.adj=0.5, 
		bicol=c("red","black"))


###################################################
### code chunk number 28: fig18
###################################################
stacked.barplot(dm, internal.margins=c(0,2,0,1), 
		skip.x=4, main.adj=0.5, 
		bicol=c("red","black"))


###################################################
### code chunk number 29: viRome.Rnw:424-427
###################################################
stacked.barplot(dm, internal.margins=c(0,2,0,1), 
		skip.x=4, main.adj=0.5, 
		samey=TRUE)


###################################################
### code chunk number 30: fig19
###################################################
stacked.barplot(dm, internal.margins=c(0,2,0,1), 
		skip.x=4, main.adj=0.5, 
		samey=TRUE)


###################################################
### code chunk number 31: viRome.Rnw:447-450
###################################################
sirna <- position.barplot(vdf = bamc, minlen = 21, 
			  maxlen = 22, reflen = 11703, 
			  samp = "SRR389184")


###################################################
### code chunk number 32: fig10
###################################################
sirna <- position.barplot(vdf = bamc, minlen = 21, 
			  maxlen = 22, reflen = 11703, 
			  samp = "SRR389184")


###################################################
### code chunk number 33: viRome.Rnw:468-469
###################################################
sirna[1:20,]


###################################################
### code chunk number 34: viRome.Rnw:476-479
###################################################
pirna <- position.barplot(vdf = bamc, minlen = 25, 
			  maxlen = 29, reflen = 11703, 
			  samp = "SRR389184")


###################################################
### code chunk number 35: fig11
###################################################
pirna <- position.barplot(vdf = bamc, minlen = 25, 
			  maxlen = 29, reflen = 11703, 
			  samp = "SRR389184")


###################################################
### code chunk number 36: viRome.Rnw:508-509
###################################################
sr <- sequence.report(bamc, minlen=1, maxlen=37)


###################################################
### code chunk number 37: viRome.Rnw:514-515
###################################################
sr[1:10,]


###################################################
### code chunk number 38: viRome.Rnw:528-530
###################################################
pwm1 <- make.pwm(bamc, minlen=25, maxlen=29, strand="neg")
pwm2 <- make.pwm(bamc, minlen=25, maxlen=29)


###################################################
### code chunk number 39: viRome.Rnw:537-538
###################################################
pwm1[,1:8]


###################################################
### code chunk number 40: viRome.Rnw:544-545
###################################################
pwm.heatmap(pwm1, col.fun=colorRampPalette(c("green","red"), space="rgb"))


###################################################
### code chunk number 41: fig12
###################################################
pwm.heatmap(pwm1, col.fun=colorRampPalette(c("green","red"), space="rgb"))


###################################################
### code chunk number 42: viRome.Rnw:555-556
###################################################
pwm.heatmap(pwm2, col.fun=colorRampPalette(c("green","red"), space="rgb"))


###################################################
### code chunk number 43: fig13
###################################################
pwm.heatmap(pwm2, col.fun=colorRampPalette(c("green","red"), space="rgb"))


###################################################
### code chunk number 44: viRome.Rnw:570-572
###################################################
library(seqLogo) # must have this installed!
seqLogo(pwm1)


###################################################
### code chunk number 45: fig14
###################################################
library(seqLogo) # must have this installed!
seqLogo(pwm1)


###################################################
### code chunk number 46: viRome.Rnw:583-585
###################################################
library(seqLogo) # must have this installed!
seqLogo(pwm2)


###################################################
### code chunk number 47: fig15
###################################################
library(seqLogo) # must have this installed!
seqLogo(pwm2)


###################################################
### code chunk number 48: viRome.Rnw:599-601
###################################################
barplot(pwm1, col=rainbow(4), legend.text=rownames(pwm1), 
	args.legend=list(x=39), xlim=c(0,45))


###################################################
### code chunk number 49: fig16
###################################################
barplot(pwm1, col=rainbow(4), legend.text=rownames(pwm1), 
	args.legend=list(x=39), xlim=c(0,45))


###################################################
### code chunk number 50: viRome.Rnw:616-617
###################################################
rdp <- read.dist.plot(sr, minlen=25, maxlen=29, method="add")


###################################################
### code chunk number 51: fig17
###################################################
rdp <- read.dist.plot(sr, minlen=25, maxlen=29, method="add")


###################################################
### code chunk number 52: viRome.Rnw:629-630
###################################################
rdp


