pwm.heatmap <- function(pwm=NULL, col.fun=colorRampPalette(c("black","red"), space="rgb"), mar=c(3,2,1,1), cex.axis=0.8) {

	# set the margins
	par(mar=mar)

	# draw the image
	image(z=t(pwm), y=1:nrow(pwm), x=1:ncol(pwm), xaxt="n", yaxt="n", col=col.fun(10), xlab="",ylab="")

	# add axes
	axis(side=2, at=1:nrow(pwm),  labels=rownames(pwm), cex.axis=cex.axis)

	# add axes
	axis(side=1, at=1:ncol(pwm), labels=colnames(pwm), cex.axis=cex.axis)
}