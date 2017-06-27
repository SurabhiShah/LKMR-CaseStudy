###################################################################
# Cross-sectional graphs for varying levels of Z1, at low/high Z2 #
###################################################################

source("PostmeanFunctions.R")
doginv  = FALSE
ngrid 	= 30
qs 		= c(0.25, 0.5, 0.75, 0.05, 0.95)
par(mfrow=c(2,3))
ylim 	= c(-6, 6)

########################################
# Top panel #
########################################

j = 1 

###############
# Timepoint 1 # 
###############

Z.Time 				= Z.1

z1.Time 			= seq(quantile(Z.Time[,1], qs[4]), quantile(Z.Time[,1], qs[5]), length=ngrid)
quants 				= quantile(Z.Time[,2], qs)
cross.sec 			= cbind(expand.grid(z1=z1.Time), z2 = quantile(Z.Time[,2], qs[j]), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]))

hgrid.cross.sec 	= newh.postmean.gfl.1(Time = 1:5, g = 1, numtime = 3, Znew=cross.sec, sel=sel)$postmean
hgrid.cross.sec 	= hgrid.cross.sec - mean(hgrid.cross.sec)
hgrid.cross.sec.sd 	= newh.postmean.gfl.1(Time = 1:5, g = 1, numtime = 3, Znew=cross.sec, sel=sel)$postvar
hgrid.sd 			= matrix(sqrt(diag(hgrid.cross.sec.sd)), length(z1.Time))
lb.grid 			= hgrid.cross.sec - 1.96*hgrid.sd
ub.grid 			= hgrid.cross.sec + 1.96*hgrid.sd

plot(z1.Time, hgrid.cross.sec, type="n", xlab="Z1", ylab="h", main=paste(names(quants)[j], "Z2, Time 1"), ylim = ylim, xlim=range(z1.Time))
polygon(c(z1.Time, rev(z1.Time)), c(lb.grid, rev(ub.grid)), border=FALSE, col="grey70")
lines(z1.Time, lb.grid, col="brown")
lines(z1.Time, ub.grid, col="brown")
lines(z1.Time, hgrid.cross.sec, lwd=2)
abline(h=0, col="gray", lty=2)

###############
# Timepoint 2 #
###############

Z.Time 				= Z.2

z1.Time 			= seq(quantile(Z.Time[,1], qs[4]), quantile(Z.Time[,1], qs[5]), length=ngrid) #
quants 				= quantile(Z.Time[,2], qs)
cross.sec 			= cbind(expand.grid(z1=z1.Time), z2 = quantile(Z.Time[,2], qs[j]), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]))

hgrid.cross.sec 	= newh.postmean.gfl.2(Time = 6:10, g = 2, numtime = 3, Znew=cross.sec, sel=sel)$postmean
hgrid.cross.sec 	= hgrid.cross.sec - mean(hgrid.cross.sec)
hgrid.cross.sec.sd 	= newh.postmean.gfl.2(Time = 6:10, g = 2, numtime = 3, Znew=cross.sec, sel=sel)$postvar
hgrid.sd 			= matrix(sqrt(diag(hgrid.cross.sec.sd)), length(z1.Time))
lb.grid 			= hgrid.cross.sec - 1.96*hgrid.sd
ub.grid 			= hgrid.cross.sec + 1.96*hgrid.sd

plot(z1.Time, hgrid.cross.sec, type="n", xlab="Z1", ylab="h", main=paste(names(quants)[j], "Z2, Time 2"), ylim = ylim, xlim=range(z1.Time))
polygon(c(z1.Time, rev(z1.Time)), c(lb.grid, rev(ub.grid)), border=FALSE, col="grey70")
lines(z1.Time, lb.grid, col="brown")
lines(z1.Time, ub.grid, col="brown")
lines(z1.Time, hgrid.cross.sec, lwd=2)
abline(h=0, col="gray", lty=2)

###############
# Timepoint 3 #
###############

Z.Time 				= Z.3

z1.Time 			= seq(quantile(Z.Time[,1], qs[4]), quantile(Z.Time[,1], qs[5]), length=ngrid)
quants 				= quantile(Z.Time[,2], qs)
cross.sec 			= cbind(expand.grid(z1=z1.Time), z2 = quantile(Z.Time[,2], qs[j]), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]))

hgrid.cross.sec 	= newh.postmean.gfl.3(Time = 11:15, g = 3, numtime = 3, Znew=cross.sec, sel=sel)$postmean
hgrid.cross.sec 	= hgrid.cross.sec - mean(hgrid.cross.sec)
hgrid.cross.sec.sd 	= newh.postmean.gfl.3(Time = 11:15, g = 3, numtime = 3, Znew=cross.sec, sel=sel)$postvar
hgrid.sd 			= matrix(sqrt(diag(hgrid.cross.sec.sd)), length(z1.Time))
lb.grid 			= hgrid.cross.sec - 1.96*hgrid.sd
ub.grid 			= hgrid.cross.sec + 1.96*hgrid.sd

plot(z1.Time, hgrid.cross.sec, type="n", xlab="Z1", ylab="h", main=paste(names(quants)[j], "Z2, Time 3"), ylim = ylim, xlim=range(z1.Time))
polygon(c(z1.Time, rev(z1.Time)), c(lb.grid, rev(ub.grid)), border=FALSE, col="grey70")
lines(z1.Time, lb.grid, col="brown")
lines(z1.Time, ub.grid, col="brown")
lines(z1.Time, hgrid.cross.sec, lwd=2)
abline(h=0, col="gray", lty=2)

########################################
# Bottom panel #
########################################

j = 3 

###############
# Timepoint 1 # 
###############

Z.Time 				= Z.1

z1.Time 			= seq(quantile(Z.Time[,1], qs[4]), quantile(Z.Time[,1], qs[5]), length=ngrid)
quants 				= quantile(Z.Time[,2], qs)
cross.sec 			= cbind(expand.grid(z1=z1.Time), z2 = quantile(Z.Time[,2], qs[j]), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]))

hgrid.cross.sec 	= newh.postmean.gfl.1(Time = 1:5, g = 1, numtime = 3, Znew=cross.sec, sel=sel)$postmean
hgrid.cross.sec 	= hgrid.cross.sec - mean(hgrid.cross.sec)
hgrid.cross.sec.sd 	= newh.postmean.gfl.1(Time = 1:5, g = 1, numtime = 3, Znew=cross.sec, sel=sel)$postvar
hgrid.sd 			= matrix(sqrt(diag(hgrid.cross.sec.sd)), length(z1.Time))
lb.grid 			= hgrid.cross.sec - 1.96*hgrid.sd
ub.grid 			= hgrid.cross.sec + 1.96*hgrid.sd

plot(z1.Time, hgrid.cross.sec, type="n", xlab="Z1", ylab="h", main=paste(names(quants)[j], "Z2, Time 1"), ylim = ylim, xlim=range(z1.Time))
polygon(c(z1.Time, rev(z1.Time)), c(lb.grid, rev(ub.grid)), border=FALSE, col="grey70")
lines(z1.Time, lb.grid, col="brown")
lines(z1.Time, ub.grid, col="brown")
lines(z1.Time, hgrid.cross.sec, lwd=2)
abline(h=0, col="gray", lty=2)

###############
# Timepoint 2 #
###############

Z.Time 				= Z.2

z1.Time 			= seq(quantile(Z.Time[,1], qs[4]), quantile(Z.Time[,1], qs[5]), length=ngrid) #
quants 				= quantile(Z.Time[,2], qs)
cross.sec 			= cbind(expand.grid(z1=z1.Time), z2 = quantile(Z.Time[,2], qs[j]), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]))

hgrid.cross.sec 	= newh.postmean.gfl.2(Time = 6:10, g = 2, numtime = 3, Znew=cross.sec, sel=sel)$postmean
hgrid.cross.sec 	= hgrid.cross.sec - mean(hgrid.cross.sec)
hgrid.cross.sec.sd 	= newh.postmean.gfl.2(Time = 6:10, g = 2, numtime = 3, Znew=cross.sec, sel=sel)$postvar
hgrid.sd 			= matrix(sqrt(diag(hgrid.cross.sec.sd)), length(z1.Time))
lb.grid 			= hgrid.cross.sec - 1.96*hgrid.sd
ub.grid 			= hgrid.cross.sec + 1.96*hgrid.sd

plot(z1.Time, hgrid.cross.sec, type="n", xlab="Z1", ylab="h", main=paste(names(quants)[j], "Z2, Time 2"), ylim = ylim, xlim=range(z1.Time))
polygon(c(z1.Time, rev(z1.Time)), c(lb.grid, rev(ub.grid)), border=FALSE, col="grey70")
lines(z1.Time, lb.grid, col="brown")
lines(z1.Time, ub.grid, col="brown")
lines(z1.Time, hgrid.cross.sec, lwd=2)
abline(h=0, col="gray", lty=2)

###############
# Timepoint 3 #
###############

Z.Time 				= Z.3

z1.Time 			= seq(quantile(Z.Time[,1], qs[4]), quantile(Z.Time[,1], qs[5]), length=ngrid)
quants 				= quantile(Z.Time[,2], qs)
cross.sec 			= cbind(expand.grid(z1=z1.Time), z2 = quantile(Z.Time[,2], qs[j]), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]))

hgrid.cross.sec 	= newh.postmean.gfl.3(Time = 11:15, g = 3, numtime = 3, Znew=cross.sec, sel=sel)$postmean
hgrid.cross.sec 	= hgrid.cross.sec - mean(hgrid.cross.sec)
hgrid.cross.sec.sd 	= newh.postmean.gfl.3(Time = 11:15, g = 3, numtime = 3, Znew=cross.sec, sel=sel)$postvar
hgrid.sd 			= matrix(sqrt(diag(hgrid.cross.sec.sd)), length(z1.Time))
lb.grid 			= hgrid.cross.sec - 1.96*hgrid.sd
ub.grid 			= hgrid.cross.sec + 1.96*hgrid.sd

plot(z1.Time, hgrid.cross.sec, type="n", xlab="Z1", ylab="h", main=paste(names(quants)[j], "Z2, Time 3"), ylim = ylim, xlim=range(z1.Time))
polygon(c(z1.Time, rev(z1.Time)), c(lb.grid, rev(ub.grid)), border=FALSE, col="grey70")
lines(z1.Time, lb.grid, col="brown")
lines(z1.Time, ub.grid, col="brown")
lines(z1.Time, hgrid.cross.sec, lwd=2)
abline(h=0, col="gray", lty=2)





