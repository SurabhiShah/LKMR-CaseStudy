###############################################################
# Relative importance of each metal at three time windows     #
# h_75 - h_25, holding the other four metals at their median  #
###############################################################

#source("OpenRData.R")
doginv=FALSE
source("PostmeanFunctions.R")

sel = resultsGFLK$sel

ylim 		= c(-3, 3)
mat.time 	= matrix(NA, nrow=numtime, ncol=grpsize)
mat.sd 		= matrix(NA, nrow=numtime, ncol=grpsize)
qs 			= c(0.25, 0.75, 0.5)

par(mfrow=c(1,3))

###############
# Timepoint 1 #
###############

Z.Time 					= Z.1
for (j in 1:grpsize) {
		cross.sec 		= rbind(apply(Z.Time, 2, median), apply(Z.Time, 2, median))
		cross.sec[,j] 	= c(quantile(Z.Time[,j], qs[2]), quantile(Z.Time[,j], qs[1]))
		hgrid.cross.sec = newh.postmean.gfl.1(Time = 1:5, g = 1, numtime = numtime, Znew=cross.sec)
		mat.time[1,j] 	= hgrid.cross.sec$postmean[1] - hgrid.cross.sec$postmean[2]
		mat.sd[1,j] 	= sqrt(hgrid.cross.sec$postvar[2,2] + hgrid.cross.sec$postvar[1,1] - 2*hgrid.cross.sec$postvar[1,2])
}
	
plot(1:5, mat.time[1,], xaxt="n",
    ylim=ylim,
    pch=15, xlab="Metal", ylab="Main effect of each metal",
    main="Time 1"
)

arrows(1:5, mat.time[1,] -1.96*mat.sd[1,], 1:5, mat.time[1,] +1.96*mat.sd[1,], length=0.05, angle=90, code=3)
axis(1, at=1:5, labels=c("Z1", "Z2", "Z3", "Z4", "Z5"))
abline(h=0)

###############
# Timepoint 2 #
###############

Z.Time 					= Z.2
for (j in 1:grpsize) {
		cross.sec 		= rbind(apply(Z.Time, 2, median), apply(Z.Time, 2, median))
		cross.sec[,j] 	= c(quantile(Z.Time[,j], qs[2]), quantile(Z.Time[,j], qs[1]))
		hgrid.cross.sec = newh.postmean.gfl.2(Time = 6:10, g = 2, numtime = numtime, Znew=cross.sec)
		mat.time[2,j] 	= hgrid.cross.sec$postmean[1] - hgrid.cross.sec$postmean[2]
		mat.sd[2,j] 	= sqrt(hgrid.cross.sec$postvar[2,2] + hgrid.cross.sec$postvar[1,1] - 2*hgrid.cross.sec$postvar[1,2])
}
	
plot(1:5, mat.time[2,], xaxt="n",
    ylim=ylim,
    pch=15, xlab="Metal", ylab="Main effect of each metal",
    main="Time 2"
)

arrows(1:5, mat.time[2,] -1.96*mat.sd[2,], 1:5, mat.time[2,] +1.96*mat.sd[2,], length=0.05, angle=90, code=3)
axis(1, at=1:5, labels=c("Z1", "Z2", "Z3", "Z4", "Z5"))
abline(h=0)

###############
# Timepoint 3 #
###############
	
Z.Time = Z.3
for (j in 1:grpsize) {
		cross.sec 		= rbind(apply(Z.Time, 2, median), apply(Z.Time, 2, median))
		cross.sec[,j] 	= c(quantile(Z.Time[,j], qs[2]), quantile(Z.Time[,j], qs[1]))
		hgrid.cross.sec = newh.postmean.gfl.3(Time = 11:15, g = 3, numtime = numtime, Znew=cross.sec, sel=sel)
		mat.time[3,j] 	= hgrid.cross.sec$postmean[1] - hgrid.cross.sec$postmean[2]
		mat.sd[3,j] 	= sqrt(hgrid.cross.sec$postvar[2,2] + hgrid.cross.sec$postvar[1,1] - 2*hgrid.cross.sec$postvar[1,2])
}
	
plot(1:5, mat.time[3,], xaxt="n",
    ylim=ylim,
    pch=15, xlab="Metal", ylab="Main effect of each metal",
    main="Time 3"
)

arrows(1:5, mat.time[3,] -1.96*mat.sd[3,], 1:5, mat.time[3,] +1.96*mat.sd[3,], length=0.05, angle=90, code=3)
axis(1, at=1:5, labels=c("Z1", "Z2", "Z3", "Z4", "Z5"))
abline(h=0)
	
