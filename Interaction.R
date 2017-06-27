##############################################
# Interaction of Mn-Zn at three time windows #
##############################################

# At Z2 75%, Z1 75% - Z1 25% (A1)
# At Z2 25%, Z1 75% - Z1 25% (A2)
# A1 - A2 is the interaction effect, while holding the other four metals at their median. 
# Interaction: (Z2_75%, Z1_75% - Z1_25%) - (Z2_25%, Z1_75% - Z1_25%) 

#######################################################################

source("PostmeanFunctions.R")

mat.time 		= NULL
mat.sd 			= NULL
qs 				= c(0.25, 0.75, 0.5)

###############
# Timepoint 1 #
###############

Z.Time 			= Z.1
cross.sec 		= rbind(apply(Z.Time, 2, median), apply(Z.Time, 2, median), apply(Z.Time, 2, median), apply(Z.Time, 2, median))
cross.sec[,1] 	= c(quantile(Z.Time[,1], qs[2]), quantile(Z.Time[,1], qs[1]), quantile(Z.Time[,1], qs[2]), quantile(Z.Time[,1], qs[1]))
cross.sec[,2] 	= c(quantile(Z.Time[,2], qs[2]), quantile(Z.Time[,2], qs[2]), quantile(Z.Time[,2], qs[1]), quantile(Z.Time[,2], qs[1]))
hgrid.cross.sec = newh.postmean.gfl.1(Time = 1:5, g = 1, numtime = numtime, Znew=cross.sec, sel=sel)
mat.time[1] 	= (hgrid.cross.sec$postmean[1] - hgrid.cross.sec$postmean[2]) - (hgrid.cross.sec$postmean[3] - hgrid.cross.sec$postmean[4])
mat.sd[1] 		= sqrt(hgrid.cross.sec$postvar[2,2] + hgrid.cross.sec$postvar[1,1] + hgrid.cross.sec$postvar[3,3] + hgrid.cross.sec$postvar[4,4] - 2*hgrid.cross.sec$postvar[1,2] - 2*hgrid.cross.sec$postvar[3,4])
		
###############
# Timepoint 2 #
###############

Z.Time 			= Z.2
cross.sec 		= rbind(apply(Z.Time, 2, median), apply(Z.Time, 2, median), apply(Z.Time, 2, median), apply(Z.Time, 2, median))
cross.sec[,1] 	= c(quantile(Z.Time[,1], qs[2]), quantile(Z.Time[,1], qs[1]), quantile(Z.Time[,1], qs[2]), quantile(Z.Time[,1], qs[1]))
cross.sec[,2] 	= c(quantile(Z.Time[,2], qs[2]), quantile(Z.Time[,2], qs[2]), quantile(Z.Time[,2], qs[1]), quantile(Z.Time[,2], qs[1]))
hgrid.cross.sec = newh.postmean.gfl.2(Time = 6:10, g = 2, numtime = numtime, Znew=cross.sec, sel=sel)
mat.time[2] 	= (hgrid.cross.sec$postmean[1] - hgrid.cross.sec$postmean[2]) - (hgrid.cross.sec$postmean[3] - hgrid.cross.sec$postmean[4])
mat.sd[2] 		= sqrt(hgrid.cross.sec$postvar[2,2] + hgrid.cross.sec$postvar[1,1] + hgrid.cross.sec$postvar[3,3] + hgrid.cross.sec$postvar[4,4] - 2*hgrid.cross.sec$postvar[1,2] - 2*hgrid.cross.sec$postvar[3,4])
	
###############
# Timepoint 3 #
###############
	
Z.Time 			= Z.3	
cross.sec 		= rbind(apply(Z.Time, 2, median), apply(Z.Time, 2, median), apply(Z.Time, 2, median), apply(Z.Time, 2, median))
cross.sec[,1] 	= c(quantile(Z.Time[,1], qs[2]), quantile(Z.Time[,1], qs[1]), quantile(Z.Time[,1], qs[2]), quantile(Z.Time[,1], qs[1]))
cross.sec[,2] 	= c(quantile(Z.Time[,2], qs[2]), quantile(Z.Time[,2], qs[2]), quantile(Z.Time[,2], qs[1]), quantile(Z.Time[,2], qs[1]))
hgrid.cross.sec = newh.postmean.gfl.3(Time = 11:15, g = 3, numtime = numtime, Znew=cross.sec, sel=sel)
mat.time[3] 	= (hgrid.cross.sec$postmean[1] - hgrid.cross.sec$postmean[2]) - (hgrid.cross.sec$postmean[3] - hgrid.cross.sec$postmean[4])
mat.sd[3] 		= sqrt(hgrid.cross.sec$postvar[2,2] + hgrid.cross.sec$postvar[1,1] + hgrid.cross.sec$postvar[3,3] + hgrid.cross.sec$postvar[4,4] - 2*hgrid.cross.sec$postvar[1,2] - 2*hgrid.cross.sec$postvar[3,4])
		
#######################################################################
	
ylim 			= range(c(mat.time -mat.sd, mat.time +mat.sd))
ylim 			= c(-3, 3)
plot(1:3, mat.time, xaxt="n",
    ylim=ylim, pch=19, xlab="Time Window", ylab="Z1-Z2 Interaction",
    main="Interaction of Z1 and Z2 at Three Critical Time Windows"
)

arrows(1:3, mat.time -1.96*mat.sd, 1:3, mat.time +1.96*mat.sd, length=0.05, angle=90, code=3)
axis(1, at=1:3, labels=c("Time 1", "Time 2", "Time 3"))
abline(h=0)




