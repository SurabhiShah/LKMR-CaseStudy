###################################
# Graphs of Z1-Z2 Heatmap Surface #
###################################

# Corresponds to LKMR-Main

doginv = FALSE
source("PostmeanFunctions.R")

ngrid 			= 30
j 				= 2 
min.plot.dist 	= 0.5 # only color level plot for points within this distance from an observed data point

ylim 			= c(-5, 5)
qs 				= c(0.1, 0.5, 0.9, 0.05, 0.95)

###############
# Timepoint 1 #
###############

Z.Time 		= Z.1

z1.1 		= seq(quantile(Z.Time[,1], qs[4]), quantile(Z.Time[,1], qs[5]), length=ngrid)
z1.2 		= seq(quantile(Z.Time[,2], qs[4]), quantile(Z.Time[,2], qs[5]), length=ngrid)

cross.sec.1 = cbind(expand.grid(z1=z1.1, z2 = z1.2), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]))

grid 		= expand.grid(z1=z1.1, z2 = z1.2)
mindists 	= rep(NA, nrow(grid))
for(i in seq_along(mindists)) {
	pt 			= as.numeric(grid[i,])
	dists 		= as.matrix(dist(rbind(pt, Z.Time[,1:2])))["pt",-1]
	mindists[i] = min(dists)
}
rm(grid, pt, dists)

hgrid.T.1 	= newh.postmean.gfl.1(Time = 1:5, g = 1, numtime = 3, Znew=cross.sec.1, sel=sel)$postmean
hgrid.T.1[mindists > min.plot.dist] = NA
hgrid.T.1 	= matrix(hgrid.T.1, nrow=ngrid)

###############
# Timepoint 2 #
###############

Z.Time 		= Z.2
z2.1 		= seq(quantile(Z.Time[,1], qs[4]), quantile(Z.Time[,1], qs[5]), length=ngrid)
z2.2 		= seq(quantile(Z.Time[,2], qs[4]), quantile(Z.Time[,2], qs[5]), length=ngrid)

cross.sec.2 = cbind(expand.grid(z1=z2.1, z2 = z2.2), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]))

grid 		= expand.grid(z1=z2.1, z2 = z2.2)
mindists 	= rep(NA, nrow(grid))
for(i in seq_along(mindists)) {
	pt 			= as.numeric(grid[i,])
	dists 		= as.matrix(dist(rbind(pt, Z.Time[,1:2])))["pt",-1]
	mindists[i] = min(dists)
}
 
hgrid.T.2 	= newh.postmean.gfl.2(Time = 6:10, g = 2, numtime = 3, Znew=cross.sec.2, sel=sel)$postmean
hgrid.T.2[mindists > min.plot.dist] = NA
hgrid.T.2 	= matrix(hgrid.T.2, nrow=ngrid)

###############
# Timepoint 3 #
###############

Z.Time 		= Z.3

z3.1 		= seq(quantile(Z.Time[,1], qs[4]), quantile(Z.Time[,1], qs[5]), length=ngrid)
z3.2 		= seq(quantile(Z.Time[,2], qs[4]), quantile(Z.Time[,2], qs[5]), length=ngrid)

cross.sec.3 = cbind(expand.grid(z1=z3.1, z2 = z3.2), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]))

grid 		= expand.grid(z1=z3.1, z2 = z3.2)
mindists 	= rep(NA, nrow(grid))
for(i in seq_along(mindists)) {
	pt 			= as.numeric(grid[i,])
	dists 		= as.matrix(dist(rbind(pt, Z.Time[,1:2])))["pt",-1]
	mindists[i] = min(dists)
}

hgrid.T.3 	= newh.postmean.gfl.3(Time = 11:15, g = 3, numtime = 3, Znew=cross.sec.3, sel=sel)$postmean
hgrid.T.3[mindists > min.plot.dist] = NA
hgrid.T.3 	= matrix(hgrid.T.3, nrow=ngrid)

################
# Heatmap plot #
################

zlim = c(-8,8)

par(oma=c( 0,0,0,4)) 
set.panel( 1,3) 

image(z1.1, z1.2, hgrid.T.1, xlab="Z1", ylab="Z2", col=tim.colors(), main="Time 1", zlim=zlim)
points(Z[,1:5])
image(z2.1, z2.2, hgrid.T.2, xlab="Z1", ylab="Z2", col=tim.colors(), main="Time 2", zlim=zlim)
points(Z[,6:10])
image(z3.1, z3.2, hgrid.T.3, xlab="Z1", ylab="Z2", col=tim.colors(), main="Time 3", zlim=zlim)
points(Z[,11:15])

par(oma=c( 0,0,0,1))
image.plot( legend.only=TRUE, zlim=zlim) 

