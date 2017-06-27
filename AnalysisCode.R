#################################
# Case Study Data Analysis Code #
#################################

set.seed(1)

library(lars); library(SuppDists); library(MCMCpack); library(magic);
library(lasso2); library(mvtnorm); library(grplasso); library(kernlab); library(MASS); library(fields)

source("DataExtraction.R")

##########
# Priors #
##########

REP 		  = 1
cofactor	= 1e-7
num.reps	= 10000
sel 		  = seq(2000, num.reps, by=5)
a.1 		  = 60 
b.1 		  = 10
a.2 		  = 45
b.2 		  = 10
sig.shape 	= 5 
sig.scale 	= 5

lambda1.sq 	= rgamma(1,shape=a.1, rate=b.1)
lambda2.sq 	= rgamma(1,shape=a.2, rate=b.2)
sig.sq0  	  = rinvgamma(1, shape = sig.shape, scale = sig.scale)
mean.bef 	  = rep(0,p)
tau.sqf  	  = rgamma(numtime, shape = (n + 1)/2, rate=lambda1.sq/2)
sig.sqf  	  = rexp(numtime-1, rate=lambda2.sq/2)

##############################################################

##################
# Initialization #
##################

#####
# h #
#####

# Construct Sigma_h inverse covariance matrix

# 1. Block diagonal matrix

list.G 						= list()
poly 						= polydot(degree=2, offset=1)
for (g in 1:numtime) {
    list.G[[g]] 			= solve(kernelMatrix(poly, Z[,(grpsize*g-(grpsize-1)):(grpsize*g)]) + cofactor*diag(n)) / tau.sqf[g]
}
cov.bf1 					= do.call(adiag, list.G)

# 2. Diagonal

list.sig.sqf 				= c()
list.sig.sqf[1] 			= 1/sig.sqf[1]
list.sig.sqf[numtime] 		= 1/sig.sqf[(numtime-1)]
for (g in 2:(numtime - 1)) {
    list.sig.sqf[g] 		= 1/sig.sqf[g] + 1/sig.sqf[g-1]
}

list.mat.sig.sqf 			= list()
for (g in 1:numtime) {
    list.mat.sig.sqf[[g]] 	= list.sig.sqf[g] * diag(n)
}

cov.bf2 					= do.call(adiag, list.mat.sig.sqf)

# 3. Off-diagonals
new.sig.sqf 				= rep(1/sig.sqf, each=n)
cov.bf3 					= offDiagonal((-1)*new.sig.sqf)

# 4. Sigma_h inverse covariance matrix
cov.bf 						= cov.bf1 + cov.bf2 + cov.bf3

# Invert cov.bf by adding a small cofactor to diagonals
SIG       					= sig.sq0*solve(cov.bf + cofactor*diag(n*numtime))
beta.fp    					= rmvnorm(1,mean=mean.bef,sigma=SIG)

#################
# for posterior #
#################

sigsq0.post = lambda1f.post = lambda2f.post = NULL
beta.f      = rbind( beta.fp,matrix(rep(NA,num.reps*p),ncol=p) )
tausqf.post = rbind( tau.sqf,matrix(rep(NA,num.reps*numtime),ncol=numtime) )
sigsqf.post = rbind( sig.sqf,matrix(rep(NA,num.reps*(numtime-1)),ncol=numtime-1) )
conf.post 	= rbind(conf, matrix(rep(NA,num.reps*c),nrow=num.reps))
conf 		= matrix(conf, nrow = c)

##################################################################

#############
# Begin job #
#############

cat(c("Job started at:",date()),fill=TRUE)
source("MCMC.R")
cat(c("Job finished at:",date()),fill=TRUE)

resultsGFLK = list(Sigsq = sigsq0.post, Lam1 = lambda1f.post, Lam2 = lambda2f.post, Beta = beta.f, Tau = tausqf.post, Omega = sigsqf.post, Z = Z, Y=Y, U= U, Conf = conf.post, sel = sel, cofactor=cofactor)

save(resultsGFLK, file = "CaseStudy3.RData")







