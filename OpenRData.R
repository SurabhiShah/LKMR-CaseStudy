# Opening RData file

library(lars); library(SuppDists); library(MCMCpack); library(magic); 
library(lasso2); library(mvtnorm); library(grplasso); library(kernlab); library(MASS); library(fields)

doginv = FALSE

offDiagonal = function(x) {
		diagBelow = diag(x)
		i = 1
		while (i <= n) {
			diagBelow=rbind(rep(0,length(x)	+i),cbind(diagBelow,rep(0,length(x) + i - 1)))
			i = i + 1
		}
		mat <- diagBelow + t(diagBelow) - diag(diag(diagBelow))
		return(mat)
}

load("CaseStudy3.RData")

Z = resultsGFLK$Z
n = nrow(Z)
grpsize = 5
numtime = 3
p		= n*numtime
Z.1 = resultsGFLK$Z[,1:grpsize]
Z.2 = resultsGFLK$Z[,(grpsize+1):(2*grpsize)]
Z.3 = resultsGFLK$Z[,(2*grpsize+1):(3*grpsize)]

sigsq0.post = resultsGFLK$Sigsq
sigsqf.post = resultsGFLK$Omega
tausqf.post = resultsGFLK$Tau
conf.post	= resultsGFLK$Conf
Y = resultsGFLK$Y
U = resultsGFLK$U
beta.f = resultsGFLK$Beta
sel = resultsGFLK$sel
cofactor = resultsGFLK$cofactor

