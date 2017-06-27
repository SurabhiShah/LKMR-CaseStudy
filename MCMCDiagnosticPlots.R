# Diagnostic plots for LKMR:

sel = 1:num.reps

#######################################
# Trace plots:
#######################################

par(mfrow = c(4,2))

plot(x = 1:length(resultsGFLK$Sigsq[sel]), y = resultsGFLK$Sigsq[sel], type = "n", main = "Trace of Sigma_Sq", xlab = "Iterations")
lines(x = 1:length(resultsGFLK$Sigsq[sel]), y = resultsGFLK$Sigsq[sel])

plot(x = 1:length(resultsGFLK$Lam1[sel]), y = resultsGFLK$Lam1[sel], type = "n", main = "Trace of Lambda_1",  xlab = "Iterations")
lines(x = 1:length(resultsGFLK$Lam1[sel]), y = resultsGFLK$Lam1[sel])

plot(x = 1:length(resultsGFLK$Lam2[sel]), y = resultsGFLK$Lam2[sel], type = "n", main = "Trace of Lambda_2",  xlab = "Iterations")
lines(x = 1:length(resultsGFLK$Lam2[sel]), y = resultsGFLK$Lam2[sel])

plot(x = 1:length(resultsGFLK$Tau[sel,1]), y = resultsGFLK$Tau[sel,1], type = "n", main = "Trace of Tau_1", xlab = "Iterations")
lines(x = 1:length(resultsGFLK$Tau[sel,1]), y = resultsGFLK$Tau[sel,1])

plot(x = 1:length(resultsGFLK$Tau[sel,2]), y = resultsGFLK$Tau[sel,2], type = "n", main = "Trace of Tau_2",  xlab = "Iterations")
lines(x = 1:length(resultsGFLK$Tau[sel,2]), y = resultsGFLK$Tau[sel,2])

plot(x = 1:length(resultsGFLK$Tau[sel,3]), y = resultsGFLK$Tau[sel,3], type = "n", main = "Trace of Tau_3",  xlab = "Iterations")
lines(x = 1:length(resultsGFLK$Tau[sel,3]), y = resultsGFLK$Tau[sel,3])

plot(x = 1:length(resultsGFLK$Omega[sel,1]), y = resultsGFLK$Omega[sel,1], type = "n", main = "Trace of Omega_1", xlab = "Iterations")
lines(x = 1:length(resultsGFLK$Omega[sel,1]), y = resultsGFLK$Omega[sel,1])

plot(x = 1:length(resultsGFLK$Omega[sel,2]), y = resultsGFLK$Omega[sel,2], type = "n", main = "Trace of Omega_2", xlab = "Iterations")
lines(x = 1:length(resultsGFLK$Omega[sel,2]), y = resultsGFLK$Omega[sel,2])

#########################################
# Trace plots for h
#########################################

# Overall h
# Average h per time window

par(mfrow = c(2,3))
h1 = apply(resultsGFLK$Beta[sel,1:n], 1, mean)
plot(density(h1), main = "Density of h_1")

h2 = apply(resultsGFLK$Beta[sel,(n+1):(2*n)], 1, mean)
plot(density(h1), main = "Density of h_2")

h3 = apply(resultsGFLK$Beta[sel,(2*n+1):(3*n)], 1, mean)
plot(density(h1), main = "Density of h_3")

plot(x = 1:length(h1), y = h1, type = "n", main = "Trace of h_1",  xlab = "Iterations")
lines(x = 1:length(h1), y = h1)

plot(x = 1:length(h2), y = h2, type = "n", main = "Trace of h_2",  xlab = "Iterations")
lines(x = 1:length(h2), y = h2)

plot(x = 1:length(h3), y = h3, type = "n", main = "Trace of h_3",  xlab = "Iterations")
lines(x = 1:length(h3), y = h3)

#########################################
# Marginal density plots
#########################################

par(mfrow = c(4,2))

plot(density(resultsGFLK$Sigsq[sel]), main = "Density of Sigma_Sq")

plot(density(resultsGFLK$Lam1[sel]), main = "Density of Lambda_1")

plot(density(resultsGFLK$Lam2[sel]), main = "Density of Lambda_2")

plot(density(resultsGFLK$Tau[sel,1]), main = "Density of Tau_1")

plot(density(resultsGFLK$Tau[sel,2]), main = "Density of Tau_2")

plot(density(resultsGFLK$Tau[sel,3]), main = "Density of Tau_3")

plot(density(resultsGFLK$Omega[sel,1]), main = "Density of Omega_1")

plot(density(resultsGFLK$Omega[sel,2]), main = "Density of Omega_2")


