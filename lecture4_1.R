##
##  Stochastic volatility model
##
##  epsilon(t) = exp(x(t)/2)*z(t), z(t) ~ N(0, 1)
##  x(t) = log(sigma(t)^2)
##  x(t+1) = omega + phi*x(t) + eta(t)
##  eta(t) ~ N(0, sigma_eta^2)
##

##--- initial setup ---##
rm(list=ls())
gc()
gc()
set.seed(1)


##--- set true values ---##

dome <- -0.1    # omega
dphi <- 0.95   # phi
dsig <- 0.2    # sigma_eta


##--- set variables ---##

ns <- 1000     # the number of observations

veps <- vx <- matrix(0, ns, 1)

dinm <- dome / (1 - dphi)      # stationary mean for initial
dimv <- dsig^2 / (1 - dphi^2)  # stationary variance for initial


##--- simulate observations ---##

ve <- rnorm(ns)           # simulation z(t)
veta <- rnorm(ns)         # simulation eta(t)

vx[1] <- dinm + sqrt(dimv) * rnorm(1);
for (t in 1:ns) {
  veps[t] <- exp(vx[t]/2) * ve[t]   # compute return
  if (t < ns) {
    vx[t+1] <- dome + dphi * vx[t] + dsig * veta[t]
  }
}


##--- draw figure ---##

dev.new()
par(mfrow = c(2, 1))
plot(exp(vx/2), type='l',
     col='2', main='Simulated volatility (sigma)',
     xlab='t', ylab='sigma(t)')
plot(veps, type='l',
     col='2', main='Simulated return (y)',
     xlab='t', ylab='epsilon(t)')


##--- output to csv file ---##

write.csv(cbind(veps, vx), "lecture4_data.csv", row.names=FALSE)

