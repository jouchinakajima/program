##
##  GARCH (1, 1) model
##
##  y(t) = a + b * y(t-1) + e(t)
##  e(t) = sig(t) * z(t)
##  sig(t)^2 = omega + beta * sig(t-1)^2 + alpha * e(t-1)^2
##  omega > 0, beta >= 0, alpha >= 0
##

##--- initial setup ---##
rm(list=ls())
gc()
gc()
set.seed(1)

##--- set true value ---##

da <- 0         # a
db <- 0         # b

domeg <- 0.5    # omega
dbeta <- 0.85   # beta
dalph <- 0.1    # alpha


##--- set variables ---##

ns <- 1500      # the number of observations

vy <- vsig <- ve <- matrix(0, ns, 1)

dini <- domeg / (1 - dbeta - dalph)   # initial volatility: stationary mean


##--- simulate observations ---##

vz <- rnorm(ns)           # simulation z(t)
vsig[1] <- sqrt(dini)     # initial volatility: sig(1)
ve[1] <- vsig[1] * vz[1]  # initial disturbance: e(1)
vy[1] <- da + ve[1]       # initial return: y(1)
vsig[2] <- sqrt(domeg + dbeta * vsig[1]^2 + dalph * ve[1]^2)
for (t in 2:ns) {
  ve[t] <- vsig[t] * vz[t]
  vy[t] <- da + db * vy[t-1] + ve[t]
  if (t < ns) {
    vsig[t+1] <- sqrt(domeg + dbeta * vsig[t]^2 + dalph * ve[t]^2)
  }
}


##--- draw figure ---##

dev.new()
par(mfrow = c(2, 1))
plot(vsig, type='l',
     col='2', main='Simulated volatility (sigma)',
     xlab='t', ylab='sigma(t)')
plot(vy, type='l',
     col='2', main='Simulated return (y)',
     xlab='t', ylab='y(t)')


##--- output to csv file ---##

vp <- matrix(0, ns, 1)
vp[1] <- 10523 + rnorm(1)
for (t in 1:(ns-1)) {
  vp[t+1] <- vp[t] * (1 + vy[t+1] / 100)
}
write.csv(vp, "psim.csv", row.names=FALSE)
