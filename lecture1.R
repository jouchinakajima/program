##
##  Stationary process: AR(1)
##
##  y(t) = mu' + phi*y(t-1) + e(t), e(t) ~ N(0, sigma^2)
##
##  Version: 2022/11/04
##

##--- initial setup ---##
rm(list=ls())
gc()
gc()

##--- set true value ---##

# dphi <- 0.95   # phi: AR(1) coefficient
sphi <- readline("phi?\n")
dphi <- as.numeric(sphi)

dmus <- 0.2   # mu': intercept

dsig <- 0.1   # sigma: error standard deviation


##--- set variables ---##

ns <- 500    # the number of observations

vy <- matrix(0, ns, 1)


##--- simulate observations ---##

if (dphi < 1) {
  vy[1] <- dmus / (1 - dphi)
} else {
  vy[1] <- 0
}

for (t in 1:(ns-1)) {
  vy[t+1] <- dmus + dphi * vy[t] + rnorm(1) * dsig
}


##--- draw figure ---##

dev.new()
plot(vy, type='l', 
     col='2', main='Simulated series',
     xlab='t', ylab='y(t)')

dev.new()
plot(vy[1:(ns-1)], vy[2:ns],
     col='2', pch=20, main='Simulated series',
     xlab='y(t-1)', ylab='y(t)')


##--- some statistics ---##

cat('\nTheoretical unconditional mean:\n ')
cat(sprintf('%.3f', dmus/(1-dphi)))
cat('\nTheoretical unconditional variance:\n ')
cat(sprintf('%.3f', dsig^2/(1-dphi^2)))
cat('\nSample mean:\n ')
cat(sprintf('%.3f', mean(vy)))
cat('\nSample variance:\n ')
cat(sprintf('%.3f', var(vy)))

