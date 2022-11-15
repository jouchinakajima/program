##
##  Stochastic volatility model
##  - Kalman filter and smoother
##
##  epsilon(t) = exp(x(t)/2)*z(t), z(t) ~ N(0, 1)
##  x(t) = log(sigma(t)^2)
##  x(t+1) = omega + phi*x(t) + eta(t)
##  eta(t) ~ N(0, sigma_eta^2)
##
##  [state space representation]
##  y(t) = log(epsilon(t)^2) + 1.27
##  u(t) = log(z(t)^2) + 1.27
##
##    y(t) = x(t) + u(t),   u(t) ~ N(0, pi^2/2)
##    x(t) = x(t-1) + eta(t)
##

##--- initial setup ---##
rm(list=ls())
gc()
gc()


##--- load data ---##

mdata <- as.matrix(read.csv("lecture4_data.csv"))

veps   <- mdata[, 1]  # epsilon(t)
vxtrue <- mdata[, 2]  # true of x(t)


##--- set parameters ---##

dome <- -0.1   # omega
dphi <- 0.95   # phi
dsig <- 0.2    # sigma_eta


##--- set variables ---##

dc <- 0.0001        # adjustment constant
vy <- log(veps^2 + dc) + 1.27   # y(t)

ns <- length(vy)    # the number of observations

vxtu <- vptu <- vn <- vF     <- matrix(0, ns, 1)
vxtf <- vptf <- vxts <- vpts <- matrix(0, ns+1, 1)

dsu2 <- pi^2 / 2  # sigma_u


##--- Kalman filter ---##

vxtf[1] <- dome / (1 - dphi)       # x(1|0)
vptf[1] <- dsig^2 / (1 - dphi^2)   # P(1|0)

for (t in 1:ns) {
  vn[t] <- vy[t] - vxtf[t]             # nu(t)
  vF[t] <- vptf[t] + dsu2              # F(t)
  
  vxtu[t] <- vxtf[t] + vptf[t] * vn[t] / vF[t] # update: x(t|t)
  vptu[t] <- vptf[t] - vptf[t]^2 / vF[t]       # update: P(t|t)
  
  vxtf[t+1] <- dome + dphi * vxtu[t]      # forecast: x(t+1|t)
  vptf[t+1] <- dphi^2 * vptu[t] + dsig^2  # forecast: P(t+1|t)
}


##--- smoother ---##

vxts[ns+1] <- vxtf[ns+1]   # x(T+1|T)
vpts[ns+1] <- vptf[ns+1]   # P(T+1|T)
for (t in ns:1) {
  dps <- dphi * vptu[t] / vptf[t+1]    # P*(t)
  
  vxts[t] <- vxtu[t] + dps * (vxts[t+1] - vxtf[t+1])   # x(t|T)
  vpts[t] <- vptu[t] + dps^2 * (vpts[t+1] - vptf[t+1]) # P(t|T)
}


##--- likelihood ---##

dllik <- -0.5 * sum(log(2*pi) + log(vF) + vn^2 / vF)  # log likelihood

cat('\nLog likelihood =', sprintf('%.4f', dllik))


##--- draw figure ---##

nd <- 200;
vyl <- c(-4, 0)
dev.new()
plot(vxtrue[1:nd], type='l', col='2', 
     xlab='t', ylab='x(t)', xlim=c(1, nd), ylim=vyl)
par(new=T)
plot(vxtu[1:nd], type='l', col='3',
     xlim=c(1, nd), ylim=vyl, ann=F)
par(new=T)
plot(vxts[1:nd], type='l', col='4',
     xlim=c(1, nd), ylim=vyl, ann=F)
legend("topleft", 
       legend = c("True", "filtered", "smoothed"), lty=1, col=2:4)
title("Kalman filter and smoother: x(t)")
