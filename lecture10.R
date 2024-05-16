##
##  MCMC algorithm: Mixture sampler
##  for stochastic volatility model
##
##  y(t) = exp(h(t)/2)*z(t), z(t) ~ N(0, 1)
##  h(t+1) = mu + phi*(h(t) - mu) + eta(t)
##  eta(t) ~ N(0, sigma^2)
##  h(1) ~ N(mu, sigma^2/(1-phi^2))
##

##--- initial setup ---##
rm(list=ls())
gc()
gc()
library(ggplot2)
library(gridExtra)
tic.time <- proc.time()
set.seed(1)


##--- load data ---##

vp <- as.matrix(read.csv("lecture10_data.csv"))      # price data

nsp <- length(vp)      # sample size for original price series
vlp <- log(vp)         # log price

vy <- (vlp[2:nsp] - vlp[1:(nsp-1)]) * 100    # compute return

ns <- length(vy)       # the number of samples

# cat(log(var(vy)))    # sample variance


##--- priors ---##

da0 <- 20     # (phi+1)/2 ~ Beta(a0, b0)
db0 <- 1.5

dm0 <- 1      # mu ~ N(mu0, v0^2)
dv0 <- 1

dn0 <- 4      # sigma^2 ~ IG(n0/2, S0/2)
dS0 <- 0.4


##--- normal mixture (Omori et al. 2007) ---##

nK <- 10
vpj <- c(  .22715,  .20674,   .18842,  .13057,    .12047,
           .05591,  .04775,   .01575,  .00609,    .00115)
vmj <- c(- .85173,  .02266, -1.97278,  .73504,  -3.46788,
         -5.55246, 1.34744, -8.68384, 1.92677, -14.65)
vvj <- c(  .62699,  .40611,   .98583,  .26768,   1.57469,
          2.54498,  .17788,  4.16591,  .11265,   7.33342)

mpj <- matrix(1, ns, 1) %*% vpj  # covert to matrix form
mmj <- matrix(1, ns, 1) %*% vmj
mvj <- matrix(1, ns, 1) %*% vvj


##--- set variables ----##

dphi <- 0.95  # initial phi
dmu  <- 1     # initial mu
dsig <- 0.2   # initial sigma

vya <- log(vy^2 + 10^-7)        # add small constant
vh <- matrix(1, ns, 1) * dmu
vs <- matrix(0, ns, 1)

vxt <- vpt <- vnt <- vFt <- vLt <- veta <- matrix(0, ns, 1)


##-- set sampling option --##

nsim  <- 100                  # the number of iterations
nburn <- 0.1 * nsim             # burn-in period
npmt  <- 3                      # the number of parameters
msamp <- matrix(0, nsim, npmt)  # store box
vsamph <- vsamphs <- matrix(0, ns, 1)

##-- MCMC sampling --##

cat('\n## Iteration:\n')

##----------------- S A M P L I N G   S T A R T ----------------##

for(m_k in (-nburn+1):nsim){

  ##--- sampling s ---##

  mlp <- log(mpj) - 0.5 * log(mvj) - 0.5 * (
           matrix((vya - vh), ns, nK) - mmj )^2 / mvj

  mp <- exp(mlp)

  vpsum <- rowSums(mp)
  mpcum <- t(apply(mp, 1, cumsum)) / matrix(vpsum, ns, nK)
  vs <- rowSums(mpcum - matrix(runif(ns), ns, nK) < 0) + 1

  vmt <- vmj[vs]
  vvt <- vvj[vs]


  ##--- sampling h ---##

  ## Kalman filter
  vxt[1] <- dmu                     # x(1|0)
  vpt[1] <- dsig^2 / (1 - dphi^2)   # P(1|0)

  for (t in 1:ns) {
    vnt[t] <- vya[t] - vmt[t] - vxt[t]     # nu(t)
    vFt[t] <- vpt[t] + vvt[t]              # F(t)
    vLt[t] <- dphi * (1 - vpt[t] / vFt[t]) # L(t)

    dxtu <- vxt[t] + vpt[t] * vnt[t] / vFt[t]  # update: x(t|t)
    dptu <- vpt[t] - vpt[t]^2 / vFt[t]         # update: P(t|t)

    if (t < ns) {
      vxt[t+1] <- (1 - dphi) * dmu + dphi * dxtu  # forecast: x(t+1|t)
      vpt[t+1] <- dphi^2 * dptu + dsig^2          # forecast: P(t+1|t)
    }
  }

  ## Simulation smoother
  dr <- dU <- 0
  for (t in ns:1) {
    dC <- dsig^2 * (1 - dsig^2 * dU)
    deps <- sqrt(dC) * rnorm(1)
    veta[t] <- dsig^2 * dr + deps
    dV <- dsig^2 * dU * vLt[t]

    dr <- vnt[t] / vFt[t] + vLt[t] * dr - dV * deps / dC
    dU <- 1 / vFt[t] + vLt[t]^2 * dU + dV^2 / dC
  }

  ## compute h
  dC <- dsig^2 / (1 - dphi^2) * (1 - dsig^2 / (1 - dphi^2) * dU)
  deps <- sqrt(dC) * rnorm(1)
  vh[1] <- dmu + dsig^2 / (1 - dphi^2) * dr + deps
  for (t in 1:(ns-1)) {
    vh[t+1] <- dmu + dphi * (vh[t] - dmu) + veta[t]
  }


  ##--- sampling phi ---##

  vhb <- vh - dmu

  dsum <- sum(vhb[1:ns-1]^2)
  dphh <- sum(vhb[1:ns-1] * vhb[2:ns]) / dsum
  dlam <- dsig / sqrt(dsum)

  dphio <- dphi
  dphin <- 1
  while (abs(dphin) >= 1) {
    dphin <- dphh + dlam * rnorm(1)
  }
  dlpn <- (da0 - 1) * log(dphin + 1) + (db0 - 1) * log(-dphin + 1) +
          0.5 * log(1 - dphin^2) -
          0.5 * (1 - dphin^2) * vhb[1]^2 / dsig^2
  dlpo <- (da0 - 1) * log(dphio + 1) + (db0 - 1) * log(-dphio + 1) +
          0.5 * log(1 - dphio^2) -
          0.5 * (1 - dphio^2) * vhb[1]^2 / dsig^2
  dfrac <- exp(dlpn - dlpo)

  if (runif(1) < dfrac) {
    dphi <- dphin
  }

  ##--- sampling mu ---##

  dvh <- 1 / sqrt(1 / dv0^2 +
                    ( (1 - dphi^2) + (ns-1) * (1 - dphi)^2) / dsig^2 )
  dmh <- dvh^2 * (dm0 / dv0^2 +
                    ( (1 - dphi^2) * vh[1] +
                      (1 - dphi) * sum(vh[2:ns] - dphi*vh[1:ns-1])) / dsig^2)

  dmu <- dmh + dvh * rnorm(1)


  ##--- sampling sigma ----##

  vhb <- vh - dmu

  dnh <- dn0 + ns
  dSh <- dS0 + (1 - dphi^2) * vhb[1]^2 +
             sum((vhb[2:ns] - dphi*vhb[1:ns-1])^2)

  dsig <- 1 / sqrt(rgamma(1, dnh/2, dSh/2))


  ##-- storing sample --##

  if(m_k >= 1){
    msamp[m_k, ] <- cbind(dphi, dmu, dsig)
    vsamph  <- vsamph  + vh
    vsamphs <- vsamphs + vh^2
  }

  if((m_k %% 100) == 0) cat(m_k, '\n')

}
##------------------- S A M P L I N G   E N D ------------------##

##-- output result --##

## MCMC sample

vmean <- colMeans(msamp)
vstdv <- sqrt(colMeans(msamp^2) - vmean^2)
vlowq <- vuppq <- matrix(0, 1, npmt)
for (i in 1:npmt) {
  vlowq[i] <- quantile(msamp[, i], .025, type=1)
  vuppq[i] <- quantile(msamp[, i], .975, type=1)
}

message('\n[MCMC sampling result]')
cat('Iteration: ', sprintf('%i', nsim))

message('\n\nPosterior mean (stdev) [2.5%, 97.5%]')
cat('phi =', sprintf('%.4f', vmean[1]))
cat(' (', sprintf('%.4f', vstdv[1]), ')')
cat(' [', sprintf('%.4f', vlowq[1]), ',', sprintf('%.4f', vuppq[1]), ']\n')
cat('mu  =', sprintf('%.4f', vmean[2]))
cat(' (', sprintf('%.4f', vstdv[2]), ')')
cat(' [', sprintf('%.4f', vlowq[2]), ',', sprintf('%.4f', vuppq[2]), ']\n')
cat('sig =', sprintf('%.4f', vmean[3]))
cat(' (', sprintf('%.4f', vstdv[3]), ')')
cat(' [', sprintf('%.4f', vlowq[3]), ',', sprintf('%.4f', vuppq[3]), ']\n')

fTsvar <- function(vx, iBm){
  ns <- length(vx)
  vz <- c(1 : ns)
  for (i in 1:ns) {
    lm.x <- lm(vx ~ vz)
    ar.x <- ar(vx, aic = TRUE)
    dsp <- ar.x$var.pred / (1 - sum(ar.x$ar))^2
  }
  return(dsp)
}

fGeweke <- function(vx, iBm) {
  ns <- length(vx)
  n1 <- floor(0.1 * ns)
  n2 <- floor(0.5 * ns)
  vx1 <- vx[1:n1]
  vx2 <- vx[(ns-n2+1):ns]
  dx1m <- mean(vx1)
  dx2m <- mean(vx2)
  dx1v <- fTsvar(vx1)
  dx2v <- fTsvar(vx2)
  dz <- (dx1m - dx2m) / sqrt(dx1v / n1 + dx2v / n2)

  return(2 * (1 - pnorm(abs(dz))))
}

vGw <- matrix(0, npmt, 1)
for (i in 1:npmt) {
  vGw[i] <- fGeweke(msamp[, i])
}

message('\nGeweke statistics (p-value)')
cat('phi:', sprintf('%.3f', vGw[1]), '\n')
cat('mu :', sprintf('%.3f', vGw[2]), '\n')
cat('sig:', sprintf('%.3f', vGw[3]), '\n')

## volatility
vsamph  <- vsamph / nsim
vsamphs <- sqrt(vsamphs / nsim - vsamph^2)
vsampl  <- vsamph - vsamphs
vsampu  <- vsamph + vsamphs

vylim <- c(floor(min(exp(vsampl/2))),
           ceiling(max(exp(vsampu/2))))

dev.new()
par(mfrow = c(2, 2))
plot(vp, type='l',
     col='2', main='Price (p)',
     xlab='t', ylab='P(t)')
plot(vy, type='l',
     col='2', main='Return (y)',
     xlab='t', ylab='y(t)')
plot(exp(vsamph/2), type='l',
     xlim=c(1, ns), ylim=vylim, 
     col='2', main='Estimated volatility (exp(h(t)/2))',
     xlab='t', ylab='sigma(t)')
par(new=T)
plot(exp(vsampl/2), type="l", col=3, lty=1,
     xlim=c(1, ns), ylim=vylim, ann=F)
par(new=T)
plot(exp(vsampu/2), type="l", col=3, lty=1,
     xlim=c(1, ns), ylim=vylim, ann=F)


## MCMC path

# function: plotline draws line plot
plotline <- function(vy, sname){

  df.out <- data.frame(n = 1:length(vy), ydata = vy)
  plot.out <-
    ggplot(df.out, aes(x=n, y=ydata, colour="red")) +
    geom_line() +
    guides(colour="none") +
    labs(x ="", y ="", title = sname) +
    theme_classic()

  return(plot.out)
}

plot.pb1 <- plotline(msamp[, 1], "phi")
plot.pb2 <- plotline(msamp[, 2], "mu")
plot.pb3 <- plotline(msamp[, 3], "sigma")

## posterior distribution

## function: plotdensity draws density plot
plotdensity <- function(vx, sname){

  df.out <- data.frame(xdata = vx)
  plot.out <-
    ggplot(df.out, aes(x=xdata)) +
    geom_line(stat = "density") +
    expand_limits(y=0) +
    theme(legend.position = 'none') +
    labs(x="", y="", title = sname) +
    theme_classic()

  return(plot.out)
}

plot.db1 <- plotdensity(msamp[, 1], "phi")
plot.db2 <- plotdensity(msamp[, 2], "mu")
plot.db3 <- plotdensity(msamp[, 3], "sigma")

dev.new()
plot.all <-
  grid.arrange(plot.pb1, plot.pb2, plot.pb3,
               plot.db1, plot.db2, plot.db3,
               ncol = 3)

##-- run time --##

total.time <- proc.time() - tic.time
message('\nTime(s): ', round(total.time[3], 1))