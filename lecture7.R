##
##  Gibbs sampler: normal mixture model
##
##  y(t)
##    = mu1 + beta*x(t) + e(t) (w.p. 1/2)
##    = mu2 + beta*x(t) + e(t) (w.p. 1/2)
##  e(t) ~ N(0, sigma^2)
##  mu1 < mu2
##
##  Prior: beta ~ N(b0, v0^2)
##         sigma^2 ~ IG(n0/2, S0/2)
##         mu1 ~ TN(mu01, w01^2)[-infty,mu2]
##         mu2 ~ TN(mu02, w02^2)[mu1,infty]
##

##--- initial setup ---##
rm(list=ls())
gc()
gc()
library(ggplot2)
library(gridExtra)
tic.time <- proc.time()
set.seed(1)


##--- set true values ---##

dbt <- 2    # beta
dst <- 0.2  # sigma
dm1t <- -1  # mu1
dm2t <- 1   # mu2   


##--- generate observation ---##

ns <- 500          # the number of observations
vx <- rnorm(ns)    # x(i)
vst <- (runif(ns) < 0.5) + 1   # true s(i)
vy <- matrix(0, ns, 1)   # y(i)
for (i in 1:ns) {
  if (vst[i] == 1) {
    vy[i] <- dm1t + dbt * vx[i] + rnorm(1) * dst
  } else {
    vy[i] <- dm2t + dbt * vx[i] + rnorm(1) * dst
  }
}

##--- priors ---##

db0 <- 0      # beta ~ N(b0, v0^2)
dv0 <- 10

dn0 <- 4      # sigma^2 ~ IG(n0/2, S0/2)
dS0 <- 0.4

dm01 <- 0     # mu1 ~ TN(mu01, w01^2)[-infty,mu2]
dw01 <- 10

dm02 <- 0     # mu2 ~ TN(mu02, w02^2)[mu1,infty]
dw02 <- 10


##--- set variables ----##

dbeta <- 0    # initial beta
dsig <- 0.5   # initial sigma
dmu1 <- -1    # initial mu1
dmu2 <- 1     # initial mu2
vs <- (runif(ns) < 0.5) + 1   # initial s

dx2 <- sum(vx^2)         # x(i)^2
vmu <- matrix(0, ns, 1)  # mu(i)


##-- set sampling option --##

nsim  <- 1000                   # the number of iterations
nburn <- 0.1 * nsim             # burn-in period
npmt  <- 4                      # the number of parameters
msamp <- matrix(0, nsim, npmt)  # store box

##-- MCMC sampling --##

cat('\n## Iteration:\n')

##----------------- S A M P L I N G   S T A R T ----------------##

for(m_k in (-nburn+1):nsim){
  
  vmu[vs == 1] <- dmu1
  vmu[vs == 2] <- dmu2
  
  ##--- sampling beta ---##
  
  dxy <- sum((vy - vmu) * vx)
  
  dvh <- (1 / dv0^2 + dx2 / dsig^2)^(-1/2)
  dbh <- dvh^2 * (db0 / dv0^2 + dxy / dsig^2)
  
  dbeta <- dbh + dvh * rnorm(1)
  
  
  ##--- sampling sigma ----##

  dnh <- dn0 + ns
  dSh <- dS0 + sum((vy - vmu - dbeta * vx)^2)

  dsig <- 1 / sqrt(rgamma(1, dnh/2, dSh/2))


  ##--- sampling mu1 ---##

  n1 <- sum(vs == 1)
  dxy <- sum(vy[vs==1] - dbeta * vx[vs==1])

  dwh <- (1 / dw01^2 + n1 / dsig^2)^(-1/2)
  dmh <- dwh^2 * (dm01 / dw01^2 + dxy / dsig^2)

  dmu1 <- dmu2
  while (dmu1 >= dmu2) {
    dmu1 <- dmh + dwh * rnorm(1)
  }

  ##--- sampling mu2 ---##

  n2 <- sum(vs == 2)
  dxy <- sum(vy[vs==2] - dbeta * vx[vs==2])

  dwh <- (1 / dw02^2 + n2 / dsig^2)^(-1/2)
  dmh <- dwh^2 * (dm02 / dw02^2 + dxy / dsig^2)

  dmu2 <- dmu1
  while (dmu1 >= dmu2) {
    dmu2 <- dmh + dwh * rnorm(1)
  }

  ##--- sampling s ---##

  vla1 <- -0.5 * (vy - dmu1 - dbeta * vx)^2 / dsig^2
  vla2 <- -0.5 * (vy - dmu2 - dbeta * vx)^2 / dsig^2

  vp <- exp(vla1) / (exp(vla1) + exp(vla2))
  vs <- (runif(ns) > vp) + 1
  
  
  ##-- storing sample --##
  
  if(m_k >= 1){
    msamp[m_k, ] <- cbind(dbeta, dsig, dmu1, dmu2)
  }
  
  if((m_k %% 100) == 0) cat(m_k, '\n')
  
}
##------------------- S A M P L I N G   E N D ------------------##

##-- output result --##

## MCMC path

## function: plotline draws line plot
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

plot.pb1 <- plotline(msamp[, 1], "beta")
plot.pb2 <- plotline(msamp[, 2], "sigma")
plot.pb3 <- plotline(msamp[, 3], "mu1")
plot.pb4 <- plotline(msamp[, 4], "mu2")

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

plot.db1 <- plotdensity(msamp[, 1], "beta")
plot.db2 <- plotdensity(msamp[, 2], "sigma")
plot.db3 <- plotdensity(msamp[, 3], "mu1")
plot.db4 <- plotdensity(msamp[, 4], "mu2")

dev.new()
plot.all1 <-
  grid.arrange(plot.pb1, plot.pb2,
               plot.db1, plot.db2,
               ncol = 2)
dev.new()
plot.all2 <-
  grid.arrange(plot.pb3, plot.pb4,
               plot.db3, plot.db4,
               ncol = 2)

##-- run time --##

total.time <- proc.time() - tic.time
message('\nTime(s): ', round(total.time[3], 1))
