##
##  Gibbs sampler: regression model
##
##  y(t) = beta*x(t) + e(t), e(t) ~ N(0, sigma^2)
##
##  Prior: beta ~ N(b0, v0^2)
##         sigma^2 ~ IG(n0/2, S0/2)
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

dbt <- 0.5   # beta
dst <- 0.2   # sigma


##--- generate observation ---##

nT <- 500          # the number of observations
vx <- rnorm(nT)    # x(t)
vy <- dbt * vx + rnorm(nT) * dst  # y(t)


##--- priors ---##

db0 <- 0      # beta ~ N(b0, v0^2)
dv0 <- 10

dn0 <- 4      # sigma^2 ~ IG(n0/2, S0/2)
dS0 <- 0.4


##--- set variables ----##

dbeta <- 0    # initial beta
dsig <- 1     # initial sigma

dx2 <- sum(vx^2)     # x(t)^2
dxy <- sum(vx * vy)  # x(t)*y(t)


##--- set sampling option ---##

nsim  <- 10000                  # the number of iterations
nburn <- 0.1 * nsim             # burn-in period
npmt  <- 2                      # the number of parameters
msamp <- matrix(0, nsim, npmt)  # store box


##--- MCMC sampling ---##

cat('\n## Iteration:\n')

##----------------- S A M P L I N G   S T A R T ----------------##

for(m_k in (-nburn+1):nsim){
  
  ##--- sampling beta ---##

  dvh <- (1 / dv0^2 + dx2 / dsig^2)^(-1/2)
  dbh <- dvh^2 * (db0^2 / dv0^2 + dxy / dsig^2)
  
  dbeta <- dbh + dvh * rnorm(1)
  
  
  ##--- sampling sigma ----##

  dnh <- dn0 + nT  
  dSh <- dS0 + sum((vy - dbeta * vx)^2)
  
  dsig <- 1 / sqrt(rgamma(1, dnh/2, dSh/2))

  
  ##-- storing sample --##
  
  if(m_k >= 1){
    msamp[m_k, ] <- cbind(dbeta, dsig)
  }
  
  if((m_k %% 100) == 0) cat(m_k, '\n')
  
}
##------------------- S A M P L I N G   E N D ------------------##

##--- output result ---##

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
## function end

plot.pb1 <- plotline(msamp[, 1], "beta")
plot.pb2 <- plotline(msamp[, 2], "sigma")

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
## function end

plot.db1 <- plotdensity(msamp[, 1], "beta")
plot.db2 <- plotdensity(msamp[, 2], "sigma")

dev.new()
plot.all <-
  grid.arrange(plot.pb1, plot.pb2, plot.db1, plot.db2,
               ncol = 2)
plot(plot.all)


##-- run time --##

total.time <- proc.time() - tic.time
message('\nTime(s): ', round(total.time[3], 1))
