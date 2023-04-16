##
##  GARCH (1, 1) model: Maximum likelihood estimation
##
##  y(t) = a + b * y(t-1) + e(t)
##  e(t) = sig(t) * z(t)
##  sig(t)^2 = omega + beta * sig(t-1)^2 + alpha * e(t-1)^2
##  omega > 0, beta >= 0, alpha >= 0
##  z(t) ~ N(0, 1)
##
##  Version: 2023/04/17
##

##--- initial setup ---##
rm(list=ls())
gc()
gc()


##--- load data ---##

vp <- as.matrix(read.csv("lecture3_data.csv"))      # price data

nsp <- length(vp)      # sample size for original price series
vlp <- log(vp)         # log price

vy <- (vlp[2:nsp] - vlp[1:(nsp-1)]) * 100    # compute return


##--- 1st step: regression for return equation ---##

ns <- length(vy)                  # sample size for return series
nsi <- ns - 1
vyi <- matrix(vy[2:ns], nsi, 1)             # dependent variable
mxi <- cbind(matrix(1, nsi, 1), vy[1:nsi])  # independent variables

nki <- ncol(mxi)
mx2 <- solve(t(mxi) %*% mxi)
vb  <- mx2 %*% (t(mxi) %*% vyi)   # OLS estimator

vei <- vyi - mxi %*% vb           # residuals
drss <- t(vei) %*% vei            # sum of squared residuals
ds2 <- drss / (nsi - nki)         # variance estimator
mS <- sum(ds2) * mx2              # covariance matrix
vbs <- sqrt(diag(mS))             # standard errors

cat('\n[1st-step result: OLS estimator (std.err.)]')
cat('\n a =', sprintf('%.4f', vb[1]))
cat(' (', sprintf('%.4f', vbs[1]), ')')
cat('\n b =', sprintf('%.4f', vb[2]))
cat(' (', sprintf('%.4f', vbs[2]), ')')


##--- 2nd step: maximum likelihood estimation for GARCH part ---##

doi <- 0.1    # initial value: omega
dbi <- 0.8    #                beta
dai <- 0.1    #                alpha

## function to compute (negative) log-likelihood of GARCH model
NLL <- function(vpar, ve) {

  nse <- length(ve)

  domeg <- vpar[1]    # omega
  dbeta <- vpar[2]    # beta
  dalph <- vpar[3]    # alpha

  if (domeg > 0 && dbeta >= 0 && dalph >= 0 && (dbeta + dalph) < 1) {

  dini <- domeg / (1 - dbeta - dalph)   # initial volatility: stationary mean

  ve2 <- ve^2
  vsig2 <- matrix(0, nse, 1)

  vsig2[1] <- dini     # initial volatility: sig(1)
  vsig2[2] <- domeg + dbeta * vsig2[1] + dalph * ve2[1]
  for (t in 2:(nse-1)) {
    vsig2[t+1] <- domeg + dbeta * vsig2[t] + dalph * ve2[t]
  }

  dlik <- -0.5 * sum(log(vsig2) + ve2 / vsig2)   # log-likelihood

  } else {
    dlik <- - 10^10
  }

  return(- dlik)        # negative log-likelihood
}
## function ends

mle <- nlm(NLL, cbind(doi, dbi, dai), ve = vei, hessian=T)    # minimize negative log-likelihood

vparm <- mle$estimate   # MLE
vpars <- sqrt(diag(solve(mle$hessian)))  # MLE standard error

cat('\n\n[2nd-step result: MLE estimator (std.err.)]')
cat('\n omega =', sprintf('%.4f', vparm[1]))
cat(' (', sprintf('%.4f', vpars[1]), ')')
cat('\n beta  =', sprintf('%.4f', vparm[2]))
cat(' (', sprintf('%.4f', vpars[2]), ')')
cat('\n alpha =', sprintf('%.4f', vparm[3]))
cat(' (', sprintf('%.4f', vpars[3]), ')')


##--- compute volatility ---##

domeg <- vparm[1]    # omega
dbeta <- vparm[2]    # beta
dalph <- vparm[3]    # alpha

dini <- domeg / (1 - dbeta - dalph)   # initial volatility: stationary mean

ve2 <- vei^2
vsig2 <- matrix(0, nsi, 1)

vsig2[1] <- dini     # initial volatility: sig(1)
vsig2[2] <- domeg + dbeta * vsig2[1] + dalph * ve2[1]
for (t in 2:(nsi-1)) {
  vsig2[t+1] <- domeg + dbeta * vsig2[t] + dalph * ve2[t]
}


##--- draw figure ---##

dev.new()
par(mfrow = c(2, 2))
plot(vp, type='l',
     col='2', main='Price (p)',
     xlab='t', ylab='P(t)',
     cex.main = 1)
plot(vy, type='l',
     col='2', main='Return (y)',
     xlab='t', ylab='y(t)',
     cex.main = 1)
plot(sqrt(vsig2), type='l',
     col='2', main='Estimated volatility (sigma)',
     xlab='t', ylab='sigma(t)',
     cex.main = 1)

