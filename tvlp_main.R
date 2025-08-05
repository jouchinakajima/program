##
##  [tvlp_main.R] Time-varying local projections
##                with stochastic volatility (SV)
##
##  Coded by: Jouchi Nakajima
##  Last update: 2025/08/05
##  Website:  http://sites.google.com/site/jnakajimaweb/
##
##  ** You may use and modify this code at your own risk
##


##--- initial setup ---##
rm(list=ls())
gc()
gc()
source("tvlp_svsamp.R")
tic.time <- proc.time()
set.seed(1)


##--- [SET] options ---##

iSV <- 1     # SV (0: off, 1: on)

nh <- 8      # maximum horizons to draw impulse responses



##--- [SET] data file ---##

mdata <- as.matrix(read.csv('tvlp_data.csv'))  # load data

nsa <- nrow(mdata)


##--------- E A C H   H O R I Z O N   S T A R T -----------##

for (h in 1:nh) {

  ##--- [SET] LP specification ---##

  ## dependent variable (y)
  vy <- mdata[(h+2):nsa, 2] - mdata[2:(nsa-h), 2]

  ## independent variables for TV coefficients (z)
  mz <- cbind(mdata[3:(nsa-h+1), 3],
              mdata[2:(nsa-h),   4],
              mdata[2:(nsa-h),   2] - mdata[1:(nsa-h-1), 2],
              matrix(1, nsa-h-1, 1))

  ## independent variable for time-invariant coefficients (x)
  mx <- cbind(mdata[(h+2):nsa,   5],
              mdata[2:(nsa-h),   5])


  ##--- set variables ---##

  ns <- length(vy)
  if (is.matrix(mx)) {
    nk <- ncol(mx)
  } else if (length(mx) > 1) {
    mx <- matrix(mx, ns, 1)
    nk <- 1
  } else nk <- 0
  if (is.matrix(mz)) {
    np <- ncol(mz)
  } else if (length(mz) > 1)  {
    mz <- matrix(mz, ns, 1)
    np <- 1
  } else message('\n### Error: No z_t! ###\n')

  vb  <- rep(0, nk)
  vv2 <- rep(1, np) * 0.1^2
  ds2 <- 0.1^2
  dw2 <- 0.1^2

  ma <- matrix(0, ns, np)
  vs2 <- rep(0.1, ns)

  vr <- rep(0, np)
  mU <- diag(np) * 0
  mT <- diag(np)
  ve <- vDi <- rep(0, ns)
  amL <- array(0, dim = c(np, np, ns))
  meta <- malpha <- matrix(0, ns, np)
  vxb <- vza <- rep(0, ns)


  ##--- set priors ---##

  vb0 <- rep(0, nk)        # beta ~ N(b0, B0)
  mB0 <- diag(nk) * 10

  vv0 <- rep(1, np) * 40   # v^2 ~ IG(v0/2, V0/2)
  vV0 <- rep(1, np) * 4

  dw0 <- 40                # w^2 ~ IG(w0/2, w0/2)
  dW0 <- 4

  ds0 <- 20                # s^2 ~ IG(s0/2, S0/2)
  dS0 <- 4

  va0 <- rep(0, np)        # a_1 ~ N(a0, H02)
  mH02 <- diag(np) * 10

  dh00 <- 3                # h_1 ~ N(h00, sig02)
  dsig02 <- 2

  if (nk > 0) {
    miB0 <- solve(mB0)
    vBb0 <- miB0 %*% vb0
  }
  vvh <- vv0 + ns - 1
  dwh <- dw0 + ns - 1
  dsh <- ds0 + ns


  ##-- set sampling option --##

  nsim  <- 20000                  # the number of iterations
  nburn <- 0.1 * nsim             # burn-in period
  npmt  <- nk + np + 1            # the number of parameters
  msamp <- matrix(0, nsim, npmt)  # store box
  msampa <- msampas <- matrix(0, ns, np)
  vsamph <- vsamphs <- rep(0, ns)
  nK <- floor(ns/30)-1            # blocks for sampling h


  ##-- MCMC sampling --##

  cat(paste('\n## [h = ', h, ']', sep=''))
  cat('\nIteration:\n')

  ##----------------- S A M P L I N G   S T A R T ----------------##

  for(k in (-nburn+1):nsim){

    if (nk > 0) {

      ##--- sampling beta ---##

      ms2 <- matrix(vs2, ns, nk)
      mBh <- solve(miB0 + t(mx) %*% (mx / ms2))
      vbh <- mBh %*% (vBb0 + t(mx) %*% ((vy - vza) / vs2))

      vb <- vbh + t(chol(mBh)) %*% rnorm(nk)
      vxb <- mx %*% vb

    }

    if (np > 0) {

      ##--- sampling alpha ---##

      vya <- vy - vxb

      va <- va0
      mP <- mH02
      if (np > 1) mH2 <- diag(vv2)
      else mH2 <- vv2
      vr <- vr * 0
      mU <- mU * 0

      for (t in 1:ns) {
        ve[t] <- vya[t] - mz[t, ] %*% va
        vDi[t] <- 1 / (mz[t, ] %*% mP %*% mz[t, ] + vs2[t])

        vK <- mP %*% mz[t, ] * vDi[t]
        amL[, , t] <- mT - vK %*% mz[t, ]

        va <- va + vK %*% ve[t]
        mP <- mP %*% t(amL[, , t]) + mH2
      }

      for (t in ns:1) {
        mC <- mH2 - mH2 %*% mU %*% mH2
        mC <- (mC + t(mC)) / 2

        mCi <- solve(mC)

        veps <- t(chol(mC)) %*% rnorm(np)
        meta[t, ] <- mH2 %*% vr + veps

        mV <- mH2 %*% mU %*% amL[, , t]
        vr <- mz[t, ] * vDi[t] * ve[t] + t(amL[, , t]) %*% vr -
          t(mV) %*% mCi %*% veps
        mU <- matrix(mz[t, ] * vDi[t], np, 1) %*% matrix(mz[t, ], 1, np) +
          t(amL[, , t]) %*% mU %*% amL[, , t] + t(mV) %*% mCi %*% mV
      }

      mC <- mH02 - mH02 %*% mU %*% mH02
      mC <- (mC + t(mC)) / 2
      veta0 <- mH02 %*% vr + chol(mC) %*% rnorm(np)

      ma[1, ] <- va0 + veta0
      for (t in 1:(ns-1)) {
        ma[t+1, ] <- ma[t, ] + meta[t, ]
      }

      vza <- rowSums(mz * ma)

      ##--- sampling v ---##

      vVh <- vV0 + colSums(meta^2)

      for (i in 1:np) {
        vv2[i] <- 1 / rgamma(1, vvh[i]/2, vVh[i]/2)
      }

    }

    if (iSV == 1) {

      ##--- sampling volatility ---##

      vya <- vy - vxb - vza
      vh <- log(vs2)

      vh <- svsamp(vya, vh, dw2, dh00, dsig02, nK)
      vs2 <- exp(vh)


      ##--- sampling w ---##

      dWh <- dW0 + sum(diff(vh)^2)

      dw2 <- 1 / rgamma(1, dwh/2, dWh/2)


    } else {

      ##--- sampling sigma ---##

      dSh <- dS0 + sum((vy - vxb - vza)^2)

      ds2 <- 1 / rgamma(1, dsh/2, dSh/2)
      vs2 <- matrix(ds2, ns, 1)
      vh  <- matrix(log(ds2), ns, 1)

    }

    if (k >= 1) {

      ##-- storing sample --##

      vsamp <- c(vb, sqrt(vv2))
      if (iSV == 1) msamp[k, ] <- c(vsamp, sqrt(dw2))
      else          msamp[k, ] <- c(vsamp, sqrt(ds2))

      if (np > 0) {
        msampa  <- msampa  + ma
        msampas <- msampas + ma^2
      }
      vsamph  <- vsamph  + vh
      vsamphs <- vsamphs + vh^2

    }

    if((k %% 2000) == 0) cat(k, '\n')

  }
  ##------------------- S A M P L I N G   E N D ------------------##

  ##--- output results ---##

  ## Parameters
  mout <- matrix(0, npmt, 4)
  for (i in 1:npmt) {
    mout[i, ] <- c(mean(msamp[, i]), sd(msamp[, i]),
                   quantile(msamp[, i], .025, type=1),
                   quantile(msamp[, i], .975, type=1))
  }
  vspar <- numeric(0)
  for (i in 1:nk) {
    vspar <- c(vspar, paste('beta-', i, sep=''))
  }
  for (i in 1:np) {
    vspar <- c(vspar, paste('v-', i, '   ', sep=''))
  }
  if (iSV == 1) {
    vspar <- c(vspar, 'w     ')
  } else {
    vspar <- c(vspar, 'sigma ')
  }

  cat('\n')
  cat('--------------------------------------------\n')
  cat('Parameter  Mean     Stdev    95%L     95%U\n')
  cat('--------------------------------------------\n')
  for (i in 1:npmt) {
    cat(paste(vspar[i], ' ',
              sprintf('%7.4f', mout[i, 1]), '',
              sprintf('%7.4f', mout[i, 2]), '',
              sprintf('%7.4f', mout[i, 3]), '',
              sprintf('%7.4f', mout[i, 4]), '\n'))
  }
  cat('--------------------------------------------\n')
  cat('TV-LP with ')
  if (iSV == 0) cat('Constant Variance\n')
  else cat('SV\n')

  ## TVP
  msampa  <- msampa / nsim
  msampas <- sqrt(msampas / nsim - msampa^2)

  ## SV
  vsamph  <- vsamph / nsim
  vsamphs <- sqrt(vsamphs / nsim - vsamph^2)

  mout <- matrix(0, ns, (np+1)*2)
  svar <- rep('', (np+1)*2)
  for (i in 1:np) {
    mout[, (i-1)*2+(1:2)] <- cbind(msampa[, i], msampas[, i])
    svar[(i-1)*2+(1:2)] <- cbind(paste('Var', i, '-Mean', sep=''),
                                 paste('Var', i, '-SD', sep=''))
  }
  mout[, i*2+(1:2)] <- cbind(vsamph, vsamphs)
  svar[i*2+(1:2)] <- c('SV-Mean', 'SV-SD')

  colnames(mout) <- svar
  write.csv(mout, paste(ifelse(iSV == 1, 'result_sv', 'result_cv'),
                        h, '.csv', sep=''), row.names = F)

}

##--- run time ---##

total.time <- proc.time() - tic.time
message('\nTime(s): ', round(total.time[3], 1))
