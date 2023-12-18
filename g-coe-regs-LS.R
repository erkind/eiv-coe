# OLS estimation of factor models for ind portfolios

# data
source(paste("http://raw.githubusercontent.com/erkind/",
             "eiv-coe/main/g-coe-data-eiv.R",
             sep = ""))
dim(data); colnames(data)[c(1:5,60:64)]
dim(R) ; colnames(R)
dim(RF); colnames(RF)
dim(X) ; colnames(X)

# packages
library(car)
library(lmtest)
library(sandwich)
library(skedastic)

# counters
nb.obs <- as.numeric(dim(R)[1])
nb.port<- as.numeric(dim(R)[2])

# time-series regressions
# 
# CAPM
M1.LS     <- list()
M1.LS.VCV <- list()
M1.LS.a   <- matrix(NA, nrow = nb.port, ncol = 1)
M1.LS.res <- matrix(NA, nrow = nb.obs, ncol = nb.port)
for(i in 1:nb.port){
   # OLS fit
   LS.fit <- lm(R[,i] ~ X[,1])
   names(LS.fit$coefficients) <- c("a", "b")
   # coefficient covariance matrix
   LS.fit.VCV <- vcovHAC(LS.fit)
   # extract summary values
   Rsq   <- summary(LS.fit)$adj.r.squared
   sigma <- summary(LS.fit)$sigma
   a     <- summary(LS.fit)$coefficients[1,1]
   a.SE  <- as.numeric(diag(LS.fit.VCV)[1])^0.5
   a.t   <- a / a.SE
   b     <- summary(LS.fit)$coefficients[2,1]
   b.SE  <- as.numeric(diag(LS.fit.VCV)[2])^0.5
   b.t   <- b / b.SE
   # excess returns time-series averages
   avg.R <- mean(R[,i])
   # factors time-series average
   avg.mkt <- mean(X[,1])
   # RF time-series average
   avg.RF <- mean(RF)
   # coe (annualized)
   coe <- (avg.RF
           + b * avg.mkt
           ) * 12
   # list to export results
   out.list <- list(
      Rsq = Rsq,
      sigma = sigma,
      a = a, a.SE = a.SE, a.t = a.t,
      b = b, b.SE = b.SE, b.t = b.t,
      avg.R = avg.R,
      avg.mkt = avg.mkt,
      avg.RF = avg.RF,
      coe = coe
      )
   # append results to current list
   M1.LS[[length(M1.LS) + 1]] <- out.list
   # coefficients covariance matrix
   M1.LS.VCV[[i]] <- vcovHAC(LS.fit)
   # intercepts
   M1.LS.a[i] <- a
   rownames(M1.LS.a) <- colnames(R)
   # residuals
   M1.LS.res[, i] <- LS.fit$residuals
   colnames(M1.LS.res) <- colnames(R)
}
# 
# FF3F
M3.LS     <- list()
M3.LS.VCV <- list()
M3.LS.a   <- matrix(NA, nrow = nb.port, ncol = 1)
M3.LS.res <- matrix(NA, nrow = nb.obs, ncol = nb.port)
for(i in 1:nb.port){
   # OLS fit
   LS.fit <- lm(R[,i] ~ X[,1] + X[,2] + X[,3])
   names(LS.fit$coefficients) <- c("a", "b", "s", "h")
   # coefficients covariance matrix
   LS.fit.VCV <- vcovHAC(LS.fit)
   # residuals
   LS.fit.res <- residuals(LS.fit)
   # extract summary values
   Rsq   <- summary(LS.fit)$adj.r.squared
   sigma <- summary(LS.fit)$sigma
   a     <- summary(LS.fit)$coefficients[1,1]
   a.SE  <- as.numeric(diag(LS.fit.VCV)[1])^0.5
   a.t   <- a / a.SE
   b     <- summary(LS.fit)$coefficients[2,1]
   b.SE  <- as.numeric(diag(LS.fit.VCV)[2])^0.5
   b.t   <- b / b.SE
   s     <- summary(LS.fit)$coefficients[3,1]
   s.SE  <- as.numeric(diag(LS.fit.VCV)[3])^0.5
   s.t   <- s / s.SE
   h     <- summary(LS.fit)$coefficients[4,1]
   h.SE  <- as.numeric(diag(LS.fit.VCV)[4])^0.5
   h.t   <- h / h.SE
   # excess returns time-series averages
   avg.R <- mean(R[,i])
   # factors time-series average
   avg.mkt <- mean(X[,1])
   avg.smb <- mean(X[,2])
   avg.hml <- mean(X[,3])
   # RF time-series average
   avg.RF  <- mean(RF)
   # coe (annualized)
   coe <- (avg.RF
           + b * avg.mkt + s * avg.smb + h * avg.hml
           ) * 12
   # list to export the results
   out.list <- list(
      Rsq = Rsq,
      sigma = sigma,
      a = a, a.SE = a.SE, a.t = a.t,
      b = b, b.SE = b.SE, b.t = b.t,
      s = s, s.SE = s.SE, s.t = s.t,
      h = h, h.SE = h.SE, h.t = h.t,
      avg.R = avg.R,
      avg.mkt = avg.mkt, avg.smb = avg.smb, avg.hml = avg.hml,
      avg.RF = avg.RF,
      coe = coe
   )
   # append results to current list
   M3.LS[[length(M3.LS) + 1]] <- out.list
   # coefficients covariance matrix
   M3.LS.VCV[[i]] <- vcovHAC(LS.fit)
   # intercepts
   M3.LS.a[i] <- a
   rownames(M3.LS.a) <- colnames(R)
   # residuals
   M3.LS.res[, i] <- LS.fit$residuals
   colnames(M3.LS.res) <- colnames(R)
}
# 
# FFC4F
M4.LS     <- list()
M4.LS.VCV <- list()
M4.LS.a   <- matrix(NA, nrow = nb.port, ncol = 1)
M4.LS.res <- matrix(NA, nrow = nb.obs, ncol = nb.port)
for(i in 1:nb.port){
   # OLS fit
   LS.fit <- lm(R[,i] ~ X[,1] + X[,2] + X[,3] + X[,6])
   names(LS.fit$coefficients) <- c("a", "b", "s", "h", "w")
   # coefficients covariance matrix
   LS.fit.VCV <- vcovHAC(LS.fit)
   # residuals
   LS.fit.res <- residuals(LS.fit)
   # extract summary values
   Rsq   <- summary(LS.fit)$adj.r.squared
   sigma <- summary(LS.fit)$sigma
   a     <- summary(LS.fit)$coefficients[1,1]
   a.SE  <- as.numeric(diag(LS.fit.VCV)[1])^0.5
   a.t   <- a / a.SE
   b     <- summary(LS.fit)$coefficients[2,1]
   b.SE  <- as.numeric(diag(LS.fit.VCV)[2])^0.5
   b.t   <- b / b.SE
   s     <- summary(LS.fit)$coefficients[3,1]
   s.SE  <- as.numeric(diag(LS.fit.VCV)[3])^0.5
   s.t   <- s / s.SE
   h     <- summary(LS.fit)$coefficients[4,1]
   h.SE  <- as.numeric(diag(LS.fit.VCV)[4])^0.5
   h.t   <- h / h.SE
   w     <- summary(LS.fit)$coefficients[5,1]
   w.SE  <- as.numeric(diag(LS.fit.VCV)[5])^0.5
   w.t   <- w / w.SE
   # excess returns time-series averages
   avg.R <- mean(R[,i])
   # factors time-series average
   avg.mkt <- mean(X[,1])
   avg.smb <- mean(X[,2])
   avg.hml <- mean(X[,3])
   avg.wml <- mean(X[,6])
   # RF time-series average
   avg.RF  <- mean(RF)
   # coe (annualized)
   coe <- (avg.RF
           + b * avg.mkt + s * avg.smb + h * avg.hml
           + w * avg.wml
           ) * 12
   # list to export the results
   out.list <- list(
      Rsq = Rsq,
      sigma = sigma,
      a = a, a.SE = a.SE, a.t = a.t,
      b = b, b.SE = b.SE, b.t = b.t,
      s = s, s.SE = s.SE, s.t = s.t,
      h = h, h.SE = h.SE, h.t = h.t,
      w = w, w.SE = w.SE, w.t = w.t,
      avg.R = avg.R,
      avg.mkt = avg.mkt, avg.smb = avg.smb, avg.hml = avg.hml,
      avg.wml = avg.wml,
      avg.RF = avg.RF,
      coe = coe
   )
   # append results to current list
   M4.LS[[length(M4.LS) + 1]] <- out.list
   # coefficients covariance matrix
   M4.LS.VCV[[i]] <- vcovHAC(LS.fit)
   # intercepts
   M4.LS.a[i] <- a
   rownames(M4.LS.a) <- colnames(R)
   # residuals
   M4.LS.res[, i] <- LS.fit$residuals
   colnames(M4.LS.res) <- colnames(R)
}
# 
# FF5F
M5.LS     <- list()
M5.LS.VCV <- list()
M5.LS.a   <- matrix(NA, nrow = nb.port, ncol = 1)
M5.LS.res <- matrix(NA, nrow = nb.obs, ncol = nb.port)
for(i in 1:nb.port){
   # OLS fit
   LS.fit <- lm(R[,i] ~ X[,1] + X[,2] + X[,3] + X[,4] + X[,5])
   names(LS.fit$coefficients) <- c("a", "b", "s", "h", "r", "c")
   # coefficients covariance matrix
   LS.fit.VCV <- vcovHAC(LS.fit)
   # residuals
   LS.fit.res <- residuals(LS.fit)
   # extract summary values
   Rsq   <- summary(LS.fit)$adj.r.squared
   sigma <- summary(LS.fit)$sigma
   a     <- summary(LS.fit)$coefficients[1,1]
   a.SE  <- as.numeric(diag(LS.fit.VCV)[1])^0.5
   a.t   <- a / a.SE
   b     <- summary(LS.fit)$coefficients[2,1]
   b.SE  <- as.numeric(diag(LS.fit.VCV)[2])^0.5
   b.t   <- b / b.SE
   s     <- summary(LS.fit)$coefficients[3,1]
   s.SE  <- as.numeric(diag(LS.fit.VCV)[3])^0.5
   s.t   <- s / s.SE
   h     <- summary(LS.fit)$coefficients[4,1]
   h.SE  <- as.numeric(diag(LS.fit.VCV)[4])^0.5
   h.t   <- h / h.SE
   r     <- summary(LS.fit)$coefficients[5,1]
   r.SE  <- as.numeric(diag(LS.fit.VCV)[5])^0.5
   r.t   <- r / r.SE
   c     <- summary(LS.fit)$coefficients[6,1]
   c.SE  <- as.numeric(diag(LS.fit.VCV)[6])^0.5
   c.t   <- c / c.SE
   # excess returns time-series averages
   avg.R <- mean(R[,i])
   # factors time-series average
   avg.mkt <- mean(X[,1])
   avg.smb <- mean(X[,2])
   avg.hml <- mean(X[,3])
   avg.rmw <- mean(X[,4])
   avg.cma <- mean(X[,5])
   # RF time-series average
   avg.RF  <- mean(RF)
   # coe (annualized)
   coe <- (avg.RF
           + b * avg.mkt + s * avg.smb + h * avg.hml
           + r * avg.rmw + c * avg.cma
           ) * 12
   # list to export the results
   out.list <- list(
      Rsq = Rsq,
      sigma = sigma,
      a = a, a.SE = a.SE, a.t = a.t,
      b = b, b.SE = b.SE, b.t = b.t,
      s = s, s.SE = s.SE, s.t = s.t,
      h = h, h.SE = h.SE, h.t = h.t,
      r = r, r.SE = r.SE, r.t = r.t,
      c = c, c.SE = c.SE, c.t = c.t,
      avg.R = avg.R,
      avg.mkt = avg.mkt, avg.smb = avg.smb, avg.hml = avg.hml,
      avg.rmw = avg.rmw, avg.cma = avg.cma,
      avg.RF = avg.RF,
      coe = coe
   )
   # append results to current list
   M5.LS[[length(M5.LS) + 1]] <- out.list
   # coefficients covariance matrix
   M5.LS.VCV[[i]] <- vcovHAC(LS.fit)
   # intercepts
   M5.LS.a[i] <- a
   rownames(M5.LS.a) <- colnames(R)
   # residuals
   M5.LS.res[, i] <- LS.fit$residuals
   colnames(M5.LS.res) <- colnames(R)
}

# collect results
M1.LS <- cbind(colnames(R), do.call(rbind, M1.LS))
M3.LS <- cbind(colnames(R), do.call(rbind, M3.LS))
M4.LS <- cbind(colnames(R), do.call(rbind, M4.LS))
M5.LS <- cbind(colnames(R), do.call(rbind, M5.LS))
View(M1.LS)
View(M3.LS)
View(M4.LS)
View(M5.LS)