# HM estimation of factor models for ind portfolios
#
# data
source(paste("http://raw.githubusercontent.com/erkind/",
             "eiv-coe/main/g-coe-data-eiv.R",
             sep = ""))
#
# packages
library(car)
library(lmtest)
library(sandwich)
library(skedastic)
#
# counters
nb.obs <- as.numeric(dim(R)[1])
nb.port<- as.numeric(dim(R)[2])
#
# function to return Z matrix
source(paste("http://raw.githubusercontent.com/erkind/",
             "eiv-coe/main/g-f-HM-Zmat.R",
             sep = ""))
#
# generating instrument matrices
Z.X1 <- generate.Z.mat(X.mat = X[,1],
                       Y.vec = R[,1],
                       n = nb.obs,
                       k = 1,
                       i = c(1, 4))
Z.X3 <- generate.Z.mat(X.mat = X[,1:3],
                       Y.vec = R[,1],
                       n = nb.obs,
                       k = 3,
                       i = c(1, 4))
Z.X4 <- generate.Z.mat(X.mat = X[,c(1:3,6)],
                       Y.vec = R[,1],
                       n = nb.obs,
                       k = 4,
                       i = c(1, 4))
Z.X5 <- generate.Z.mat(X.mat = X[,1:5],
                       Y.vec = R[,1],
                       n = nb.obs,
                       k = 5,
                       i = c(1, 4))
#
# check Z matices
head(Z.X1, n = 3); dim(Z.X1)
head(Z.X3, n = 3); dim(Z.X3)
head(Z.X4, n = 3); dim(Z.X4)
head(Z.X5, n = 3); dim(Z.X5)
#
# regs X on Z
# for each X[i] run X[i] ~ all Z's and save residuals
# CAPM
w.hat.X11 <- resid(lm(X[,1] ~ Z.X1))
# FF3F
w.hat.X31 <- resid(lm(X[,1] ~ Z.X3))
w.hat.X32 <- resid(lm(X[,2] ~ Z.X3))
w.hat.X33 <- resid(lm(X[,3] ~ Z.X3))
# FFC4F
w.hat.X41 <- resid(lm(X[,1] ~ Z.X4))
w.hat.X42 <- resid(lm(X[,2] ~ Z.X4))
w.hat.X43 <- resid(lm(X[,3] ~ Z.X4))
w.hat.X44 <- resid(lm(X[,6] ~ Z.X4))
# FF5F
w.hat.X51 <- resid(lm(X[,1] ~ Z.X5))
w.hat.X52 <- resid(lm(X[,2] ~ Z.X5))
w.hat.X53 <- resid(lm(X[,3] ~ Z.X5))
w.hat.X54 <- resid(lm(X[,4] ~ Z.X5))
w.hat.X55 <- resid(lm(X[,5] ~ Z.X5))
#
# time-series regressions
#
# CAPM
M1.HM <- list()
M1.HM.VCV <- list()
M1.HM.a   <- matrix(NA, nrow = nb.port, ncol = 1)
M1.HM.res <- matrix(NA, nrow = nb.obs, ncol = nb.port)
for(i in 1:nb.port){
   # estimate HM regression
   HM.fit <- lm(R[,i] ~ X[,1] + w.hat.X11)
   names(HM.fit$coefficients) <- c("a", "b",
                                   "psi.b")
   # coefficient covariance matrix
   HM.fit.VCV <- vcovHAC(HM.fit)
   # calculate IV statistics
   DWH <- as.numeric(
      linearHypothesis(HM.fit,
                       c("psi.b = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$Chisq[2])
   DWH.pval <- as.numeric(
      linearHypothesis(HM.fit,
                       c("psi.b = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$`Pr(>Chisq)`[2])
   OR <- nb.obs * summary(lm(resid(HM.fit) ~ Z.X1))$r.squared
   OR.pval <- 1 - pchisq(q = OR, df = dim(Z.X1)[2] - 1)
   # extract summary values
   Rsq   <- summary(HM.fit)$adj.r.squared
   sigma <- summary(HM.fit)$sigma
   a     <- summary(HM.fit)$coefficients[1,1]
   a.SE  <- as.numeric(diag(HM.fit.VCV)[1])^0.5
   a.t   <- a / a.SE
   b     <- summary(HM.fit)$coefficients[2,1]
   b.SE  <- as.numeric(diag(HM.fit.VCV)[2])^0.5
   b.t   <- b / b.SE
   psi.b   <- summary(HM.fit)$coefficients[3,1]
   psi.b.SE<- as.numeric(diag(HM.fit.VCV)[3])^0.5
   psi.b.t <- psi.b / psi.b.SE
   # excess returns time-series averages
   avg.R <- mean(R[,i])
   # factors time-series average
   avg.mkt <- mean(X[,1])
   # RF time-series average
   avg.RF  <- mean(RF)
   # coe (annualized)
   coe <- (avg.RF
           + b * avg.mkt
           ) * 12
   # list to export the results
   out.list <- list(
      Rsq = Rsq,
      sigma = sigma,
      DWH = DWH, DWH.pval = DWH.pval,
      J = OR, J.pval = OR.pval,
      a = a, a.SE = a.SE, a.t = a.t,
      b = b, b.SE = b.SE, b.t = b.t,
      psi.b = psi.b, psi.b.SE = psi.b.SE, psi.b.t = psi.b.t,
      avg.R = avg.R,
      avg.mkt = avg.mkt,
      avg.RF = avg.RF,
      coe = coe
   )
   # append results to current list
   M1.HM[[i]] <- out.list
   # coefficients covariance matrix
   M1.HM.VCV[[i]] <- vcovHAC(HM.fit)
   # intercepts
   M1.HM.a[i] <- a
   rownames(M1.HM.a) <- colnames(R)
   # residuals
   M1.HM.res[,i] <- HM.fit$residuals
   colnames(M1.HM.res) <- colnames(R)
}
#
# FF3F
M3.HM <- list()
M3.HM.VCV <- list()
M3.HM.a   <- matrix(NA, nrow = nb.port, ncol = 1)
M3.HM.res <- matrix(NA, nrow = nb.obs, ncol = nb.port)
for(i in 1:nb.port){
   # estimate HM regression
   HM.fit <- lm(R[,i] ~ X[,1] + X[,2] + X[,3] + 
                        w.hat.X31 + w.hat.X32 + w.hat.X33)
   names(HM.fit$coefficients) <- c("a", "b", "s", "h",
                                   "psi.b", "psi.s", "psi.h")
   # coefficient covariance matrix
   HM.fit.VCV <- vcovHAC(HM.fit)
   # calculate IV statistics
   DWH <- as.numeric(
      linearHypothesis(HM.fit,
                       c("psi.b = 0", "psi.s = 0", "psi.h = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$Chisq[2])
   DWH.pval <- as.numeric(
      linearHypothesis(HM.fit,
                       c("psi.b = 0", "psi.s = 0", "psi.h = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$`Pr(>Chisq)`[2])
   OR <- nb.obs * summary(lm(resid(HM.fit) ~ Z.X3))$r.squared
   OR.pval <- 1 - pchisq(q = OR, df = dim(Z.X3)[2] - 1)
   # extract summary values
   Rsq   <- summary(HM.fit)$adj.r.squared
   sigma <- summary(HM.fit)$sigma
   a     <- summary(HM.fit)$coefficients[1,1]
   a.SE  <- as.numeric(diag(HM.fit.VCV)[1])^0.5
   a.t   <- a / a.SE
   b     <- summary(HM.fit)$coefficients[2,1]
   b.SE  <- as.numeric(diag(HM.fit.VCV)[2])^0.5
   b.t   <- b / b.SE
   s     <- summary(HM.fit)$coefficients[3,1]
   s.SE  <- as.numeric(diag(HM.fit.VCV)[3])^0.5
   s.t   <- s / s.SE
   h     <- summary(HM.fit)$coefficients[4,1]
   h.SE  <- as.numeric(diag(HM.fit.VCV)[4])^0.5
   h.t   <- h / h.SE
   psi.b   <- summary(HM.fit)$coefficients[5,1]
   psi.b.SE<- as.numeric(diag(HM.fit.VCV)[5])^0.5
   psi.b.t <- psi.b / psi.b.SE
   psi.s   <- summary(HM.fit)$coefficients[6,1]
   psi.s.SE<- as.numeric(diag(HM.fit.VCV)[6])^0.5
   psi.s.t <- psi.s / psi.s.SE
   psi.h   <- summary(HM.fit)$coefficients[7,1]
   psi.h.SE<- as.numeric(diag(HM.fit.VCV)[7])^0.5
   psi.h.t <- psi.h / psi.h.SE
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
      DWH = DWH, DWH.pval = DWH.pval,
      J = OR, J.pval = OR.pval,
      a = a, a.SE = a.SE, a.t = a.t,
      b = b, b.SE = b.SE, b.t = b.t,
      s = s, s.SE = s.SE, s.t = s.t,
      h = h, h.SE = h.SE, h.t = h.t,
      psi.b = psi.b, psi.b.SE = psi.b.SE, psi.b.t = psi.b.t,
      psi.s = psi.s, psi.s.SE = psi.s.SE, psi.s.t = psi.s.t,
      psi.h = psi.h, psi.h.SE = psi.h.SE, psi.h.t = psi.h.t,
      avg.R = avg.R,
      avg.mkt = avg.mkt, avg.smb = avg.smb, avg.hml = avg.hml,
      avg.RF = avg.RF,
      coe = coe
   )
   # append results to current list
   M3.HM[[i]] <- out.list
   # coefficients covariance matrix
   M3.HM.VCV[[i]] <- vcovHAC(HM.fit)
   # intercepts
   M3.HM.a[i] <- a
   rownames(M3.HM.a) <- colnames(R)
   # residuals
   M3.HM.res[,i] <- HM.fit$residuals
   colnames(M3.HM.res) <- colnames(R)
}
#
# FFC4F
M4.HM <- list()
M4.HM.VCV <- list()
M4.HM.a   <- matrix(NA, nrow = nb.port, ncol = 1)
M4.HM.res <- matrix(NA, nrow = nb.obs, ncol = nb.port)
for(i in 1:nb.port){
   # estimate HM regression
   HM.fit <- lm(R[,i] ~ X[,1] + X[,2] + X[,3] + X[,6] + 
                        w.hat.X41 + w.hat.X42 + w.hat.X43 + w.hat.X44)
   names(HM.fit$coefficients) <- c("a", "b", "s", "h", "w",
                                   "psi.b", "psi.s", "psi.h", "psi.w")
   # coefficient covariance matrix
   HM.fit.VCV <- vcovHAC(HM.fit)
   # calculate IV statistics
   DWH <- as.numeric(
      linearHypothesis(HM.fit,
                       c("psi.b = 0", "psi.s = 0", "psi.h = 0",
                         "psi.w = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$Chisq[2])
   DWH.pval <- as.numeric(
      linearHypothesis(HM.fit,
                       c("psi.b = 0", "psi.s = 0", "psi.h = 0",
                         "psi.w = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$`Pr(>Chisq)`[2])
   OR <- nb.obs * summary(lm(resid(HM.fit) ~ Z.X4))$r.squared
   OR.pval <- 1 - pchisq(q = OR, df = dim(Z.X4)[2] - 1)
   # extract summary values
   Rsq   <- summary(HM.fit)$adj.r.squared
   sigma <- summary(HM.fit)$sigma
   a     <- summary(HM.fit)$coefficients[1,1]
   a.SE  <- as.numeric(diag(HM.fit.VCV)[1])^0.5
   a.t   <- a / a.SE
   b     <- summary(HM.fit)$coefficients[2,1]
   b.SE  <- as.numeric(diag(HM.fit.VCV)[2])^0.5
   b.t   <- b / b.SE
   s     <- summary(HM.fit)$coefficients[3,1]
   s.SE  <- as.numeric(diag(HM.fit.VCV)[3])^0.5
   s.t   <- s / s.SE
   h     <- summary(HM.fit)$coefficients[4,1]
   h.SE  <- as.numeric(diag(HM.fit.VCV)[4])^0.5
   h.t   <- h / h.SE
   w     <- summary(HM.fit)$coefficients[5,1]
   w.SE  <- as.numeric(diag(HM.fit.VCV)[5])^0.5
   w.t   <- w / w.SE
   psi.b   <- summary(HM.fit)$coefficients[6,1]
   psi.b.SE<- as.numeric(diag(HM.fit.VCV)[6])^0.5
   psi.b.t <- psi.b / psi.b.SE
   psi.s   <- summary(HM.fit)$coefficients[7,1]
   psi.s.SE<- as.numeric(diag(HM.fit.VCV)[7])^0.5
   psi.s.t <- psi.s / psi.s.SE
   psi.h   <- summary(HM.fit)$coefficients[8,1]
   psi.h.SE<- as.numeric(diag(HM.fit.VCV)[8])^0.5
   psi.h.t <- psi.h / psi.h.SE
   psi.w   <- summary(HM.fit)$coefficients[9,1]
   psi.w.SE<- as.numeric(diag(HM.fit.VCV)[9])^0.5
   psi.w.t <- psi.w / psi.w.SE
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
      DWH = DWH, DWH.pval = DWH.pval,
      J = OR, J.pval = OR.pval,
      a = a, a.SE = a.SE, a.t = a.t,
      b = b, b.SE = b.SE, b.t = b.t,
      s = s, s.SE = s.SE, s.t = s.t,
      h = h, h.SE = h.SE, h.t = h.t,
      w = w, w.SE = w.SE, w.t = w.t,
      psi.b = psi.b, psi.b.SE = psi.b.SE, psi.b.t = psi.b.t,
      psi.s = psi.s, psi.s.SE = psi.s.SE, psi.s.t = psi.s.t,
      psi.h = psi.h, psi.h.SE = psi.h.SE, psi.h.t = psi.h.t,
      psi.w = psi.w, psi.w.SE = psi.w.SE, psi.w.t = psi.w.t,
      avg.R = avg.R,
      avg.mkt = avg.mkt, avg.smb = avg.smb, avg.hml = avg.hml,
      avg.wml = avg.wml,
      avg.RF = avg.RF,
      coe = coe
   )
   # append results to current list
   M4.HM[[i]] <- out.list
   # coefficients covariance matrix
   M4.HM.VCV[[i]] <- vcovHAC(HM.fit)
   # intercepts
   M4.HM.a[i] <- a
   rownames(M4.HM.a) <- colnames(R)
   # residuals
   M4.HM.res[,i] <- HM.fit$residuals
   colnames(M4.HM.res) <- colnames(R)
}
#
# FF5F
M5.HM <- list()
M5.HM.VCV <- list()
M5.HM.a   <- matrix(NA, nrow = nb.port, ncol = 1)
M5.HM.res <- matrix(NA, nrow = nb.obs, ncol = nb.port)
for(i in 1:nb.port){
   # estimate HM regression
   HM.fit <- lm(R[,i] ~ X[,1] + X[,2] + X[,3] + X[,4] + X[,5] +
                        w.hat.X51 + w.hat.X52 + w.hat.X53 +
                        w.hat.X54 + w.hat.X55)
   names(HM.fit$coefficients) <- c("a", "b", "s", "h", "r", "c",
                                   "psi.b", "psi.s", "psi.h",
                                   "psi.r", "psi.c")
   # coefficient covariance matrix
   HM.fit.VCV <- vcovHAC(HM.fit)
   # calculate IV statistics
   DWH <- as.numeric(
      linearHypothesis(HM.fit,
                       c("psi.b = 0", "psi.s = 0", "psi.h = 0",
                         "psi.r = 0", "psi.c = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$Chisq[2])
   DWH.pval <- as.numeric(
      linearHypothesis(HM.fit,
                       c("psi.b = 0", "psi.s = 0", "psi.h = 0",
                         "psi.r = 0", "psi.c = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$`Pr(>Chisq)`[2])
   OR <- nb.obs * summary(lm(resid(HM.fit) ~ Z.X5))$r.squared
   OR.pval <- 1 - pchisq(q = OR, df = dim(Z.X5)[2] - 1)
   # extract summary values
   Rsq   <- summary(HM.fit)$adj.r.squared
   sigma <- summary(HM.fit)$sigma
   a     <- summary(HM.fit)$coefficients[1,1]
   a.SE  <- as.numeric(diag(HM.fit.VCV)[1])^0.5
   a.t   <- a / a.SE
   b     <- summary(HM.fit)$coefficients[2,1]
   b.SE  <- as.numeric(diag(HM.fit.VCV)[2])^0.5
   b.t   <- b / b.SE
   s     <- summary(HM.fit)$coefficients[3,1]
   s.SE  <- as.numeric(diag(HM.fit.VCV)[3])^0.5
   s.t   <- s / s.SE
   h     <- summary(HM.fit)$coefficients[4,1]
   h.SE  <- as.numeric(diag(HM.fit.VCV)[4])^0.5
   h.t   <- h / h.SE
   r     <- summary(HM.fit)$coefficients[5,1]
   r.SE  <- as.numeric(diag(HM.fit.VCV)[5])^0.5
   r.t   <- r / r.SE
   c     <- summary(HM.fit)$coefficients[6,1]
   c.SE  <- as.numeric(diag(HM.fit.VCV)[6])^0.5
   c.t   <- c / c.SE
   psi.b   <- summary(HM.fit)$coefficients[7,1]
   psi.b.SE<- as.numeric(diag(HM.fit.VCV)[7])^0.5
   psi.b.t <- psi.b / psi.b.SE
   psi.s   <- summary(HM.fit)$coefficients[8,1]
   psi.s.SE<- as.numeric(diag(HM.fit.VCV)[8])^0.5
   psi.s.t <- psi.s / psi.s.SE
   psi.h   <- summary(HM.fit)$coefficients[9,1]
   psi.h.SE<- as.numeric(diag(HM.fit.VCV)[9])^0.5
   psi.h.t <- psi.h / psi.h.SE
   psi.r   <- summary(HM.fit)$coefficients[10,1]
   psi.r.SE<- as.numeric(diag(HM.fit.VCV)[10])^0.5
   psi.r.t <- psi.r / psi.r.SE
   psi.c   <- summary(HM.fit)$coefficients[11,1]
   psi.c.SE<- as.numeric(diag(HM.fit.VCV)[11])^0.5
   psi.c.t <- psi.c / psi.c.SE
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
      DWH = DWH, DWH.pval = DWH.pval,
      J = OR, J.pval = OR.pval,
      a = a, a.SE = a.SE, a.t = a.t,
      b = b, b.SE = b.SE, b.t = b.t,
      s = s, s.SE = s.SE, s.t = s.t,
      h = h, h.SE = h.SE, h.t = h.t,
      r = r, r.SE = r.SE, r.t = r.t,
      c = c, c.SE = c.SE, c.t = c.t,
      psi.b = psi.b, psi.b.SE = psi.b.SE, psi.b.t = psi.b.t,
      psi.s = psi.s, psi.s.SE = psi.s.SE, psi.s.t = psi.s.t,
      psi.h = psi.h, psi.h.SE = psi.h.SE, psi.h.t = psi.h.t,
      psi.r = psi.r, psi.r.SE = psi.r.SE, psi.r.t = psi.r.t,
      psi.c = psi.c, psi.c.SE = psi.c.SE, psi.c.t = psi.c.t,
      avg.R = avg.R,
      avg.mkt = avg.mkt, avg.smb = avg.smb, avg.hml = avg.hml,
      avg.rmw = avg.rmw, avg.cma = avg.cma,
      avg.RF = avg.RF,
      coe = coe
   )
   # append results to current list
   M5.HM[[i]] <- out.list
   # coefficients covariance matrix
   M5.HM.VCV[[i]] <- vcovHAC(HM.fit)
   # intercepts
   M5.HM.a[i] <- a
   rownames(M5.HM.a) <- colnames(R)
   # residuals
   M5.HM.res[,i] <- HM.fit$residuals
   colnames(M5.HM.res) <- colnames(R)
}
#
# collect output
M1.HM <- cbind(colnames(R), do.call(rbind, M1.HM))
M3.HM <- cbind(colnames(R), do.call(rbind, M3.HM))
M4.HM <- cbind(colnames(R), do.call(rbind, M4.HM))
M5.HM <- cbind(colnames(R), do.call(rbind, M5.HM))
# 
# visualize output
View(M1.HM)
View(M3.HM)
View(M4.HM)
View(M5.HM)
#
# testing for weak instruments
# run x[i] ~ Z.mat, check F-stat > 10 (Stock, Yogo, 2005) or F-stat > 24
# (Olea, Pflueger, 2013) -> If yes, Z do not suffer from weak ins problem
# CAPM
summary(lm(X[,1] ~ Z.X1))$fstatistic
# FF3F
summary(lm(X[,1] ~ Z.X3))$fstatistic
summary(lm(X[,2] ~ Z.X3))$fstatistic
summary(lm(X[,3] ~ Z.X3))$fstatistic
# FFC4F
summary(lm(X[,1] ~ Z.X4))$fstatistic
summary(lm(X[,2] ~ Z.X4))$fstatistic
summary(lm(X[,3] ~ Z.X4))$fstatistic
summary(lm(X[,4] ~ Z.X4))$fstatistic
# FF5F
summary(lm(X[,1] ~ Z.X5))$fstatistic
summary(lm(X[,2] ~ Z.X5))$fstatistic
summary(lm(X[,3] ~ Z.X5))$fstatistic
summary(lm(X[,4] ~ Z.X5))$fstatistic
summary(lm(X[,5] ~ Z.X5))$fstatistic
#
# testing for exogenous instruments
# CAPM
ex.iv.M1 <- lapply(1:dim(R)[2],
                   function(i)
                      lm(R[,i] ~ X[,1] +
                                 w.hat.X11
                         )
                   )
ex.iv.M1 <- lapply(1:dim(R)[2],
                   function(i)
                      lm(ex.iv.M1[[i]]$resid ~ -1 + Z.X1
                         )
                   )
ex.iv.M1 <- lapply(ex.iv.M1, summary)
# FF3F
ex.iv.M3 <- lapply(1:dim(R)[2],
                   function(i)
                      lm(R[,i] ~ X[,1] + X[,2] + X[,3] + 
                                 w.hat.X31 +
                                 w.hat.X32 +
                                 w.hat.X33
                         )
                   )
ex.iv.M3 <- lapply(1:dim(R)[2],
                   function(i)
                      lm(ex.iv.M3[[i]]$resid ~ -1 + Z.X3
                         )
                   )
ex.iv.M3 <- lapply(ex.iv.M3, summary)
# FFC4F
ex.iv.M4 <- lapply(1:dim(R)[2],
                   function(i)
                      lm(R[,i] ~ X[,1] + X[,2] + X[,3] + X[,6] + 
                            w.hat.X41 +
                            w.hat.X42 +
                            w.hat.X43 +
                            w.hat.X44
                         )
                   )
ex.iv.M4 <- lapply(1:dim(R)[2],
                   function(i)
                      lm(ex.iv.M4[[i]]$resid ~ -1 + Z.X4
                         )
                   )
ex.iv.M4 <- lapply(ex.iv.M4, summary)
# FF5F
ex.iv.M5 <- lapply(1:dim(R)[2],
                   function(i)
                      lm(R[,i] ~ X[,1] + X[,2] + X[,3] + X[,4] + X[,5] +
                            w.hat.X51 +
                            w.hat.X52 +
                            w.hat.X53 +
                            w.hat.X54 +
                            w.hat.X55
                         )
                   )
ex.iv.M5 <- lapply(1:dim(R)[2],
                   function(i)
                      lm(ex.iv.M5[[i]]$resid ~ -1 + Z.X5
                         )
                   )
ex.iv.M5 <- lapply(ex.iv.M5, summary)
# generate data frames
exo.test.M1 <- data.frame(
   industry = character(0),
   R.sq     = numeric(0),
   F.stat   = numeric(0),
   avg.p.val= numeric(0)
   )
for(i in 1:length(ex.iv.M1)){
   model.summary <- ex.iv.M1[[i]]
   industry      <- colnames(R)[i]
   R.sq          <- model.summary$r.squared
   F.stat        <- model.summary$fstatistic[1]
   avg.p.val     <- mean(model.summary$coefficients[,4])
   #
   exo.test.M1 <- rbind(exo.test.M1,
                        data.frame(
                           industry = industry,
                           R.sq     = R.sq,
                           F.stat   = F.stat,
                           avg.p.val= avg.p.val)
                        )
}
exo.test.M3 <- data.frame(
   industry = character(0),
   R.sq     = numeric(0),
   F.stat   = numeric(0),
   avg.p.val= numeric(0)
)
for(i in 1:length(ex.iv.M3)){
   model.summary <- ex.iv.M3[[i]]
   industry      <- colnames(R)[i]
   R.sq          <- model.summary$r.squared
   F.stat        <- model.summary$fstatistic[1]
   avg.p.val     <- mean(model.summary$coefficients[,4])
   #
   exo.test.M3 <- rbind(exo.test.M3,
                        data.frame(
                           industry = industry,
                           R.sq     = R.sq,
                           F.stat   = F.stat,
                           avg.p.val= avg.p.val)
   )
}
exo.test.M4 <- data.frame(
   industry = character(0),
   R.sq     = numeric(0),
   F.stat   = numeric(0),
   avg.p.val= numeric(0)
)
for(i in 1:length(ex.iv.M4)){
   model.summary <- ex.iv.M4[[i]]
   industry      <- colnames(R)[i]
   R.sq          <- model.summary$r.squared
   F.stat        <- model.summary$fstatistic[1]
   avg.p.val     <- mean(model.summary$coefficients[,4])
   #
   exo.test.M4 <- rbind(exo.test.M4,
                        data.frame(
                           industry = industry,
                           R.sq     = R.sq,
                           F.stat   = F.stat,
                           avg.p.val= avg.p.val)
   )
}
exo.test.M5 <- data.frame(
   industry = character(0),
   R.sq     = numeric(0),
   F.stat   = numeric(0),
   avg.p.val= numeric(0)
)
for(i in 1:length(ex.iv.M5)){
   model.summary <- ex.iv.M5[[i]]
   industry      <- colnames(R)[i]
   R.sq          <- model.summary$r.squared
   F.stat        <- model.summary$fstatistic[1]
   avg.p.val     <- mean(model.summary$coefficients[,4])
   #
   exo.test.M5 <- rbind(exo.test.M5,
                        data.frame(
                           industry = industry,
                           R.sq     = R.sq,
                           F.stat   = F.stat,
                           avg.p.val= avg.p.val)
   )
}
# merge data frames 
library(dplyr)
# add column for model specification
exo.test.M1 <- exo.test.M1 %>% mutate(Specification = "M1")
exo.test.M3 <- exo.test.M3 %>% mutate(Specification = "M3")
exo.test.M4 <- exo.test.M4 %>% mutate(Specification = "M4")
exo.test.M5 <- exo.test.M5 %>% mutate(Specification = "M5")
exo.test <- bind_rows(exo.test.M1, exo.test.M3, exo.test.M4, exo.test.M5)
#
View(exo.test)



