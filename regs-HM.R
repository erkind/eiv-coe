dim(data); colnames(data)[c(1:5,60:64)]
dim(R); colnames(R)
dim(F1); colnames(F1)
dim(F3); colnames(F3)
dim(F4); colnames(F4)
dim(F5); colnames(F5)

# packages --------------------------------------------------------------------
library(car)
library(lmtest)
library(sandwich)
library(skedastic)

# lists -----------------------------------------------------------------------
# counters
nb.obs <- as.numeric(dim(R)[1])
nb.port<- as.numeric(dim(R)[2])
# CAPM
M1.HM    <- list()
M1.HM.VCV<- list()
M1.HM.res<- matrix(NA, nrow = nb.obs, ncol = nb.port)
M1.HM.a  <- matrix(NA, nrow = nb.port, ncol = 1)
M1.HM.DWH<- matrix(NA, nrow = nb.port, ncol = 1)
M1.HM.OR <- matrix(NA, nrow = nb.port, ncol = 1)
M1.HM.DWH.pval<- matrix(NA, nrow = nb.port, ncol = 1)
M1.HM.OR.pval <- matrix(NA, nrow = nb.port, ncol = 1)
# FF3F
M3.HM    <- list()
M3.HM.VCV<- list()
M3.HM.res<- matrix(NA, nrow = nb.obs, ncol = nb.port)
M3.HM.a  <- matrix(NA, nrow = nb.port, ncol = 1)
M3.HM.DWH<- matrix(NA, nrow = nb.port, ncol = 1)
M3.HM.OR <- matrix(NA, nrow = nb.port, ncol = 1)
M3.HM.DWH.pval<- matrix(NA, nrow = nb.port, ncol = 1)
M3.HM.OR.pval <- matrix(NA, nrow = nb.port, ncol = 1)
# FFC4F
M4.HM    <- list()
M4.HM.VCV<- list()
M4.HM.res<- matrix(NA, nrow = nb.obs, ncol = nb.port)
M4.HM.a  <- matrix(NA, nrow = nb.port, ncol = 1)
M4.HM.DWH<- matrix(NA, nrow = nb.port, ncol = 1)
M4.HM.OR <- matrix(NA, nrow = nb.port, ncol = 1)
M4.HM.DWH.pval<- matrix(NA, nrow = nb.port, ncol = 1)
M4.HM.OR.pval <- matrix(NA, nrow = nb.port, ncol = 1)
# FF5F
M5.HM    <- list()
M5.HM.VCV<- list()
M5.HM.res<- matrix(NA, nrow = nb.obs, ncol = nb.port)
M5.HM.a  <- matrix(NA, nrow = nb.port, ncol = 1)
M5.HM.DWH<- matrix(NA, nrow = nb.port, ncol = 1)
M5.HM.OR <- matrix(NA, nrow = nb.port, ncol = 1)
M5.HM.DWH.pval<- matrix(NA, nrow = nb.port, ncol = 1)
M5.HM.OR.pval <- matrix(NA, nrow = nb.port, ncol = 1)

# Z matrix --------------------------------------------------------------------
# functions
source(paste("C:/Users/68596/Documents/docs/",
             "projects/R-scripts/HME-Zmat.R",
             sep = ""))
source(paste("C:/Users/68596/Documents/docs/",
             "projects/R-scripts/HME-Zmat-labels.R",
             sep = ""))
# instrument matrix
Z.F1 <- generate_Zmat(n = nb.obs, k = 1,
                      i = c(1, 4), X.mat = F1)
Z.F3 <- generate_Zmat(n = nb.obs, k = 3,
                      i = c(1, 4), X.mat = F3)
Z.F4 <- generate_Zmat(n = nb.obs, k = 4,
                      i = c(1, 4), X.mat = F4)
Z.F5 <- generate_Zmat(n = nb.obs, k = 5,
                      i = c(1, 4), X.mat = F5)
# instrument labels
colnames(Z.F1) <- generate_Zmat_labels(i = c(1, 4), k = 1)
colnames(Z.F3) <- generate_Zmat_labels(i = c(1, 4), k = 3)
colnames(Z.F4) <- generate_Zmat_labels(i = c(1, 4), k = 4)
colnames(Z.F5) <- generate_Zmat_labels(i = c(1, 4), k = 5)
# check the number of columns
head(Z.F1, n = 3); dim(Z.F1)
head(Z.F3, n = 3); dim(Z.F3)
head(Z.F4, n = 3); dim(Z.F4)
head(Z.F5, n = 3); dim(Z.F5)

# reg X on Z ------------------------------------------------------------------
# for each X[i] run X[i] ~ all Z's and save residuals
# CAPM
w.hat.F1.MKT <- resid(lm(F1[,1] ~ Z.F1))
# FF3F
w.hat.F3.MKT <- resid(lm(F3[,1] ~ Z.F3))
w.hat.F3.SMB <- resid(lm(F3[,2] ~ Z.F3))
w.hat.F3.HML <- resid(lm(F3[,3] ~ Z.F3))
# FFC4F
w.hat.F4.MKT <- resid(lm(F4[,1] ~ Z.F4))
w.hat.F4.SMB <- resid(lm(F4[,2] ~ Z.F4))
w.hat.F4.HML <- resid(lm(F4[,3] ~ Z.F4))
w.hat.F4.WML <- resid(lm(F4[,4] ~ Z.F4))
# FF5F
w.hat.F5.MKT <- resid(lm(F5[,1] ~ Z.F5))
w.hat.F5.SMB <- resid(lm(F5[,2] ~ Z.F5))
w.hat.F5.HML <- resid(lm(F5[,3] ~ Z.F5))
w.hat.F5.RMW <- resid(lm(F5[,4] ~ Z.F5))
w.hat.F5.CMA <- resid(lm(F5[,5] ~ Z.F5))

# estimations -----------------------------------------------------------------
# CAPM
for(i in 1:nb.port){
   #
   M1.HM[[i]] <- lm(R[,i] ~ F1[,1] + w.hat.F1.MKT)
   names(M1.HM[[i]]$coefficients) <- c("a", "b", "psi.b")
   M1.HM.VCV[[i]] <- vcovHAC(M1.HM[[i]], type = "HC3")
   M1.HM.res[,i] <- M1.HM[[i]]$residuals
   M1.HM.a[i] <- summary(M1.HM[[i]])$coefficients[1,1]
   #
   M1.HM.DWH[i] <- as.numeric(
      linearHypothesis(M1.HM[[i]],
                       c("psi.b = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$Chisq[2])
   M1.HM.DWH.pval[i] <- as.numeric(
      linearHypothesis(M1.HM[[i]],
                       c("psi.b = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$`Pr(>Chisq)`[2])
   M1.HM.OR[i] <- nb.obs * summary(lm(resid(M1.HM[[i]]) ~ Z.F1))$r.squared
   M1.HM.OR.pval[i] <- 1 - pchisq(q = M1.HM.OR[i], df = dim(Z.F1)[2] - 1)
}
# FF3F
for(i in 1:nb.port){
   #
   M3.HM[[i]] <- lm(R[,i] ~ F3[,1] + F3[,2] + F3[,3] + 
                       w.hat.F3.MKT + w.hat.F3.SMB + w.hat.F3.HML)
   names(M3.HM[[i]]$coefficients) <- c("a", "b", "s", "h",
                                       "psi.b", "psi.s", "psi.h")
   M3.HM.VCV[[i]] <- vcovHAC(M3.HM[[i]], type = "HC3")
   M3.HM.res[,i] <- M3.HM[[i]]$residuals
   M3.HM.a[i] <- summary(M3.HM[[i]])$coefficients[1,1]
   #
   M3.HM.DWH[i] <- as.numeric(
      linearHypothesis(M3.HM[[i]],
                       c("psi.b = 0", "psi.s = 0", "psi.h = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$Chisq[2])
   M3.HM.DWH.pval[i] <- as.numeric(
      linearHypothesis(M3.HM[[i]],
                       c("psi.b = 0", "psi.s = 0", "psi.h = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$`Pr(>Chisq)`[2])
   M3.HM.OR[i] <- nb.obs * summary(lm(resid(M3.HM[[i]]) ~ Z.F3))$r.squared
   M3.HM.OR.pval[i] <- 1 - pchisq(q = M3.HM.OR[i], df = dim(Z.F3)[2] - 3)
}
# FF4F
for(i in 1:nb.port){
   #
   M4.HM[[i]] <- lm(R[,i] ~ F4[,1] + F4[,2] + F4[,3] + F4[,4] +  
                       w.hat.F3.MKT + w.hat.F3.SMB +
                       w.hat.F3.HML + w.hat.F4.WML)
   names(M4.HM[[i]]$coefficients) <- c("a", "b", "s", "h", "w",
                                       "psi.b", "psi.s", "psi.h", "psi.w")
   M4.HM.VCV[[i]] <- vcovHAC(M4.HM[[i]], type = "HC3")
   M4.HM.res[,i] <- M4.HM[[i]]$residuals
   M4.HM.a[i] <- summary(M4.HM[[i]])$coefficients[1,1]
   #
   M4.HM.DWH[i] <- as.numeric(
      linearHypothesis(M4.HM[[i]],
                       c("psi.b = 0", "psi.s = 0", "psi.h = 0",
                         "psi.w = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$Chisq[2])
   M4.HM.DWH.pval[i] <- as.numeric(
      linearHypothesis(M4.HM[[i]],
                       c("psi.b = 0", "psi.s = 0", "psi.h = 0",
                         "psi.w = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$`Pr(>Chisq)`[2])
   M4.HM.OR[i] <- nb.obs * summary(lm(resid(M4.HM[[i]]) ~ Z.F4))$r.squared
   M4.HM.OR.pval[i] <- 1 - pchisq(q = M4.HM.OR[i], df = dim(Z.F4)[2] - 4)
}
# FF5F
for(i in 1:nb.port){
   #
   M5.HM[[i]] <- lm(R[,i] ~ F5[,1] + F5[,2] + F5[,3] +
                    F5[,4] + F5[,5] +
                    w.hat.F5.MKT + w.hat.F5.SMB + w.hat.F5.HML + 
                    w.hat.F5.RMW + w.hat.F5.CMA)
   names(M5.HM[[i]]$coefficients) <- c("a", "b", "s", "h", "r", "c",
                                       "psi.b", "psi.s", "psi.h",
                                       "psi.r", "psi.c")
   M5.HM.VCV[[i]] <- vcovHAC(M5.HM[[i]], type = "HC3")
   M5.HM.res[,i] <- M5.HM[[i]]$residuals
   M5.HM.a[i] <- summary(M5.HM[[i]])$coefficients[1,1]
   #
   M5.HM.DWH[i] <- as.numeric(
      linearHypothesis(M5.HM[[i]],
                       c("psi.b = 0", "psi.s = 0", "psi.h = 0",
                         "psi.r = 0", "psi.c = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$Chisq[2])
   M5.HM.DWH.pval[i] <- as.numeric(
      linearHypothesis(M5.HM[[i]],
                       c("psi.b = 0", "psi.s = 0", "psi.h = 0",
                         "psi.r = 0", "psi.c = 0"),
                       test = "Chisq",
                       white.adjust = c("hc0"))$`Pr(>Chisq)`[2])
   M5.HM.OR[i] <- nb.obs * summary(lm(resid(M5.HM[[i]]) ~ Z.F5))$r.squared
   M5.HM.OR.pval[i] <- 1 - pchisq(q = M5.HM.OR[i], df = dim(Z.F5)[2] - 5)
}

# results ---------------------------------------------------------------------
# empty tables
M1.HM.output <- matrix(NA, nrow = nb.port, ncol = 15)
M3.HM.output <- matrix(NA, nrow = nb.port, ncol = 27)
M4.HM.output <- matrix(NA, nrow = nb.port, ncol = 33)
M5.HM.output <- matrix(NA, nrow = nb.port, ncol = 39)
# CAPM
for(i in 1:nb.port){
   # industry
   M1.HM.output[i,1] <- colnames(R)[i]
   # adj R-sq
   M1.HM.output[i,2]<- summary(M1.HM[[i]])$adj.r.squared
   # DWH stat, p-value
   M1.HM.output[i,3]<- M1.HM.DWH[i]
   M1.HM.output[i,4]<- M1.HM.DWH.pval[i]
   # J-stat, p-value
   M1.HM.output[i,5]<- M1.HM.OR[i]
   M1.HM.output[i,6]<- M1.HM.OR.pval[i]
   # alpha, SE, t-stat
   M1.HM.output[i,7] <- summary(M1.HM[[i]])$coefficients[1,1]
   M1.HM.output[i,8] <- summary(M1.HM[[i]])$coefficients[1,2]
   M1.HM.output[i,9] <- summary(M1.HM[[i]])$coefficients[1,3]
   # beta, SE, t-stat
   M1.HM.output[i,10]<- summary(M1.HM[[i]])$coefficients[2,1]
   M1.HM.output[i,11]<- summary(M1.HM[[i]])$coefficients[2,2]
   M1.HM.output[i,12]<- summary(M1.HM[[i]])$coefficients[2,3]
   # psi.beta, SE, t-stat
   M1.HM.output[i,13]<- summary(M1.HM[[i]])$coefficients[3,1]
   M1.HM.output[i,14]<- summary(M1.HM[[i]])$coefficients[3,2]
   M1.HM.output[i,15]<- summary(M1.HM[[i]])$coefficients[3,3]
}
colnames(M1.HM.output) <- c("industry", "Rsq",
                            "DWH", "DWH-p-val",
                            "J", "J-p-val",
                            "a", "SE(a)", "t(a)",
                            "b", "SE(b)", "t(b)",
                            "psi.b", "SE(psi.b)", "t(psi.b)"
                            )
# FF3F
for(i in 1:nb.port){
   # industry
   M3.HM.output[i,1] <- colnames(R)[i]
   # adj R-sq
   M3.HM.output[i,2]<- summary(M3.HM[[i]])$adj.r.squared
   # DWH stat, p-value
   M3.HM.output[i,3]<- M3.HM.DWH[i]
   M3.HM.output[i,4]<- M3.HM.DWH.pval[i]
   # J-stat, p-value
   M3.HM.output[i,5]<- M3.HM.OR[i]
   M3.HM.output[i,6]<- M3.HM.OR.pval[i]
   # alpha, SE, t-stat
   M3.HM.output[i,7] <- summary(M3.HM[[i]])$coefficients[1,1]
   M3.HM.output[i,8] <- summary(M3.HM[[i]])$coefficients[1,2]
   M3.HM.output[i,9] <- summary(M3.HM[[i]])$coefficients[1,3]
   # beta, SE, t-stat
   M3.HM.output[i,10]<- summary(M3.HM[[i]])$coefficients[2,1]
   M3.HM.output[i,11]<- summary(M3.HM[[i]])$coefficients[2,2]
   M3.HM.output[i,12]<- summary(M3.HM[[i]])$coefficients[2,3]
   # smb, SE, t-stat
   M3.HM.output[i,13]<- summary(M3.HM[[i]])$coefficients[3,1]
   M3.HM.output[i,14]<- summary(M3.HM[[i]])$coefficients[3,2]
   M3.HM.output[i,15]<- summary(M3.HM[[i]])$coefficients[3,3]
   # hml, SE, t-stat
   M3.HM.output[i,16]<- summary(M3.HM[[i]])$coefficients[4,1]
   M3.HM.output[i,17]<- summary(M3.HM[[i]])$coefficients[4,2]
   M3.HM.output[i,18]<- summary(M3.HM[[i]])$coefficients[4,3]
   # psi.beta, SE, t-stat
   M3.HM.output[i,19]<- summary(M3.HM[[i]])$coefficients[5,1]
   M3.HM.output[i,20]<- summary(M3.HM[[i]])$coefficients[5,2]
   M3.HM.output[i,21]<- summary(M3.HM[[i]])$coefficients[5,3]
   # psi.smb, SE, t-stat
   M3.HM.output[i,22]<- summary(M3.HM[[i]])$coefficients[6,1]
   M3.HM.output[i,23]<- summary(M3.HM[[i]])$coefficients[6,2]
   M3.HM.output[i,24]<- summary(M3.HM[[i]])$coefficients[6,3]
   # psi.hml, SE, t-stat
   M3.HM.output[i,25]<- summary(M3.HM[[i]])$coefficients[7,1]
   M3.HM.output[i,26]<- summary(M3.HM[[i]])$coefficients[7,2]
   M3.HM.output[i,27]<- summary(M3.HM[[i]])$coefficients[7,3]
}
colnames(M3.HM.output) <- c("industry", "Rsq",
                            "DWH", "DWH-p-val",
                            "J", "J-p-val",
                            "a", "SE(a)", "t(a)",
                            "b", "SE(b)", "t(b)",
                            "s", "SE(s)", "t(s)",
                            "h", "SE(h)", "t(h)",
                            "psi.b", "SE(psi.b)", "t(psi.b)",
                            "psi.s", "SE(psi.s)", "t(psi.s)",
                            "psi.h", "SE(psi.h)", "t(psi.h)"
                            )
# FFC4F
for(i in 1:nb.port){
   # industry
   M4.HM.output[i,1] <- colnames(R)[i]
   # adj R-sq
   M4.HM.output[i,2]<- summary(M4.HM[[i]])$adj.r.squared
   # DWH stat, p-value
   M4.HM.output[i,3]<- M4.HM.DWH[i]
   M4.HM.output[i,4]<- M4.HM.DWH.pval[i]
   # J-stat, p-value
   M4.HM.output[i,5]<- M4.HM.OR[i]
   M4.HM.output[i,6]<- M4.HM.OR.pval[i]
   # alpha, SE, t-stat
   M4.HM.output[i,7] <- summary(M4.HM[[i]])$coefficients[1,1]
   M4.HM.output[i,8] <- summary(M4.HM[[i]])$coefficients[1,2]
   M4.HM.output[i,9] <- summary(M4.HM[[i]])$coefficients[1,3]
   # beta, SE, t-stat
   M4.HM.output[i,10]<- summary(M4.HM[[i]])$coefficients[2,1]
   M4.HM.output[i,11]<- summary(M4.HM[[i]])$coefficients[2,2]
   M4.HM.output[i,12]<- summary(M4.HM[[i]])$coefficients[2,3]
   # smb, SE, t-stat
   M4.HM.output[i,13]<- summary(M4.HM[[i]])$coefficients[3,1]
   M4.HM.output[i,14]<- summary(M4.HM[[i]])$coefficients[3,2]
   M4.HM.output[i,15]<- summary(M4.HM[[i]])$coefficients[3,3]
   # hml, SE, t-stat
   M4.HM.output[i,16]<- summary(M4.HM[[i]])$coefficients[4,1]
   M4.HM.output[i,17]<- summary(M4.HM[[i]])$coefficients[4,2]
   M4.HM.output[i,18]<- summary(M4.HM[[i]])$coefficients[4,3]
   # wml, SE, t-stat
   M4.HM.output[i,19]<- summary(M4.HM[[i]])$coefficients[5,1]
   M4.HM.output[i,20]<- summary(M4.HM[[i]])$coefficients[5,2]
   M4.HM.output[i,21]<- summary(M4.HM[[i]])$coefficients[5,3]
   # psi.beta, SE, t-stat
   M4.HM.output[i,22]<- summary(M4.HM[[i]])$coefficients[6,1]
   M4.HM.output[i,23]<- summary(M4.HM[[i]])$coefficients[6,2]
   M4.HM.output[i,24]<- summary(M4.HM[[i]])$coefficients[6,3]
   # psi.smb, SE, t-stat
   M4.HM.output[i,25]<- summary(M4.HM[[i]])$coefficients[7,1]
   M4.HM.output[i,26]<- summary(M4.HM[[i]])$coefficients[7,2]
   M4.HM.output[i,27]<- summary(M4.HM[[i]])$coefficients[7,3]
   # psi.hml, SE, t-stat
   M4.HM.output[i,28]<- summary(M4.HM[[i]])$coefficients[8,1]
   M4.HM.output[i,29]<- summary(M4.HM[[i]])$coefficients[8,2]
   M4.HM.output[i,30]<- summary(M4.HM[[i]])$coefficients[8,3]
   # psi.wml, SE, t-stat
   M4.HM.output[i,31]<- summary(M4.HM[[i]])$coefficients[9,1]
   M4.HM.output[i,32]<- summary(M4.HM[[i]])$coefficients[9,2]
   M4.HM.output[i,33]<- summary(M4.HM[[i]])$coefficients[9,3]
}
colnames(M4.HM.output) <- c("industry", "Rsq",
                            "DWH", "DWH-p-val",
                            "J", "J-p-val",
                            "a", "SE(a)", "t(a)",
                            "b", "SE(b)", "t(b)",
                            "s", "SE(s)", "t(s)",
                            "h", "SE(h)", "t(h)",
                            "w", "SE(w)", "t(w)",
                            "psi.b", "SE(psi.b)", "t(psi.b)",
                            "psi.s", "SE(psi.s)", "t(psi.s)",
                            "psi.h", "SE(psi.h)", "t(psi.h)",
                            "psi.w", "SE(psi.w)", "t(psi.w)"
                            )
# FF5F
for(i in 1:nb.port){
   # industry
   M5.HM.output[i,1] <- colnames(R)[i]
   # adj R-sq
   M5.HM.output[i,2]<- summary(M5.HM[[i]])$adj.r.squared
   # DWH stat, p-value
   M5.HM.output[i,3]<- M5.HM.DWH[i]
   M5.HM.output[i,4]<- M5.HM.DWH.pval[i]
   # J-stat, p-value
   M5.HM.output[i,5]<- M5.HM.OR[i]
   M5.HM.output[i,6]<- M5.HM.OR.pval[i]
   # alpha, SE, t-stat
   M5.HM.output[i,7] <- summary(M5.HM[[i]])$coefficients[1,1]
   M5.HM.output[i,8] <- summary(M5.HM[[i]])$coefficients[1,2]
   M5.HM.output[i,9] <- summary(M5.HM[[i]])$coefficients[1,3]
   # beta, SE, t-stat
   M5.HM.output[i,10]<- summary(M5.HM[[i]])$coefficients[2,1]
   M5.HM.output[i,11]<- summary(M5.HM[[i]])$coefficients[2,2]
   M5.HM.output[i,12]<- summary(M5.HM[[i]])$coefficients[2,3]
   # smb, SE, t-stat
   M5.HM.output[i,13]<- summary(M5.HM[[i]])$coefficients[3,1]
   M5.HM.output[i,14]<- summary(M5.HM[[i]])$coefficients[3,2]
   M5.HM.output[i,15]<- summary(M5.HM[[i]])$coefficients[3,3]
   # hml, SE, t-stat
   M5.HM.output[i,16]<- summary(M5.HM[[i]])$coefficients[4,1]
   M5.HM.output[i,17]<- summary(M5.HM[[i]])$coefficients[4,2]
   M5.HM.output[i,18]<- summary(M5.HM[[i]])$coefficients[4,3]
   # rmw, SE, t-stat
   M5.HM.output[i,19]<- summary(M5.HM[[i]])$coefficients[5,1]
   M5.HM.output[i,20]<- summary(M5.HM[[i]])$coefficients[5,2]
   M5.HM.output[i,21]<- summary(M5.HM[[i]])$coefficients[5,3]
   # cma, SE, t-stat
   M5.HM.output[i,22]<- summary(M5.HM[[i]])$coefficients[6,1]
   M5.HM.output[i,23]<- summary(M5.HM[[i]])$coefficients[6,2]
   M5.HM.output[i,24]<- summary(M5.HM[[i]])$coefficients[6,3]
   # psi.beta, SE, t-stat
   M5.HM.output[i,25]<- summary(M5.HM[[i]])$coefficients[7,1]
   M5.HM.output[i,26]<- summary(M5.HM[[i]])$coefficients[7,2]
   M5.HM.output[i,27]<- summary(M5.HM[[i]])$coefficients[7,3]
   # psi.smb, SE, t-stat
   M5.HM.output[i,28]<- summary(M5.HM[[i]])$coefficients[8,1]
   M5.HM.output[i,29]<- summary(M5.HM[[i]])$coefficients[8,2]
   M5.HM.output[i,30]<- summary(M5.HM[[i]])$coefficients[8,3]
   # psi.hml, SE, t-stat
   M5.HM.output[i,31]<- summary(M5.HM[[i]])$coefficients[9,1]
   M5.HM.output[i,32]<- summary(M5.HM[[i]])$coefficients[9,2]
   M5.HM.output[i,33]<- summary(M5.HM[[i]])$coefficients[9,3]
   # psi.rmw, SE, t-stat
   M5.HM.output[i,34]<- summary(M5.HM[[i]])$coefficients[10,1]
   M5.HM.output[i,35]<- summary(M5.HM[[i]])$coefficients[10,2]
   M5.HM.output[i,36]<- summary(M5.HM[[i]])$coefficients[10,3]
   # psi.cma, SE, t-stat
   M5.HM.output[i,37]<- summary(M5.HM[[i]])$coefficients[11,1]
   M5.HM.output[i,38]<- summary(M5.HM[[i]])$coefficients[11,2]
   M5.HM.output[i,39]<- summary(M5.HM[[i]])$coefficients[11,3]
}
colnames(M5.HM.output) <- c("industry", "Rsq",
                            "DWH", "DWH-p-val",
                            "J", "J-p-val",
                            "a", "SE(a)", "t(a)",
                            "b", "SE(b)", "t(b)",
                            "s", "SE(s)", "t(s)",
                            "h", "SE(h)", "t(h)",
                            "r", "SE(r)", "t(r)",
                            "c", "SE(c)", "t(c)",
                            "psi.b", "SE(psi.b)", "t(psi.b)",
                            "psi.s", "SE(psi.s)", "t(psi.s)",
                            "psi.h", "SE(psi.h)", "t(psi.h)",
                            "psi.r", "SE(psi.r)", "t(psi.r)",
                            "psi.c", "SE(psi.c)", "t(psi.c)"
                            )

# exports ---------------------------------------------------------------------
# set path
my.path <- paste("C:/Users/68596/Dropbox/eiv-coe/",
                 "scripts-output/", sep = "")
# export results
write.csv(M1.HM.output,
          file = paste(my.path,
                       "M1.HM.output.csv",
                       sep = ""))
write.csv(M3.HM.output,
          file = paste(my.path,
                       "M3.HM.output.csv",
                       sep = ""))
write.csv(M4.HM.output,
          file = paste(my.path,
                       "M4.HM.output.csv",
                       sep = ""))
write.csv(M5.HM.output,
          file = paste(my.path,
                       "M5.HM.output.csv",
                       sep = ""))
