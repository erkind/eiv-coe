# HM estimation of factor models for ind portfolios

# data ------------------------------------------------------------------------
# # loading from hard drive, mind customizing the path
# source(
#    paste(
#       "C:/Users/68596/Dropbox/eiv-coe/",
#       "scripts/data-eiv-coe.R",
#       sep = "")
#    )
# loading from GitHub
source(
   paste(
      "http://raw.githubusercontent.com/erkind/",
      "eiv-coe/main/data-eiv-coe.R",
      sep = "")
   )
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
# choose the industry between 1 and 44
my.ind <- 1
# fitting different models with OLS vs. HM for the 1st industry
# CAPM
M1.LS <- lm(R[,my.ind] ~ F1[,1])

M1.HM <- lm(R[,my.ind] ~ F1[,1] + w.hat.F1.MKT)
# FF3F
M3.LS <- lm(R[,my.ind] ~ F3[,1] + F3[,2] + F3[,3])
M3.HM <- lm(R[,my.ind] ~ F3[,1] + F3[,2] + F3[,3] + 
               w.hat.F3.MKT + w.hat.F3.SMB + w.hat.F3.HML)
# FFC4F
M4.LS <- lm(R[,my.ind] ~ F4[,1] + F4[,2] + F4[,3] + F4[,4])
M4.HM <- lm(R[,my.ind] ~ F4[,1] + F4[,2] + F4[,3] + F4[,4] +  
               w.hat.F3.MKT + w.hat.F3.SMB +
               w.hat.F3.HML + w.hat.F4.WML)
# FF5F
M5.LS <- lm(R[,my.ind] ~ F5[,1] + F5[,2] + F5[,3] + F5[,4] + F5[,5])
M5.HM <- lm(R[,my.ind] ~ F5[,1] + F5[,2] + F5[,3] + F5[,4] + F5[,5] +
               w.hat.F5.MKT + w.hat.F5.SMB + w.hat.F5.HML + 
               w.hat.F5.RMW + w.hat.F5.CMA)
#
names(M1.LS$coefficients) <- c("a", "b")
names(M3.LS$coefficients) <- c("a", "b", "s", "h")
names(M4.LS$coefficients) <- c("a", "b", "s", "h", "w")
names(M5.LS$coefficients) <- c("a", "b", "s", "r", "c")
names(M1.HM$coefficients) <- c("a", "b",
                               "psi.b")
names(M3.HM$coefficients) <- c("a", "b", "s", "h",
                               "psi.b", "psi.s", "psi.h")
names(M4.HM$coefficients) <- c("a", "b", "s", "h", "w",
                               "psi.b", "psi.s", "psi.h", "psi.w")
names(M5.HM$coefficients) <- c("a", "b", "s", "h", "r", "c",
                               "psi.b", "psi.s", "psi.h", "psi.r", "psi.c")

# output ----------------------------------------------------------------------
summary(M1.LS)
summary(M1.HM)
summary(M3.LS)
summary(M3.HM)
summary(M4.LS)
summary(M4.HM)
summary(M5.LS)
summary(M5.HM)


