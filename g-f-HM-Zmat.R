# calculating Z-matrix of higher-moments instruments
#
generate.Z.mat <- function(
      X.mat,		      # n by k matrix of independent variables
      Y.vec = NULL,     # n by 1 vector of dependent variable
      n,                # sample size 
      k,                # nb of predictors
      iv.list = c(1, 4) # list of higher-moments
)
{
   #
   # nb of cols in X.mat must be equal to k
   if(k != dim(as.matrix(X.mat))[2]){
      stop("The number of columns in X.mat must be equal to k.")
   }
   #
   # calculate Y.vec and X.mat in mean deviations form
   y.vec <- scale(Y.vec, scale = FALSE)
   x.mat <- scale(X.mat, scale = FALSE)
   #
   # calculate higher-moment instruments
   #
   # z0
   # calculate z0 = vector of ones
   Z.mat <- matrix(rep(1, n), nrow = n, ncol = 1)
   # label z0
   colnames(Z.mat) <- "z0"
   #
   # column labels index for additional instruments
   col.index <- 1
   #
   # z1
   if(is.element(1, iv.list)){
      # calculate z1
      z1 <- x.mat * x.mat
      # append Z.mat
      Z.mat <- cbind(Z.mat, z1)
      # label z1
      new.cols <- paste0("z1", 1:k)
      colnames(Z.mat)[(col.index + 1):(col.index + k)] <- new.cols
      col.index <- col.index + k
   }
   # z2
   if(is.element(2, iv.list)){
      # calculate z2
      z2 <- matrix(0, nrow = n, ncol = k)
      for(j in 1:k){
         z2[,j] <- x.mat[,j] * y.vec
      }
      # append Z.mat
      Z.mat <- cbind(Z.mat, z2)
      # label z2
      new.cols <- paste0("z2", 1:k)
      colnames(Z.mat)[(col.index + 1):(col.index + k)] <- new.cols
      col.index <- col.index + k
   }
   # z3
   if((is.element(3, iv.list))){
      # calculate z3
      z3 <- y.vec * y.vec
      # append Z.mat
      Z.mat <- cbind(Z.mat, z3)
      # label z3
      colnames(Z.mat)[col.index + 1] <- "z31"
      col.index <- col.index + 1
   }
   # z4
   if((is.element(4, iv.list))){
      # calculate z4
      # z4.1 and z4.2
      z4.1 <- x.mat * x.mat * x.mat
      z4.2 <- 3 * x.mat %*% (((t(x.mat) %*% x.mat) / n) * diag(k))
      #
      z4 <- z4.1 - z4.2
      # append Z.mat
      Z.mat <- cbind(Z.mat, z4)
      # label z4
      new.cols <- paste0("z4", 1:k)
      colnames(Z.mat)[(col.index + 1):(col.index + k)] <- new.cols
      col.index <- col.index + k
   }
   # z5
   if((is.element(5, iv.list))){
      # calculate z5
      # z5.1
      z5.1 <- matrix(0, nrow = n, ncol = k)
      for(j in 1:k){
         z5.1[,j] <- x.mat[,j] * x.mat[,j] * y.vec
      }
      # z5.2
      z5.2 <- 2 * (x.mat %*% 
                      diag(
                         as.vector(
                            (1/n) * (t(x.mat) %*% y.vec)
                         ),
                         nrow = k, ncol = k
                      )
      )
      # z5.3
      # z5.3.1 and z5.3.2
      z5.3.1 <- matrix(1, nrow = k, ncol = 1)
      z5.3.2 <- ((t(x.mat) %*% x.mat) / n) * diag(k)
      # 
      z5.3 <- y.vec %*% (t(z5.3.1) %*% z5.3.2)
      # 
      z5 <- z5.1 - z5.2 - z5.3
      # append Z.mat
      Z.mat <- cbind(Z.mat, z5)
      # label z5
      new.cols <- paste0("z5", 1:k)
      colnames(Z.mat)[(col.index + 1):(col.index + k)] <- new.cols
      col.index <- col.index + k
   }
   # z6
   if((is.element(6, iv.list))){
      # calculate z6
      # z6.1
      z6.1 <- matrix(0, nrow = n, ncol = k)
      for(j in 1:k){
         z6.1[,j] <- x.mat[,j] * y.vec * y.vec
      }
      # z6.2
      scalar <- (1/n) * (t(y.vec) %*% y.vec)
      z6.2 <- x.mat * as.numeric(scalar)
      # z6.3
      z6.3 <- 2 * (
         y.vec %*% ((1/n) * (t(y.vec) %*% x.mat))
      )
      #
      z6 <- z6.1 - z6.2 - z6.3
      # append Z.mat
      Z.mat <- cbind(Z.mat, z6)
      # label z6
      new.cols <- paste0("z6", 1:k)
      colnames(Z.mat)[(col.index + 1):(col.index + k)] <- new.cols
      col.index <- col.index + k
   }
   # z7
   if((is.element(7, iv.list))){
      # calculate z7
      z7 <- matrix(0, nrow = n, ncol = 1)
      #
      z7.1 <- y.vec * y.vec * y.vec
      #
      scalar <- (1/n) * (t(y.vec) %*% y.vec)
      z7.2 <- 3 * y.vec * as.numeric(scalar)
      #
      z7 <- z7.1 - z7.2
      # append Z.mat
      Z.mat <- cbind(Z.mat, z7)
      # label z7
      colnames(Z.mat)[col.index + 1] <- "z71"
      col.index <- col.index + 1
   }
   #
   return(Z.mat)
}