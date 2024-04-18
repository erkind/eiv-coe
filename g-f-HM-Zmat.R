# calculating Z-matrix of higher-moments instruments
#
generate.Z.mat <- function(
      X.mat,		  # n by k matrix of independent variables
      Y.vec = NULL, # n by 1 vector of dependent variable
      n,            # sample size 
      k,            # nb of predictors
      i = c(1, 4)   # list of higher-moments
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
   # z0 = unit vector
   Z.mat <- rep(1, n, 1)
   # z1
   if(is.element(1, i)){
      Z.mat <- cbind(Z.mat, x.mat * x.mat)
   }
   # z2
   if(is.element(2, i)){
      my.value <- matrix(0, nrow = n, ncol = ncol(x.mat))
      for(j in 1:ncol(x.mat)){
         my.value[,j] <- x.mat[,j] * y.vec
      }
      Z.mat <- cbind(Z.mat, my.value)
   }
   # z3
   if((is.element(3, i))){
      Z.mat <- cbind(Z.mat, y.vec * y.vec)
   }
   # z4
   if((is.element(4, i))){
      part1 <- x.mat * x.mat * x.mat
      part2 <- 3 * x.mat %*% (mean((t(x.mat) %*% x.mat) / n) * diag(k))
      my.value <- part1 - part2
      Z.mat <- cbind(Z.mat, my.value)
   }
   # z5
   if((is.element(5, i))){
      #
      part1 <- matrix(0, nrow = n, ncol = ncol(x.mat))
      for(j in 1:ncol(x.mat)){
         part1[,j] <- x.mat[,j] * x.mat[,j] * y.vec
      }
      part2 <- matrix(0, nrow = n, ncol = ncol(x.mat))
      part2 <- 2 * x.mat %*% (mean(t(x.mat) %*% y.vec) * diag(ncol(x.mat)))
      # 
      part3 <- matrix(0, nrow = n, ncol = ncol(x.mat))
      part3 <- y.vec %*% t(rep(1, ncol(x.mat))) %*%
         (mean((t(x.mat) %*% x.mat) / n) * diag(k))
      # 
      my.value <- part1 - part2 - part3
      Z.mat <- cbind(Z.mat, my.value)
   }
   # z6
   if((is.element(6, i))){
      part1 <- x.mat * y.vec * y.vec
      part2 <- x.mat * mean((t(y.vec) %*% y.vec / n))
      part3 <- 2 * y.vec * mean((t(y.vec) %*% x.mat / n))
      my.value <- part1 - part2 - part3
      Z.mat <- cbind(Z.mat, my.value)
   }
   # z7
   if((is.element(7, i))){
      my.value <- y.vec * y.vec * y.vec - 3 * y.vec * mean(t(y.vec) %*% y.vec / n)
      Z.mat <- cbind(Z.mat, my.value)
   }
   #
   # Z.mat column labels
   #
   # define valid range for i
   i.range <- 1:7
   # i values consistent with the nb of columns in Z.mat, i in 1:7
   # 
   # define valid range for k
   k.range <- 1:5
   # k values consistent with the max nb of regressors among different
   # specifications. the default value set as 5 to handle the FF 5-factor
   # model.
   # if k > 5, the user should add additional lines after the command line
   # beginning with " if(k == 5){... " below.
   #
   # check if the values of i and k are within the valid ranges
   if(any(!i %in% i.range)){
      stop("Invalid value(s) for i
           This should be an integer between 1 and 7.")
   }
   if (!(k %in% k.range)) {
      stop("Invalid value for k.
           It should be an integer between 1 and 5.")
   }
   #
   # initialize empty vector
   labels <- c()
   #
   # always include z0 in the labels vector
   labels <- c(labels, "z0")
   #
   # if k is 1, add additional labels for 1:1 for each value of i
   if(k == 1){
      for (ii in i){
         labels <- c(labels, paste0("z", ii, 1:1))
      }
   }
   # if k is 2, add additional labels for 1:2 for each value of i
   if(k == 2){
      for (ii in i){
         labels <- c(labels, paste0("z", ii, 1:2))
      }
   }
   # if k is 3, add additional labels for 1:3 for each value of i
   if(k == 3){
      for (ii in i){
         labels <- c(labels, paste0("z", ii, 1:3))
      }
   }
   # if k is 4, add additional labels for 1:4 for each value of i
   if(k == 4){
      for (ii in i){
         labels <- c(labels, paste0("z", ii, 1:4))
      }
   }
   # if k is 5, add additional labels for 1:5 for each value of i
   if(k == 5){
      for (ii in i){
         labels <- c(labels, paste0("z", ii, 1:5))
      }
   }
   #
   # rename Z.mat columns
   colnames(Z.mat) <- labels
   #
   # return Z.mat
   return(Z.mat)
}