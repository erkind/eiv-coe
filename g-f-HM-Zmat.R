# Calculating Z-matrix of higher-moments instruments
#
generate.Z.mat <- function(
      X.mat,		# n by k matrix of independent variables
      Y.mat,		# n by 1 vector of dependent variable
      n,          # sample size 
      k,          # nb of predictors
      i = c(1, 4) # list of higher-moments
      )
   {
   #
   # # set the sample size and the number of predictors
   # n <- dim(X.mat)[1]
   # k <- dim(X.mat)[2]
   #
   # calculate X.mat in mean deviations form
   x.mat <- apply(as.data.frame(X.mat), 2, function(a){return(a - mean(a))})
   #
   # calculate Y.mat in mean deviations form
   y.mat <- apply(as.data.frame(Y.mat), 2, function(a){return(a - mean(a))})
   #
   # calculating higher-moments instruments
   Z.mat <- rep(1, n, 1)
   # Z.mat <- rep(1, n)
   if(is.element(1, i)){
      Z.mat <- cbind(Z.mat, x.mat * x.mat)
   }
   if(is.element(2, i)){
      Z.mat <- cbind(Z.mat, x.mat * y.mat)
   }
   if((is.element(3, i))){
      Z.mat <- cbind(Z.mat, y.mat * y.mat)
   }
   if((is.element(4, i))){
      part1 <- x.mat * x.mat * x.mat
      part2 <- 3 * x.mat %*% (mean((t(x.mat) %*% x.mat) / n) * diag(k))
      my.value <- part1 - part2
      Z.mat <- cbind(Z.mat, my.value)
   }
   if((is.element(5, i))){
      part1 <- x.mat * x.mat * y.mat
      part2 <- 2 * x.mat %*% (mean(t(x.mat) %*% y.mat) * diag(k))
      part3 <- y.mat * sum(mean(t(x.mat) %*% x.mat / n) * diag(k))
      my.value <- part1 - part2 - part3
      Z.mat <- cbind(Z.mat, my.value)
   }
   if((is.element(6, i))){
      part1 <- x.mat * y.mat * y.mat
      part2 <- x.mat * mean((t(y.mat) %*% y.mat / n))
      part3 <- 2 * y.mat * mean((t(y.mat) %*% x.mat / n))
      my.value <- part1 - part2 - part3
      Z.mat <- cbind(Z.mat, my.value)
   }
   if((is.element(7, i))){
      my.value <- y.mat * y.mat * y.mat - 3 * y.mat * mean(t(y.mat) %*% y.mat / n)
      Z.mat <- cbind(Z.mat, my.value)
   }
   #
   # define valid ranges for i and j
   valid_i <- 1:7
   valid_k <- 1:5
   # check if the values of i and k are within the valid ranges
   if(any(!i %in% valid_i)){
      stop("Invalid value(s) for i
           This should be an integer between 1 and 7.")
   }
   if (!(k %in% valid_k)) {
      stop("Invalid value for k.
           It should be an integer between 1 and 5.")
   }
   # initialize an empty vector to store labels
   labels <- c()
   # always include z0 in the labels vector
   labels <- c(labels, "z0")
   # # loop through each value of i
   # for (ii in i){
   #    # loop through all numbers up to j for each value of i
   #    for (kk in 1:k) {
   #       # generate labels and concatenate to the labels vector
   #       labels <- c(labels, paste0("z", ii, kk))
   #    }
   # }
   # if k is 1, add additional labels for 1:1 for each value of i
   if (k == 1){
      for (ii in i){
         labels <- c(labels, paste0("z", ii, 1:1))
      }
   }
   # if k is 2, add additional labels for 1:2 for each value of i
   if (k == 2){
      for (ii in i){
         labels <- c(labels, paste0("z", ii, 1:2))
      }
   }
   # if k is 3, add additional labels for 1:3 for each value of i
   if (k == 3){
      for (ii in i){
         labels <- c(labels, paste0("z", ii, 1:3))
      }
   }
   # if k is 4, add additional labels for 1:4 for each value of i
   if (k == 4){
      for (ii in i){
         labels <- c(labels, paste0("z", ii, 1:4))
      }
   }
   # if k is 5, add additional labels for 1:5 for each value of i
   if (k == 5){
      for (ii in i){
         labels <- c(labels, paste0("z", ii, 1:5))
      }
   }
   colnames(Z.mat) <- labels
   # return Z.mat
   return(Z.mat)
}


