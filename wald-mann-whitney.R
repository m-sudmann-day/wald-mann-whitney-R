# Matthew Sudmann-Day
# Barcelona GSE Data Science
# Implementations of Wald and Mann-Whitney tests

# Reject if we are confident that mean of x is less than mean of y.
wald.test <- function(Y, X, alpha) {

  # Calculate the Wald test statistic, stat.
  stat <- (mean(Y) - mean(X)) / sqrt(var(Y)/length(Y) + var(X)/length(X))

  # Determine the area under the curve to the left of the point @ stat.
  pval <- pnorm(stat, lower.tail=FALSE)

  # Reject H0 if the p-value is less than the provided alpha.
  reject <- pval < alpha
  return(list(stat=stat, pval=pval, reject=reject))
}

# Reject if we are confident that mean of x is less than mean of y.
mann.whitney.test <- function(Y, X, alpha) {
  
  # Calculate the Mann Whitney test statistic, stat.
  n1 <- length(Y)
  n2 <- length(X)
  R <- sum(rank(c(Y,X))[1:n1])
  U <- n1*n2 + (n1*(n1+1))/2 - R
  mean_U <- n1*n2/2
  var_U <- (n1*n2*(n1+n2+1))/12
  stat <- (U - mean_U)/sqrt(var_U)
  
  # Determine the area under the curve to the left of the point @ stat.
  pval <- pnorm(stat, lower.tail=TRUE)
  
  # Reject H0 if the p-value is less than the provided alpha.
  reject <- pval < alpha
  
  return(list(stat=stat, pval=pval, reject=reject))
}

# Return the rate of rejections to test size of passed-in function.
mc.test.size <- function(test, alpha, S, N) {
  
  N1 <- N[,1]
  N2 <- N[,2]
  rejects <- rep(0, length(N1))
  
  # Loop through the list of experiments we've been handed.
  for (repl in 1:length(N1)) {
    
    # Reset our number of rejects for this experiment.
    r <- 0
    
    # Loop through the number of simulations in this experiment.
    # For each one, generate samples from the student t-distribution.
    # If rejected, increment r, our local rejection counter.
    for (sim in 1:S) {
      x <- rt(N1[repl], 2.1)
      y <- rt(N2[repl], 2.1)
      r <- r + test(x, y, alpha)[[3]]
    }
    
    # Store the number of rejections for this experiment.
    rejects[repl] <- r
  }
  
  # Return a vector of probabilities (the number of rejections divided
  # by the number of simulations) for each experiment.
  return(rejects/S)
}

# Return the rate of rejections to test size of passed-in function.
mc.test.power <- function(test, alpha, S, N, delta) {
  
  rejects <- rep(0, length(delta))
  n1 <- N[1]
  n2 <- N[2]
  
  # Loop through the list of experiments we've been handed.
  for (repl in 1:length(delta)) {
    
    d <- delta[repl]

    # Reset our number of rejects for this experiment.
    r <- 0

    # Loop through the number of simulations in this experiment.
    # For each one, generate samples from the student t-distribution.
    # If rejected, increment r, our local rejection counter.
    for (sim in 1:S) {
      y <- rt(n1, 2.1) + d
      x <- rt(n2, 2.1)
      r <- r + test(y, x, alpha)[[3]]
    }
    
    # Store the number of rejections for this experiment.
    rejects[repl] <- r
  }
  
  # Return a vector of probabilities (the number of rejections divided
  # by the number of simulations) for each experiment.
  return(rejects/S)
}


mc.test.power.dumb <- function(test, alpha, S, N, delta) {
  
  rejects <- rep(0, S, length(delta))

  
  # Loop through the list of experiments we've been handed.
  for (repl in 1:length(delta)) {
    
    # Reset our number of rejects for this experiment.

    # Loop through the number of simulations in this experiment.
    # For each one, generate samples from the student t-distribution.
    # If rejected, increment r, our local rejection counter.
    for (sim in 1:S) {
      y <- rt(N[1], 2.1) + delta[repl]
      x <- rt(N[2], 2.1)
      if (test(y, x, alpha)$rejects == TRUE)
        rejects[sim,repl] <- 1
    }
    
  }
  
  return(colMeans(rejects))
}

delta <- seq(0.0,500,0.1)
system.time(mc.test.power(mann.whitney.test , 0.10 , 200 , c(10,13) , delta ))
system.time(mc.test.power.dumb(mann.whitney.test , 0.10 , 200 , c(10,13) , delta ))


N <- matrix( c(4,5,6,7,8,9,10,11,20,21,500,520) , 6,2)
N
wt.size <- mc.test.size( wald.test , 0.1 , 200 , N )
mwt.size <- mc.test.size( mann.whitney.test , 0.1 , 200 , N )



plot( wt.size , t="l", col="darkred" , lwd=2 , ylim=c(0,0.20))
grid()
lines( mwt.size , t="l", col="darkblue" , lwd=2 )
lines( rep( 0.1 , nrow(N) ) , col="black" , lwd=1 )
legend( "topright", c("wald","mann-whitney") , col=c("darkred","darkblue") , lty=1 )

delta <- seq(0.0,5,0.05)
wt.power <- mc.test.power( wald.test , 0.10 , 200 , c(10,13) , delta )
mwt.power <- mc.test.power( mann.whitney.test , 0.10 , 200 , c(10,13) , delta )

plot( delta , wt.power , t="l", col="darkred" , lwd=2 , ylim=c(0,1) )
grid()
lines( delta , mwt.power , t="l", col="darkblue" , lwd=2 )
legend( "bottomright", c("wald","mann-whitney") , col=c("darkred","darkblue") , lty=1 )
