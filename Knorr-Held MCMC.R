
## KNORR-HELD AND RASSER - DISEASE CLUSTERING ALGORITHM ##
## CODE WRITTEN BY SAM LEONARD ##



## INSTRUCTIONS: ##

  ##  1.  Download the 'Knorr Held Functions and Data.R' file, if you haven't already.  
  ##  2.  Install and load the packages below.  You may have them already if you use R a lot.    
  ##  3.  Read the functions file by running the two lines under '# Functions'.  Select 'Knorr Held Functions and Data.R' whe prompted.
  ##  4.  You can now run the algorithm with the German dataset.  It's the code under the MCMC heading.  I have set 10,000 iterations
  ##      after 10,000 burn-in.  The algorithm will take about 3 minutes to run.
  ##  5.  You can use code under the RESULTS heading to view the some of the MCMC results.  This is a small run, so 
  ##      you will see variation in the acceptance rates and distribution of k (the number of clusters).  


# Packages:

install.packages('expm')
install.packages('spam')
install.packages('nimble')
install.packages('stats')
library('expm')
library('spam')
library('nimble')
library('stats')

# Functions:

FUNCTIONS_FILE <- file.choose()
source(FUNCTIONS_FILE)

    


#####################################################
##################### MCMC ##########################
#####################################################

for(j in 1:bbb){
  
  if (length(Gk) == 1){
    step <- sample(c(1, 3, 4, 5, 6), 1, prob = c(.8, .05, .05, .05, .05))
  }else if(length(Gk) == n){
    step <- sample(c(2, 3, 4, 5, 6), 1, prob = c(.8, .05, .05, .05, .05))
  }else{
    step <- sample(1:6, 1, prob = c(.4, .4, .05, .05, .05, .05))
  }
  
  y.sums.star <- c()
  y.sums <- c()
  e.sums.star <- c()
  e.sums <- c()
  
  ## Birth step
  
  if (step == 1){
    
    not.centers <- regions[! regions %in% Gk]
    birth.candidate <- sample(not.centers, 1)
    birth.location <- sample(k+1, 1)
    
    Gk.star <- birth.step(k, Gk, regions, birth.candidate, birth.location)
    C.star <- set.clusters(k+1, distances, Gk.star, n, regions)
    
    y.sums.star <- sapply(C.star, function(x)sum(y[x]))
    e.sums.star <- sapply(C.star, function(x)sum(e[x]))
    y.sums <- sapply(C, function(x)sum(y[x]))
    e.sums <- sapply(C, function(x)sum(e[x]))
    
    H.star <- birth.add.height(H, mu.tilde, sigma2.tilde, y.sums.star[birth.location], 
                               e.sums.star[birth.location])
    
    if (log(runif(1, 0, 1)) < birth.alpha.simple(c, H, H.star, mu, sigma2, mu.tilde, sigma2.tilde, C, C.star, birth.location, y.sums, y.sums.star, 
                                          e.sums, e.sums.star)){
      Gk <- Gk.star
      C <- C.star
      H <- H.star
    }
  }
  
  # Death step
  
  if (step == 2){
    
    death.location <- sample(k, 1)
    Gk.star <- Gk[-death.location]
    H.star <- H[-death.location]
    C.star <- set.clusters(k-1, distances, Gk.star, n, regions)
  
    y.sums.star <- sapply(C.star, function(x)sum(y[x]))
    e.sums.star <- sapply(C.star, function(x)sum(e[x]))
    y.sums <- sapply(C, function(x)sum(y[x]))
    e.sums <- sapply(C, function(x)sum(e[x]))

    
    if (log(runif(1, 0, 1)) < death.alpha.simple(c, H, H.star, mu, sigma2, mu.tilde, sigma2.tilde, C, C.star, 
                                          death.location, y.sums, y.sums.star, e.sums, e.sums.star)){
      Gk <- Gk.star
      C <- C.star
      H <- H.star
    }
  }
  
  
  # Shift step
  
  if (step == 3){
  
    Gk.star <- shift.step(Gk, adj, regions)
    C.star <- set.clusters(k, distances, Gk.star, n, regions)
    
    if (log(runif(1, 0, 1)) < shift.alpha(k, C, C.star, Gk, Gk.star, H, adj)){
      C <- C.star
      Gk <- Gk.star
    }
  }
  
  # Switch step
  
  if (step == 4){
    
    switch.candidates <- sample(k, 2)
    Gk.star <- switch.step(switch.candidates, Gk, n, regions)
    C.star <- set.clusters(k, distances, Gk.star, n, regions)
    H.star <- H
    H.star[c(switch.candidates)] <- H[c(switch.candidates[2], switch.candidates[1])]
    
    if (log(runif(1, 0, 1)) < switch.alpha(k, C, C.star, H, H.star)){
      Gk <- Gk.star
      C <- C.star
      H <- H.star
    }
  }
  
  # Height step
  
  if (step == 5){
  
    y.sums <- sapply(C, function(x)sum(y[x]))
    e.sums <- sapply(C, function(x)sum(e[x]))
    
    H.star <- mapply(rgamma, 1, y.sums+(mu.tilde**2/sigma2.tilde), e.sums+(mu.tilde/sigma2.tilde))
    H.star.alpha <- height.alpha(H, H.star, y, y.sums, e.sums, mu.tilde, sigma2.tilde)
    
    for(i in 1:k){
      if(log(runif(1, 0, 1)) < H.star.alpha[i]){
        H[i] <- H.star[i]
      }
    }
  }
  
  # Hyperparameter step
  
  if (step == 6){
    mu.mean <-(1/k)*sum(log(H))
    mu.sigma <- sqrt((1/k)*sigma2)
    sigma2.a <- a + (k/2)
    sigma2.b <- b + .5 * sum((log(H)-mu)**2)
    
    mu <- rnorm(1, mu.mean, mu.sigma)
    sigma2 <- 1/rgamma(1, sigma2.a, sigma2.b)
    mu.tilde <- exp(mu+.5*sigma2)
    sigma2.tilde <- exp(sigma2)*(exp(sigma2)-1)*exp(2*mu)
  }
  
  k <- length(Gk)
}




for(j in 1:rrr){
  
  if (length(Gk) == 1){
    step <- sample(c(1, 3, 4, 5, 6), 1, prob = c(.8, .05, .05, .05, .05))
  }else if(length(Gk) == n){
    step <- sample(c(2, 3, 4, 5, 6), 1, prob = c(.8, .05, .05, .05, .05))
  }else{
  step <- sample(1:6, 1, prob = c(.4, .4, .05, .05, .05, .05))
  }
  
  y.sums.star <- c()
  y.sums <- c()
  e.sums.star <- c()
  e.sums <- c()
  
  ## Birth step
  
  if (step == 1){
    birth.attempts = birth.attempts + 1
    
    not.centers <- regions[! regions %in% Gk]
    birth.candidate <- sample(not.centers, 1)
    birth.location <- sample(k+1, 1)
    
    Gk.star <- birth.step(k, Gk, regions, birth.candidate, birth.location)
    C.star <- set.clusters(k+1, distances, Gk.star, n, regions)
    
    y.sums.star <- sapply(C.star, function(x)sum(y[x]))
    e.sums.star <- sapply(C.star, function(x)sum(e[x]))
    y.sums <- sapply(C, function(x)sum(y[x]))
    e.sums <- sapply(C, function(x)sum(e[x]))
    
    H.star <- birth.add.height(H, mu.tilde, sigma2.tilde, y.sums.star[birth.location], 
                               e.sums.star[birth.location])
    
    mu.tilde <- exp(mu+.5*sigma2)
    sigma2.tilde <- exp(sigma2)*(exp(sigma2)-1)*exp(2*mu)
    
    if (log(runif(1, 0, 1)) < birth.alpha.simple(c, H, H.star, mu, sigma2, mu.tilde, sigma2.tilde, C, C.star, birth.location, y.sums, y.sums.star, 
                                         e.sums, e.sums.star)){
      birth.accept = birth.accept + 1
      Gk <- Gk.star
      C <- C.star
      H <- H.star
    }
  }
  
  # Death step
  
  if (step == 2){
    
    death.attempts = death.attempts + 1
    death.location <- sample(k, 1)
    Gk.star <- Gk[-death.location]
    H.star <- H[-death.location]
    C.star <- set.clusters(k-1, distances, Gk.star, n, regions)
    
    y.sums.star <- sapply(C.star, function(x)sum(y[x]))
    e.sums.star <- sapply(C.star, function(x)sum(e[x]))
    y.sums <- sapply(C, function(x)sum(y[x]))
    e.sums <- sapply(C, function(x)sum(e[x]))
    
    mu.tilde <- exp(mu+.5*sigma2)
    sigma2.tilde <- exp(sigma2)*(exp(sigma2)-1)*exp(2*mu)
    
    death.alpha.simple(c, H, H.star, mu, sigma2, mu.tilde, sigma2.tilde, C, C.star, death.location, 
                       y.sums, y.sums.star, e.sums, e.sums.star)
      
    if (log(runif(1, 0, 1)) < death.alpha.simple(c, H, H.star, mu, sigma2, mu.tilde, sigma2.tilde, C, C.star, 
                                         death.location, y.sums, y.sums.star, e.sums, e.sums.star)){
      death.accept = death.accept + 1
      Gk <- Gk.star
      C <- C.star
      H <- H.star
    }
  }
  

  # Shift step
  
  if (step == 3){
    shift.attempts = shift.attempts + 1
    Gk.star <- shift.step(Gk, adj, regions)
    C.star <- set.clusters(k, distances, Gk.star, n, regions)
  
    if (log(runif(1, 0, 1)) < shift.alpha(k, C, C.star, Gk, Gk.star, H, adj)){
      C <- C.star
      Gk <- Gk.star
      shift.accept = shift.accept + 1
    }
  }
  
  # Switch step
  
  if (step == 4){
    switch.attempts = switch.attempts + 1
    switch.candidates <- sample(k, 2)
    Gk.star <- switch.step(switch.candidates, Gk, n, regions)
    C.star <- set.clusters(k, distances, Gk.star, n, regions)
    H.star <- H
    H.star[c(switch.candidates)] <- H[c(switch.candidates[2], switch.candidates[1])]
    
   
      if (log(runif(1, 0, 1)) < switch.alpha(k, C, C.star, H, H.star)){
        Gk <- Gk.star
        C <- C.star
        switch.accept = switch.accept + 1
        H <- H.star
      }
  }
  
  # Height step

  if (step == 5){
    height.attempts = height.attempts + 1
    accept = FALSE
    
    y.sums <- sapply(C, function(x)sum(y[x]))
    e.sums <- sapply(C, function(x)sum(e[x]))
    
    mu.tilde <- exp(mu+.5*sigma2)
    sigma2.tilde <- exp(sigma2)*(exp(sigma2)-1)*exp(2*mu)
    H.star <- mapply(rgamma, 1, y.sums+(mu.tilde**2/sigma2.tilde), e.sums+(mu.tilde/sigma2.tilde))
    H.star.alpha <- height.alpha(H, H.star, y, y.sums, e.sums, mu.tilde, sigma2.tilde)
    height.attempts1 = height.attempts1 + k
  
    for(i in 1:k){
      if(log(runif(1, 0, 1)) < H.star.alpha[i]){
        H[i] <- H.star[i]
        height.accept1 = height.accept1 + 1
        accept = TRUE
      }
    }
    if (accept){
      height.accept = height.accept + 1
    }
  }
  
  # Hyperparameter step

  if (step == 6){
    mu.mean <- (1/k)*sum(log(H))
    mu.sigma <- sqrt((1/k)*sigma2)
    sigma2.a <- a + (k/2)
    sigma2.b <- b + .5 * sum((log(H)-mu)**2)
    
    mu <- rnorm(1, mu.mean, mu.sigma)
    sigma2 <- 1/rgamma(1, sigma2.a, sigma2.b)
  }
  
  k <- length(Gk)
  
  if (j %in% saved.values){
    hyper.out[j/skip, ] <- c(mu, sigma2)
    Gk.out[j/skip, ] <- c(Gk, rep(0, (n - length(Gk))))
    k.out[j/skip] <- k
  
    for (i in 1:length(C)){
      H.out[j/skip, C[[i]]] = H[i]
    }
  }  
} 

#####################################################
##################  RESULTS #########################
#####################################################


### Colored map of median relative risk at each region ###

H.medians <- apply(H.out, 2, median)
map.landkreis(H.medians, col = NULL, zlim = range(H.medians), add = FALSE, 
              legendpos = c(.88, .9, .05, .4))

### Acceptance rates ###

print('birth:')
print(birth.accept/birth.attempts)

print('death:')
print(death.accept/death.attempts)

print('height:')
print(height.accept/height.attempts)

print('switch')
print(switch.accept/switch.attempts)

print('shift')
print(shift.accept/shift.attempts)

### Posterior distribution of k ###

k.count <- rep(0, 544)
for (i in 1:544){
  k.count[i] <- length(which(k.out == i))
}
k.prob <- k.count/(length(saved.values))
plot(k.prob, pch = 20, xlab = 'k', ylab = '', type = 'l')




