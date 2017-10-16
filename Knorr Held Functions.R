#####################################################
##################### FUNCTIONS #####################
#####################################################


##### GENERATES DISTANCE MATRIX ######

boundary.count <- function(adj){
  num = sqrt(length(adj))
  adj.powers = array(0, c(num, num, num))
  distance = array(0, c(num, num))
  d = 1
  dummy = 0
  
  for (i in 1:num){
    adj.powers[ , , i] = adj %^% i
  }
  
  for (j in 1:(num-1)){
    for (k in (j+1):num){
      while (dummy == 0){
        if (adj.powers[k, j, d] > 0){
          distance[k, j] = distance[j, k] = d
          dummy = 1
        }
        d = d + 1
      }
      dummy = 0
      d = 1
    }
  }
  return(distance)
}

##### GENERATES CLUSTERS FROM CLUSTER CENTER LIST ######

set.clusters <- function(k, distances, Gk, n, regions){
  
  not.centers <- regions[! regions %in% Gk]
  
  if(k == 1){
    return(list(append(Gk, not.centers)))
  }
  
  if (k == n){
    return(sapply(1:n, list))
  }
  
  if (k == n - 1){
    distances.to.centers <- matrix(distances[not.centers, Gk], nrow = 1)
  }else{
    distances.to.centers <- distances[not.centers, Gk]
  }
  
  
  pairing.vector <- apply(distances.to.centers, 1, which.min)
  
  center.matrix <- matrix(c(Gk, Gk), ncol = 2)
  not.center.matrix <- matrix(c(not.centers, Gk[pairing.vector]), ncol = 2)
  pairing.matrix <- rbind(center.matrix, not.center.matrix)
  
  C <- list()
  
  for(i in 1:length(Gk)){
    C[[i]] <- pairing.matrix[(pairing.matrix[ , 2] == Gk[i]), 1]
  }
  
  return(C)
}

#####  GENERATES SWITCHED CLUSTER CENTER LIST FOR SWITCH STEP ##### 

switch.step <- function(switch.candidates, Gk, n, regions){
  not.centers <- regions[! regions %in% Gk]
  Gk.star <- replace(Gk, switch.candidates, Gk[rev(switch.candidates)])
  return(Gk.star)
}

##### GENERATES SHIFTED CLUSTER CENTER LIST FOR SHIFT STEP #####

shift.step <- function(Gk, adj, regions){
  not.centers <- regions[! regions %in% Gk]
  centers.w.neighbors <- Gk[which(apply(adj[not.centers, Gk], 2, sum) > 0)]
  
  if(length(Gk[which(apply(adj[not.centers, Gk], 2, sum) > 0)]) == 1){
    g <- Gk[which(apply(adj[not.centers, Gk], 2, sum) > 0)]
  }else{
    g <- sample(centers.w.neighbors, 1)
  }
  if(length(not.centers[which(adj[g, not.centers] > 0)])==1){
    g.star <- not.centers[which(adj[g, not.centers] > 0)]
  } else{
    g.star <- sample(not.centers[which(adj[g, not.centers] > 0)], 1)
  }
  return(replace(Gk, which(Gk == g), g.star))
}

##### GENERATES CLUSTER CENTER LIST WITH A BIRTH #####

birth.step <- function(k, Gk, regions, birth.candidate, birth.location){
  Gk.star <- rep(0, k+1)
  Gk.star <- replace(Gk.star, birth.location, birth.candidate)
  Gk.star <- replace(Gk.star, which(Gk.star == 0), Gk)
  return(Gk.star)
}

##### GENERATES NEW HEIGHT LIST WITH A BIRTH IN SAME LOCATION #####

birth.add.height <- function(H, mu.tilde, sigma2.tilde, y.sums.birth, e.sums.birth){
  H.star <- rep(1, length(H)+1)
  new.h <- rgamma(1, y.sums.birth+(mu.tilde**2/sigma2.tilde), 
                  e.sums.birth+(mu.tilde/sigma2.tilde))
  H.star <- replace(H.star, birth.location, new.h)
  H.star <- replace(H.star, which(H.star == 1), H)
  return(H.star)
}

##### CALCULATES ACCEPTANCE FOR SHIFT STEP #####    

shift.alpha <- function(k, C, C.star, Gk, Gk.star, H, adj){
  likelihood.C = 0
  likelihood.C.star = 0
  
  for (i in 1:k)
  {
    likelihood.C = likelihood.C + sum(mapply(dpois, y[C[[i]]], H[i]*e[C[[i]]], log = TRUE))
  }
  
  for (i in 1:k)
  {
    likelihood.C.star = likelihood.C.star + sum(mapply(dpois, y[C.star[[i]]], H[i]*e[C.star[[i]]], log = TRUE))
  }
  
  not.centers.Gk <- regions[! regions %in% Gk]
  centers.w.neighbors.Gk <- Gk[which(apply(adj[not.centers.Gk, Gk], 2, sum) > 0)]
  
  not.centers.Gk.star <- regions[! regions %in% Gk.star]
  centers.w.neighbors.Gk.star <- Gk.star[which(apply(adj[not.centers.Gk.star, Gk.star], 2, sum) > 0)]
  
  n <- length(centers.w.neighbors.Gk)
  n.star <- length(centers.w.neighbors.Gk.star)
  
  m <- sum(adj[Gk[which(Gk - Gk.star != 0)], not.centers.Gk])
  m.star <- sum(adj[Gk.star[which(Gk - Gk.star != 0)], not.centers.Gk.star])
  
  return(likelihood.C.star - likelihood.C + log(n) - log(n.star) + log(m) - log(m.star))
}  

##### CALCULATES ACCEPTANCE FOR HEIGHT STEP #####

height.alpha <- function(H, H.star, y, y.sums, e.sums, mu.tilde, sigma2.tilde){
  
  likelihood.H = sum(y.sums*log(H) - e.sums*H)
  likelihood.H.star = sum(y.sums*log(H.star) - e.sums*H.star)
  
  p.H <- -log(H)-(1/(2*sigma2))*(log(H)-mu)**2
  p.H.star <- -log(H.star)-(1/(2*sigma2))*(log(H.star)-mu)**2
  
  q.H <- dgamma(H, y.sums + mu.tilde**2/sigma2.tilde, e.sums + 
                  mu.tilde/sigma2.tilde, log = TRUE)
  
  q.H.star <- dgamma(H.star, y.sums + mu.tilde**2/sigma2.tilde, e.sums + 
                       mu.tilde/sigma2.tilde, log = TRUE)
  
  return(likelihood.H.star - likelihood.H + p.H.star - p.H + q.H - q.H.star)
}

##### CALCULATE ACCEPTANCE FOR SWITCH STEP #####

switch.alpha <- function(k, C, C.star, H, H.star){
  likelihood.C = 0
  likelihood.C.star = 0
  
  for (i in 1:k)
  {
    likelihood.C = likelihood.C + sum(mapply(dpois, y[C[[i]]], H[i]*e[C[[i]]], log = TRUE))
  }
  
  for (i in 1:k)
  {
    likelihood.C.star = likelihood.C.star + sum(mapply(dpois, y[C.star[[i]]], 
                                                       H.star[i]*e[C.star[[i]]], log = TRUE))
  }
  
  return(likelihood.C.star - likelihood.C)
}

##### CALCULATE ACCEPTANCE FOR BIRTH STEP #####

birth.alpha.simple <- function(c, H, H.star, mu, sigma2, mu.tilde, sigma2.tilde, C, C.star, birth.location, y.sums, 
                               y.sums.star, e.sums, e.sums.star){
  
  likelihood.C <- 0
  likelihood.C.star <- 0
  k <- length(H)
  
  p.h.star <- -log(H.star[birth.location])  - log(sqrt(2*pi*sigma2)) - 
    (1/(2*sigma2))*(log(H.star[birth.location]) - mu)**2
  
  q.h.star <- dgamma(H.star[birth.location], y.sums.star[birth.location] + mu.tilde**2/sigma2.tilde, 
                     e.sums.star[birth.location] + mu.tilde/sigma2.tilde, log = TRUE)
  
  likelihood.C.star = sum(dpois(y[C.star[[birth.location]]], e[C.star[[birth.location]]]*
                                  H.star[birth.location], log = TRUE))
  
  H.full <- rep(0, n)
  for (i in 1:length(C)){
    H.full[C[[i]]] = H[i]
  }
  
  
  likelihood.C = sum(dpois(y[C.star[[birth.location]]], e[C.star[[birth.location]]]*
                             H.full[C.star[[birth.location]]], log = TRUE))
  
  return(log(1-c) + p.h.star - q.h.star + likelihood.C.star - likelihood.C)
}

##### CALCULATE ACCEPTANCE FOR DEATH STEP #####

death.alpha.simple <- function(c, H, H.star, mu, sigma2, mu.tilde, sigma2.tilde, C, C.star, death.location, y.sums, 
                               y.sums.star, e.sums, e.sums.star){
  likelihood.C <- 0
  likelihood.C.star <- 0
  k <- length(H)
  
  p.h <- -log(H[death.location]) - log(sqrt(2*pi*sigma2)) -
    (1/(2*sigma2))*(log(H[death.location]) - mu)**2
  
  q.h <- dgamma(H[death.location], y.sums[death.location] + mu.tilde**2/sigma2.tilde, 
                e.sums[death.location] + mu.tilde/sigma2.tilde, log = TRUE)
  
  likelihood.C <- sum(dpois(y[C[[death.location]]], e[C[[death.location]]]*H[death.location], log = TRUE))
  
  H.star.full <- rep(0, n)
  for (i in 1:length(C.star)){
    H.star.full[C.star[[i]]] = H.star[i]
  }
  
  likelihood.C.star = sum(dpois(y[C[[death.location]]], e[C[[death.location]]]*
                                  H.star.full[C[[death.location]]], log = TRUE))
  
  return(-log(1-c) - p.h + q.h + likelihood.C.star - likelihood.C)
}


##### GENERATE DATA SUMS FOR EACH CLUSTER #####

data.sums <- function(C, data, sum.vector){
  for (i in 1:length(C.star)){
    sum.vector[i] = sum(data[C[[i]]])
  }
  return(sum.vector)
}

###################################################################
################# TO LOAD GERMANY DATASET #########################
###################################################################

wege <- scan("~/Downloads/wege.txt", sep = "")
wege <- matrix(wege, ncol = 544, byrow = TRUE)

adjacent <- scan("~/Downloads/adjacent.txt", sep = "")
adjacent <- matrix(adjacent, nrow = 544, byrow = TRUE)
adjacent <- adjacent[, -1]
adj <- matrix(0, nrow = 544, ncol = 544)

for (i in 1:544){
  adj[i, adjacent[i, which(adjacent[i, ]>0)]] = 1
}

diag(adj) = 1

library(readr)

mund <- read_delim("~/Downloads/mund.dat", 
                   " ", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE)

mund <- data.matrix(mund)

distances <- wege
e <- mund[, 2]
y <- mund[, 1]

###################################################################
##################### DATA AND PARAMETERS #########################
###################################################################

bbb <- 10000
rrr <- 10000
skip <- rrr/10000
a <- 1
b <- .01
c <- .02
saved.values <- seq(from = skip, to = rrr, by = skip)
n <- length(y)
k <- sample(1:544, 1)
mu <- 2
sigma2 <- 2
mu.tilde <- exp(mu+.5*sigma2)
sigma2.tilde <- exp(sigma2)*(exp(sigma2)-1)*exp(2*mu)
regions <- seq(1:length(y))
Gk <- sample(regions, k)
C <- set.clusters(k, distances, Gk, n, regions)
H <- rep(1, k)
hyper.out <- matrix(rep(0, length(saved.values)*2), ncol = 2)
H.out <- matrix(0, nrow = length(saved.values), ncol = n)
Gk.out <- matrix(0, nrow = length(saved.values), ncol = n)
k.out <- rep(0, length(saved.values))
shift.accept <- 0
shift.attempts <- 0
switch.accept <- 0
switch.attempts <- 0
height.accept <- 0
height.attempts <- 0
birth.attempts <- 0
birth.accept <- 0
death.attempts <- 0
death.accept <- 0
height.accept1 <- 0 
height.attempts1 <- 0
