
## Given a cluster configuration 'C', a proposed cluster configuration under the birth step 'C.star', and corresponding
## vectors of heights 'H' and 'H.star', this function returns the probability of accepting the new cluster configuration.

birth.alpha <- function(c, H, H.star, mu, sigma2, mu.tilde, sigma2.tilde, C, C.star, birth.location, y.sums, 
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
