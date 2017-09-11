## Given a cluster configuration 'C', a proposed cluster under the death step 'C.star', and corresponding heights, this function
## returns the probability of accepting the new cluster configuration 'C.star'.

death.alpha <- function(c, H, H.star, mu, sigma2, mu.tilde, sigma2.tilde, C, C.star, death.location, y.sums, 
                        y.sums.star, e.sums, e.sums.star){
  likelihood.C <- 0
  likelihood.C.star <- 0
  k <- length(H)
  
  p.h <- -log(H[death.location]) + log(sqrt(2*pi*sigma2)) -
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
