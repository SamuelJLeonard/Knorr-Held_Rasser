## Given a cluster configuration 'C' and a configuration generated under the switch step 'C.star', this function returns the
## probability of accepting 'C.star' as the new cluster configuration.

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
