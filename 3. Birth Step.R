 ## This function generates a new list of cluster centers 'Gk.star' under the conditions of the birth step.
 
 birth.step <- function(k, Gk, regions, birth.candidate, birth.location){
      Gk.star <- rep(0, k+1)
      Gk.star <- replace(Gk.star, birth.location, birth.candidate)
      Gk.star <- replace(Gk.star, which(Gk.star == 0), Gk)
      return(Gk.star)
    }
