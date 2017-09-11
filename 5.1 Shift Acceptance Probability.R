## Returns the probability of accepting a cluster configuration 'C.star' generated under the shift step.

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
