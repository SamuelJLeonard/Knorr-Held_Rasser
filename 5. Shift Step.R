## This function returns a new cluster configuration under the shift step.  If 'n' is a positive integer, 'sample(n, 1)' R will
## return a random sample from 1 to n.  This causes problems if there is only one eligible cluster for the shift step.  
## The if statements below are a quick solution to this problem.

shift.step <- function(Gk, adj){
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
