# Given a distance matrix and a vector of cluster centers 'Gk', this function returns the corresponding clustering 
# configuration as a list.

set.clusters <- function(k, distances, Gk, n){
  regions <- seq(1:n)
  not.centers <- regions[! regions %in% Gk]

  distances.to.centers <- distances[not.centers, Gk]
  pairing.vector <- apply(distances.to.centers, 1, which.min)

  center.matrix <- matrix(c(Gk, Gk), ncol = 2)
  not.center.matrix <- matrix(c(not.centers, Gk[pairing.vector]), ncol = 2)
  pairing.matrix <- rbind(center.matrix, not.center.matrix)

  C.not.centers <- list()
  
  for (i in 1:length(Gk)){
    C.not.centers[[i]] <- not.centers[which(pairing.vector == i)]
  } 
  
  C <- list()
  
  for (i in 1:length(Gk)){
    C[[i]] <- Gk[i]
  }

  for (i in 1:length(Gk)){
    C[[i]] <- append(C[[i]], C.not.centers[[i]])
  }
  
  return(C)
}
 


