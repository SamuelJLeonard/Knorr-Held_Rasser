# Given an adjacency matrix 'adj', this function returns a matrix in which entry i, j is the distance between regions i and j.
# Note: The operation %^% is matrix exponentiation.  This operation is part of the expm package.

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

    
