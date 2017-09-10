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
