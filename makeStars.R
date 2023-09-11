makeStars = function(x){
  stars = c("****", "***", "**", "*", "ns")
  vec = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i = findInterval(x, vec)
  stars[i]
}