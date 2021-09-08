samplesize <- function(x, min, max) {
  mean <- (min + max)/2
  if (x < mean){
    print('Sample size is not enough')
  } else if ( mean <= x){
    print('This sample can be used')
  }
}
