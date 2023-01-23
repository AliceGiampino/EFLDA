corpus_DTM <- function(corpus){
  # from a matrix with column Doc and Word to DTM
  D <- length(corpus)
  V <- ncol(corpus[[1]]$doc)

  DTM <- matrix(NA, ncol=V, nrow=D)
  #colnames(DTM) <- colnames(corpus[[1]]$doc)

  for(d in 1:D) DTM[d,] <- colSums(corpus[[d]]$doc)

  return(DTM)
}
