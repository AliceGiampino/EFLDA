#' Corpus from DTM
#'
#' @param corpus matrix with column Doc and Word
#'
#' @return DTM
#' @export
corpus_DTM <- function(corpus){
  D <- length(corpus)
  V <- ncol(corpus[[1]]$doc)

  DTM <- matrix(NA, ncol=V, nrow=D)
  #colnames(DTM) <- colnames(corpus[[1]]$doc)

  for(d in 1:D) DTM[d,] <- colSums(corpus[[d]]$doc)

  class(DTM) <- "DTM"

  return(DTM)
}
