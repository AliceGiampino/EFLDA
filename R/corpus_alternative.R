#' Corpus alternative
#'
#' @param corpus a list with: x$doc and x$topic_words
#'
#' @return a matrix with column Word, Doc and Topic
#' @export
Corpus_alternative <- function(corpus){
  N <- sum(unlist(lapply(corpus, function(x) nrow(x$doc))))

  Corpus_alternative <- matrix(NA, ncol=3, nrow=N)
  colnames(Corpus_alternative) <- c("Word", "Doc", "Topic")

  Corpus_alternative[,1] <- unlist(lapply(corpus, function(x) apply(x$doc,1, function(rr) which(rr==1))))
  Corpus_alternative[,2] <- unlist(lapply(1:length(corpus), function(x) rep(x, length(corpus[[x]]$topic_words))))
  Corpus_alternative[,3] <- unlist(lapply(corpus, function(cc) cc$topic_words))

  Corpus_alternative <- as.data.frame(Corpus_alternative)
  return(Corpus_alternative)
}
