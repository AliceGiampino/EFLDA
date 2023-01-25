#' DTM to matrix
#'
#' @param corpus_dtm
#'
#' @return a matrix
#' @export
#'
#' @examples
DTM_to_matrix <- function(corpus_dtm){
  # from DTM to a matrix with columns Word and Doc
  N_words <- sum(corpus_dtm)

  data_matrix <- matrix(NA, ncol=2, nrow=N_words)
  colnames(data_matrix) <- c("Word", "Doc")
  count <- 0

  for(d in 1:nrow(corpus_dtm)){
    for(w in 1:ncol(corpus_dtm)){

      if(corpus_dtm[d,w] >0){
        data_matrix[(count+1):(count+corpus_dtm[d,w]), 1] <- w
        data_matrix[(count+1):(count+corpus_dtm[d,w]), 2] <- d

        count <- count + corpus_dtm[d,w]
      }

    }
  }

  return(as.data.frame(data_matrix))
}
