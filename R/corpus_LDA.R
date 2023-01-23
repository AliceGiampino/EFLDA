#' Corpus created by LDA
#'
#' @param D number of documents
#' @param V number of words in the vocabulary
#' @param K number of topics
#' @param alpha parameter of the Dirichlet
#' @param beta parameter of the Dirichlet
#' @param eps vector of (non-negative) means for the Poisson
#'
#' @return a corpus
#' @export
#'
#' @examples

corpus_LDA <- function(D, V, K, alpha, beta, eps){
  Phi <- LearnBayes::rdirichlet(K, beta)
  corpus <-lapply(1:D, function(x){
    N_d <- rpois(1, eps)

    theta_d <- as.vector(LearnBayes::rdirichlet(1, alpha))

    # doc_d is a matrix with 0 and 1
    doc_d <- matrix(0, ncol=V, nrow=N_d)
    #colnames(doc_d) <- as.character(1:V)

    # extract topic for every word in a document:
    z_d <- sample(1:K, N_d, prob=theta_d, replace=T)

    for(k in unique(z_d)){
      ww_k <- sample(1:V, sum(z_d==k), prob=Phi[k,], replace=T)
      doc_d[cbind(which(z_d==k), ww_k)] <- 1
    }

    return(list(doc=doc_d, topic_words=z_d, theta_d=theta_d, Phi=Phi, N_d=N_d))
  })
  return(corpus)
}
