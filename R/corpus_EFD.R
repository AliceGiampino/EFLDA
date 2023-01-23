#' Corpus from Extended Flexible Dirichlet
#'
#' @param D number of documents
#' @param V number of words in a vocabulary
#' @param K number of topics
#' @param alpha parameter of the Dirichlet
#' @param tau parameter that changes the mean of the clusters
#' @param p vector of probabilities
#' @param beta parameter of the Dirichlet
#' @param eps vector of (non-negative) means for the Poisson
#'
#' @return
#' @export
#'
#' @examples
corpus_EFD <- function(D, V, K, alpha, tau, p, beta, eps){
  Phi <- LearnBayes::rdirichlet(K, beta)
  corpus <- lapply(1:D, function(x){
    N_d <- rpois(1, eps)

    theta_d <- as.vector(rEFDir(1, alpha=alpha, tau=tau, p=p,label=F))

    # doc_d e' una matrice con 0 e 1
    doc_d <- matrix(0, ncol=V, nrow=N_d)
    #colnames(doc_d) <- as.character(1:V)

    # estraggo i topic di ogni parola del documento:
    z_d <- sample(1:K, N_d, prob=theta_d, replace=T)

    for(k in unique(z_d)){
      ww_k <- sample(1:V, sum(z_d==k), prob=Phi[k,], replace=T)
      doc_d[cbind(which(z_d==k), ww_k)] <- 1
    }

    return(list(doc=doc_d, topic_words=z_d, theta_d=theta_d, Phi=Phi, N_d=N_d))
  })
  #lapply(corpus, function(x) x$doc)
  return(corpus)
}
