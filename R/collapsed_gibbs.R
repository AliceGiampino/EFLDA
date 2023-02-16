#' Collapsed Gibbs sampling
#'
#' @param data_DTM data in format DTM, it can be also a matrix with columns Word and Doc, or a data.frame
#' @param K number of topic
#' @param alpha parameter of the Dirichlet
#' @param beta parameter of the Dirichlet
#' @param tau parameter of the EFD changes mean of clusters
#' @param p probabilities vector
#' @param type LDA or EFD
#' @param thin thinning period
#' @param niter number of iteration
#' @param warmup percentage of the warmup
#' @param seed of the analysis
#' @param init.z initial values of the topic
#' @param verbose if you want print of the iteration
#' @param all.post True if you want all the niter samples from the posterior
#'
#' @return list of the results from the gibbs sampling
#' @export
collapsed_gibbs <- function(data_DTM,
                            # data_DTM is a DTM, it can be also a matrix with columns Word and Doc or data.frame
                             K, alpha=NULL, beta=NULL, tau=NULL, p=NULL,type="LDA",
                             thin=1, niter=5000, warmup=0.5, seed=42, init.z=NULL, verbose=T,
                             all.post = F
){
  #if(class(data)[1]!="matrix"& class(data)[1]!="data.frame") data <- DTM_to_matrix(data_DTM)
  if(class(data)[1]=="DTM") data <- DTM_to_matrix(data_DTM)
  if(class(data)[1]=="matrix") data <- as.data.frame(data)

  if(!(type %in% c("LDA","EFD"))) stop("type must be either LDA or EFD.")
  if(type=="EFD" & (is.null(tau)|is.null(p))) stop("You must specify both tau and p.")
  if(sum(p) != 1 & any(p<=0)) stop("Elements of p must be in the (0,1) interval and must suming to 1.")

  if(verbose==T){vb <- 1}else{vb <- 0}

  colnames(data) <- c("Word", "Doc")
  #init.z = NULL
  seed = seed
  N = nrow(data_DTM)
  V = length(unique(data_DTM[,1]))

  if(is.null(init.z)){
    set.seed(seed)
    data$Init_Topic <- as.factor(sample(1:K, nrow(data), T))
  } else {
    data$Init_Topic <- as.factor(init.z)
  }
  if(all.post) {
    keep_index <- seq(1:niter)
  } else {
    keep_index <- seq(niter*warmup+1, niter, by=thin)
  }

  if(is.null(alpha)) alpha <- rep(1, K)
  if(is.null(beta))  beta <- rep(1, V)

  data <- data.matrix(data)
  if(type=="LDA"){
    res <- collapsed_lda_cpp(data=data, K=K, alpha=alpha, beta=beta,
                             thin=thin, niter=niter, warmup=warmup, keep_index = keep_index, verbose = vb)
  } else if(type=="EFD"){
    res <- collapsed_efd_cpp(data=data, K=K, alpha=alpha, beta=beta, tau=tau, p=p,
                             thin=thin, niter=niter, warmup=warmup, keep_index = keep_index, verbose = vb)
  }

  return(res)
}
