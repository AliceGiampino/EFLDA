#' Collapsed Gibbs sampling
#'
#' @param data_DTM data in format DTM
#' @param K number of topic
#' @param alpha parameter of the Dirichlet
#' @param beta parameter of the Dirichlet
#' @param tau parameter of the EFD changes mean of clusters
#' @param p probabilities vector
#' @param type LDA or EFD or EFD_Phi
#' @param thin thinning period
#' @param niter number of iteration
#' @param warmup percentage of the warmup
#' @param seed of the analysis
#' @param init.z initial values of the topic, it is a list with D elements and each has length Nd
#' @param verbose if you want print of the iteration
#' @param all.post True if you want all the niter samples from the posterior
#' @param data_output True if you want the training data in the output
#'
#' @return list of the results from the gibbs sampling
#' @export
collapsed_gibbs <- function(data_DTM,
                            # data_DTM is a DTM, it can be also a matrix with columns Word and Doc or data.frame
                             K, alpha=NULL, beta=NULL, tau=NULL, p=NULL,type="LDA",
                             thin=1, niter=5000, warmup=0.5, seed=42, init.z=NULL, verbose=T,
                             all.post = F, data_output=T
){
  # data_DTM is a DTM, it can be also a matrix with columns different words
  # and rows different documents, the values are the counts of each word in each document
  data <- transform_data(data_DTM)

  if(!(type %in% c("LDA","EFD", "EFD_Phi"))) stop("type must be either LDA or EFD or EFD_Phi.")
  if(type=="EFD" & (is.null(tau)|is.null(p))) stop("You must specify both tau and p.")
  if(sum(p) != 1 & any(p<=0)) stop("Elements of p must be in the (0,1) interval and must suming to 1.")

  if(verbose==T){vb <- 1; nupd <- round(niter/10)}else{vb <- 0; nupd=0}

  seed = seed
  N = sum(data)
  V = ncol(data)
  D = nrow(data)

  if(is.null(init.z)){
    set.seed(seed)
    z_init <- list()
    for(d in 1:D){

      z_init[[d]] <- sample(1:K, sum(data[d,]), T)

    }
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
                             z_init = z_init,
                             niter=niter, keep_index = keep_index,
                             verbose = vb, nupd = nupd)
  } else if(type=="EFD"){
    res <- collapsed_efd_cpp(data=data, K=K, alpha=alpha, beta=beta, tau=tau, p=p,
                             z_init = z_init, niter=niter, keep_index = keep_index,
                             verbose = vb, nupd = nupd)
  } else if(type=="EFD_Phi"){
    res <- collapsed_efd_cpp_Phi(data=data, K=K, alpha=alpha, beta=beta, tau=tau, p=p,
                             z_init = z_init, niter=niter, keep_index = keep_index,
                             verbose = vb, nupd = nupd)
  }
  if(data_output==T){

    res <- rlist::list.append(res,data=data)

  }

  return(res)
}
