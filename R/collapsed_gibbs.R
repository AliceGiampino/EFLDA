#' Collapsed Gibbs sampling
#'
#' @param data_DTM data in format DTM, it can be also a matrix with columns "Word" and "Doc" or "i" and "j", or a data.frame
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
                             all.post = F, data_output=F
){
  # data_DTM is a DTM, it can be also a matrix with columns different words
  # and rows different documents, the values are the counts of each word in each document
  data <- transform_data(data_DTM)

  if(!(type %in% c("LDA","EFD"))) stop("type must be either LDA or EFD.")
  if(type=="EFD" & (is.null(tau)|is.null(p))) stop("You must specify both tau and p.")
  if(sum(p) != 1 & any(p<=0)) stop("Elements of p must be in the (0,1) interval and must suming to 1.")

  if(verbose==T){vb <- 1; nupd <- round(niter/10)}else{vb <- 0; nupd=0}

  seed = seed
  N = sum(data)
  V = ncol(data)
  D = nrow(data)

  if(is.null(init.z)){
    set.seed(seed)
    Init_Topic <- as.factor(sample(1:K, nrow(data), T))
  } else {
    Init_Topic <- as.factor(init.z)
  }

  if(dim(data)[2]==2){data <- cbind(data, Init_Topic)}
  if(dim(data)[2]==3){data[,3] <- Init_Topic}

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
                             niter=niter, keep_index = keep_index,
                             verbose = vb, nupd = nupd)
  } else if(type=="EFD"){
    res <- collapsed_efd_cpp(data=data, K=K, alpha=alpha, beta=beta, tau=tau, p=p,
                             niter=niter, keep_index = keep_index,
                             verbose = vb, nupd = nupd)
  }
  if(data_output==T){

    res <- rlist::list.append(res,data=data)

  }

  return(res)
}
