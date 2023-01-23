collapsed_gibbs <- function(data_DTM, # data_DTM is a DTM
                             K, alpha=NULL, beta=NULL, tau=NULL, p=NULL,type="LDA",
                             thin=1, niter=5000, warmup=0.5, seed=42, init.z=NULL, verbose=T,
                             all.post = F# Do you want all the niter samples from the posterior?
){
  data <- DTM_to_matrix(data_DTM)

  if(!(type %in% c("LDA","EFD"))) stop("type must be either LDA or EFD.")
  if(type=="EFD" & (is.null(tau)|is.null(p))) stop("You must specify both tau and p.")
  if(sum(p) != 1 & any(p<=0)) stop("Elements of p must be in the (0,1) interval and must suming to 1.")

  if(verbose==T){vb <- 1}else{vb <- 0}

  colnames(data) <- c("Word", "Doc")
  init.z = NULL
  seed = seed
  N = nrow(data_DTM)
  V = length(unique(data_DTM[,1]))

  if(is.null(init.z)){
    set.seed(seed)
    data$"Init_Topic" <- as.factor(sample(1:K, N, T))
  } else {
    data$"Init_Topic" <- as.factor(init.z)
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
