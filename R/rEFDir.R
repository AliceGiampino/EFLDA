#' Random realization from extended flexible Dirichlet
#'
#' @describeIn The function generates random values from the extended flexible Dirichlet (EFD) distribution.
#'
#' @param n the number of values to generate.
#' @param alpha a vector with positive elements.
#' @param p a vector with positive elements summing to 1.
#' @param tau a vector with positive elements.
#' @param label logical. If TRUE, an additional column pointing the component of the mixture generating the observation is returned.
#'
#' @return A matrix with \code{n} rows, where each row represents a sample from the EFD.
#'
#' @example
#'
#' n <- 4
#' alpha <- c(4, 6, 2)
#' p <- c(.2, .3, .5)
#' tau <- c(7, 2, 6)
#' rEFDir(n,alpha,p,tau,label=FALSE)
#'
#' @references{
#'  Ongaro, A., Migliorati, S., Ascari, R. (2020). A new mixture model on the simplex. Statistics and Computing, \bold{30}, 749-770, doi: 10.1007/s11222-019-09920-x.
#' }
#'
#' @import stats
#'
#' @export
#'
rEFDir <- function(n,alpha,p,tau,label=FALSE){
  D<-length(alpha)
  multin<-matrix(stats::rmultinom(n,1,p), ncol=D,byrow=TRUE)
  x<-matrix(NA, nrow=n, ncol=(D+1),byrow=TRUE)
  if(n==1){
    for (j in 1: D){
      if (multin[1,j]==1) {x[1,j]<-stats::rgamma(1,alpha[j]+tau[j]);x[1,D+1]<-j} else x[1,j]<-stats::rgamma(1,alpha[j])
    }
    somma<-sum(x[,1:D])
    if(label) c(x[1,D+1],x[1,1:D]/(somma)) else cbind(x[1,1:D]/(somma))

  } else{
    for (i in 1:n){ for (j in 1: D){
      if (multin[i,j]==1) {x[i,j]<-stats::rgamma(1,alpha[j]+tau[j]);x[i,D+1]<-j} else x[i,j]<-stats::rgamma(1,alpha[j])
    }}
    somma<-apply(x[,1:D],1,sum)
    if(label) cbind(x[,D+1],x[,1:D]/as.vector(somma)) else cbind(x[,1:D]/as.vector(somma))
  }

}
