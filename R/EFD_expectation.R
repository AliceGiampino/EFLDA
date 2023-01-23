#' Expected value of Extended Flexible Dirichlet
#'
#' @param a paramater of the Dirichlet
#' @param p_vect vector of probabilities
#' @param t parameter that changes the mean of the clusters
#'
#' @return
#' @export
#'
#' @examples
E.EFD<-function(a,p_vect,t){
  aplus<-sum(a)
  a*sum(p_vect/(aplus+t))+t*(p_vect/(aplus+t))
}
