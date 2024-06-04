#' cluster allocation
#'
#' @param model
#'
#' @return the cluster allocation matrix
#' @export
#'
#' @examples
cluster_allocation <- function(model){

  if(model$type == "LDA"){stop("It can be calculated only with EFD prior.")}

  theta_post <- model$theta_post
  D <- dim(theta_post)[1]
  K <- dim(theta_post)[2]
  n_post <- dim(theta_post)[3]
  alpha <- model$alpha
  tau <- model$tau
  p <- model$p

  return(cluster_allocation(arma::cube theta_post,
                     int D,
                     int n_post,
                     int K,
                     arma::colvec alpha,
                     arma::colvec tau,
                     arma::colvec p))


}
