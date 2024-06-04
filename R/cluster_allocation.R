#' cluster allocation
#'
#' @param model is the EFLDA fitted model
#'
#' @return the cluster allocation matrix
#' @export
cluster_allocation <- function(model){

  if(model$type_model == "LDA"){stop("It can be calculated only with EFD prior.")}

  theta_post <- model$theta_post
  D <- dim(theta_post)[1]
  K <- dim(theta_post)[2]
  n_post <- dim(theta_post)[3]
  alpha <- model$alpha
  tau <- model$tau
  p <- model$p

  return(cluster_allocation_cpp(theta_post,
                     D,
                     n_post,
                     K,
                     alpha,
                    tau,
                     p))


}
