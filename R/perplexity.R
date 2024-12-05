#' perplexity
#'
#' @param model the fitted model
#' @param newdata test dataset
#' @param posterior_mean use posterior mean estimates
#'
#' @return the value of the perplexity
#' @export
perplexity <- function(model, newdata=NULL, posterior_mean = FALSE){

  if(is.null(newdata)){

    if(posterior_mean==FALSE){

      # use chain estimates from CGS:

      loglik <- model$loglik

      den <- sum(model$data)

    }else{ # posterior_mean == TRUE

      # use posterior mean for phi and theta:
      phi_post_mean <- as.data.frame(apply(model$phi_post, c(1,2) , mean))
      theta_post_mean <- as.data.frame(apply(model$theta_post, c(1,2) , mean))

      loglik <- c()

      D <- nrow(theta_post_mean)
      K = ncol(theta_post_mean)

      for(d in 1:D){

        temp_w_list <- unlist(sapply(1:ncol(model$data), function(x) rep(x, model$data[d,x])))

        subtotal = 0

        for(w in 1:length(model$z_post[[d]])){

          wid <- temp_w_list[w]

          subtotal = subtotal + log(sum(phi_post_mean[wid,]*theta_post_mean[d,]))

        }
        loglik[d] = subtotal
      }

      den = sum(model$data)

    }

  }else{

    # use newdata:

    newdata <- transform_data(newdata)

    phi_post_mean <- as.data.frame(apply(model$phi_post, c(1,2) , mean))
    theta_post_mean <- as.data.frame(apply(model$theta_post, c(1,2) , mean))
    K = ncol(theta_post_mean)
    z_init_test <- vector(mode="list", length=nrow(newdata))

    for(d in 1:length(z_init_test)){

      temp_w_list <- unlist(sapply(1:ncol(newdata), function(x) rep(x, newdata[d,x])))

      for(w in 1:sum(newdata[d,])){

        wid <- temp_w_list[w]

        z_init_test[[d]][w]  = sample(1:K, 1, replace=F, prob=phi_post_mean[wid,])

        # TO DO: if a new word is passed, casually sample z_init_test from 1:K with prob= 1/K
        # the corresponding phi_w becomes 0

      }

    }

    K = ncol(theta_post_mean)

    niter <- model$niter
    keep_index <- model$keep_index
    alpha <- model$alpha
    V_test <- ncol(newdata)
    beta <- rep(1, V_test)

    method <- model$type_model

    if(method=="LDA"){
      K = ncol(theta_post_mean)
      model <- collapsed_lda_cpp_pred(newdata,
                                      phi_post_mean = as.matrix(phi_post_mean),
                                      alpha=alpha,
                                      beta=beta,
                                      K = K,
                                      niter=niter,
                                      keep_index = keep_index,
                                      z_init=z_init_test,
                                      verbose=0)

      model <- rlist::list.append(model,data=newdata)

    }
    if(method=="EFLDA"){

      tau = model$tau
      p = model$p
      K = ncol(theta_post_mean)

      model <- collapsed_efd_cpp_pred(newdata,
                                        phi_post_mean = as.matrix(phi_post_mean),
                                        alpha=alpha,
                                        beta=beta,
                                        tau = tau,
                                        p = p,
                                        K = K,
                                        niter=niter,
                                        keep_index = keep_index,
                                        z_init=z_init_test,
                                        verbose=0)

      model <- rlist::list.append(model,data=newdata)

    }
    if(method=="EFLDA_Phi"){

      tau = model$tau
      p = model$p
      K = ncol(theta_post_mean)

      model <- collapsed_efd_cpp_Phi_pred(newdata,
                                      phi_post_mean = as.matrix(phi_post_mean),
                                      alpha=alpha,
                                      beta=beta,
                                      tau = tau,
                                      p = p,
                                      K = K,
                                      niter=niter,
                                      keep_index = keep_index,
                                      z_init=z_init_test,
                                      verbose=0)

      model <- rlist::list.append(model,data=newdata)

    }
    if(posterior_mean == TRUE){

      # phi_mean train and theta_mean test

      # use posterior mean for phi and theta:
      theta_post_mean <- as.data.frame(apply(model$theta_post, c(1,2) , mean))
      K = ncol(theta_post_mean)
      loglik <- c()

      D <- nrow(theta_post_mean)

      for(d in 1:D){

        temp_w_list <- unlist(sapply(1:ncol(model$data), function(x) rep(x, model$data[d,x])))

        subtotal = 0

        for(w in 1:length(model$z_post[[d]])){

          wid <- temp_w_list[w]

          subtotal = subtotal + log(sum(phi_post_mean[wid,]*theta_post_mean[d,]))

        }
        loglik[d] = subtotal
      }

    }else{

      # phi_post_mean train and theta_post_(b)_Gibbs test
      loglik <- model$loglik_mean

    }

    den <- sum(newdata)

  }


  return(exp(-sum(loglik)/den))

}
