// I only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List collapsed_lda_cpp(arma::mat data,
                             arma::colvec alpha,
                             arma::colvec beta,
                             int K,
                             int thin,
                             int niter,
                             double warmup,
                             arma::colvec keep_index,
                             int verbose){

  // data is a matrix with columns Word, Doc, Init_Topic (0, 1, 2)
  // alpha is the vector of alpha for the Dirichlet distr over topics
  // beta is the vector of mass for the Dirichlet distr over documents
  // K number of topics
  // thin is the thinning period
  // niter is the number of iterations
  // warmup is the percentage of burning
  // keep_index is the vector of indexes to use for the chains

  arma::colvec wuni = unique(data.col(0));
  arma::colvec duni = unique(data.col(1));
  int V = wuni.size();
  int D = duni.size();
  arma::colvec num = data.col(0);
  int N = num.size();

  if(alpha.has_nan()){for(int k=0; k<K; k++) alpha(k) = 1;}
  if(beta.has_nan()){for(int v=0; v<V; v++) beta(v) = 1;}

  //double aplus = accu(alpha);
  double bplus = accu(beta);

  // Count matrices:
  // number of words assigned to topic k (col) in document d (row)
  arma::mat n_d_k(D, K);
  // number of times word w is assigned to topic k
  arma::mat n_k_w(K, V);
  // total number of times any word is assigned to topic k
  arma::colvec n_k(K);

  // number of words assigned to topic k in document d
  for(int d=0; d<D; d++){
    arma::uvec loc_d = find(data.col(1) == d+1);
    for(int k=0; k<K; k++){
      arma::colvec temp = data.col(2);
      n_d_k(d,k) = accu(temp(loc_d)==k+1);
    }
  }

  // number of times word w is assigned to topic k
  for(int w=0; w<V; w++){
    arma::uvec sub = find(data.col(0) == w+1);
    for(int k=0; k<K; k++){
      arma::colvec tempw = data.col(2);
      n_k_w(k,w) = accu(tempw(sub)==k+1);
      n_k(k) = accu(tempw==k+1);
    }
  }

  arma::colvec z = data.col(2);
  int n_post = keep_index.size();
  arma::mat z_post(n_post, N);

  for(int iter=0; iter<niter; iter++){

    if (verbose > 0) Rprintf("Iteration number %d over %i \n", iter+1, niter);

    for(int i=0; i<N; i++){
      int word = data.col(0)(i);
      int topic = z(i);
      int d = data.col(1)(i);

      // Decrease:
      n_d_k(d-1,topic-1) = n_d_k(d-1,topic-1) - 1;
      n_k_w(topic-1,word-1) = n_k_w(topic-1,word-1) - 1;
      n_k(topic-1) = n_k(topic-1) - 1;

      arma::rowvec ndk = n_d_k.row(d-1);
      arma::colvec nkw = n_k_w.col(word-1);
      arma::colvec full_cond(K);

      for(int k=0; k<K; k++){
        full_cond(k)  = (ndk(k)+alpha(k))*(nkw(k) + beta(word-1))/(n_k(k)+bplus);
      }
      full_cond = full_cond/accu(full_cond);
      arma::uvec opts = arma::linspace<arma::uvec>(1, K, K);
      topic = Rcpp::RcppArmadillo::sample(opts, 1, 1, full_cond)[0];
      z(i) = topic;

      // Increase:
      n_d_k(d-1,topic-1) = n_d_k(d-1,topic-1) + 1;
      n_k_w(topic-1,word-1) = n_k_w(topic-1,word-1) + 1;
      n_k(topic-1) = n_k(topic-1) + 1;
    }


    arma::rowvec rowz = arma::conv_to<arma::rowvec>::from(z);
    for (int j=0; j<n_post; j++) {
      if (keep_index(j) == iter+1) {

        z_post.row(j) = rowz;
      }

    }
  }
  // phi & theta estimation:
  arma::cube phi_post(V,K,n_post);
  arma::cube theta_post(D,K,n_post);

  arma::mat temp_phi(V,K);
  arma::mat temp_theta(D,K);
  arma::mat n_d_k_Post(D,K);
  arma::mat n_k_w_Post(K,V);
  arma::colvec n_k_Post(K);

  for(int j=0; j<n_post; j++){

    arma::rowvec zp = z_post.row(j);

    // number of words assigned to topic k in document d
    for(int d=0; d<D; d++){
      arma::uvec sub = find(data.col(1) == d+1);
      for(int k=0; k<K; k++){
        n_d_k_Post(d,k) = accu(zp(sub)==k+1);
      }
    }
    // number of times word w is assigned to topic k
    for(int w=0; w<V; w++){
      arma::uvec sub = find(data.col(0) == w+1);
      for(int k=0; k<K; k++){
        n_k_w_Post(k,w) = accu(zp(sub)==k+1);
      }
    }
    // total number of times any word is assigned to topic k
    for(int k=0; k<K; k++){
      n_k_Post(k) = accu(z==k+1);
    }

    for(int d=0; d<D; d++){

      for(int k=0; k<K; k++){
        theta_post(d,k,j) = (n_d_k_Post(d,k) + alpha(k))/accu(n_d_k_Post.row(d)+alpha(k));
      }
    }

    for(int k=0; k<K; k++){
      for(int w=0; w<V; w++){
        phi_post(w,k,j) = (n_k_w_Post(k,w) + beta(w))/accu(n_k_w_Post.row(k)+beta(w));
      }
    }
  }


  return Rcpp::List::create(Rcpp::Named("z_post")=z_post,
                            Rcpp::Named("theta_post")=theta_post,
                            Rcpp::Named("phi_post")=phi_post,
                            Rcpp::Named("data")=data,
                            Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("type_model")="EFD");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List collapsed_efd_cpp(arma::mat data,
                             arma::colvec alpha,
                             arma::colvec beta,
                             arma::colvec tau,
                             arma::colvec p,
                             int K,
                             int thin,
                             int niter,
                             double warmup,
                             arma::colvec keep_index,
                             int verbose){

  // data is a matrix with columns Word, Doc, Init_Topic (0, 1, 2)
  // alpha is the vector of alpha for the Dirichlet distr over topics
  // beta is the vector of mass for the Dirichlet distr over documents
  // tau and p are the parameter of the extended dirichlet
  // K number of topics
  // thin is the thinning period
  // niter is the number of iterations
  // warmup is the percentage of burning
  // keep_index is the vector of indexes to use for the chains

  arma::colvec wuni = unique(data.col(0));
  arma::colvec duni = unique(data.col(1));
  int V = wuni.size();
  int D = duni.size();
  arma::colvec num = data.col(0);
  int N = num.size();

  if(alpha.has_nan()){for(int k=0; k<K; k++) alpha(k) = 1;}
  if(beta.has_nan()){for(int v=0; v<V; v++) beta(v) = 1;}

  double aplus = accu(alpha);
  double bplus = accu(beta);

  // Count matrices:
  // number of words assigned to topic k (col) in document d (row)
  arma::mat n_d_k(D, K);
  // number of times word w is assigned to topic k
  arma::mat n_k_w(K, V);
  // total number of times any word is assigned to topic k
  arma::colvec n_k(K);

  // number of words assigned to topic k in document d
  for(int d=0; d<D; d++){
    arma::uvec loc_d = find(data.col(1) == d+1);
    for(int k=0; k<K; k++){
      arma::colvec temp = data.col(2);
      n_d_k(d,k) = accu(temp(loc_d)==k+1);
    }
  }

  // number of times word w is assigned to topic k
  for(int w=0; w<V; w++){
    arma::uvec sub = find(data.col(0) == w+1);
    for(int k=0; k<K; k++){
      arma::colvec tempw = data.col(2);
      n_k_w(k,w) = accu(tempw(sub)==k+1);
      n_k(k) = accu(tempw==k+1);
    }
  }

  arma::colvec z = data.col(2);
  int n_post = keep_index.size();
  arma::mat z_post(n_post, N);

  arma::mat A(K,3);

  A.col(0) = alpha;
  A.col(1) = alpha+tau;
  A.col(2) = tau+aplus;

  arma::mat Agamma = lgamma(A);
  arma::colvec logp = log(p);

  arma::mat B(K,3);
  arma::colvec log_t_long(K);
  arma::colvec logsum(2);
  arma::colvec t1(K);

  for(int iter=0; iter<niter; iter++){

    if (verbose > 0) Rprintf("Iteration number %d over %i \n", iter+1, niter);

    for(int i=0; i<N; i++){
      int word = data.col(0)(i);
      int topic = z(i);
      int d = data.col(1)(i);

      // Decrease:
      n_d_k(d-1,topic-1) = n_d_k(d-1,topic-1) - 1;
      n_k_w(topic-1,word-1) = n_k_w(topic-1,word-1) - 1;
      n_k(topic-1) = n_k(topic-1) - 1;

      arma::rowvec ndk = n_d_k.row(d-1);
      arma::colvec nkw = n_k_w.col(word-1);
      arma::colvec full_cond(K);

      arma::colvec rowd = arma::conv_to<arma::colvec>::from(n_d_k.row(d-1));
      B.col(0) = alpha+tau+rowd;
      B.col(1) = alpha+rowd;
      B.col(2) = aplus+tau+1+accu(rowd);

      arma::mat Bgamma = lgamma(B);

      arma::colvec fact_log_temp = logp+Agamma.col(0)-Agamma.col(1)+Agamma.col(2)+Bgamma.col(0)-Bgamma.col(1)-Bgamma.col(2);

      arma::colvec term_tau = (tau/(alpha+rowd));

      for(int k=0; k<K; k++){
        t1(k) = ((ndk(k)+alpha(k))*(nkw(k) + beta(word-1))/(n_k(k)+bplus));
      }

      double logfact = log(sum(exp(fact_log_temp - max(fact_log_temp)))) + max(fact_log_temp);

      logsum(0) = logfact;

      for(int k=0; k<K; k++){

        double second = fact_log_temp(k) + log(term_tau(k));

        logsum(1) = second;

        log_t_long(k) = log(sum(exp(logsum - max(logsum)))) + max(logsum);
      }

      arma::colvec sumlog = log(t1) + log_t_long;

      full_cond = exp(log(t1) + log_t_long - (log(sum(exp(sumlog - max(sumlog)))) + max(sumlog)));
      arma::uvec opts = arma::linspace<arma::uvec>(1, K, K);
      topic = Rcpp::RcppArmadillo::sample(opts, 1, 1, full_cond)[0];
      z(i) = topic;

      // Increase:
      n_d_k(d-1,topic-1) = n_d_k(d-1,topic-1) + 1;
      n_k_w(topic-1,word-1) = n_k_w(topic-1,word-1) + 1;
      n_k(topic-1) = n_k(topic-1) + 1;
    }

    arma::rowvec rowz = arma::conv_to<arma::rowvec>::from(z);
    for (int j=0; j<n_post; j++) {
      if (keep_index(j) == iter+1) {

        z_post.row(j) = rowz;
      }

    }
  }
  // phi & theta estimation:
  arma::cube phi_post(V,K,n_post);
  arma::cube theta_post(D,K,n_post);

  arma::mat temp_phi(V,K);
  arma::mat temp_theta(D,K);
  arma::mat n_d_k_Post(D,K);
  arma::mat n_k_w_Post(K,V);
  arma::colvec n_k_Post(K);

  arma::mat Bpost(K, 3);
  arma::colvec p_star(K);

  for(int j=0; j<n_post; j++){

    arma::rowvec zp = z_post.row(j);

    // number of words assigned to topic k in document d
    for(int d=0; d<D; d++){
      arma::uvec sub = find(data.col(1) == d+1);
      for(int k=0; k<K; k++){
        n_d_k_Post(d,k) = accu(zp(sub)==k+1);
      }
    }
    // number of times word w is assigned to topic k
    for(int w=0; w<V; w++){
      arma::uvec sub = find(data.col(0) == w+1);
      for(int k=0; k<K; k++){
        n_k_w_Post(k,w) = accu(zp(sub)==k+1);
      }
    }
    // total number of times any word is assigned to topic k
    for(int k=0; k<K; k++){
      n_k_Post(k) = accu(z==k+1);
    }

    for(int d=0; d<D; d++){

      arma::colvec rowdp = arma::conv_to<arma::colvec>::from(n_d_k_Post.row(d));
      Bpost.col(0) = alpha+tau+rowdp;
      Bpost.col(1) = alpha+rowdp;
      Bpost.col(2) = aplus+tau+1+accu(rowdp);

      arma::mat BgammaP = lgamma(Bpost);

      p_star = logp+Agamma.col(0)+BgammaP.col(0)+Agamma.col(2)-BgammaP.col(1)-Agamma.col(1)-BgammaP.col(2);
      double ppp = log(sum(exp(p_star - max(p_star)))) + max(p_star);
      arma::colvec p_starxp = exp(p_star - ppp);

      for(int k=0; k<K; k++){

        double one = (n_d_k_Post(d,k) + alpha(k))*accu(p_starxp/(accu(rowdp + alpha)+tau));
        double two = tau(k)*(p_starxp(k)/(accu(rowdp + alpha)+tau(k)));
        theta_post(d,k,j) = one+two;
      }
    }

    for(int k=0; k<K; k++){
      arma::colvec rowkw = arma::conv_to<arma::colvec>::from(n_k_w_Post.row(k));
      for(int w=0; w<V; w++){
        phi_post(w,k,j) = (n_k_w_Post(k,w) + beta(w))/accu(rowkw+beta);
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("z_post")=z_post,
                            Rcpp::Named("theta_post")=theta_post,
                            Rcpp::Named("phi_post")=phi_post,
                            Rcpp::Named("data")=data,
                            Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("tau")=tau,
                            Rcpp::Named("p")=p,
                            Rcpp::Named("type_model")="EFD");
}
