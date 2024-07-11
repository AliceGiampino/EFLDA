// I only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <string>
#include <vector>
using namespace Rcpp;

// whichOne
//
// Gives the index of the element equal to 1.
// [[Rcpp::export]]
int whichOne(IntegerVector x) {
  int index = -1;
  int n = x.size();
  for (int i = 0; i < n; ++i) {
    if (x[i] == 1) {
      if (index == -1) {
        index = i;
      } else {
        Rcpp::Rcout << "There is more than one value equal to 1 in the vector." << std::endl;
      }
    }
  }

  // If the element does not apper:
  if (index == -1) { Rcpp::Rcout << "Error: There is no corresponding element. Vector is " << x << std::endl;}
  // if the element appears return the index:
  return index;
}

// whichIndex
//
// Gives the index of the element if present in the vector.
// [[Rcpp::export]]
int whichIndex(NumericVector x, double val) {
  int index = -1;
  int n = x.size();

  for (int i = 0; i < n; ++i) {
    if (x[i] == val) {
      if (index == -1) {
        index = i;
      } else {
        Rcpp::Rcout << "Error: The value appears multiple times." << std::endl;
      }
    }
  }

  // If the element appears, return the index.
  if (index == -1) { Rcpp::Rcout << "Error: The value does not appear in this vector." << std::endl;}
  return index;
}

// whichMax
//
// Gives the index of the max element of the vector, correspond to which.max() in R
// [[Rcpp::export]]
int whichMax(arma::colvec vec){

  int K = vec.size();
  double max_val = vec[0];
  int max_ind = 0;

  for(int i = 0; i < K; i++){
    if(vec[i] > max_val) {
      max_val = vec[i];
      max_ind = i;
    }
  }

  return (max_ind + 1); // Adding 1 to convert from 0-based to 1-based indexing
}

// oneSampleMultinom
//
// Samples a single vector from a single draw of the Multinomial
// distribution with parameters specified by probabilities.

// [[Rcpp::export]]
IntegerVector oneSampleMultinom(NumericVector probs) {
  int K = probs.size();
  IntegerVector ans(K);
  rmultinom(1, probs.begin(), K, ans.begin());
  return(ans);
}


// whichMultinom
//
// Returns the index of the sampled multinomial vector which is 1.
// [[Rcpp::export]]
int whichMultinom(NumericVector probs, int d = 0, int w = 0) {
  int K = probs.size();
  int out = -1;
  IntegerVector vec = oneSampleMultinom(probs);
  out = whichOne(vec);
  if (out == -1) {
    Rcpp::Rcout << "Error: Internal error. "
                << "(Row, Column) is ("
                << d+1 << ", " << w+1 << ")"<< std::endl;

    Rcpp::Rcout << "Probabilities vector: (";
    for (int i = 0; i < (K - 1); ++i) {
      Rcpp::Rcout << probs[i] << ", ";
    }
    Rcpp::Rcout << probs[(K-1)] << ")" << std::endl;


  }
  return out;
}

// [[Rcpp::export]]
bool isInVector(int value, NumericVector vec) {
  for(int i = 0; i < vec.size(); ++i) {
    if(vec[i] == value) {
      return true;
    }
  }
  return false;
}

// Function to calculate the density of a Dirichlet distribution
// for a given value x and parameter alpha and tau
// [[Rcpp::export]]
double fDir(arma::colvec x, arma::colvec alpha, arma::colvec tau, int position) {

  // Extract dimension of the Dirichlet distribution
  int K = x.size();

  arma::colvec base(K); // Initialize vector with zeros
  base[position-1] = 1.0; // Set value at position 'position' to 1
  arma::colvec tau_k(K);

  for(int k=0; k<K; k++){
    tau_k[k] = tau[k]*base[k];
  }
  // Calculate the density of the Dirichlet distribution
  double density = lgamma(arma::sum(alpha+tau_k)) - arma::sum(lgamma(alpha+tau_k));

  double somma = 0.0;

  for(int k=0; k<K; k++){

    somma = somma+ (alpha[k]+tau_k[k]-1)*log(x[k]);
  }
  density = density+somma;
  return density;
}

// [[Rcpp::export]]
arma::mat cluster_allocation_cpp(arma::cube theta_post,
                             int D,
                             int n_post,
                             int K,
                             arma::colvec alpha,
                             arma::colvec tau,
                             arma::colvec p){

  // cluster allocation:
  arma::mat cl_alloc(n_post, D);

  for(int j = 0; j < n_post; j++){

    arma::mat theta_post_mat = theta_post.slice(j); // DxK

    for(int d = 0; d < D; d++){

      arma::colvec l_q(K);
      arma::rowvec theta_row_d = theta_post_mat.row(d);
      arma::colvec theta_col_d = theta_row_d.t();

      for(int k = 0; k < K; k++){

        l_q[k] = log(p[k]) + fDir(theta_col_d, alpha, tau, k+1);

      } // close k
      int max_ind = whichMax(l_q);

      cl_alloc(j, d) = max_ind; // Transpose the indices
    } // close d

  } // close j

  return cl_alloc;

}


// [[Rcpp::export]]
Rcpp::List collapsed_lda_cpp(NumericMatrix& data,
                             arma::colvec& alpha,
                             arma::colvec& beta,
                             int K,
                             int niter,
                             NumericVector& keep_index,
                             Rcpp::List& z_init,
                             int verbose,
                             int nupd = 0){

  // data is a matrix with columns Word, Doc, Init_Topic (0, 1, 2)
  // alpha is the vector of alpha for the Dirichlet distr over topics
  // beta is the vector of mass for the Dirichlet distr over documents
  // K number of topics
  // thin is the thinning period
  // niter is the number of iterations
  // warmup is the percentage of burning
  // keep_index is the vector of indexes to use for the chains

  int V = data.ncol();
  int D = data.nrow();

  // number of posterior iterations to
  int n_post = keep_index.size();
  // phi & theta estimation:
  arma::cube phi_post(V,K,n_post);
  arma::cube theta_post(D,K,n_post);

  //double aplus = accu(alpha);
  double bplus = accu(beta);

  // COUNT MATRICES INITIALIZATION:
  // number of words assigned to topic k (col) in document d (row)
  arma::mat n_d_k(D, K);
  // number of times word w is assigned to topic k
  arma::mat n_k_w(K, V);
  // total number of times any word is assigned to topic k
  arma::colvec n_k(K);

  for(int d=0; d<D; d++){

    // topic assignments in document d
    arma::vec z_d = z_init[d];

    int count = 0;

    for(int w=0; w<V; w++){

      // counts in cell (d,w)
      int ndw = data(d,w);

      if(ndw >0){

        for(int i=count; i<(ndw+count); i++){

          int initZ = 0;
          // i-th element of z_d, index from 0 to (K-1)
          initZ = z_d(i)-1;
          n_d_k(d, initZ ) += 1;
          n_k_w(initZ, w) += 1;
          n_k(initZ) += 1;

        } // close all ndw
        count += ndw;

      } // close ndw > 0
    } // close word

  } // finish count matrices initialization

  Rcpp::List z_post(D);
  //Rcpp::List z_post_m1(D);

  // initialize print iteration
  if(nupd == 0){
    nupd = (int) (niter / 10);
  }

  Rcpp::List z_init_up = clone(z_init);

  // Gibbs
  for(int iter=0; iter<niter; iter++){

    // print the current completed work
    if(iter == 1){
      Rcpp::Rcout << "Completed:\t" << (iter) << "/" << niter << std::endl;
    }
    if((iter + 1) % nupd == 0){
      Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << niter << std::endl;
    }

    for(int d=0; d<D; d++){

      // topic assignments in document d
      arma::vec z_d = z_init_up[d];
      int count = 0;

      for(int w=0; w<V; w++){

        // Only assign a new topic to word 'w' if it occurs in document 'd'.
        int ndw = data(d,w);

        if (ndw > 0) {

          for(int i=count; i<(ndw+count); i++){

            // topic are indexed 1-k
            // This means that we will have to re-index from 0-(k-1)
            int topic = z_d[i]-1;

            // Decrease:

            n_d_k(d, topic) -= 1;
            n_k_w(topic, w) -= 1;
            n_k(topic) -= 1;

            // 'full_cond' will hold the parameters for the multinomial distribution.
            // It is zero-indexed.
            NumericVector full_cond(K);

            // 'sum' will hold the sum of the full_cond so that we can normalize
            // the multinomial probabilities to sum to 1.
            double sum = 0;

            for (int k = 0; k < K; ++k) {
              double prob = (n_d_k(d, k) + alpha(k));
              prob = prob * (n_k_w(k, w) + beta(w));
              prob = prob / (n_k(k) + bplus);

              full_cond[k] = prob;
              sum += prob;
            }

            // Normalize params so that the probabilities sum to 1.
            for (int k = 0; k < K; ++k) {
              full_cond[k] = full_cond[k] / sum;
            }

            // Receive a zero-indexed topic from whichMultinom
            topic = whichMultinom(full_cond, d, w);
            // Add 1 to it to reindex in 1-k the topic
            z_d[i] = topic+1;

            // Increase:
            n_d_k(d, topic)  += 1;
            n_k_w(topic, w)  += 1;
            n_k(topic) += 1;
          } // close ndw for

          count += ndw;

        } // close data(d,w) > 0

        if (isInVector(iter+1, keep_index)) {

          int j = whichIndex(keep_index, iter+1);

          for(int k=0; k<K; k++){
            phi_post(w,k,j) = (n_k_w(k,w) + beta(w))/accu(n_k_w.row(k)+beta(w));
          }

        } // close keep index

      } // close word for

      if((iter+1)==niter){
        z_post[d] = z_d;
      }

      //if((iter+1)==(niter-1)){
      //  z_post_m1[d] = z_d;
      //}

      if (isInVector(iter+1, keep_index)) {

        int j = whichIndex(keep_index, iter+1);

        for(int k=0; k<K; k++){
          theta_post(d,k,j) = (n_d_k(d,k) + alpha(k))/accu(n_d_k.row(d)+alpha(k));

        }


      } // close keep index

      z_init_up[d] = z_d;

    } // close document for

  } // close iter gibbs

  // loglikelihood:
  arma::colvec loglik(D);

  for(int d=0; d<D; d++){

    float total = 0.0;

    for(int iter=0; iter<niter; iter++){
      if (isInVector(iter+1, keep_index)) {

        int j = whichIndex(keep_index, iter+1);
        int count = 0;

        float sub_tot = 0.0;

        for(int w=0; w<V; w++){

          // Only assign a new topic to word 'w' if it occurs in document 'd'.
          int ndw = data(d,w);

          if (ndw > 0) {

            for(int i=count; i<(ndw+count); i++){

              float somma_k = 0.0;

              for(int k=0; k<K; k++){

                somma_k += phi_post(w,k,j)*theta_post(d,k,j);
              } // close k, sum_k is the sum over all topic

              sub_tot += log(somma_k);
            } //close Nd

          } // close ndw > 0

        } // close word

        total += sub_tot;

      } // close keep_index
    } // closer iter

    loglik[d] = total/keep_index.size();

  } // close d

  // close loglikelihood
  return Rcpp::List::create(Rcpp::Named("z_post")=z_post,
                            Rcpp::Named("theta_post")=theta_post,
                            Rcpp::Named("phi_post")=phi_post,
                            Rcpp::Named("loglik")=loglik,
                            Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("niter")=niter,
                            Rcpp::Named("keep_index")=keep_index,
                            Rcpp::Named("type_model")="LDA");
}


// [[Rcpp::export]]
Rcpp::List collapsed_efd_cpp(NumericMatrix data,
                               arma::colvec& alpha,
                               arma::colvec& beta,
                               arma::colvec& tau,
                               arma::colvec& p,
                               int K,
                               int niter,
                               NumericVector& keep_index,
                               Rcpp::List z_init,
                               int verbose,
                               int nupd = 0){

  // data is a matrix with columns Word, Doc, Init_Topic (0, 1, 2)
  // alpha is the vector of alpha for the Dirichlet distr over topics
  // beta is the vector of mass for the Dirichlet distr over documents
  // K number of topics
  // niter is the number of iterations
  // keep_index is the vector of indexes to use for the chains

  int V = data.ncol();
  int D = data.nrow();

  // number of posterior iterations to
  int n_post = keep_index.size();
  // phi & theta estimation:
  arma::cube phi_post(V,K,n_post);
  arma::cube theta_post(D,K,n_post);

  arma::mat A(K,3);
  double aplus = accu(alpha);

  A.col(0) = alpha;
  A.col(1) = alpha+tau;
  A.col(2) = tau+aplus;

  arma::mat Agamma = lgamma(A);
  arma::colvec logp = log(p);

  arma::mat B(K,3);
  arma::colvec log_t_long(K);
  arma::colvec logsum(2);
  arma::colvec t1(K);

  arma::mat Bpost(K, 3);
  arma::colvec p_star(K);

  double bplus = accu(beta);

  // COUNT MATRICES INITIALIZATION:
  // number of words assigned to topic k (col) in document d (row)
  arma::mat n_d_k(D, K);
  // number of times word w is assigned to topic k
  arma::mat n_k_w(K, V);
  // total number of times any word is assigned to topic k
  arma::colvec n_k(K);

  for(int d=0; d<D; d++){

    // topic assignments in document d
    arma::vec z_d = z_init[d];
    int count = 0;
    for(int w=0; w<V; w++){
      int initZ = 0;

      if(data(d,w)>0){
        // counts in cell (d,w)
        int ndw = data(d,w);

        for(int i=count; i<(ndw+count); i++){

          // i-th element of z_d
          initZ = z_d(i);

          n_d_k(d, (initZ - 1) ) += 1;
          n_k_w((initZ - 1),w) += 1;
          n_k((initZ - 1)) += 1;

        }
        count += ndw;
      }
    }
  } // finish count matrices initialization

  Rcpp::List z_post(D);
  //Rcpp::List z_post_m1(D);

  // initialize print iteration
  if(nupd == 0){
    nupd = (int) (niter / 10);
  }

  Rcpp::List z_init_up = clone(z_init);

  arma::cube ndk_mat(D,K,niter);
  arma::cube nkw_mat(K,V,niter);
  arma::cube nk_mat(K,1,niter);

  // Gibbs
  for(int iter=0; iter<niter; iter++){

    if(verbose > 0){
      // print the current completed work
      if(iter == 1){
        Rcpp::Rcout << "Completed:\t" << (iter) << "/" << niter << std::endl;
      }
      if((iter + 1) % nupd == 0){
        Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << niter << std::endl;
      }
    }

    for(int d=0; d<D; d++){

      // topic assignments in document d
      arma::vec z_d = z_init_up[d];
      int count = 0;

      for(int w=0; w<V; w++){

        // Only assign a new topic to word 'w' if it occurs in document 'd'.
        int ndw = data(d,w);

        if (ndw > 0) {

          for(int i=count; i<(ndw+count); i++){

            int topic = z_d(i)-1;

            // Decrease:
            // topic are indexed 1-k
            // This means that we will have to re-index from 0-(k-1)
            n_d_k(d, topic) -= 1;
            n_k_w(topic, w) -= 1;
            n_k(topic) -= 1;

            NumericVector full_cond(K);
            arma::colvec rowd = arma::conv_to<arma::colvec>::from(n_d_k.row(d));
            B.col(0) = alpha+tau+rowd;
            B.col(1) = alpha+rowd;
            B.col(2) = aplus+tau+1+accu(rowd);

            arma::mat Bgamma = lgamma(B);

            arma::colvec fact_log_temp = logp+Agamma.col(0)-Agamma.col(1)+Agamma.col(2)+Bgamma.col(0)-Bgamma.col(1)-Bgamma.col(2);

            arma::colvec term_tau = (tau/(alpha+rowd));

            for(int k=0; k<K; k++){
              t1(k) = ((n_d_k(d,k)+alpha(k))*(n_k_w(k,w) + beta(w))/(n_k(k)+bplus));
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
            // Receive a zero-indexed topic from whichMultinom
            topic = whichMultinom(full_cond, d, w);
            // Add 1 to it to reindex in 1-k the topic
            z_d[i] = topic+1;

            // Increase:
            n_d_k(d, topic)  += 1;
            n_k_w(topic, w)  += 1;
            n_k(topic) += 1;
          } // close ndw for

          count += ndw;

        } // close data(d,w) > 0

        if (isInVector(iter+1, keep_index)) {

          int j = whichIndex(keep_index, iter+1);

          for(int k=0; k<K; k++){
            phi_post(w,k,j) = (n_k_w(k,w) + beta(w))/accu(n_k_w.row(k)+beta(w));
          }

        } // close keep index

      } // close word for

      if((iter+1)==niter){
        z_post[d] = z_d;
      }

      if (isInVector(iter+1, keep_index)) {

        int j = whichIndex(keep_index, iter+1);

        arma::colvec rowdp = arma::conv_to<arma::colvec>::from(n_d_k.row(d));
        Bpost.col(0) = alpha+tau+rowdp;
        Bpost.col(1) = alpha+rowdp;
        Bpost.col(2) = aplus+tau+1+accu(rowdp);

        arma::mat BgammaP = lgamma(Bpost);

        p_star = logp+Agamma.col(0)+BgammaP.col(0)+Agamma.col(2)-BgammaP.col(1)-Agamma.col(1)-BgammaP.col(2);
        double ppp = log(sum(exp(p_star - max(p_star)))) + max(p_star);
        arma::colvec p_starxp = exp(p_star - ppp);

        for(int k=0; k<K; k++){

          double one = (n_d_k(d,k) + alpha(k))*accu(p_starxp/(accu(rowdp + alpha)+tau));
          double two = tau(k)*(p_starxp(k)/(accu(rowdp + alpha)+tau(k)));
          theta_post(d,k,j) = one+two;

        }

      } // close keep index

      z_init_up[d] = z_d;

      if((iter+1)==niter){
        z_post[d] = z_d;
      }

    } // close document for

    ndk_mat.slice(iter) = n_d_k;
    nkw_mat.slice(iter) = n_k_w;
    nk_mat.slice(iter) = n_k;

  } // close iter for

  // loglikelihood:
  arma::colvec loglik(D);

  for(int d=0; d<D; d++){

    float total = 0.0;

    for(int iter=0; iter<niter; iter++){
      if (isInVector(iter+1, keep_index)) {

        int j = whichIndex(keep_index, iter+1);
        int count = 0;

        float sub_tot = 0.0;

        for(int w=0; w<V; w++){

          // Only assign a new topic to word 'w' if it occurs in document 'd'.
          int ndw = data(d,w);

          if (ndw > 0) {

            for(int i=count; i<(ndw+count); i++){

              float somma_k = 0.0;

              for(int k=0; k<K; k++){

                somma_k += phi_post(w,k,j)*theta_post(d,k,j);
              } // close k, sum_k is the sum over all topic

              sub_tot += log(somma_k);
            } //close Nd

          } // close ndw > 0

        } // close word

        total += sub_tot;

      } // close keep_index
    } // closer iter

    loglik[d] = total/keep_index.size();

  } // close d

  // close loglikelihood

  return Rcpp::List::create(Rcpp::Named("z_post")=z_post,
                            Rcpp::Named("theta_post")=theta_post,
                            Rcpp::Named("phi_post")=phi_post,
                            Rcpp::Named("loglik")=loglik,
                            Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("tau")=tau,
                            Rcpp::Named("p")=p,
                            Rcpp::Named("n_d_k")=ndk_mat,
                            Rcpp::Named("n_k_w")=nkw_mat,
                            Rcpp::Named("n_k")=nk_mat,
                            Rcpp::Named("niter")=niter,
                            Rcpp::Named("keep_index")=keep_index,
                            Rcpp::Named("type_model")="EFLDA");
}


// [[Rcpp::export]]
Rcpp::List collapsed_lda_cpp_pred(NumericMatrix& data,
                                  arma::colvec& alpha,
                                  arma::colvec& beta,
                                  arma::mat phi_post_mean,
                                  int K,
                                  int niter,
                                  NumericVector& keep_index,
                                  Rcpp::List& z_init,
                                  int verbose,
                                  int nupd = 0){

  // data is a matrix with columns Word, Doc, Init_Topic (0, 1, 2)
  // alpha is the vector of alpha for the Dirichlet distr over topics
  // beta is the vector of mass for the Dirichlet distr over documents
  // K number of topics
  // thin is the thinning period
  // niter is the number of iterations
  // warmup is the percentage of burning
  // keep_index is the vector of indexes to use for the chains

  int V = data.ncol();
  int D = data.nrow();

  // number of posterior iterations to
  int n_post = keep_index.size();
  // phi & theta estimation:
  arma::cube phi_post(V,K,n_post);
  arma::cube theta_post(D,K,n_post);

  //double aplus = accu(alpha);
  double bplus = accu(beta);

  // COUNT MATRICES INITIALIZATION:
  // number of words assigned to topic k (col) in document d (row)
  arma::mat n_d_k(D, K);
  // number of times word w is assigned to topic k
  arma::mat n_k_w(K, V);
  // total number of times any word is assigned to topic k
  arma::colvec n_k(K);

  for(int d=0; d<D; d++){

    // topic assignments in document d
    arma::vec z_d = z_init[d];

    int count = 0;

    for(int w=0; w<V; w++){

      // counts in cell (d,w)
      int ndw = data(d,w);

      if(ndw >0){

        for(int i=count; i<(ndw+count); i++){

          int initZ = 0;
          // i-th element of z_d, index from 0 to (K-1)
          initZ = z_d(i)-1;
          n_d_k(d, initZ ) += 1;
          n_k_w(initZ, w) += 1;
          n_k(initZ) += 1;

        } // close all ndw
        count += ndw;

      } // close ndw > 0
    } // close word

  } // finish count matrices initialization

  Rcpp::List z_post(D);

  // initialize print iteration
  if(nupd == 0){
    nupd = (int) (niter / 10);
  }

  Rcpp::List z_init_up = clone(z_init);

  // Gibbs
  for(int iter=0; iter<niter; iter++){

    // print the current completed work
    if(verbose >0){

      if(iter == 1){
        Rcpp::Rcout << "Completed:\t" << (iter) << "/" << niter << std::endl;
      }
      if((iter + 1) % nupd == 0){
        Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << niter << std::endl;
      }

    }

    for(int d=0; d<D; d++){

      // topic assignments in document d
      arma::vec z_d = z_init_up[d];
      int count = 0;

      for(int w=0; w<V; w++){

        // Only assign a new topic to word 'w' if it occurs in document 'd'.
        int ndw = data(d,w);

        if (ndw > 0) {

          for(int i=count; i<(ndw+count); i++){

            // topic are indexed 1-k
            // This means that we will have to re-index from 0-(k-1)
            int topic = z_d[i]-1;

            // Decrease:

            n_d_k(d, topic) -= 1;
            n_k_w(topic, w) -= 1;
            n_k(topic) -= 1;

            // 'full_cond' will hold the parameters for the multinomial distribution.
            // It is zero-indexed.
            NumericVector full_cond(K);

            // 'sum' will hold the sum of the full_cond so that we can normalize
            // the multinomial probabilities to sum to 1.
            double sum = 0;

            for (int k = 0; k < K; ++k) {
              double prob = (n_d_k(d, k) + alpha(k));
              prob = prob * (n_k_w(k, w) + beta(w));
              prob = prob / (n_k(k) + bplus);

              full_cond[k] = prob;
              sum += prob;
            }

            // Normalize params so that the probabilities sum to 1.
            for (int k = 0; k < K; ++k) {
              full_cond[k] = full_cond[k] / sum;
            }

            // Receive a zero-indexed topic from whichMultinom
            topic = whichMultinom(full_cond, d, w);
            // Add 1 to it to reindex in 1-k the topic
            z_d[i] = topic+1;

            // Increase:
            n_d_k(d, topic)  += 1;
            n_k_w(topic, w)  += 1;
            n_k(topic) += 1;
          } // close ndw for

          count += ndw;

        } // close data(d,w) > 0

        if (isInVector(iter+1, keep_index)) {

          int j = whichIndex(keep_index, iter+1);

          for(int k=0; k<K; k++){
            phi_post(w,k,j) = (n_k_w(k,w) + beta(w))/accu(n_k_w.row(k)+beta(w));
          }

        } // close keep index

      } // close word for

      if((iter+1)==niter){
        z_post[d] = z_d;
      }

      if (isInVector(iter+1, keep_index)) {

        int j = whichIndex(keep_index, iter+1);

        for(int k=0; k<K; k++){
          theta_post(d,k,j) = (n_d_k(d,k) + alpha(k))/accu(n_d_k.row(d)+alpha(k));
        }

      } // close keep index

      z_init_up[d] = z_d;

    } // close document for

  } // close iter gibbs


  // loglikelihood:
  arma::colvec loglik(D);
  arma::colvec loglik_mean(D);

  for(int d=0; d<D; d++){

    float total = 0.0;
    float tot_mean = 0.0;

    for(int iter=0; iter<niter; iter++){
      if (isInVector(iter+1, keep_index)) {

        int j = whichIndex(keep_index, iter+1);
        int count = 0;

        float sub_tot = 0.0;
        float sub_tot_mean = 0.0;

        for(int w=0; w<V; w++){

          // Only assign a new topic to word 'w' if it occurs in document 'd'.
          int ndw = data(d,w);

          if (ndw > 0) {

            for(int i=count; i<(ndw+count); i++){

              float somma_k = 0.0;
              float somma_k_mean = 0.0;

              for(int k=0; k<K; k++){

                somma_k += phi_post(w,k,j)*theta_post(d,k,j);
                somma_k_mean += phi_post_mean(w,k)*theta_post(d,k,j);
              } // close k, sum_k is the sum over all topic

              sub_tot += log(somma_k);
              sub_tot_mean += log(somma_k_mean);
            } //close Nd

          } // close ndw > 0

        } // close word

        total += sub_tot;
        tot_mean += sub_tot_mean;

      } // close keep_index
    } // closer iter

    loglik[d] = total/keep_index.size();
    loglik_mean[d] = tot_mean/keep_index.size();

  } // close d

  // close loglikelihood

  return Rcpp::List::create(Rcpp::Named("z_post")=z_post,
                            Rcpp::Named("theta_post")=theta_post,
                            Rcpp::Named("phi_post_mean")=phi_post_mean,
                            Rcpp::Named("phi_post")=phi_post,
                            Rcpp::Named("loglik")=loglik,
                            Rcpp::Named("loglik_mean")=loglik_mean,
                            Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("niter")=niter,
                            Rcpp::Named("type_model")="LDA");
}


// [[Rcpp::export]]
Rcpp::List collapsed_efd_cpp_pred(NumericMatrix data,
                                    arma::colvec& alpha,
                                    arma::colvec& beta,
                                    arma::colvec& tau,
                                    arma::colvec& p,
                                    arma::mat phi_post_mean,
                                    int K,
                                    int niter,
                                    NumericVector& keep_index,
                                    Rcpp::List z_init,
                                    int verbose,
                                    int nupd = 0){

  // data is a matrix with columns Word, Doc, Init_Topic (0, 1, 2)
  // alpha is the vector of alpha for the Dirichlet distr over topics
  // beta is the vector of mass for the Dirichlet distr over documents
  // K number of topics
  // niter is the number of iterations
  // keep_index is the vector of indexes to use for the chains

  int V = data.ncol();
  int D = data.nrow();

  // number of posterior iterations to
  int n_post = keep_index.size();
  // phi & theta estimation:
  arma::cube phi_post(V,K,n_post);
  arma::cube theta_post(D,K,n_post);

  arma::mat A(K,3);
  double aplus = accu(alpha);

  A.col(0) = alpha;
  A.col(1) = alpha+tau;
  A.col(2) = tau+aplus;

  arma::mat Agamma = lgamma(A);
  arma::colvec logp = log(p);

  arma::mat B(K,3);
  arma::colvec log_t_long(K);
  arma::colvec logsum(2);
  arma::colvec t1(K);

  arma::mat Bpost(K, 3);
  arma::colvec p_star(K);

  double bplus = accu(beta);

  // COUNT MATRICES INITIALIZATION:
  // number of words assigned to topic k (col) in document d (row)
  arma::mat n_d_k(D, K);
  // number of times word w is assigned to topic k
  arma::mat n_k_w(K, V);
  // total number of times any word is assigned to topic k
  arma::colvec n_k(K);

  for(int d=0; d<D; d++){

    // topic assignments in document d
    arma::vec z_d = z_init[d];
    int count = 0;
    for(int w=0; w<V; w++){
      int initZ = 0;

      if(data(d,w)>0){
        // counts in cell (d,w)
        int ndw = data(d,w);

        for(int i=count; i<(ndw+count); i++){

          // i-th element of z_d
          initZ = z_d(i);

          n_d_k(d, (initZ - 1) ) += 1;
          n_k_w((initZ - 1),w) += 1;
          n_k((initZ - 1)) += 1;

        }
        count += ndw;
      }
    }
  } // finish count matrices initialization

  Rcpp::List z_post(D);
  //Rcpp::List z_post_m1(D);

  // initialize print iteration
  if(nupd == 0){
    nupd = (int) (niter / 10);
  }

  Rcpp::List z_init_up = clone(z_init);
  arma::cube ndk_mat(D,K,niter);
  arma::cube nkw_mat(K,V,niter);
  arma::cube nk_mat(K,1,niter);
  // Gibbs
  for(int iter=0; iter<niter; iter++){

    if(verbose > 0){
      // print the current completed work
      if(iter == 1){
        Rcpp::Rcout << "Completed:\t" << (iter) << "/" << niter << std::endl;
      }
      if((iter + 1) % nupd == 0){
        Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << niter << std::endl;
      }
    }

    for(int d=0; d<D; d++){

      // topic assignments in document d
      arma::vec z_d = z_init_up[d];
      int count = 0;

      for(int w=0; w<V; w++){

        // Only assign a new topic to word 'w' if it occurs in document 'd'.
        int ndw = data(d,w);

        if (ndw > 0) {

          for(int i=count; i<(ndw+count); i++){

            int topic = z_d(i)-1;

            // Decrease:
            // topic are indexed 1-k
            // This means that we will have to re-index from 0-(k-1)
            n_d_k(d, topic) -= 1;
            n_k_w(topic, w) -= 1;
            n_k(topic) -= 1;

            NumericVector full_cond(K);
            arma::colvec rowd = arma::conv_to<arma::colvec>::from(n_d_k.row(d));
            B.col(0) = alpha+tau+rowd;
            B.col(1) = alpha+rowd;
            B.col(2) = aplus+tau+1+accu(rowd);

            arma::mat Bgamma = lgamma(B);

            arma::colvec fact_log_temp = logp+Agamma.col(0)-Agamma.col(1)+Agamma.col(2)+Bgamma.col(0)-Bgamma.col(1)-Bgamma.col(2);

            arma::colvec term_tau = (tau/(alpha+rowd));

            for(int k=0; k<K; k++){
              t1(k) = ((n_d_k(d,k)+alpha(k))*(n_k_w(k,w) + beta(w))/(n_k(k)+bplus));
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
            // Receive a zero-indexed topic from whichMultinom
            topic = whichMultinom(full_cond, d, w);
            // Add 1 to it to reindex in 1-k the topic
            z_d[i] = topic+1;

            // Increase:
            n_d_k(d, topic)  += 1;
            n_k_w(topic, w)  += 1;
            n_k(topic) += 1;
          } // close ndw for

          count += ndw;

        } // close data(d,w) > 0

        if (isInVector(iter+1, keep_index)) {

          int j = whichIndex(keep_index, iter+1);

          for(int k=0; k<K; k++){
            phi_post(w,k,j) = (n_k_w(k,w) + beta(w))/accu(n_k_w.row(k)+beta(w));
          }

        } // close keep index

      } // close word for

      if((iter+1)==niter){
        z_post[d] = z_d;
      }

      if (isInVector(iter+1, keep_index)) {

        int j = whichIndex(keep_index, iter+1);

        arma::colvec rowdp = arma::conv_to<arma::colvec>::from(n_d_k.row(d));
        Bpost.col(0) = alpha+tau+rowdp;
        Bpost.col(1) = alpha+rowdp;
        Bpost.col(2) = aplus+tau+1+accu(rowdp);

        arma::mat BgammaP = lgamma(Bpost);

        p_star = logp+Agamma.col(0)+BgammaP.col(0)+Agamma.col(2)-BgammaP.col(1)-Agamma.col(1)-BgammaP.col(2);
        double ppp = log(sum(exp(p_star - max(p_star)))) + max(p_star);
        arma::colvec p_starxp = exp(p_star - ppp);

        for(int k=0; k<K; k++){

          double one = (n_d_k(d,k) + alpha(k))*accu(p_starxp/(accu(rowdp + alpha)+tau));
          double two = tau(k)*(p_starxp(k)/(accu(rowdp + alpha)+tau(k)));
          theta_post(d,k,j) = one+two;

        }

      } // close keep index

      z_init_up[d] = z_d;

      if((iter+1)==niter){
        z_post[d] = z_d;
      }

    } // close document for
    ndk_mat.slice(iter) = n_d_k;
    nkw_mat.slice(iter) = n_k_w;
    nk_mat.slice(iter) = n_k;
  } // close iter for


  // loglikelihood:
  arma::colvec loglik(D);
  arma::colvec loglik_mean(D);

  for(int d=0; d<D; d++){

    float total = 0.0;
    float tot_mean = 0.0;

    for(int iter=0; iter<niter; iter++){
      if (isInVector(iter+1, keep_index)) {

        int j = whichIndex(keep_index, iter+1);
        int count = 0;

        float sub_tot = 0.0;
        float sub_tot_mean = 0.0;

        for(int w=0; w<V; w++){

          // Only assign a new topic to word 'w' if it occurs in document 'd'.
          int ndw = data(d,w);

          if (ndw > 0) {

            for(int i=count; i<(ndw+count); i++){

              float somma_k = 0.0;
              float somma_k_mean = 0.0;

              for(int k=0; k<K; k++){

                somma_k += phi_post(w,k,j)*theta_post(d,k,j);
                somma_k_mean += phi_post_mean(w,k)*theta_post(d,k,j);
              } // close k, sum_k is the sum over all topic

              sub_tot += log(somma_k);
              sub_tot_mean += log(somma_k_mean);
            } //close Nd

          } // close ndw > 0

        } // close word

        total += sub_tot;
        tot_mean += sub_tot_mean;

      } // close keep_index
    } // closer iter

    loglik[d] = total/keep_index.size();
    loglik_mean[d] = tot_mean/keep_index.size();

  } // close d

  // close loglikelihood

  return Rcpp::List::create(Rcpp::Named("z_post")=z_post,
                            Rcpp::Named("theta_post")=theta_post,
                            Rcpp::Named("phi_post")=phi_post,
                            Rcpp::Named("phi_post_mean")=phi_post_mean,
                            Rcpp::Named("loglik")=loglik,
                            Rcpp::Named("loglik_mean")=loglik_mean,
                            Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("tau")=tau,
                            Rcpp::Named("p")=p,
                            Rcpp::Named("n_d_k")=ndk_mat,
                            Rcpp::Named("n_k_w")=nkw_mat,
                            Rcpp::Named("n_k")=nk_mat,
                            Rcpp::Named("niter")=niter,
                            Rcpp::Named("keep_index")=keep_index,
                            Rcpp::Named("type_model")="EFLDA");
}
