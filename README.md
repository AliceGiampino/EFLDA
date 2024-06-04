# EFLDA

<!-- badges: start -->

<!-- badges: end -->

The goal of EFLDA is to provide a package for topic modelling tecnique exploiting LDA and the Extended Flexible LDA.

## Installation

You can install the development version of EFLDA from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AliceGiampino/EFLDA")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# load the package:
library(EFLDA)

# create data:

set.seed(123)
K <- 3 # Number of Topic
V <- 10 # Length of vocabulary
eps <- 5000 # average of words per document; it influences N_d
D <- 25 # Number of documents

beta <- rep(1,V)

# Corpus EFD:

alpha <- c(10,25,15)
tau <- rep(30,3)
p <- c(.5,.3,.2)

if(length(alpha) != K) stop("Invalid dimension for alpha")
if(length(p) != K) stop("Invalid dimension for p")
if(length(tau) != K) stop("Invalid dimension for tau")

Phi <- LearnBayes::rdirichlet(K, beta)
N_d_vec <- rpois(D, eps)
N <- sum(N_d_vec)
theta_d <- matrix(NA, ncol=K, nrow=D)
alloc_real <- c()
for(d in 1:D){

  # EFLDA (FLDA: taus are all equal)
  all <- as.vector(rEFDir(1, alpha=alpha, tau=tau, p=p,label=T))
  theta_d[d,] <- all[-1]
  alloc_real[d] <- all[1]
}

# TRAINING
corpus <-
  lapply(1:D, function(x) {
    
    # doc_d is a matrix with 0 and 1:
    doc_d <-Matrix::Matrix(matrix(0, ncol=V, nrow=N_d_vec[x]), sparse=T)
    
    # extract real topics:
    z_d <- sample(1:K, N_d_vec[x], prob=theta_d[x,], replace=T)
    
    for(k in unique(z_d)){
      ww_k <- sample(1:V, sum(z_d==k), prob=Phi[k,], replace=T)
      doc_d[cbind(which(z_d==k), ww_k)] <- 1
    }
    
    return(list(doc=doc_d, topic_words=z_d))
  })

data_list <- (Corpus_alternative(corpus))[,1:2]

z_init <- list()
for(d in 1:D){
  
  z_init[[d]] <- sample(1:K, sum(data_list$Doc==d), T)
  
}

data_doc_train <- as.data.frame(table(data_list$Doc, data_list$Word))

data_doc_train <- as.matrix(reshape(data_doc_train, idvar = "Var1",
                                    timevar = "Var2", direction = "wide")[,-1])

# Test set
set.seed(123)
D_test <- 5

N_d_vec <- rpois(D_test, eps)
N <- sum(N_d_vec)
theta_d_test <- matrix(NA, ncol=K, nrow=D_test)

for(d in 1:D_test){
  
  # EFLDA (FLDA: taus are all equal)
  theta_d_test[d,] <- as.vector(rEFDir(1, alpha=alpha, tau=tau, p=p,label=F))
  
}

corpus_test <-
  lapply(1:D_test, function(x) {
    # doc_d e' una matrice con 0 e 1
    # doc_d <- matrix(0, ncol=V, nrow=N_d_vec[x])
    doc_d <-Matrix::Matrix(matrix(0, ncol=V, nrow=N_d_vec[x]), sparse=T)
    #colnames(doc_d) <- as.character(1:V)
    
    # estraggo i topic di ogni parola del documento:
    z_d <- sample(1:K, N_d_vec[x], prob=theta_d_test[x,], replace=T)
    
    for(k in unique(z_d)){
      ww_k <- sample(1:V, sum(z_d==k), prob=Phi[k,], replace=T)
      doc_d[cbind(which(z_d==k), ww_k)] <- 1
    }
    
    return(list(doc=doc_d, topic_words=z_d))
  })

data_list_test <- (Corpus_alternative(corpus_test))[,1:2]

data_doc_test <- as.data.frame(table(data_list_test$Doc, data_list_test$Word))

data_doc_test <- as.matrix(reshape(data_doc_test, idvar = "Var1",
                                   timevar = "Var2", direction = "wide")[,-1])

niter = 1000

lda <- collapsed_gibbs(data_doc_train, K = K,
                       alpha=alpha, beta=beta,
                       type="LDA", thin=1, niter=niter,
                       warmup=0.5, seed=42, init.z=z_init,
                       verbose=T,
                       all.post = F, data_output=T)

flda <- collapsed_gibbs(data_doc_train, K = K,
                        alpha=alpha, beta=beta, tau=tau, p=p,
                        type="EFD", thin=1, niter=niter,
                        warmup=0.5, seed=42, init.z=NULL,
                        verbose=T,
                        all.post = F, data_output=T)
                        
                        
# Find the cluster allocation from the FD mixture:

cl_alloc = cluster_allocation(flda)

cl_allocation <- as.numeric(apply(cl_alloc, 2, function(x) names(which.max(table(x)))))

cbind(Document = 1:D, Cluster = cl_allocation, alloc_real = alloc_real)

table(cl_allocation, alloc_real)

# calculate the optimal permutation for label switching:
library(combinat)
# # optional code for L_p distances
#set p
p_norm = 2

perms = permn(1:K)
L_p_vec = rep(0,length(perms))

theta_mean_post <- as.data.frame(apply(flda$theta_post, c(1,2) , mean))

for(i in 1:length(perms)){

  L_p = 0
  for(k in 1:K){
    L_p = L_p + sum((theta_d[,k] - theta_mean_post[,perms[[i]][k]])^p_norm)
  }
  L_p_vec[i] = L_p^(1/p_norm)
}
#find the permutaiton which minimizes label switching
perm = perms[[which(L_p_vec == min(L_p_vec))]]

table(cl_allocation, alloc_real)

# Create a mapping between original labels and new labels
mapping <- setNames(perm, unique(cl_allocation))

# Apply the mapping to relabel cl_allocation
relabeled_cl_allocation <- mapping[as.character(cl_allocation)]

table(relabeled_cl_allocation, alloc_real)                        


# Perplexity on train data:

perplexity(lda) 
perplexity(flda) 

# Perplexity for train data using the posterior means:

perplexity(lda, posterior_mean=TRUE) 
perplexity(flda, posterior_mean=TRUE) 


# Perplexity for new data:

perplexity(lda, newdata=data_doc_test)
perplexity(flda, newdata=data_doc_test)

# Perplexity for new data using the posterior means:

perplexity(lda, newdata=data_doc_test, posterior_mean=TRUE)
perplexity(flda, newdata=data_doc_test, posterior_mean=TRUE)

```
