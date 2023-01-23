rEFDir <- function(n,alpha,p,tau,label=FALSE){
  D<-length(alpha)
  multin<-matrix(rmultinom(n,1,p), ncol=D,byrow=TRUE)
  x<-matrix(NA, nrow=n, ncol=(D+1),byrow=TRUE)
  if(n==1){
    for (j in 1: D){
      if (multin[1,j]==1) {x[1,j]<-rgamma(1,alpha[j]+tau[j]);x[1,D+1]<-j} else x[1,j]<-rgamma(1,alpha[j])
    }
    somma<-sum(x[,1:D])
    if(label) c(x[1,D+1],x[1,1:D]/(somma)) else cbind(x[1,1:D]/(somma))

  } else{
    for (i in 1:n){ for (j in 1: D){
      if (multin[i,j]==1) {x[i,j]<-rgamma(1,alpha[j]+tau[j]);x[i,D+1]<-j} else x[i,j]<-rgamma(1,alpha[j])
    }}
    somma<-apply(x[,1:D],1,sum)
    if(label) cbind(x[,D+1],x[,1:D]/as.vector(somma)) else cbind(x[,1:D]/as.vector(somma))
  }

}
