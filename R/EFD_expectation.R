E.EFD<-function(a,p_vect,t){
  aplus<-sum(a)
  a*sum(p_vect/(aplus+t))+t*(p_vect/(aplus+t))
}
