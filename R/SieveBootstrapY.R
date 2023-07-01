SieveBootstrapY <-
function(Ys,B=1000,stepsback=2,offset,display=TRUE,norm="Frob",usecor=FALSE){
  n=ncol(Ys)
  p=nrow(Ys)
  eps_mat = matrix(,ncol=n,nrow=p-stepsback)
  exp_mat = matrix(,ncol=n,nrow=p-stepsback)
  rho_mat = matrix(,nrow=n,ncol=stepsback) # last column is one step back.
  INDs=rep(1:stepsback,p-stepsback)+ rep(0:(p-stepsback-1),each=stepsback)
  for(i in 1:n){
    A = matrix(Ys[INDs,i],nrow=p-stepsback,ncol=stepsback,byrow=T)    
    c = Ys[(stepsback+1):p,i]
    rho_mat[i,]=solve(t(A)%*%A)%*%t(A)%*%c
    exp_mat[,i]=A%*%rho_mat[i,]
    eps_mat[,i]=c-exp_mat[,i]
  }
  dboot = matrix(,nrow=p,ncol=B)
  for(b in 1:B){
    INDboot=sample(1:(p-stepsback),replace=TRUE)
    Ysboot=rbind(Ys[1:stepsback,],exp_mat + eps_mat[INDboot,])
    dboot[(offset):(p-offset-1),b]=d_all(Ysboot,offset,p-offset-1,norm,usecor)
    if(display){ProgressBar(B,b)}
  }
  return(dboot[(offset):(p-offset-1),])
}
