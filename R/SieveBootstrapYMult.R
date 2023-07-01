SieveBootstrapYMult <-
function(Ys,B=1000,stepsback=2,offset,display=TRUE,norm="Frob",usecor=FALSE){
  n=ncol(Ys[[1]])
  p=nrow(Ys[[1]])
  eps_mat_ls = list()
  exp_mat_ls = list()
  for(ii in 1:length(Ys)){
    ele=Ys[[ii]]
    eps_mat = matrix(,ncol=n,nrow=p-stepsback)
    exp_mat = matrix(,ncol=n,nrow=p-stepsback)
    rho_mat = matrix(,nrow=n,ncol=stepsback) # last column is one step back.
    INDs=rep(1:stepsback,p-stepsback)+ rep(0:(p-stepsback-1),each=stepsback)
    for(i in 1:n){
      A = matrix(ele[INDs,i],nrow=p-stepsback,ncol=stepsback,byrow=T)    
      c = ele[(stepsback+1):p,i]
      rho_mat[i,]=solve(t(A)%*%A)%*%t(A)%*%c
      exp_mat[,i]=A%*%rho_mat[i,]
      eps_mat[,i]=c-exp_mat[,i]
    }
    eps_mat_ls[[ii]]=eps_mat
    exp_mat_ls[[ii]]=exp_mat
  }
  dboot = matrix(,nrow=p,ncol=B)
  for(b in 1:B){
    INDboot=sample(1:(p-stepsback),replace=TRUE)
    Ysboot_ls = list()
    for(ii in 1:length(Ys)){
      ele=Ys[[ii]]
      exp_mat = exp_mat_ls[[ii]]
      eps_mat = eps_mat_ls[[ii]]
      Ysboot_ls[[ii]]=rbind(ele[1:stepsback,],exp_mat + eps_mat[INDboot,])
    }
    dboot[(offset):(p-offset-1),b]=d_allMult(Ysboot_ls,offset,p-offset-1,norm,usecor)
    if(display){ProgressBar(B,b)}
  }
  return(dboot[(offset):(p-offset-1),])
}
