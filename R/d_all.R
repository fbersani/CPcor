d_all <-
function(Y,str,stp,norm="Frob",usecor=FALSE){
  p=nrow(Y)
  n=ncol(Y)
  ret_v = rep(0,stp-str+1)
  if(!usecor){
    S1 = t(Y[1:str,])%*%Y[1:str,]*(1/str)
    S2 = (t(Y[(str+1):p,])%*%Y[(str+1):p,])*(1/(p-str))
  }
  for(i in 1:(stp-str+1)){
    if(usecor){
      S1=cor(Y[1:(str+i),])
      S2=cor(Y[(str+i+1):p,])
    }else{
      changemat=matrix(t(Y[str+i,]),ncol=1)%*%matrix(Y[str+i,],nrow=1)
      S1=(S1*(str+i-1)+changemat)/(str+i)
      S2=(S2*(p-str-i+1)-changemat)/(p-str-i)      
    }
    difmat=S1-S2
    if(norm=="Frob"){
      ret_v[i]=sum(diag((difmat)%*%difmat))      
    }else if(norm=="Max"){
      diag(difmat)=0
      ret_v[i]=max(abs(difmat))
    }
  }
  return(ret_v)
}
