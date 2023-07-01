d_allMult <-
function(Y,str,stp,norm="Frob",usecor=FALSE){
  p=nrow(Y[[1]])
  n=ncol(Y[[1]])
  ret_v = rep(0,stp-str+1)
  S1_v = list()
  S2_v = list()
  if(!usecor){
    for(ii in 1:length(Y)){
      ele=Y[[ii]]
      S1_v[[ii]] = t(ele[1:str,])%*%ele[1:str,]*(1/str)
      S2_v[[ii]] = (t(ele[(str+1):p,])%*%ele[(str+1):p,])*(1/(p-str))    
    }
  }
  for(i in 1:(stp-str+1)){
    S1 = matrix(0,nrow=n,ncol=n)
    S2 = matrix(0,nrow=n,ncol=n)
    for(ii in 1:length(Y)){
      ele = Y[[ii]]
      if(usecor){
        S1_v[[ii]]=cor(ele[1:(str+i),])
        S2_v[[ii]]=cor(ele[(str+i+1):p,])
      }else{
        changemat=matrix(t(ele[str+i,]),ncol=1)%*%matrix(ele[str+i,],nrow=1)
        S1_v[[ii]]=(S1_v[[ii]]*(str+i-1)+changemat)/(str+i)
        S2_v[[ii]]=(S2_v[[ii]]*(p-str-i+1)-changemat)/(p-str-i)              
      }
      S1 = S1+S1_v[[ii]]/length(Y)
      S2 = S2+S2_v[[ii]]/length(Y)
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
