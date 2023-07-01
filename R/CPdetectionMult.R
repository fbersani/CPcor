CPdetectionMult <-
function(Y,B=500,offset=NA,stepsback=2,norm="Frob",usecor=FALSE,display=TRUE){
  n=ncol(Y[[1]])
  p=nrow(Y[[1]])
  N=length(Y)
  if(is.na(offset)){offset=max(2,ceiling(n/N))}
  if(2*offset>=p-2){stop("Offset too large. Not enough data.")}
  Ys=list()
  for(i in 1:length(Y)){
    Ys[[i]]=scale(Y[[i]])/sqrt((p-1)/p)  
  }
  dboot = SieveBootstrapYMult(Ys,B=B,stepsback=stepsback,offset=offset,display=display,norm=norm,usecor=usecor)
  b_mean = apply(dboot,1,mean)
  b_sd = apply(dboot,1,sd)
  dboot = (dboot-b_mean)/b_sd
  maxd_boot=apply(dboot,2,max)
  dreal=d_allMult(Ys,offset,p-offset-1,norm=norm,usecor)
  dreal = (dreal-b_mean)/b_sd
  realmax = max(dreal)
  INDmax = which(dreal==realmax)  
  pval = length(which(maxd_boot>realmax))/B
  return(list("pvalue"=pval,"CPind"=INDmax+offset,"dscaled"=dreal))
  
}
