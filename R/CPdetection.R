CPdetection <-
function(Y,B=500,offset=NA,stepsback=1,norm="Frob",usecor=FALSE,display=TRUE){
  if(norm!="Frob" && norm!="Max"){stop("You must use 'Frob' or 'Max' as the norm")}
  if(is.list(Y)){
    return(CPdetectionMult(Y,B,offset,stepsback,norm,usecor,display))
  }
  n=ncol(Y)
  p=nrow(Y)
  if(is.na(offset)){offset=n}
  if(2*offset>=p-2){stop("Offset too large. Not enough data.")}
  Ys = scale(Y)/sqrt((p-1)/p)  
  dboot = SieveBootstrapY(Ys,B=B,stepsback=stepsback,offset=offset,display=display,norm=norm,usecor=usecor)
  b_mean = apply(dboot,1,mean)
  b_sd = apply(dboot,1,sd)
  dboot = (dboot-b_mean)/b_sd
  maxd_boot=apply(dboot,2,max)
  dreal=d_all(Ys,offset,p-offset-1,norm=norm,usecor=usecor)
  dreal = (dreal-b_mean)/b_sd
  realmax = max(dreal)
  INDmax = which(dreal==realmax)  
  pval = length(which(maxd_boot>realmax))/B
  return(list("pvalue"=pval,"CPind"=INDmax+offset,"dscaled"=dreal))
}
