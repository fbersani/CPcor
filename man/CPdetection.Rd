\name{CPdetection}
\alias{CPdetection}
\title{
Change point detection for correlation networks observed over time.
}
\description{
When correlation networks change over time, this is manifested by changes in the correlations between nodes. The sieve bootstrap is used to resample and construct the null distribution for the matrix norm of the difference between the empirical correlation matrix before and after a potential change point. 
}
\usage{
CPdetection(Y, B = 500, offset = NA, stepsback = 1, norm = "Frob", usecor = FALSE, display = TRUE)
}
\arguments{
  \item{Y}{
  The n by T matrix of n node values observed over T time points.
  }
  \item{B}{
Number of bootstrap iterations.
}
  \item{offset}{
Defaults to n. Change point detection is only attempted if 2*offset<T-2
  }
  \item{stepsback}{
Defaults to 1. Number of autoregression terms in the sieve bootstrap.
}
  \item{norm}{
Defaults to "Frob". Choice of matrix norm to use. Must use either the Frobenius norm ("Frob") or the maximum norm ("Max").
}
  \item{usecor}{
Boolean value. Defaults to FALSE. Indicates whether to use the correlation matrix (diagonals are set to 1) or not.
}
  \item{display}{
Boolean valued. Defaults to TRUE. Indicates whether to display the progess bar or not.
}
}
\details{

}
\value{
\item{pvalue}{The p-value for the test of change point detection.}
\item{CPind}{The location (column index) of the most significant change point.}
\item{dscaled}{The list of scaled distances for each column index tested as a change point.}
}
\references{
"Change Point Detection in Correlation Networks." Ian Barnett and J.P. Onnela. Scientific Reports. 2016.
}
\author{
Ian Barnett}

\examples{
n = 10
p = 500
rho=0.5
Sigma1 = diag(n)*(1-rho)+matrix(rho,nrow=n,ncol=n) 
Sigma2t = cbind(diag(floor(n/2))*(1-rho)+matrix(rho,nrow=floor(n/2),ncol=floor(n/2)),matrix(0,nrow=floor(n/2),ncol=n-floor(n/2)))
Sigma2b = cbind(matrix(0,nrow=floor(n/2),ncol=n-floor(n/2)),diag(floor(n/2))*(1-rho)+matrix(rho,nrow=floor(n/2),ncol=floor(n/2)))
Sigma2 = rbind(Sigma2t,Sigma2b)

Y = rbind(rmvnorm(floor(p/2),mean=rep(0,n),sigma=Sigma1),rmvnorm(floor(p/2),mean=rep(0,n),sigma=Sigma2))
cp.null.out=CPdetection(Y,B=200,offset=n,display=FALSE)
cp.null.out$CPind #the location of the change point should be close to p/2
plot(n:(p-n-1),cp.null.out$dscaled,xlab="column index",ylab="change point detection statistic")
}
\keyword{ distribution }