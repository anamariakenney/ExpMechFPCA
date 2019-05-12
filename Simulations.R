################################################################################################################
#---------------------------------------------------------------------------------------------------------------
# R code for:
# simulating functional data
# running exponential mechanism for FPCA
# calculating performance measurements
#---------------------------------------------------------------------------------------------------------------
################################################################################################################

#---------------------------------------------------------------------------------------------------------------
# Load necessary packages
#---------------------------------------------------------------------------------------------------------------
library(fda)
library(rstiefel)
library(fdapace)
library(MASS)
library(stats)

#---------------------------------------------------------------------------------------------------------------
# Main Functions
#---------------------------------------------------------------------------------------------------------------

#----------------------------------------------------
# 2 norm
#----------------------------------------------------
norm_2 <- function(x) sqrt(sum(x^2))

#----------------------------------------------------
# get projection matrix
#----------------------------------------------------
proj <- function(V){
  return(V%*%ginv(t(V)%*%V)%*%t(V))
}

#----------------------------------------------------
# Frobenius norm
#----------------------------------------------------
frobenius <- function(V){
  return(sqrt(sum(diag(t(V)%*%V))))
}

#----------------------------------------------------
# subspace norm measurement
#----------------------------------------------------
subspace <- function(Vtrue,Vest){
  measure = .5*(frobenius(proj(Vtrue)-proj(Vest)))^2
  return(measure)
}

#----------------------------------------------------
# variance ratio measurement
#----------------------------------------------------
propvar <- function(Vest, Vtrue, X){
  measure = (frobenius(proj(Vest)%*%t(X)))^2/(frobenius(proj(Vtrue)%*%t(X)))^2
  return(measure)
}

#----------------------------------------------------
# generate curves
#----------------------------------------------------
gen.data <- function(n,p,t,m,errvar){
  Btrueeigf = create.fourier.basis(rangeval=c(0,1),nbasis=p)
  time = seq(0,1,length.out=t)
  Btrueeigf = eval.basis(time,Btrueeigf)
  
  # Generate coefficients matrix X
  Xtruec = matrix(,nrow=n,ncol=p)
  for(i in 1:n){
    Xtruec[i,] = rnorm(p,1,.1)
    for(j in 1:5){
      Xtruec[i,j] = (1/j^2)%*%Xtruec[i,j]
    }
    for(j in 6:p){
      Xtruec[i,j] = (1/j^2)%*%Xtruec[i,j]
    }
  }
  # Check over variance explained in components
  Xtruec = Xtruec/max(apply(Xtruec,1,norm_2))
  Xtruec[,1] = rep(0,dim(Xtruec)[1])
  Xtruec[,4] = Xtruec[,4]+0.1
  # generate functions using basis and coefficients
  error = rnorm(n*t,0,errvar)
  error = matrix(error,nrow=n,ncol=t)
  Xtruef = Xtruec%*%t(Btrueeigf) + error
  
  # Using fourier as our basis, smooth discrete data, create functional object
  FB = create.fourier.basis(rangeval=c(0,1),nbasis=m)
  mypar = fdPar(FB,2,.000001) # lowest gcv
  FB = smooth.basis(time,t(Xtruef),mypar)
  Xestc = t(coef(FB$fd))
  
  # Generate Lambda inverse
  values = seq(1,m)
  evals = sapply(values, FUN=function(x)  1/x^3)
  eye = diag(evals)
  Lambda_inv = ginv(eye)
  
  # return relevant terms
  ans<-list(Xestc=Xestc, Lambda_inv=Lambda_inv, Xtruef = Xtruef,Xtruec=Xtruec,Btrueeigf=Btrueeigf)
  return(ans)
}

#----------------------------------------------------
# Gibbs sampler (also calls on data generation function) procedure to draw V matrix
#----------------------------------------------------
sample.V <- function(n,p,t,m,k,eps,itermax,errvar){
  
  # Generate data
  Data = gen.data(n,p,t,m,errvar)
  p1=Data$p1
  Xestc = Data$Xestc
  Lambda_inv = Data$Lambda_inv
  Xtruef = Data$Xtruef
  Xtruec = Data$Xtruec
  Btrueeigf = Data$Btrueeigf
  
  # Initialize V
  if(k == 1){
    V0 = rnorm(m)
    V0 = V0/sqrt(sum(V0^2))
  } else{
    B = diag(k)
    V0 <- diag(m)
    V0 <- V0[,1:k]
  }
  F = vector(,length=itermax)
  VT = V0
  ptm = proc.time()
  VNORM=vector(,length=itermax)
  Vncollect = vector()
  for(i in 1:itermax){
    # Update A for BMF
    A =(1/2)*((eps)*(t(Xestc)%*%Xestc)-Lambda_inv)
    
    # Sample V from BMF
    if(k == 1){
      Vn = rbing.vector.gibbs(A = A, x = V0)
    } else{
      Vn = rbing.matrix.gibbs(A=A, B = B, X=V0)
    }
    if(i > itermax-100){
      Vncollect = cbind(Vncollect,Vn)
    }
    # Measurement to check for gibbs convergence: used in Chaudhuri (2013)
    expres = (1/i)*(Vn+VT)
    
    VT = Vn+VT
    V0 = Vn
    F[i] = (1/sqrt(k))*norm(expres,"F")
  }
  comptime = proc.time() - ptm
  print(paste("Comp time: ",comptime))
  # Check if we drew norm 1 columns
  if(k ==1){
    check = norm_2(Vn)
  } else{
    check = apply(Vn, 2, norm_2)
  }
  
  # Compute true components if applying PCA on estimated X coefficient matrix
  Vsvd.est = svd(Xestc)$v
  
  # Compute true components if applying PCA on true X coefficient matrix
  Vsvd.true = svd(Xtruec)$v
  
  # Plots for convergence
  par(mfrow=c(1,1)
  plot(1:itermax,log(F),main = "Chaudhuri measurement",type="l")
  # Return relevant terms
  ans<-list(Xestc=Xestc, Lambda_inv=Lambda_inv, Xtruef = Xtruef,Xtruec=Xtruec,Btrueeigf=Btrueeigf,
            V = Vn, Vsvd.est = Vsvd.est, Vsvd.true = Vsvd.true, F=F, comptime=comptime, check=check,
            Vncollect=Vncollect)
  return(ans)
}

#---------------------------------------------------------------------------------------------------------------
# Parameter specification
#---------------------------------------------------------------------------------------------------------------
n = 500             # sample size
t = 100             # no. of time points
p = 21              # no. of true components
m = 41              # no. of orthogonal basis functions for estimation
k = 1               # no. of components desired
eps = 1             # privacy budget
itermax = 20000     # no. of iterations for Gibbs
errvar = 1          # variance of error added when generating data

#---------------------------------------------------------------------------------------------------------------
# Run sampler and compute performance measurements
#---------------------------------------------------------------------------------------------------------------
Data <- sample.V(n,p,t,m,k,eps,itermax,errvar)
Xtruef = Data$Xtruef
Vnp.est = Data$Vsvd.est[,1:k]
Vnp.true = Data$Vsvd.true
Btrueeigf = Data$Btrueeigf
X = Data$Xestc
Vest = Data$V
Vncollect= Data$Vncollect
PROP = propvar(Vest,Vnp.est,X)
SUB = subspace(Vnp.est,Vest)