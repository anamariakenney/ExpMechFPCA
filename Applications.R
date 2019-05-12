################################################################################################################
#---------------------------------------------------------------------------------------------------------------
# R code for:
# Running FPCA on Berkeley growth and DTI data sets with the exponential mechanism
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
library(kernlab)
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
# project for using Gaussian kernel
#----------------------------------------------------
projection_basis <- function(y, eigenvect, M_integ)
{
    N <- dim(y)[1]
    y_mat_inprods <- matrix(NA, dim(eigenvect)[2],N)
    for (i in 1:N)
    {
        for(j in 1:dim(eigenvect)[2])
        {
            y_mat_inprods[j,i] <- sum(eigenvect[,j]*y[i,])/M_integ
        }
    }
    return(y_mat_inprods)
}

#----------------------------------------------------
# Gibbs sampler procedure to draw V matrix
#----------------------------------------------------
sample.V <- function(Xestc,Lambda_inv,k,eps,itermax){
  m = dim(Xestc)[2]

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
  
  # Plots for convergence
  par(mfrow=c(1,2))
  plot(1:itermax,log(F),main = "Chaudhuri measurement",type="l")
  ans<-list(Xestc=Xestc, Lambda_inv=Lambda_inv,
            V = Vn, Vsvd.est = Vsvd.est, F=F, comptime=comptime, check=check,
            Vncollect=Vncollect)
  return(ans)
}

#---------------------------------------------------------------------------------------------------------------
# Berkeley: load growth data
#---------------------------------------------------------------------------------------------------------------
data(growth)
names(growth)

#---------------------------------------------------------------------------------------------------------------
# create Y
#---------------------------------------------------------------------------------------------------------------
Y = rbind(t(growth$hgtm),t(growth$hgtf))
#---------------------------------------------------------------------------------------------------------------
# Create time grid
#---------------------------------------------------------------------------------------------------------------
T_domain = growth$age

#---------------------------------------------------------------------------------------------------------------
# DTI: load in data
#---------------------------------------------------------------------------------------------------------------
#data("DTI")
#names(DTI)

#---------------------------------------------------------------------------------------------------------------
# create Y
#---------------------------------------------------------------------------------------------------------------
#Y = DTI$cca
#idx = 0
#for(i in 1:dim(Y)[1]){
#  if(length(which(is.na(Y[i,])))>0){
#    idx = c(idx,i)
#  }
#}
#idx
#Y = Y[-idx[-1],]
#length(which(is.na(Y)))

#---------------------------------------------------------------------------------------------------------------
# Create time grid
#---------------------------------------------------------------------------------------------------------------
#T_domain = 1:93

#---------------------------------------------------------------------------------------------------------------
# specify parameters
#---------------------------------------------------------------------------------------------------------------
k = 1              # no. of components to release
eps = 1            # privacy budget  
itermax = 20000    # no. of gibbs iterations
smoothpar = 0.03   # for kernel, select value that requires ~5 eigenvalues to explain >99% variation
thres = 0.99       # threshold for variation explained by basis functions
n = dim(Y)[1]      # sample size
kernel = "gaussian" # kernel to use
#---------------------------------------------------------------------------------------------------------------
# setup kernel, get eigenvalues and vectors for basis representation
#---------------------------------------------------------------------------------------------------------------
M_integ <- length(T_domain)/diff(range(T_domain))
rbfkernel <- rbfdot(sigma = smoothpar)
mat <- kernelMatrix(rbfkernel, T_domain)
kernel_def <- list(vectors = eigen(mat)$vectors*sqrt(M_integ),
                   values = eigen(mat)$values/M_integ)

# get number of significant eigen values
num_eigen <- which((cumsum(kernel_def$values)/sum(kernel_def$values))>thres)[1]
eigen_chosen <- 1:num_eigen

# defintion of the eigenvalues and eigenvectors of the kernel
autoval <- kernel_def$values[eigen_chosen]
autovett <- kernel_def$vectors[,eigen_chosen]

kernelfinal = list(eigenvect = autovett, eigenval = autoval)

#---------------------------------------------------------------------------------------------------------------
# construct lambda and lambda_inv matrix
#---------------------------------------------------------------------------------------------------------------
Lambda_inv = diag(1/kernelfinal$eigenval)

#---------------------------------------------------------------------------------------------------------------
# re-scale Y as needed
#---------------------------------------------------------------------------------------------------------------
Y_imp_c = scale(Y,scale=FALSE)
Y_imp_c = Y_imp_c/max(apply(Y_imp_c,1,norm_2))

#---------------------------------------------------------------------------------------------------------------
# Get X coefficient matrix for Gaussian kernel basis
#---------------------------------------------------------------------------------------------------------------
Xestc=t(projection_basis(Y_imp_c,kernelfinal$eigenvect,M_integ))
max(apply(Xestc,1,norm_2))

Beigf = kernelfinal$eigenvect
Xestf = Xestc%*%t(Beigf)

#---------------------------------------------------------------------------------------------------------------
# start algorithm
#---------------------------------------------------------------------------------------------------------------
data = sample.V(Xestc,Lambda_inv,k,eps,itermax)

Vest = data$V
Vnp.est = data$Vsvd.est[,1:k]
Vfest = t(Vest)%*%t(Beigf)
Vfnp.est = t(Vnp.est)%*%t(Beigf)
Vncollect= data$Vncollect

PROP = propvar(Vest,Vnp.est,Xestc)
SUB = subspace(Vnp.est,Vest)