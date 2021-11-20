library(plgp)
library(randtoolbox)

# function k(X1, X2)
#       X1 and X2 are two sets of data
#       each one can have two dimensions, of shape (N,2)
#       this function will use X1 to compute a meshgrid for 
cov_ortho_gen <- function(X1, X2 = NULL, theta, f.sim, df.sim = NULL, cpara, #theta->gamma
                          inputBounds = NULL, MC.num = 500, nugget = 1e-8){
  
  if(is.null(inputBounds)) inputBounds <- apply(X1, 2, range)#get the maximum or minimum value of x1 in rows
  
  # get k0(X1, X2), which is the covariance matrix with lengthscale d
  # when X2=NULL, it has shape (N,N)
  # when X2 has shape (L,2), it has shape (N,L)
  cov.ori <- covar.sep(X1, X2, d=-1/(4*log(theta)), g=0)
  
  # sampling a grid for integration, MC.num is the number of grid points, denoted as M here
  # xi has the shape of (M,2)
  if(ncol(X1) == 1){
    xi <- matrix(seq(0, 1, length.out = MC.num), ncol = ncol(X1))#seq ->linespace
  }else{
    xi <- sobol(MC.num, ncol(X1))
  }
  
  # resize the grid to the bounds of input
  xi <- xi * matrix(rep(inputBounds[2,] - inputBounds[1,], each = nrow(xi)), ncol = ncol(xi)) + 
    matrix(rep(inputBounds[1,], each = nrow(xi)), ncol = ncol(xi)) 
  
  # get the covariance matrix of xi, k0(xi,xi), of shape (M,M)
  Wm <- covar.sep(xi, d=-1/(4*log(theta)), g=0)
  
  # get the covariance matrix between xi and X1, k0(xi,X1), of shape (M,N)
  # actually, this equals 
  #     wx1 <- covar.sep(xi, X1, d=-1/(4*log(theta)), g=0))
  wx1 <- apply(X1, 1, function(x) covar.sep(matrix(x, ncol = ncol(X1)), xi, d=-1/(4*log(theta)), g=0))
  
  if(is.null(X2)) {
    wx2 <- wx1
  }else{
    # get the covariance matrix between xi and X2, k0(xi,X2), of shape (M,L)
    # actually, this equals 
    #     wx2 <- covar.sep(xi, X2, d=-1/(4*log(theta)), g=0))
    wx2 <- apply(X2, 1, function(x) covar.sep(matrix(x, ncol = ncol(X2)), xi, d=-1/(4*log(theta)), g=0))
  }
  
  # the gradient of f: df/d\theta, of shape (M,3)
  Fm <- df.sim(xi, cpara)
  if(length(cpara) == 1) Fm <- matrix(Fm, ncol = 1)
  
  # find useless data and delete them
  rm.index <- apply(Fm, 2, function(x) all(x==0))
  if(all(rm.index)){
    cov.ortho <- cov.ori
  }else{
    if(length(cpara) == 1){
      FWF_inverse <- 1/(t(Fm) %*% Wm %*% Fm) # t is transform, change row and col
    }else{
      if(any(rm.index)) {
        Fm <- Fm[,!rm.index]
      }
      # H_\theta = grad f(xi) * k0(xi,xi) * grad^T f(xi), of shape (3,3)
      FWF <- t(Fm) %*% Wm %*% Fm                            
      FWF_inverse <- solve(FWF + diag(nugget * diag(FWF)))  # (H_\theta)^-1
    }
    # the final k(X1,X2), of shape (N,N) when X2=NULL and (N,L) when X2 has shape (L,2)
    cov.ortho <- cov.ori - t(wx1) %*% Fm %*% FWF_inverse %*% t(Fm) %*% wx2
  }
  
  return(cov.ortho)
}

# the loss function for optimization
#  L(X,Y;\theta) = [Y-f(X;\theta)]^T * K^-1 * [Y-f(X;\theta)] - 1/2 det(K)
nl <- function(para, X, Y) 
{
  d <- ncol(X)
  theta <- para[1:d]                             
  cpara <- para[(length(para)-d):length(para)]
  n <- length(Y)
  K <- cov_ortho_gen(X, NULL, theta, f.sim, df.sim, cpara) 
  Ki <- solve(K+diag(1e-4,nrow(K)))
  ldetK <- determinant(K+diag(1e-4,nrow(K)), logarithm=TRUE)$modulus
  ll <- - (n/2)*log(t(Y-f.sim(X,cpara)) %*% Ki %*% (Y-f.sim(X,cpara))) - (1/2)*ldetK
  return(-ll)
}

# predict the output at new X, and give the uncertainty as sigma_p
predict.func <- function(out, newdata){
  d <- ncol(X)
  print(d)
  para <- out$par
  theta <- para[1:d]                             
  cpara <- para[(length(para)-d):length(para)]
  print(cpara)
  
  KX <- cov_ortho_gen(newdata, X, theta, f.sim, df.sim, cpara) 
  KXX <- cov_ortho_gen(newdata, NULL, theta, f.sim, df.sim, cpara) 
  K <- cov_ortho_gen(X, NULL, theta, f.sim, df.sim, cpara) 
  Ki <- solve(K+diag(1e-4,nrow(K)))
  mup <- f.sim(newdata,cpara) + KX %*% Ki %*% (Z - f.sim(X,cpara))
  tau2hat <- drop(t(Z-f.sim(X,cpara)) %*% Ki %*% (Z-f.sim(X,cpara)) / nrow(X))
  Sigmap <- diag(tau2hat*(KXX + diag(1e-4, nrow(newdata)) - KX %*% Ki %*% t(KX)))
  
  return(list(mean=mup, sd=sqrt(pmax(Sigmap,0)), model=f.sim(newdata,cpara)))
}

# function f(x,\theta)
f.sim <- function(x, cpara) {
  if(is.null(dim(x))){
    x1 <- x[1]
    x2 <- x[2]
  }else{
    x1 <- x[,1]
    x2 <- x[,2]
  }  
  
  out <- cpara[1] + x1*cpara[2] + x2*cpara[3]
  return(c(out))
}

# function df/d\theta
df.sim <- function(x, cpara) {
  if(is.null(dim(x))){
    return(c(1,x[,1],x[,2]))
  }else{
    return(cbind(1,x))
  }
}

# read data
data.df <- read.csv("raw alpha for exp_data.csv")
X <- as.matrix(data.df[,1:2])
Z <- data.df[,3]

# estimate parameters
out <- optim(c(0.5, 0.5, 0, 0, 0), nl, method="L-BFGS-B", lower=c(1e-6,1e-6,-5,-5,-5), upper=c(1-1e-6, 1-1e-6,5,5,5), X=X, Y=Z) 

# prediction
newdata.df <- read.csv("alpha.csv")
newdata <- as.matrix(newdata.df[,1:2])

predictions <- predict.func(out, newdata)

# output predictions
output.df <- data.frame(newdata, mean=predictions$mean, true = newdata.df[,3], model=predictions$model, upper=predictions$mean+qnorm(0.975)*predictions$sd,
                        lower=predictions$mean+qnorm(0.025)*predictions$sd)

write.csv(output.df, file = "predictions.csv", row.names = FALSE)


