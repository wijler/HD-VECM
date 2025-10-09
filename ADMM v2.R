rm(list=ls())
setwd("C:/Github projects/HD-VECM")
library(urca)

#Miscellaneous functions
choose_r = function(vecm_fit,alpha=0.05){
  alphas = c(0.1,0.05,0.01)
  if(!(alpha %in% alphas)){
    stop("invalid significance level chosen")
  }
  ind = which(alphas==alpha)
  vecm_test_stats = vecm_fit@teststat
  vecm_cvs = vecm_fit@cval[,ind]
  r_max = length(vecm_cvs)
  r = c(r_max:1)[max(which(vecm_test_stats>vecm_cvs))]
  r
}
softt_matrix = function(X,lambda){
  X_sign = sign(X)
  X_tmp = abs(X)-lambda
  X_tmp[X_tmp<0] = 0
  X_sign*X_tmp
}
prox_L2 = function(x,lambda){
  c_mult = 1 - lambda/sqrt(sum(x^2))
  if(c_mult>0){
    x_new = c_mult*x
  }else{
    x_new = rep(0,length(x))
  }
  x_new
}
loss_function = function(Y_D,Y_lag,A,B,lambda_L1,lambda_L2,omegas){
  t = nrow(Y_D)+1
  n = ncol(A)
  l2_loss = sum(sapply(1:n,function(j){
    omegas[j]*sqrt(sum(A[,j:n]^2) + sum(B[,j:n]^2))
  }))
  loss = sum((Y_D - Y_lag%*%B%*%t(A)))^2/t +
    lambda_L1*(sum(abs(A))+sum(abs(B))) +
    lambda_L2*l2_loss
  loss
}
C_prox = function(C,M,AB,YY,DYY,rho,step_size,lambda,max_iter = 1000,thresh=1e-6){
  
  C_new = C
  step_lambda = step_size*lambda
  n_par = nrow(C)^2
  iterations = 0
  while(TRUE){
    C_old = C_new
    nabla_C = 2*(C_old%*%YY - DYY) + M + rho*(C_old - AB)
    C_tmp = C_old - step_size*nabla_C
    C_tmp_svd = svd(C)
    d_softt = softt_matrix(C_tmp_svd$d,step_lambda)
    C_new = C_tmp_svd$u%*%diag(d_softt)%*%t(C_tmp_svd$v)
    
    #Check convergence
    iterations = iterations + 1
    if(norm(C_new-C_old,"F")/nrow(C_old)^2 < thresh){
      convergence=1
      break
    }
    if(iterations > max_iter){
      convergence=0
      break
    }
  }
  list(C = C_new,iterations=iterations,convergence=convergence,d=d_softt)
}
AB_prox = function(A,B,S,rho,lambda,omegas,step_size,max_iter=1000,thresh=1e-6){
  n = ncol(A)
  A_new = A
  B_new = B
  AB_new = A%*%t(B)
  iterations=0
  while(TRUE){
    A_old = A_new
    B_old = B_new
    AB_old = AB_new
    nabla_A = -rho*(S-AB_old)%*%B_old
    nabla_B = -rho*t(S-AB_old)%*%A_old
    for(j in 1:ncol(A_old)){
      c_tmp = c(A_old[,j],B_old[,j]) - step_size*c(nabla_A[,j],nabla_B[,j])
      c_prox = prox_L2(c_tmp,step_size*lambda*omegas[j])
      A_new[,j] = head(c_prox,n)
      B_new[,j] = tail(c_prox,n)
    }
    AB_new = A_new%*%t(B_new)
    
    #Check convergence
    iterations = iterations + 1
    if(norm(AB_new-AB_old,"F")/length(AB_old)^2 < thresh){
      convergence=1
      break
    }
    if(iterations > max_iter){
      convergence=0
      break
    }
  }
  list(A=A_new,B=B_new,AB = AB_new,iterations=iterations,convergence=convergence)
}

#Generate data
set.seed(7406)
n = 10; t = 500
burnin = 50
t_total = t+burnin
alpha = c(-0.5,rep(0,n-1))
beta = c(1,rep(-1,n-1))
A_true = cbind(alpha,matrix(0,n,n-1))
B_true = cbind(beta,matrix(0,n,n-1))
Phi = diag(0,n)
Pi = alpha%*%t(beta)
Phi_1 = diag(n) + Pi + Phi
Phi_2 = -Phi
eps = matrix(rnorm(t_total*n),t_total,n)
Y = matrix(0,t_total,n)
for(s in 3:t_total){
  Y[s,] = Phi_1%*%Y[s-1,] + Phi_2%*%Y[s-2,] + eps[s,]
}
Y = Y[-c(1:burnin),]
matplot(Y,type="l")

#Initial parameters
colnames(Y) = paste("x",1:n,sep="")
vecm_fit = ca.jo(Y,K=2)
r = choose_r(vecm_fit,0.05)
B_new = vecm_fit@V
A_new = vecm_fit@W
AB_new = C_new = A_new%*%t(B_new)
M_new = matrix(0,nrow(A_new),nrow(A_new))

#Pre-calculate quantities
Y_lag = Y[-t,]
Y_diff = diff(Y)
Y_lag_cross = crossprod(Y_lag)/t
DY_Y_lag = crossprod(Y_diff,Y_lag)/t
omegas = 1/eigen(crossprod(Y)/t^2)$values

#ADMM settings
rho = 1; lambda_nuclear = 1e-2; lambda_L2 = 1e-2
counter = 0
while(counter < 1001){
  
  counter = counter+1
  C_old = C_new
  A_old = A_new
  B_old = B_new
  AB_old = AB_new
  M_old = M_new
  
  #Update C
  C_update = C_prox(C_old,M_old,AB_old,YY = Y_lag_cross,DYY = DY_Y_lag,
                 rho = rho,step_size = 1,lambda = lambda_nuclear)
  C_new = C_update$C
  
  #Update A,B
  S = C_new + M_old/rho
  
  A=A_old;B=B_old;lambda=lambda_L2;step_size=1;max_iter=1000;thresh=1e-6
  AB_update = AB_prox(A_old,B_old,S,rho,lambda_L2,omegas,step_size = 1)
}
