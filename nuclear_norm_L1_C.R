rm(list=ls())
setwd("C:/Github projects/HD-VECM")
library(urca)
source("source.R")

#Miscellaneous functions
generate_VECM = function(n,t,burnin,A,B){
  
  t_total = t+burnin
  Pi = A%*%t(B)
  Phi_1 = diag(n) + Pi
  eps = matrix(rnorm(t_total*n),t_total,n)
  Y = matrix(0,t_total,n)
  for(s in 2:t_total){
    Y[s,] = Phi_1%*%Y[s-1,] + eps[s,]
  }
  Y = head(Y,t)
  list(Y=Y,Pi=Pi)
}
choose_r = function(vecm_fit,alpha=0.05){
  alphas = c(0.1,0.05,0.01)
  if(!(alpha %in% alphas)){
    stop("invalid significance level chosen")
  }
  ind = which(alphas==alpha)
  n_stats = length(vecm_fit@teststat)
  vecm_test_stats = vecm_fit@teststat[n_stats:1]
  vecm_cvs = vecm_fit@cval[n_stats:1,ind]
  r = min(which(vecm_test_stats<vecm_cvs)) - 1
  r
}
BIC_value = function(Y,A,B){
  log(det(crossprod(diff(Y) - Y[-nrow(Y),]%*%B%*%t(A))/nrow(Y)))
   + log(nrow(Y))*(sum(A!=0)+sum(B!=0))/nrow(Y) 
}
AIC_value = function(Y,A,B){
  log(det(crossprod(diff(Y) - Y[-nrow(Y),]%*%B%*%t(A))/nrow(Y)))
  + 2*(sum(A!=0)+sum(B!=0))
}
trace_R2 = function(A_hat,A){
  
  zero_cols = apply(A_hat,2,function(x) all(x==0))
  if(all(!zero_cols)){
    max_r = ncol(A_hat)
  }else{
    max_r = min(which(zero_cols))-1
  }
  A_hat_sub = A_hat[,1:max_r]
  P_A_hat = A_hat_sub%*%solve(crossprod(A_hat_sub))%*%t(A_hat_sub)
  stat = sum((P_A_hat%*%A)^2)/sum(A^2)
  stat
}

#Generate data
set.seed(7406)
t = 500
my_data = draw_VECM(t)
A = cbind(my_data$A,matrix(0,11,7))
B = cbind(my_data$B,matrix(0,11,7))
Pi = A%*%t(B)
Y = my_data$Y
r = my_data$r
matplot(Y,type="l")
n = nrow(A)

#Initial parameters
colnames(Y) = paste("x",1:n,sep="")
vecm_fit = ca.jo(Y,K=2)
r = choose_r(vecm_fit,0.05)
A_Johan = vecm_fit@W[,1:r]
B_Johan = vecm_fit@V[,1:r]
C_init = C_Johan = A_Johan%*%t(B_Johan)
M_init = matrix(0,n,n)
rho = 1
lambda_NN = 1; lambda_L1 = 0.1

#Functions to update C
prox_nuclear = function(X,lambda){
  X_svd = svd(X)
  d_softt = softt_matrix(X_svd$d,lambda)
  X_new = X_svd$u%*%diag(d_softt)%*%t(X_svd$v)
  list(X_new=X_new,d=d_softt)
}

ADMM_NN_L1 = function(Y,C_init,M_init,lambda_NN,lambda_L1,rho,
                      adaptive_NN = FALSE,
                      eps_abs=1e-3,eps_rel=1e-3,max_iter=1000,
                      print_resid = FALSE){
  
  #Transform data
  t=nrow(Y); n = ncol(Y)
  Y_lag = Y[-t,]
  Y_diff = diff(Y)
  Y_lag_cross = crossprod(Y_lag)/t
  DY_Y_lag = crossprod(Y_diff,Y_lag)/t
  C1_den = solve(Y_lag_cross + diag(rho/2,n))
  C1_cst = C1_den%*%DY_Y_lag
  
  if(adaptive_NN){
    omegas = 1/eigen(crossprod(Y)/t^2)$values
  }else{
    omegas=NULL
  }
  
  #Transform desired accuracies
  eps_abs_scaled = eps_abs*n
  
  #Initialize objects to update
  C1_new = C_init
  C2_new = C_init
  C_new = C_init
  M1_new = M_init
  M2_new = M_init
  
  #Run ADMM algorithm
  counter = 0
  while(TRUE){
    
    counter = counter+1
    C1_old = C_new
    C2_old = C_new
    C_old = C_new
    M1_old = M1_new
    M2_old = M2_new
    
    #Update C1
    U1 = C_old - M1_old/rho
    C1_new = C1_cst + (rho/2)*C1_den%*%U1
    
    #Update C2
    U2 = C_old - M2_old/rho
    U2_SVD = svd(U2)
    d_softt = softt_matrix(U2_SVD$d,lambda_NN/rho,omegas)
    C2_new = U2_SVD$u%*%diag(d_softt)%*%t(U2_SVD$v)
    
    #Update C3
    W_tilde = (C1_new + M1_old/rho + C2_new + M2_old/rho)/2
    C_new = softt_matrix(W_tilde,lambda_L1/(2*rho))
    
    #Update M
    M1_new = M1_old + rho*(C1_new - C1_old)
    M2_new = M2_old + rho*(C2_new - C2_old)
    
    #Compute primal and dual residuals
    dual_resid = rho*(C_old-C_new)
    dual_resid_norm = sqrt(sum(dual_resid^2))
    primal_resid = 2*C_new - C1_new - C2_new
    primal_resid_norm = sqrt(sum(primal_resid^2))
    if(print_resid){
      print(c(primal_resid_norm,dual_resid_norm))
    }
    
    #Check convergence
    C_norms = c(norm(C1_new,"F"),norm(C2_new,"F"),norm(C_new,"F"))
    eps_primal = eps_abs_scaled + eps_rel*max(C_norms)
    M_norms = c(norm(M1_new,"F"),norm(M2_new,"F"))
    eps_dual = eps_abs_scaled + eps_rel*max(M_norms)
    if(primal_resid_norm < eps_primal & dual_resid_norm < eps_dual){
      conv=0
      break
    }
    if(counter == max_iter){
      conv=1
      break
    }
  }
  rank = sum(d_softt!=0)
  output = list(C1=C1_new,C2=C2_new,C=C_new,M1=M1_new,M2=M2_new,
                d=d_softt,rank=rank,
                conv=conv,iterations=counter)
  output
}

test = ADMM_NN_L1(Y,C_init,M_init,lambda_NN = 2,lambda_L1 = 0.0001,
                  rho = 1,adaptive_NN = FALSE,max_iter = 1e5,
                  eps_abs = 1e-5,eps_rel = 1e-5,print_resid = T)
test$rank
c(norm(Pi - C_init,"F"),norm(Pi - test$C1,"F"),norm(Pi - test$C2,"F"),norm(Pi - test$C,"F"))
test$iterations
test$C
test$conv

test_ssvd = rssvd(Y_diff,Y_lag,11)
U = test_ssvd$U
D = test_ssvd$D
V = test_ssvd$V
C_hat = U%*%diag(D)%*%t(V)
norm(Pi - C_hat,"F")
