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
A_init = A_Johan = vecm_fit@W
B_init = B_Johan = vecm_fit@V
C_Johan = A_Johan%*%t(B_Johan)
MA_init = MB_init = matrix(0,n,n)
rho = 1
lambda_L2 = 1; lambda_L1 = 0.1
omegas = NULL

#Functions to update C
ADMM_L2_L1 = function(Y,A_init,B_init,MA_init,MB_init,
                      lambda_L2,lambda_L1,
                      rho,omegas=NULL,
                      eps_abs=1e-3,eps_rel=1e-3,max_iter=1000,
                      print_resid = FALSE,adaptive_rho=T){
  
  #Transform data
  t=nrow(Y); n = ncol(Y)
  Y_lag = Y[-t,]
  Y_diff = diff(Y)
  Y_lag_cross = crossprod(Y_lag)/t
  DY_Y_lag = crossprod(Y_diff,Y_lag)/t

  #Set omegas if needed
  if(is.null(omegas)){
    omegas = 1/eigen(crossprod(Y)/t^2)$values
  }
  
  #Transform desired accuracies
  eps_abs_scaled = eps_abs*n
  
  #Initialize objects to update
  A1_new = A_new = A_init
  B1_new = B_new = B_init
  C_new = rbind(A_init,B_init)
  MA_new = MA_init
  MB_new = MB_init

  #Run ADMM algorithm
  counter = 0
  while(TRUE){
    counter = counter+1
    
    #rho quantities
    lambda_L1_div_rho = lambda_L1/rho
    rho_div_2 = diag(rho/2,n)
    rho_div_2_kron = diag(rho/2,n^2)
    
    #Update old parameters
    A1_old = A1_new
    B1_old = B1_new
    A_old = A_new
    B_old = B_new
    C_old = C_new
    MA_old = MA_new
    MB_old = MB_new

    #Update A1
    UA = A_old - MA_old/rho
    A1_new = (DY_Y_lag%*%B1_old + (rho/2)*UA)%*%
      solve(t(B1_old)%*%Y_lag_cross%*%B1_old +rho_div_2)
    A1_new[abs(A1_new)<1e-10] = 0
    
    #Update B1
    UB = B_old - MB_old/rho
    B1_new = matrix(solve(kronecker(crossprod(A1_new),Y_lag_cross) + rho_div_2_kron)%*%
                      c(t(DY_Y_lag)%*%A1_new + (rho/2)*UB),n,n)
    B1_new[abs(B1_new)<1e-10] = 0
    
    #Update A and B jointly
    C1_new = rbind(A1_new,B1_new)
    M_old = rbind(MA_old,MB_old)
    WC = C1_new + M_old/rho
    for(j in 1:n){
      w_j = WC[,j]
      c_softt_j = softt_matrix(w_j,lambda_L1_div_rho)
      c_mult = (1-(lambda_L2*omegas[j]/rho)/sqrt(sum(c_softt_j^2)))
      if(c_mult>0){
        C_new[,j] = c_mult*c_softt_j
      }else{
        C_new[,j] = 0
      }
    }
    A_new = C_new[1:n,]
    B_new = C_new[(n+1):(2*n),]
    
    #Update M
    MA_new = MA_old + rho*(A1_new - A_new)
    MB_new = MB_old + rho*(B1_new - B_new)
    M_new = rbind(MA_new,MB_new)
    
    #Compute primal and dual residuals
    #dual_resid = rho*(C_old-C_new)
    dual_resid = rho*(C_new-C_old)
    dual_resid_norm = sqrt(sum(dual_resid^2))
    primal_resid = C1_new - C_new
    primal_resid_norm = sqrt(sum(primal_resid^2))
    if(print_resid){
      print(c(primal_resid_norm,dual_resid_norm))
    }
    
    #Check convergence
    C_norms = c(norm(C1_new,"F"),norm(C_new,"F"))
    eps_primal = eps_abs_scaled + eps_rel*max(C_norms)
    M_norm = norm(M_new,"F")
    eps_dual = eps_abs_scaled + eps_rel*M_norm
    
    if(primal_resid_norm < eps_primal & dual_resid_norm < eps_dual){
      conv=0
      break
    }
    if(counter == max_iter){
      conv=1
      break
    }
    
    #Update rho
    if(adaptive_rho){
      if(primal_resid_norm > 10*dual_resid_norm){
        rho = 2*rho
      }else if(dual_resid_norm > 10*primal_resid_norm){
        rho = rho/2
      }
    }
  }
  
  rank = sum(apply(C_new,2,function(x) !all(x==0)))
  output = list(A=A_new,B=B_new,
                A1=A1_new,B1 = B1_new,rank=rank,
                M=M_new,conv=conv,iterations=counter)
  output
}

test = ADMM_L2_L1(Y,A_init,B_init,MA_init,MB_init,
                  lambda_L2=0.01,lambda_L1=0.01,
                  rho=10,omegas=NULL,
                  eps_abs=1e-3,eps_rel=1e-3,max_iter=1e5,
                  print_resid = T,adaptive_rho=T)

test$conv
test$rank
C_hat = test$A%*%t(test$B)
matplot(cbind(c(Pi),c(C_hat),c(C_Johan)),type="l")
c(norm(Pi - C_Johan,"F"),norm(Pi - C_hat,"F"))
test$iterations

test_ssvd = rssvd(Y_diff,Y_lag,11)
U = test_ssvd$U
D = test_ssvd$D
V = test_ssvd$V
C_hat = U%*%diag(D)%*%t(V)
norm(Pi - C_hat,"F")
