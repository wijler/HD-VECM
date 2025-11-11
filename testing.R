rm(list=ls())
setwd("C:/Github projects/HD-VECM")
library(urca)
library(rrpack)
source("source.R")
Rcpp::sourceCpp("source_Rcpp.cpp")

#Generate data
set.seed(2935)
t = 300
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
r_Johan = choose_r(vecm_fit,0.05)
A_init = A_Johan = vecm_fit@W
B_init = B_Johan = vecm_fit@V
C_init = C_Johan = A_Johan%*%t(B_Johan)
M_init = MA_init = MB_init = matrix(0,n,n)

#Testing
lambda_L2_grid = exp(seq(log(1e-3),0,length.out=10))
lambda_L1_grid = lambda_L2_grid
lambda_grid = as.matrix(expand.grid(lambda_L1_grid,lambda_L2_grid)[,c(2,1)])


rho="auto";rho_init=0.1;rho_mult=0.5;mu=10
step_size = "auto";step_init=1;step_mult=0.5
step_max_iter=100;max_iter = 1000
ADMM_thresh=1e-4;C_thresh=1e-4;AB_thresh=1e-4
ADMM_NN_L2_L1 = AB_NN_L2_L1_tuned(C_init,A_init,B_init,M_init,Y,0.01,lambda_grid,"BIC",
                                  rho="auto",rho_init=0.1,rho_mult=0.5,mu=10,
                                  step_size = "auto",step_init=1,step_mult=0.5,
                                  step_max_iter=100,max_iter = 1000,
                                  ADMM_thresh=1e-4,C_thresh=1e-4,AB_thresh=1e-4)
                            

ADMM_L2_L1_rcpp = AB_L2_L1_tuned_Rcpp(Y,A_init,B_init,MA_init,MB_init,lambda_grid,rho=10,omegas=-999,
                                eps_abs=1e-3,eps_rel=1e-3,max_iter=1e5,
                                adaptive_rho=T,crit="AIC")
                      
SSVD_BIC = rssvd(diff(Y),Y[-nrow(Y),],nrank = 11,ic.type="BIC")
SSVD_AIC = rssvd(diff(Y),Y[-nrow(Y),],nrank = 11,ic.type="AIC")

C_ADMM_AIC = ADMM_L2_L1_rcpp$A%*%t(ADMM_L2_L1_rcpp$B)
ind_BIC = which.min(ADMM_L2_L1_rcpp$BICs)
C_ADMM_BIC = ADMM_L2_L1_rcpp$coefs[[ind_BIC]][,,1]%*%t(ADMM_L2_L1_rcpp$coefs[[ind_BIC]][,,2])
C_SVD_BIC = SSVD_BIC$U%*%diag(SSVD_BIC$D)%*%t(SSVD_BIC$V)
C_SVD_AIC = SSVD_AIC$U%*%diag(SSVD_AIC$D)%*%t(SSVD_AIC$V)
MSEs = c(norm(Pi - C_Johan,"F"),
                norm(Pi - C_ADMM_AIC,"F"),norm(Pi - C_ADMM_BIC,"F"),
                norm(Pi - C_SVD_AIC,"F"),norm(Pi - C_SVD_BIC,"F"))
ranks = c(r_Johan,ADMM_L2_L1_rcpp$r,ADMM_L2_L1_rcpp$ranks[ind_BIC],
          SSVD_AIC$rank,SSVD_BIC$rank)
results = rbind(MSEs,ranks)
colnames(results) = c("Johansen","ADMM_AIC","ADMM_BIC","RSVD_AIC","RSVD_BIC")
rownames(results) = c("MSE","rank")
results

#Plot estimates
matplot(cbind(c(Pi),c(C_ADMM),c(C_SVD)),type="l")


                  