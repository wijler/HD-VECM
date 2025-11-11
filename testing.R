rm(list=ls())
setwd("C:/Github projects/HD-VECM")
library(urca)
source("source.R")
Rcpp::sourceCpp("source_Rcpp.cpp")

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

#Testing
lambda_L2_grid = exp(seq(log(1e-3),0,length.out=10))
lambda_L1_grid = lambda_L2_grid
lambda_grid = as.matrix(expand.grid(lambda_L1_grid,lambda_L2_grid)[,c(2,1)])


ADMM_L2_L1 = AB_L2_L1_tuned(Y,A_init,B_init,MA_init,MB_init,
                      lambda_grid=lambda_grid,rho=10,omegas=NULL,
                      eps_abs=1e-3,eps_rel=1e-3,max_iter=1e5,
                      adaptive_rho=T,crit="AIC")

ADMM_L2_L1_rcpp = AB_L2_L1_tuned_Rcpp(Y,A_init,B_init,MA_init,MB_init,lambda_grid,rho=10,omegas=-999,
                                eps_abs=1e-3,eps_rel=1e-3,max_iter=1e5,
                                adaptive_rho=T,crit="AIC")
                      



C_ADMM = ADMM_L2_L1$A%*%t(ADMM_L2_L1$B)
C_ADMM_rcpp = ADMM_L2_L1_rcpp$A%*%t(ADMM_L2_L1_rcpp$B)
est_metrics = c(norm(Pi - C_Johan,"F"),norm(Pi - C_ADMM,"F"),norm(Pi - C_ADMM_rcpp,"F"))
est_metrics

# opt_metrics = c(ADMM_L2_L1$conv,ADMM_L2_L1$iterations,ADMM_L2_L1$rank)
# names(opt_metrics) = c("convergence","iterations","rank")
# opt_metrics
                      
# PGD_L2_L1 = VECM_SG_tuned_Rcpp(A=A_init,B=B_init,
#                                      Y=Y,lambda_grid,-999,lambda_R=1e-2,
#                                      step_init=1e-3,step_mult=0.5,
#                                      step_max_iter=100,
#                                      max_iter = 1e7,thresh=1e-5,
#                                      crit="BIC")



#Plot estimates
matplot(cbind(c(Pi),c(C_hat),c(C_Johan)),type="l")


                  