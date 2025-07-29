rm(list=ls())
setwd("C:/Github projects/HD-VECM")
library(urca)
source("source.R")

#Generate data
my_data = draw_VECM(t=500,spec = "HD_sparse",a=-0.6)
Y = my_data$Y
A_true = my_data$A
B_true = my_data$B
Pi_true = A_true%*%t(B_true)
r_true = my_data$r
n = ncol(Y)

#Obtain initial values
colnames(Y) = paste("x",1:n,sep="")
vecm_fit = ca.jo(Y,K=2)
B_init = vecm_fit@V
A_init = vecm_fit@W
r_johan = choose_r(vecm_fit)

#Fit VECM
# rho=100;lambda_L1 = 1e-6;lambda_L2 = 1;lambda_F = 1e-4
# A_init = A_init;B_init = B_init;adaptive = T;max_iter = 1000
rho = 1
my_vecm = VECM_ADMM(Y,rho=rho,lambda_L1 = 1e-2*rho,
                    lambda_L2 = 0.005*rho,lambda_F = 1e-5,
                    A_init = A_init,B_init = B_init,adaptive = T,
                    max_iter = 10000,thresh=1e-3)
r_ADMM = my_vecm$r
c(my_vecm$r,my_vecm$iterations)
Pi_ADMM = my_vecm$A%*%t(my_vecm$B)
Pi_ADMM_restr = my_vecm$A[,1:r_ADMM]%*%t(my_vecm$B[,1:r_ADMM])
Pi_Johan = A_init[,1:r_johan]%*%t(B_init[,1:r_johan])
MSEs = c(sqrt(mean((Pi_ADMM - Pi_true)^2)),
         sqrt(mean((Pi_ADMM_restr - Pi_true)^2)),
         sqrt(mean((Pi_Johan - Pi_true)^2)))
names(MSEs) = c("ADMM","ADMM_restr","Johansen")
MSEs