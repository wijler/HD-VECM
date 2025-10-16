rm(list=ls())
setwd("C:/Github projects/HD-VECM")
library(urca)
source("source.R")

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

#Generate data
set.seed(746)
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
B_hat = vecm_fit@V
A_hat = vecm_fit@W
Pi_hat = A_hat%*%t(B_hat)

start_indices = c(rep(1,4),2:n)
local_names = c("LS","L2",paste("G",rep(1:n),sep=""))
As_k = MAs_k = Bs_k = MBs_k = list()
Zeros = matrix(0,n,n)
for(k in 1:length(start_indices)){
  j = start_indices[k]
  As_k[[k]] = A_hat[,j:n]
  MAs_k[[k]] = Zeros[,j:n]
  Bs_k[[k]] = B_hat[,j:n]
  MBs_k[[k]] = Zeros[,j:n]
}
names(As_k) = c("A",paste("A",local_names,sep="_"))
names(MAs_k) = c("MA",paste("MA",local_names,sep="_"))
names(Bs_k) = c("B",paste("B",local_names,sep="_"))
names(MBs_k) = c("MB",paste("MB",local_names,sep="_"))

As_kp1 = As_k
MAs_kp1 = MAs_k
Bs_kp1 = Bs_k
MBs_kp1 = MBs_k

#Compute MSEs for Johansen's method and OLS
A_johansen = cbind(A_hat[,1],matrix(0,n,n-1))
B_johansen = cbind(B_hat[,1],matrix(0,n,n-1))
Pi_johansen = A_johansen%*%t(B_johansen)
MSEs_johansen = c(sqrt(sum((A-A_johansen)^2)),
                  sqrt(sum((B-B_johansen)^2)),
                  sqrt(sum((Pi-Pi_johansen)^2)))
MSEs_OLS = c(sqrt(sum((A-A_hat)^2)),
             sqrt(sum((B-B_hat)^2)),
             sqrt(sum((Pi-Pi_hat)^2)))

#ADMM algorithm
Y_lag = Y[-t,]
Y_diff = diff(Y)
t = t-1
Y_lag_cross = crossprod(Y_lag)/t
DY_Y_lag = crossprod(Y_diff,Y_lag)/t

omegas = 1/eigen(crossprod(Y)/t^2)$values
N_j = c(1,1:n)
rho = 100; lambda_L1 = 1e-8*rho; lambda_L2 = 0.01*rho; lambda_F = 1e-4
I_rhod2 = diag(rho/2,n)
I_rhod2_kron = diag(rho/2,n^2)
MSEs = matrix(0,5,3)
rownames(MSEs) = c("ADMM","ADMM_restr","ADMM_norm","Johansen","OLS")
colnames(MSEs) = c("A","B","Pi")
MSEs[c("Johansen","OLS"),] = rbind(MSEs_johansen,MSEs_OLS)
counter = 0
loss = loss_function(Y_diff,Y_lag,A_hat,B_hat,lambda_L1,lambda_L2,omegas)
while(counter < 10001){
  counter = counter + 1
  
  loss_old = loss
  
  #Update A_LS
  As_kp1$A_LS = (DY_Y_lag%*%Bs_k$B_LS + (rho/2)*As_k$A - (1/2)*MAs_k$MA_LS)%*%
    solve(t(Bs_k$B_LS)%*%Y_lag_cross%*%Bs_k$B_LS + I_rhod2)
  
  #Update B_LS
  Bs_kp1$B_LS = matrix(solve(kronecker(crossprod(As_kp1$A_LS),Y_lag_cross) + I_rhod2_kron,
                           c(t(DY_Y_lag)%*%As_kp1$A_LS + 0.5*(rho*Bs_k$B - MBs_k$MB_LS))),n,n)
  
  #Update A_L2
  As_kp1$A_L2 = (rho*As_kp1$A - MAs_k$MA_L2)/(2*lambda_F + rho)
  
  
  #Update B_L2
  Bs_kp1$B_L2 = (rho*Bs_kp1$B - MBs_k$MB_L2)/(2*lambda_F + rho)
  
  #Update (Aj,Bj), j=1,...,N
  for(j in 1:n){
    omega_j = omegas[j] #Change later
    diff_j = c(cbind(As_k$A[,j:n],Bs_k$B[,j:n]) - cbind(MAs_k[[j+3]],MBs_k[[j+3]])/rho)
    shrink_mult = 1-lambda_L2*omega_j/(rho*sqrt(sum(diff_j^2)))
    if(shrink_mult>0){
      AB_j = matrix(shrink_mult*diff_j,nrow=n)
    }else{
      AB_j = matrix(0,n,2*(n+1-j))
    }
    As_kp1[[j+3]] = AB_j[,1:(n+1-j)]
    Bs_kp1[[j+3]] = AB_j[,(n+2-j):(2*(n+1-j))]
  }
  
  #Update A and B
  
  A_kp1_sums = As_kp1$A_LS + As_kp1$A_L2
  B_kp1_sums = Bs_kp1$B_LS + Bs_kp1$B_L2
  W_j = diag(1/(1:n+2))
  for(j in 1:n){
    A_kp1_sums[,j:n] = A_kp1_sums[,j:n] + As_kp1[[j+3]]
    B_kp1_sums[,j:n] = B_kp1_sums[,j:n] + Bs_kp1[[j+3]]
  }
  As_kp1$A = softt_matrix(A_kp1_sums,lambda_L1/rho)%*%W_j
  Bs_kp1$B = softt_matrix(B_kp1_sums,lambda_L1/rho)%*%W_j
  # As_kp1$A = softt_matrix(A_kp1_sums%*%W_j,lambda_L1/rho)
  # Bs_kp1$B = softt_matrix(B_kp1_sums%*%W_j,lambda_L1/rho)
  
  #As_kp1$A; Bs_kp1$B
  
  #Update MA and MB
  MAs_kp1$MA_LS = MAs_k$MA_LS + rho*(As_kp1$A_LS - As_kp1$A)
  MAs_kp1$MA_L2 = MAs_k$MA_L2 + rho*(As_kp1$A_L2 - As_kp1$A)
  MBs_kp1$MB_LS = MBs_k$MB_LS + rho*(Bs_kp1$B_LS - Bs_kp1$B)
  MBs_kp1$MB_L2 = MBs_k$MB_L2 + rho*(Bs_kp1$B_L2 - Bs_kp1$B)
  for(j in 1:n+3){
    MAs_kp1[[j]] = MAs_k[[j]] + rho*(As_kp1[[j]] - As_kp1$A[,(j-3):n])
    MBs_kp1[[j]] = MBs_k[[j]] + rho*(Bs_kp1[[j]] - Bs_kp1$B[,(j-3):n])
  }
  
  #Update values at k
  As_k = As_kp1; Bs_k = Bs_kp1
  MAs_k = MAs_kp1; MBs_k = MBs_kp1
  
  #Compute ADMM MSEs
  A_ADMM = As_kp1$A
  B_ADMM = Bs_kp1$B
  Pi_ADMM = A_ADMM%*%t(B_ADMM)
  MSEs["ADMM",] = c(sqrt(sum((A-A_ADMM)^2)),
                    sqrt(sum((B-B_ADMM)^2)),
                    sqrt(sum((Pi-Pi_ADMM)^2)))
  zero_js = which(sapply(1:n,function(j){
    As_kp1[[j+1]][1]==0
  }))
  if(length(zero_js)==0){
    r_hat = n
  }else if(length(zero_js)==n){
    r_hat = 1
  }else{
    r_hat = max(c(min(zero_js)-1,1))
  }
  
  #Compute A,B with rank restriction
  A_ADMM_restr = cbind(A_ADMM[,1:r_hat],matrix(0,n,n-r_hat))
  B_ADMM_restr = cbind(B_ADMM[,1:r_hat],matrix(0,n,n-r_hat))
  Pi_ADMM_restr = A_ADMM_restr%*%t(B_ADMM_restr)
  MSEs["ADMM_restr",] = c(sqrt(sum((A-A_ADMM_restr)^2)),
                    sqrt(sum((B-B_ADMM_restr)^2)),
                    sqrt(sum((Pi-Pi_ADMM_restr)^2)))
  
  #Compute A,B with rank restriction and normalization
  A_ADMM_norm = A_ADMM_restr*-0.5/A_ADMM_restr[1,1]
  B_ADMM_norm = B_ADMM_restr*A_ADMM_restr[1,1]/(-0.5)
  MSEs["ADMM_norm",] = c(sqrt(sum((A-A_ADMM_norm)^2)),
                          sqrt(sum((B-B_ADMM_norm)^2)),
                          MSEs["ADMM_restr","Pi"])
  
  
  #print(MSEs[1,])
  
  #Compute loss
  loss = loss_function(Y_diff,Y_lag,A_ADMM,B_ADMM,lambda_L1,lambda_L2,omegas)
  loss_diff = (loss-loss_old)/loss_old
  #print(loss)
  #if(loss_diff>-0.001){break}
  
  print(sum(cbind(A_ADMM,B_ADMM)^2))
  
  #Compute primal residual
  # r_norm = 0
  # for(k in 1:(n+1)){
  #   j = N_j[k]
  #   r_norm = r_norm + sum((cbind(As_kp1[[k]],Bs_kp1[[k]]) - cbind(As_kp1$A[,j:n],Bs_kp1$B[,j:n]))^2)
  # }
  # r_norm = sqrt(r_norm)
  # print(r_norm)
}
rank_hat = min(which(sapply(1:n,function(j){
  As_kp1[[j+1]][1]==0
})))-1
rank_hat
MSEs

loss_function(Y_diff,Y_lag,A_johansen,B_johansen,lambda_L1,lambda_L2,omegas)
