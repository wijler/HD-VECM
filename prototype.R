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
loss_function = function(Y_D,Y_lag,A,B,lambda_1,lambda_g,omegas){
  t = nrow(Y_D)+1
  n = ncol(A)
  l2_loss = sum(sapply(1:n,function(j){
    omegas[j]*sqrt(sum(A[,j:n]^2) + sum(B[,j:n]^2))
  }))
  loss = sum((Y_D - Y_lag%*%B%*%t(A)))^2/t +
    lambda_1*(sum(abs(A))+sum(abs(B))) +
    lambda_g*l2_loss
  loss
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
B_hat = vecm_fit@V
A_hat = vecm_fit@W
Pi_hat = A_hat%*%t(B_hat)

As_k = MAs_k = Bs_k = MBs_k = list()
Zeros = matrix(0,n,n)
As_k[[1]] = As_k[[n+2]] = A_hat
MAs_k[[1]] = MAs_k[[n+2]] = Zeros
Bs_k[[1]] = Bs_k[[n+2]] = B_hat
MBs_k[[1]] = MBs_k[[n+2]] = Zeros
A_names = c(paste("A",rep(0:n),sep=""),"A")
MA_names = c(paste("MA",rep(0:n),sep=""),"MA")
B_names = c(paste("B",rep(0:n),sep=""),"B")
MB_names = c(paste("MB",rep(0:n),sep=""),"MB")
for(j in 2:(n+1)){
  As_k[[j]] = A_hat[,(j-1):n]
  MAs_k[[j]] = Zeros[,(j-1):n]
  Bs_k[[j]] = B_hat[,(j-1):n]
  MBs_k[[j]] = Zeros[,(j-1):n]
}
names(As_k) = A_names
names(MAs_k) = MA_names
names(Bs_k) = B_names
names(MBs_k) = MB_names

As_kp1 = As_k
MAs_kp1 = MAs_k
Bs_kp1 = Bs_k
MBs_kp1 = MBs_k

#Compute MSEs for Johansen's method and OLS
A_johansen = cbind(A_hat[,1],matrix(0,n,n-1))
B_johansen = cbind(B_hat[,1],matrix(0,n,n-1))
Pi_johansen = A_johansen%*%t(B_johansen)
MSEs_johansen = c(sqrt(sum((A_true-A_johansen)^2)),
                  sqrt(sum((B_true-B_johansen)^2)),
                  sqrt(sum((Pi-Pi_johansen)^2)))
MSEs_OLS = c(sqrt(sum((A_true-A_hat)^2)),
             sqrt(sum((B_true-B_hat)^2)),
             sqrt(sum((Pi-Pi_hat)^2)))

#ADMM algorithm
Y_lag = Y[-t,]
Y_diff = diff(Y)
t = t-1
Y_lag_cross = crossprod(Y_lag)/t
DY_Y_lag = crossprod(Y_diff,Y_lag)/t

omegas = 1/eigen(crossprod(Y)/t^2)$values
N_j = c(1,1:n)
rho = 100; lambda_1 = 1e-8*rho; lambda_g = 0.01*rho
I_rhod2 = diag(rho/2,n)
I_rhod2_kron = diag(rho/2,n^2)
MSEs = matrix(0,5,3)
rownames(MSEs) = c("ADMM","ADMM_restr","ADMM_norm","Johansen","OLS")
colnames(MSEs) = c("A","B","Pi")
MSEs[c("Johansen","OLS"),] = rbind(MSEs_johansen,MSEs_OLS)
loss = loss_function(Y_diff,Y_lag,A_hat,B_hat,lambda_1,lambda_g,omegas)
max_iter = 10000
loss_vec = rep(NA,max_iter)
counter = 0
while(counter < max_iter + 1){
  counter = counter + 1
  
  loss_old = loss
  
  #Update A0
  As_kp1$A0 = (DY_Y_lag%*%Bs_k$B0 + (rho/2)*As_k$A - (1/2)*MAs_k$MA0)%*%
    solve(t(Bs_k$B0)%*%Y_lag_cross%*%Bs_k$B0 + I_rhod2)
  
  #Update B0
  Bs_kp1$B0 = matrix(solve(kronecker(crossprod(As_kp1$A0),Y_lag_cross) + I_rhod2_kron,
                           c(t(DY_Y_lag)%*%As_kp1$A0 + 0.5*(rho*Bs_k$B - MBs_k$MB0))),n,n)
  
  #Update (Aj,Bj), j=1,...,N
  for(j in 1:n){
    omega_j = omegas[j] #Change later
    diff_j = c(cbind(As_k$A[,j:n],Bs_k$B[,j:n]) - cbind(MAs_k[[j+1]],MBs_k[[j+1]])/rho)
    shrink_mult = 1-lambda_g*omega_j/(rho*sqrt(sum(diff_j^2)))
    if(shrink_mult>0){
      AB_j = matrix(shrink_mult*diff_j,nrow=n)
    }else{
      AB_j = matrix(0,n,2*(n+1-j))
    }
    As_kp1[[j+1]] = AB_j[,1:(n+1-j)]
    Bs_kp1[[j+1]] = AB_j[,(n+2-j):(2*(n+1-j))]
  }
  
  #Update A and B
  
  A_kp1_sums = As_kp1$A0
  B_kp1_sums = Bs_kp1$B0
  W_j = diag(1/(1:n+1))
  for(j in 1:n){
    A_kp1_sums[,j:n] = A_kp1_sums[,j:n] + As_kp1[[j+1]]
    B_kp1_sums[,j:n] = B_kp1_sums[,j:n] + Bs_kp1[[j+1]]
  }
  As_kp1$A = softt_matrix(A_kp1_sums,lambda_1/rho)%*%W_j
  Bs_kp1$B = softt_matrix(B_kp1_sums,lambda_1/rho)%*%W_j
  # As_kp1$A = softt_matrix(A_kp1_sums%*%W_j,lambda_1/rho)
  # Bs_kp1$B = softt_matrix(B_kp1_sums%*%W_j,lambda_1/rho)
  
  #As_kp1$A; Bs_kp1$B
  
  #Update MA and MB
  MAs_kp1$MA0 = MAs_k$MA0 + rho*(As_kp1$A0 - As_kp1$A)
  MBs_kp1$MB0 = MBs_k$MB0 + rho*(Bs_kp1$B0 - Bs_kp1$B)
  for(j in 2:(n+1)){
    MAs_kp1[[j]] = MAs_k[[j]] + rho*(As_kp1[[j]] - As_kp1$A[,(j-1):n])
    MBs_kp1[[j]] = MBs_k[[j]] + rho*(Bs_kp1[[j]] - Bs_kp1$B[,(j-1):n])
  }
  
  #Update values at k
  As_k = As_kp1; Bs_k = Bs_kp1
  MAs_k = MAs_kp1; MBs_k = MBs_kp1
  
  #Compute ADMM MSEs
  A_ADMM = As_kp1$A
  B_ADMM = Bs_kp1$B
  Pi_ADMM = A_ADMM%*%t(B_ADMM)
  MSEs["ADMM",] = c(sqrt(sum((A_true-A_ADMM)^2)),
                    sqrt(sum((B_true-B_ADMM)^2)),
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
  MSEs["ADMM_restr",] = c(sqrt(sum((A_true-A_ADMM_restr)^2)),
                    sqrt(sum((B_true-B_ADMM_restr)^2)),
                    sqrt(sum((Pi-Pi_ADMM_restr)^2)))
  
  #Compute A,B with rank restriction and normalization
  A_ADMM_norm = A_ADMM_restr*-0.5/A_ADMM_restr[1,1]
  B_ADMM_norm = B_ADMM_restr*A_ADMM_restr[1,1]/(-0.5)
  MSEs["ADMM_norm",] = c(sqrt(sum((A_true-A_ADMM_norm)^2)),
                          sqrt(sum((B_true-B_ADMM_norm)^2)),
                          MSEs["ADMM_restr","Pi"])
  
  
  #print(MSEs[1,3])
  
  #Compute loss
  loss = loss_function(Y_diff,Y_lag,A_ADMM,B_ADMM,lambda_1,lambda_g,omegas)
  loss_diff = (loss-loss_old)/loss_old
  loss_vec[counter] = loss
  print(loss)
  #if(loss_diff>-0.001){break}
  
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

loss_function(Y_diff,Y_lag,A_johansen,B_johansen,lambda_1,lambda_g,omegas)
