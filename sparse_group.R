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

#Functions to update A and B
softt_matrix = function(X,lambda){
  X_sign = sign(X)
  X_tmp = abs(X)-lambda
  X_tmp[X_tmp<0] = 0
  X_sign*X_tmp
}
prox_L2 = function(x,n,lambda,omegas){
  n2 = 2*n
  for(j in 1:n){
    ind_j = 1:n2 + (j-1)*n2
    c_mult = 1 - omegas[j]*lambda/sqrt(sum(x[ind_j]^2))
    if(c_mult>0){
      x[ind_j] = c_mult*x[ind_j]
    }else{
      x[ind_j] = 0
    }
  }
  x
}
prox_L2_L1 = function(x,n,lambda_L2,lambda_L1,omegas){
  n2 = 2*n
  for(j in 1:n){
    ind_j = 1:n2 + (j-1)*n2
    x_tmp = softt_matrix(x[ind_j],lambda_L1)
    c_mult = 1 - omegas[j]*lambda_L2/sqrt(sum(x_tmp^2))
    if(c_mult>0){
      x[ind_j] = c_mult*x_tmp
    }else{
      x[ind_j] = 0
    }
  }
  x
}
AB_loss = function(DY,Y_lag,AB){
  sum((DY-Y_lag%*%t(AB))^2)/nrow(Y_lag)
}
AB_step_size = function(DY,Y_lag,A,B,nabla_A,nabla_B,
                        lambda_L2,lambda_L1=NULL,
                        omegas,step_init=1,step_mult=0.5,
                        max_iter = 1000){
  
  n = nrow(A); n2 = 2*n
  AB = A%*%t(B)
  c = c(rbind(A,B))
  ind_A = c(sapply(1:n,function(j) 1:n + (j-1)*n2))
  ind_B = ind_A + n
  nabla_c = c(rbind(nabla_A,nabla_B))
  g_old = AB_loss(DY,Y_lag,AB)
  step_size = step_init/step_mult
  count = 0
  while(TRUE){
    count = count+1
    step_size = step_mult*step_size
    step_lambda_L2 = step_size*lambda_L2
    c_tmp = c - step_size*nabla_c
    
    #proximal update
    if(is.null(lambda_L1)){
      c_new = prox_L2(c_tmp,n,step_lambda_L2,omegas)
    }else{
      step_lambda_L1 = step_size*lambda_L1
      c_new = prox_L2_L1(c_tmp,n,step_lambda_L2,step_lambda_L1,omegas)
    }
    
    #Update parameters
    A_new = matrix(c_new[ind_A],n,n)
    B_new =  matrix(c_new[ind_B],n,n)
    AB_new = A_new%*%t(B_new)
    
    #Backtracking line search
    G_new = (c - c_new)/step_size
    g_new = AB_loss(DY,Y_lag,AB_new)
    LB = g_old - step_size*t(nabla_c)%*%G_new + step_size*sum(G_new^2)/2
    if(!(g_new > LB)){
      step_convergence = 0
      break
    }else if(count>max_iter){
      step_convergence = 1
      break
    }
  }
  list(step_size=step_size,A_new=A_new,B_new=B_new,
       AB_new = AB_new,step_iterations=count,
       step_convergence = step_convergence)
}

VECM_SG = function(A,B,Y,
                   lambda_L2,lambda_L1=NULL,omegas,
                   step_size="auto",step_init=1,step_mult=0.5,
                   step_max_iter=100,
                   max_iter=1000,thresh=1e-5,
                   print_dist=TRUE){
  
  #Transform data
  n = ncol(Y)
  Y_lag = Y[-nrow(Y),]
  DY = diff(Y)
  t = nrow(Y_lag)
  DYYx2 = 2*crossprod(DY,Y_lag)/t
  YYx2 = 2*crossprod(Y_lag)/t
  omegas = 1/eigen(crossprod(Y)/(nrow(Y)^2))$values
  
  #Conduct PGD
  step_convergence = rep(NA,max_iter)
  A_new = A
  B_new = B
  AB_new = A%*%t(B)
  iterations=0
  while(TRUE){
    A_old = A_new
    B_old = B_new
    AB_old = AB_new
    nabla_A = (AB_old%*%YYx2 - DYYx2)%*%B_old
    nabla_B = (YYx2%*%t(AB_old) - t(DYYx2))%*%A_old
    if(step_size == "auto"){
      
      AB_step_size_obj = AB_step_size(DY,Y_lag,
                                      A_old,B_old,nabla_A,nabla_B,
                                      lambda_L2,lambda_L1,omegas,
                                      step_init,step_mult,step_max_iter)
                   
      A_new = AB_step_size_obj$A_new
      B_new = AB_step_size_obj$B_new
      AB_new = AB_step_size_obj$AB_new
      step_convergence[iterations+1] = AB_step_size_obj$step_convergence
    }else if(is.numeric(step_size)){
      
      # for(j in 1:ncol(A_old)){
      #   c_tmp = c(A_old[,j],B_old[,j]) - step_size*c(nabla_A[,j],nabla_B[,j])
      #   c_prox = prox_L2(c_tmp,step_size*lambda*omegas[j])
      #   A_new[,j] = head(c_prox,n)
      #   B_new[,j] = tail(c_prox,n)
      # }
      # AB_new = A_new%*%t(B_new)
      
    }
    
    #Check convergence
    iterations = iterations + 1
    AB_dist = sqrt(sum((AB_new-AB_old)^2))/n
    if(AB_dist < thresh){
      convergence=0
      break
    }
    if(iterations > max_iter){
      convergence=1
      break
    }
    if(print_dist){
      print(AB_dist)
    }
  }
  step_convergence = step_convergence[!is.na(step_convergence)]
  list(A=A_new,B=B_new,AB = AB_new,iterations=iterations,
       convergence=convergence,step_convergence=step_convergence)
}

VECM_nuclear_tuned = function(C_init,A_init,B_init,M_init,Y,
                              lambda_nuclear,lambda_grid,
                              rho="auto",rho_init=0.1,rho_mult=0.5,mu=10,
                              step_size = "auto",step_init=1,step_mult=0.5,
                              step_max_iter=100,max_iter = 1000,
                              ADMM_thresh=1e-4,C_thresh=1e-4,AB_thresh=1e-4){
  
  n_lambdas = nrow(lambda_grid)
  lambda_L2_changes = which(diff(lambda_grid[,1])!=0) #Used for updating initialization
  if(ncol(lambda_grid)==1){lambda_L1 = NULL}
  C_init_old = C_init; A_init_old = A_init; 
  B_init_old = B_init; M_init_old = M_init;
  coefs = array(NA,dim=c(n,n,4,n_lambdas),
                dimnames = list(row=1:n,col=1:n,
                                mat=c("C","A","B","M"),
                                lambda = 1:n_lambdas))
  BICs = PI_MSEs = r_C = r_AB = rep(0,n_lambdas,2);
  for(j in 1:n_lambdas){
    lambda_L2 = lambda_grid[j,1]
    if(ncol(lambda_grid)==2){ lambda_L1 = lambda_grid[j,2] }
    
    VECM_obj = VECM_nuclear(C_init = C_init,A_init = A_init,B_init = B_init,
                            M_init = M_init,Y = Y,lambda_nuclear = 0.001,
                            lambda_L2 = lambda_L2,lambda_L1 = lambda_L1,
                            rho = rho,rho_init=rho_init,rho_mult = rho_mut,
                            mu=mu,step_size = step_size,step_init = step_init,
                            step_mult = step_mult,step_max_iter = step_max_iter,
                            max_iter = max_iter,ADMM_thresh = ADMM_thresh,
                            C_thresh = C_thresh,AB_thresh = AB_thresh)
    coefs[,,1,j] = VECM_obj$C
    coefs[,,2,j] = VECM_obj$A
    coefs[,,3,j] = VECM_obj$B
    coefs[,,4,j] = VECM_obj$M
    
    #Update initializers for warm start
    if(j %in% lambda_L2_changes){
      last_lambda_index = min(which(lambda_grid[,1]==lambda_grid[j,1]))
      C_init = coefs[,,1,last_lambda_index]
      A_init = coefs[,,2,last_lambda_index]
      B_init = coefs[,,3,last_lambda_index]
      M_init = coefs[,,4,last_lambda_index]
    }else{
      C_init = VECM_obj$C; A_init = VECM_obj$A; 
      B_init = VECM_obj$B; M_init = VECM_obj$M;
    }
    
    #Compute BICs
    r_C[j] = sum(VECM_obj$d!=0)
    r_AB[j] = sum(apply(A_init,2,function(x) !all(x==0)))
    BICs[j] = BIC_value(Y,A_init,B_init)
    PI_MSEs[j] = sqrt(sum((Pi-A_init%*%t(B_init))^2))
    print(j)
  }
  lambdas_opt = lambda_grid[which.min(PI_MSEs),] #Cheating a bit for now
  ind_opt = which.min(PI_MSEs)
  A_opt = coefs[,,"A",ind_opt]
  B_opt = coefs[,,"B",ind_opt]
  C_opt = coefs[,,"C",ind_opt]
  M_opt = coefs[,,"M",ind_opt]
  
  #return results
  list(lambda_opt = lambdas_opt,A = A_opt, B = B_opt, 
       C = C_opt, M = M_opt, BICs = BICs, PI_MSEs = PI_MSEs)
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

# n = 10; t = 500; burnin = 50
# t_total = t+burnin
# alpha = c(-0.5,rep(0,n-1))
# beta = c(1,rep(-1,n-1))
# A = cbind(alpha,matrix(0,n,n-1))
# B = cbind(beta,matrix(0,n,n-1))
# VECM_sim = generate_VECM(n,t,burnin,A,B)
# Pi = VECM_sim$Pi
# Y = VECM_sim$Y
# matplot(Y,type="l")

#Initial parameters
colnames(Y) = paste("x",1:n,sep="")
vecm_fit = ca.jo(Y,K=2)
r = choose_r(vecm_fit,0.05)
A_Johan = vecm_fit@W
B_Johan = vecm_fit@V
AB_Johan = C_Johan = A_Johan%*%t(B_Johan)


#Testing PGD
A_init=A_Johan;B_init=B_Johan
lambda_L2=0.2;lambda_L1=NULL;
step_size = "auto";step_init=1;step_mult=0.5
step_max_iter=100;max_iter = 1000;thresh=1e-4
test = VECM_SG(A=A_Johan,B=B_Johan,
               Y=Y,lambda_L2=0.1,lambda_L1=0.1,
               step_size = "auto",step_init=1,step_mult=0.5,
               step_max_iter=100,max_iter = 1000,thresh=1e-5,
               print_dist=T)
test$convergence
mean(test$step_convergence)
c(norm(Pi - AB_Johan,"F"),norm(Pi - test$AB,"F"))
test$A








#Tuned results
rho=0.001
lambda_L2_grid = rho*exp(seq(log(1e-3),0,length.out=))
lambda_grid = matrix(lambda_L2_grid,ncol=1)
VECM_L2 = VECM_nuclear_tuned(C_init = C_Johan,A_init = A_Johan,B_init = B_Johan,
                             M_init = M_init,Y = Y,
                             lambda_nuclear = 0.001,lambda_grid,
                             rho = rho,rho_init=0.1,rho_mult = 0.5,mu=10,
                             step_size = "auto",step_init = 1,step_mult = 0.5,
                             step_max_iter = 100,max_iter = 1000,
                             ADMM_thresh = 1e-5,C_thresh = 1e-6,AB_thresh = 1e-5)

lambda_L1_grid = lambda_L2_grid
lambda_grid = expand.grid(lambda_L1_grid,lambda_L2_grid)[,c(2,1)]
VECM_L2_L1 = VECM_nuclear_tuned(C_init = C_Johan,A_init = A_Johan,B_init = B_Johan,
                             M_init = M_init,Y = Y,
                             lambda_nuclear = 0.001,lambda_grid,
                             rho = rho,rho_init=0.1,rho_mult = 0.5,mu=10,
                             step_size = "auto",step_init = 1,step_mult = 0.5,
                             step_max_iter = 100,max_iter = 1000,
                             ADMM_thresh = 1e-5,C_thresh = 1e-6,AB_thresh = 1e-5)
#Compute losses

#VECM-L2
A_L2 = VECM_L2$A
B_L2 = VECM_L2$B
AB_L2 = A_L2%*%t(B_L2)
A_L2_norm = A_L2*A[1,1]/A_L2[1,1]
B_L2_norm = B_L2*A_L2[1,1]/A[1,1]
r_L2 = sum(apply(A_L2,2,function(x) !all(x==0)))

#VECM L2 + L1
A_L2_L1 = VECM_L2_L1$A
B_L2_L1 = VECM_L2_L1$B
AB_L2_L1 = A_L2_L1%*%t(B_L2_L1)
A_L2_L1_norm = A_L2_L1*A[1,1]/A_L2_L1[1,1]
B_L2_L1_norm = B_L2_L1*A_L2_L1[1,1]/A[1,1]
r_L2_L1 = sum(apply(A_L2_L1,2,function(x) !all(x==0)))

#VECM RR
A_RR = A_Johan; A_RR[,-c(1:r)] = 0
B_RR = B_Johan; B_RR[,-c(1:r)] = 0
AB_RR = A_RR%*%t(B_RR)

fits_Pi = c(norm(Pi - AB_Johan,"F"),norm(Pi - AB_RR,"F"),
            norm(Pi - AB_L2,"F"),norm(Pi - AB_L2_L1,"F"))

fits_A = c(norm(A - A_Johan,"F"),norm(A - A_RR,"F"),
           norm(A - A_L2,"F"),norm(A - A_L2_L1,"F"))

fits_A_norm = c(norm(A - A_Johan*A[1,1]/A_Johan[1,1],"F"),
                norm(A - A_RR*A[1,1]/A_RR[1,1],"F"),
                norm(A - A_L2_norm,"F"),norm(A - A_L2_L1_norm,"F"))

fits_B = c(norm(B - B_Johan,"F"),norm(B - B_RR,"F"),
           norm(B - B_L2,"F"),norm(B - B_L2_L1,"F"))
fits_B_norm = c(norm(B - B_Johan*A_Johan[1,1]/A[1,1],"F"),
                norm(B - B_RR*A_RR[1,1]/A[1,1],"F"),
                norm(B - B_L2_norm,"F"),norm(B - B_L2_L1_norm,"F"))
ranks = c(ncol(A_Johan),r,r_L2,r_L2_L1)
fits = rbind(fits_Pi,fits_A,fits_A_norm,fits_B,fits_B_norm,ranks)
rownames(fits) = c("Pi","A","A_norm","B","B_norm","rank")
colnames(fits) = c("OLS","Johansen","VECM_L2","VECM_L2_L1")
fits





# #Compute VECM based on ADMM with Nuclear norm, L2 and L1 penalty
# VECM_NN_L2_L1 = VECM_nuclear(C_init = C_Johan,A_init = A_Johan,B_init = B_Johan,
#                              M_init = M_init,Y = Y,lambda_nuclear = 0.01,
#                              lambda_L2 = 0.0001,lambda_L1 = 0.0001,
#                              rho = 0.001,rho_init=0.1,rho_mult = 0.5,mu=10,
#                              step_size = "auto",step_init = 1,step_mult = 0.5,
#                              step_max_iter = 100,max_iter = 1000,
#                              ADMM_thresh = 1e-5,C_thresh = 1e-6,AB_thresh = 1e-5)
#                          
# VECM_NN_L2 = VECM_nuclear(C_init = C_Johan,A_init = A_Johan,B_init = B_Johan,
#                           M_init = M_init,Y = Y,lambda_nuclear = 0.01,
#                           lambda_L2 = 0.0001,lambda_L1 = NULL,
#                           rho = "auto",rho_init=0.1,rho_mult = 0.5,mu=10,
#                           step_size = "auto",step_init = 1,step_mult = 0.5,
#                           step_max_iter = 100,max_iter = 1000,
#                           ADMM_thresh = 1e-5,C_thresh = 1e-6,AB_thresh = 1e-4)