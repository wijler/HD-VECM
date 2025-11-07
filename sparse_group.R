rm(list=ls())
setwd("C:/Github projects/HD-VECM")
library(urca)
source("source.R")
Rcpp::sourceCpp("sparse_group_fast.cpp")

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
  + 2*(sum(A!=0)+sum(B!=0))/nrow(Y) 
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
AB_loss = function(DY,Y_lag,AB,lambda_R,A=NULL,B=NULL){
  loss = sum((DY-Y_lag%*%t(AB))^2)/nrow(Y_lag)
  if(lambda_R>0){
     loss = loss + lambda_R*(sum(A^2)+sum(B^2))
  }
  loss
}
AB_step_size = function(DY,Y_lag,A,B,nabla_A,nabla_B,
                        lambda_L2,lambda_L1=0,lambda_R=0,
                        omegas,step_init=1,step_mult=0.5,
                        max_iter = 1000){
  
  n = nrow(A); n2 = 2*n
  AB = A%*%t(B)
  c = c(rbind(A,B))
  ind_A = c(sapply(1:n,function(j) 1:n + (j-1)*n2))
  ind_B = ind_A + n
  nabla_c = c(rbind(nabla_A,nabla_B))
  g_old = AB_loss(DY,Y_lag,AB,lambda_R,A,B)
  step_size = step_init/step_mult
  count = 0
  while(TRUE){
    count = count+1
    step_size = step_mult*step_size
    step_lambda_L2 = step_size*lambda_L2
    c_tmp = c - step_size*nabla_c
    
    #proximal update
    if(lambda_L1==0){
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
    g_new = AB_loss(DY,Y_lag,AB_new,lambda_R,A_new,B_new)
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
                   lambda_L2,lambda_L1=0,lambda_R=0,omegas=NULL,
                   step_init=1,step_mult=0.5,
                   step_max_iter=100,
                   max_iter=1000,thresh=1e-5,
                   print_dist=TRUE,print_loss=FALSE){
  
  #Transform data
  n = ncol(Y)
  Y_lag = Y[-nrow(Y),]
  DY = diff(Y)
  t = nrow(Y_lag)
  DYYx2 = 2*crossprod(DY,Y_lag)/t
  YYx2 = 2*crossprod(Y_lag)/t
  if(is.null(omegas)){
    omegas = 1/eigen(crossprod(Y)/(nrow(Y)^2))$values
  }
 
  #Conduct PGD
  step_convergence = rep(NA,max_iter)
  A_new = A_old = A
  B_new = B_old = B
  AB_new = A%*%t(B)
  iterations=0
  while(TRUE){
    
    #Update old coefficients
    A_old = A_new
    B_old = B_new
    AB_old = AB_new
    
    #Update gradient
    if(lambda_R > 0){
      nabla_A = (AB_old%*%YYx2 - DYYx2)%*%B_old + 2*lambda_R*A_old
      nabla_B = (YYx2%*%t(AB_old) - t(DYYx2))%*%A_old + 2*lambda_R*B_old
    }else{
      nabla_A = (AB_old%*%YYx2 - DYYx2)%*%B_old
      nabla_B = (YYx2%*%t(AB_old) - t(DYYx2))%*%A_old
    }
    
    #Perform proximal update with backtracking
    AB_step_size_obj = AB_step_size(DY,Y_lag,
                                    A_old,B_old,nabla_A,nabla_B,
                                    lambda_L2,lambda_L1,lambda_R,omegas,
                                    step_init,step_mult,step_max_iter)
    
    A_new = AB_step_size_obj$A_new
    B_new = AB_step_size_obj$B_new
    AB_new = AB_step_size_obj$AB_new
    step_convergence[iterations+1] = AB_step_size_obj$step_convergence
    
    
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
    if(print_dist & (iterations %% 100 == 0)){
      print(AB_dist)
    }
    if(print_loss & (iterations %% 100 == 0)){
      loss_new = AB_loss(DY,Y_lag,AB_new,lambda_R,A_new,B_new)
      print(loss_new)
    }
  }
  step_convergence = step_convergence[!is.na(step_convergence)]
  list(A=A_new,B=B_new,AB = AB_new,iterations=iterations,
       convergence=convergence,step_convergence=step_convergence)
}

VECM_SG_tuned = function(A,B,Y,lambda_grid,lambda_R,omegas=NULL,
                         step_init=1,step_mult=0.5,
                         step_max_iter=100,max_iter=1000,
                         thresh=1e-5,crit="BIC"){
  
  n_lambdas = nrow(lambda_grid)
  lambda_L2_changes = which(diff(lambda_grid[,1])!=0) #Used for updating initialization
  if(ncol(lambda_grid)==1){lambda_L1 = 0}
  if(is.null(omegas)){
    omegas = 1/eigen(crossprod(Y)/(nrow(Y)^2))$values
    }
  
  A_old = A;
  B_old = B;
  coefs = array(NA,dim=c(n,n,2,n_lambdas),
                dimnames = list(row=1:n,col=1:n,
                                mat=c("A","B"),
                                lambda = 1:n_lambdas))
  BICs = AICs = PI_MSEs = r_AB = rep(0,n_lambdas,2);
  for(j in 1:n_lambdas){
    lambda_L2 = lambda_grid[j,1]
    if(ncol(lambda_grid)==2){ lambda_L1 = lambda_grid[j,2] }
    
    VECM_obj = VECM_SG(A_old,B_old,Y,
                       lambda_L2,lambda_L1,lambda_R,omegas,
                       step_init,step_mult,
                       step_max_iter,max_iter,thresh,
                       print_dist=FALSE,print_loss=FALSE)
    coefs[,,1,j] = VECM_obj$A
    coefs[,,2,j] = VECM_obj$B

    #Update initializers for warm start
    if(j %in% lambda_L2_changes){
      last_lambda_index = min(which(lambda_grid[,1]==lambda_grid[j,1]))
      A_old = coefs[,,1,last_lambda_index]
      B_old = coefs[,,2,last_lambda_index]
    }else{
      A_old = VECM_obj$A;
      B_old = VECM_obj$B;
    }
    
    #Compute BICs
    r_AB[j] = sum(apply(A_init,2,function(x) !all(x==0)))
    BICs[j] = BIC_value(Y,VECM_obj$A,VECM_obj$B)
    AICs[j] = AIC_value(Y,VECM_obj$A,VECM_obj$B)
    PI_MSEs[j] = sqrt(sum((Pi-VECM_obj$AB)^2))
    print(j)
  }
  
  if(crit=="MSE"){
    lambdas_opt = lambda_grid[which.min(PI_MSEs),]
    ind_opt = which.min(PI_MSEs)
  }else if(crit=="AIC"){
    lambdas_opt = lambda_grid[which.min(AICs),]
    ind_opt = which.min(AICs)
  }else if(crit=="BIC"){
    lambdas_opt = lambda_grid[which.min(BICs),]
    ind_opt = which.min(BICs)
  }
  
  A_opt = coefs[,,1,ind_opt]
  B_opt = coefs[,,2,ind_opt]
  
  #return results
  list(lambda_opt = lambdas_opt,ind_opt=ind_opt,
       A = A_opt, B = B_opt,
       AICs = AICs,BICs = BICs, PI_MSEs = PI_MSEs)
  
}


#Nested version
prox_L2_nested = function(X,lambda,omegas){
  n = ncol(X)
  for(j in 1:n){
    group_norm = sqrt(sum(X[,(n+1-j):n]^2))
    threshold = group_norm*lambda*omegas[n+1-j]
    X_j = X[,(n+1-j):n]
    if(group_norm > threshold){
      X[,(n+1-j):n] = (1-threshold/group_norm)*X[,(n+1-j):n]
    }else{
      X[,(n+1-j):n] = 0
    }
    
    #Get rid of numerical noise
    X[,(n+1-j):n][abs(X[,(n+1-j):n])<1e-10] = 0
  }
  X
}
AB_step_size_nested = function(DY,Y_lag,A,B,nabla_A,nabla_B,
                               lambda_L2,lambda_R=0,
                               omegas,step_init=1,step_mult=0.5,
                               max_iter = 1000){
  
  n = nrow(A); n2 = 2*n
  AB = A%*%t(B)
  C = rbind(A,B)
  nabla_C = rbind(nabla_A,nabla_B)
  g_old = AB_loss(DY,Y_lag,AB,lambda_R,A,B)
  step_size = step_init/step_mult
  count = 0
  while(TRUE){
    count = count+1
    step_size = step_mult*step_size
    step_lambda_L2 = step_size*lambda_L2
    C_tmp = C - step_size*nabla_C
    
    #proximal update
    C_new = prox_L2_nested(C_tmp,step_lambda_L2,omegas)
    
    #Update parameters
    A_new = C_new[1:n,]
    B_new =  C_new[(n+1):n2,]
    AB_new = A_new%*%t(B_new)
    
    #Backtracking line search
    G_new = (C - C_new)/step_size
    g_new = AB_loss(DY,Y_lag,AB_new,lambda_R,A_new,B_new)
    LB = g_old - step_size*sum(nabla_C*G_new) + step_size*sum(G_new^2)/2
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
VECM_SG_nested = function(A,B,Y,
                          lambda_L2,lambda_R=0,omegas=NULL,
                          step_init=1,step_mult=0.5,
                          step_max_iter=100,
                          max_iter=1000,thresh=1e-5,
                          print_dist=TRUE,print_loss=FALSE){
  
  #Transform data
  n = ncol(Y)
  Y_lag = Y[-nrow(Y),]
  DY = diff(Y)
  t = nrow(Y_lag)
  DYYx2 = 2*crossprod(DY,Y_lag)/t
  YYx2 = 2*crossprod(Y_lag)/t
  if(is.null(omegas)){
    #omegas = 1/cumsum(eigen(crossprod(Y)/(nrow(Y)^2))$values[n:1])[n:1]
    omegas = 1/eigen(crossprod(Y)/(nrow(Y)^2))$values
  }
  
  #Conduct PGD
  step_convergence = rep(NA,max_iter)
  A_new = A_old = A
  B_new = B_old = B
  AB_new = A%*%t(B)
  iterations=0
  while(TRUE){
    
    #Update old coefficients
    A_old = A_new
    B_old = B_new
    AB_old = AB_new
    
    #Update gradient
    if(lambda_R > 0){
      nabla_A = (AB_old%*%YYx2 - DYYx2)%*%B_old + 2*lambda_R*A_old
      nabla_B = (YYx2%*%t(AB_old) - t(DYYx2))%*%A_old + 2*lambda_R*B_old
    }else{
      nabla_A = (AB_old%*%YYx2 - DYYx2)%*%B_old
      nabla_B = (YYx2%*%t(AB_old) - t(DYYx2))%*%A_old
    }
    
    #Perform proximal update with backtracking
    AB_step_size_obj = AB_step_size_nested(DY,Y_lag,
                                           A_old,B_old,nabla_A,nabla_B,
                                           lambda_L2,lambda_R,omegas,
                                           step_init,step_mult,step_max_iter)
                                    
    
    A_new = AB_step_size_obj$A_new
    B_new = AB_step_size_obj$B_new
    AB_new = AB_step_size_obj$AB_new
    step_convergence[iterations+1] = AB_step_size_obj$step_convergence
    
    
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
    if(print_dist & (iterations %% 100 == 0)){
      print(AB_dist)
    }
    if(print_loss & (iterations %% 100 == 0)){
      loss_new = AB_loss(DY,Y_lag,AB_new,lambda_R,A_new,B_new)
      print(loss_new)
    }
  }
  step_convergence = step_convergence[!is.na(step_convergence)]
  list(A=A_new,B=B_new,AB = AB_new,iterations=iterations,
       convergence=convergence,step_convergence=step_convergence,omegas=omegas)
}
                   


#Accelerated version
AB_step_size_accelerated = function(DY,Y_lag,A,B,nabla_A,nabla_B,
                        lambda_L2,lambda_L1=0,lambda_R=0,
                        omegas,step_init=1,step_mult=0.5,
                        max_iter = 1000){
  
  n = nrow(A); n2 = 2*n
  AB = A%*%t(B)
  c = c(rbind(A,B))
  ind_A = c(sapply(1:n,function(j) 1:n + (j-1)*n2))
  ind_B = ind_A + n
  nabla_c = c(rbind(nabla_A,nabla_B))
  g_old = AB_loss(DY,Y_lag,AB,lambda_R,A,B)
  step_size = step_init/step_mult
  count = 0
  while(TRUE){
    count = count+1
    step_size = step_mult*step_size
    step_lambda_L2 = step_size*lambda_L2
    c_tmp = c - step_size*nabla_c
    
    #proximal update
    if(lambda_L1==0){
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
    c_diff = c_new - c
    g_new = AB_loss(DY,Y_lag,AB_new,lambda_R,A_new,B_new)
    LB = g_old + t(nabla_c)%*%c_diff + sum(c_diff^2)/(2*step_size)
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
VECM_SG_accelerated = function(A,B,Y,
                               lambda_L2,lambda_L1=0,lambda_R=0,omegas,
                               step_size="auto",step_init=1,step_mult=0.5,
                               step_max_iter=100,
                               max_iter=1000,thresh=1e-5,
                               print_dist=TRUE,print_loss=FALSE){
  
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
  A_new = A_old = A
  B_new = B_old = B
  AB_new = A%*%t(B)
  step_new = step_old = step_init
  iterations=0
  while(TRUE){
    
    #Update old coefficients
    k_mult = (iterations-1)/(iterations+2)
    A_old = A_new + k_mult*(A_new - A_old)
    B_old = B_new + k_mult*(A_new - A_old)
    AB_old = A_old%*%t(B_old)
    step_old = step_new
    
    #Update gradient
    if(lambda_R > 0){
      nabla_A = (AB_old%*%YYx2 - DYYx2)%*%B_old + 2*lambda_R*A_old
      nabla_B = (YYx2%*%t(AB_old) - t(DYYx2))%*%A_old + 2*lambda_R*B_old
    }else{
      nabla_A = (AB_old%*%YYx2 - DYYx2)%*%B_old
      nabla_B = (YYx2%*%t(AB_old) - t(DYYx2))%*%A_old
    }
    
    #Perform proximal update with backtracking
    AB_step_size_obj = AB_step_size_accelerated(DY,Y_lag,
                                    A_old,B_old,nabla_A,nabla_B,
                                    lambda_L2,lambda_L1,lambda_R,omegas,
                                    step_old,step_mult,step_max_iter)
    
    A_new = AB_step_size_obj$A_new
    B_new = AB_step_size_obj$B_new
    AB_new = AB_step_size_obj$AB_new
    step_new = AB_step_size_obj$step_size
    step_convergence[iterations+1] = AB_step_size_obj$step_convergence
    
    
    #Check convergence
    iterations = iterations + 1
    AB_dist = sqrt(sum((A_new - A_old)^2) + sum((B_new - B_old)^2))/step_new
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
    if(print_loss){
      loss_new = AB_loss(DY,Y_lag,AB_new,lambda_R,A_new,B_new)
      print(loss_new)
    }
  }
  step_convergence = step_convergence[!is.na(step_convergence)]
  list(A=A_new,B=B_new,AB = AB_new,iterations=iterations,
       convergence=convergence,step_convergence=step_convergence)
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
A_Johan = vecm_fit@W
B_Johan = vecm_fit@V
C_init = AB_Johan = C_Johan = A_Johan%*%t(B_Johan)


#Testing PGD nested
A_init=A_Johan;B_init=B_Johan
lambda_L2=1;
step_init=1;step_mult=0.5
step_max_iter=100;max_iter = 1000;thresh=1e-4
test = VECM_SG_nested(A=A_Johan,B=B_Johan,
               Y=Y,lambda_L2=0.01,lambda_R=1e-2,
               step_init=1e-3,step_mult=0.5,
               step_max_iter=100,
               max_iter = 1e7,thresh=1e-5,
               print_dist=F,print_loss = T)
test$convergence
mean(test$step_convergence)
c(norm(Pi - AB_Johan,"F"),norm(Pi - test$AB,"F"))
test$A



#Testing RCPP implementation
test_Rcpp = VECM_SG_Rcpp(A=A_Johan,B=B_Johan,
                         Y=Y,lambda_L2=0.01,omegas=matrix(-999,1,1),lambda_L1=0.01,lambda_R=1e-2,
                         step_init=1e-3,step_mult=0.5,
                         step_max_iter=100,
                         max_iter = 1e7,thresh=1e-5,
                         print_dist=F,print_loss = T)
                    
test_Rcpp$convergence
mean(test_Rcpp$step_convergence)
c(norm(Pi - AB_Johan,"F"),norm(Pi - test_Rcpp$AB,"F"))
test_Rcpp$B


#Tuned results
lambda_L2_grid = exp(seq(log(0.001),log(1),length.out=10))
lambda_grid = matrix(lambda_L2_grid,ncol=1)
VECM_L2 = VECM_SG_tuned_Rcpp(A=A_Johan,B=B_Johan,
                                Y=Y,lambda_grid,-999,lambda_R=1e-2,
                                step_init=1e-3,step_mult=0.5,
                                step_max_iter=100,
                                max_iter = 1e7,thresh=1e-5,
                                crit="BIC")
                           

lambda_L1_grid = lambda_L2_grid
lambda_grid = as.matrix(expand.grid(lambda_L1_grid,lambda_L2_grid)[,c(2,1)])
VECM_L2_L1 = VECM_SG_tuned_Rcpp(A=A_Johan,B=B_Johan,
                                      Y=Y,lambda_grid,-999,lambda_R=1e-2,
                                      step_init=1e-3,step_mult=0.5,
                                      step_max_iter=100,
                                      max_iter = 1e7,thresh=1e-5,
                                      crit="BIC")
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