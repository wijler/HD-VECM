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
  r = c(r_max:1)[min(which(vecm_test_stats>vecm_cvs))]
  r
}
softt_matrix = function(X,lambda,omegas=NULL){
  X_sign = sign(X)
  if(is.null(omegas)){
    X_tmp = abs(X)-lambda
  }else{
    X_tmp = abs(X)-lambda*omegas
  }
  X_tmp[X_tmp<0] = 0
  X_sign*X_tmp
}
loss_function = function(Y_D,Y_lag,A,B,lambda_L1,lambda_L2,lambda_F,omegas){
  t = nrow(Y_D)+1
  n = ncol(A)
  l2_loss = sum(sapply(1:n,function(j){
    omegas[j]*sqrt(sum(A[,j:n]^2) + sum(B[,j:n]^2))
  }))
  loss = sum((Y_D - Y_lag%*%B%*%t(A)))^2/t +
    lambda_L1*(sum(abs(A))+sum(abs(B))) +
    lambda_L2*l2_loss + lambda_F*(sum(A^2)+sum(B^2))
  loss
}
BIC_value = function(Y,A,B){
  log(det(crossprod(diff(Y) - Y[-nrow(Y),]%*%B%*%t(A))/nrow(Y))) + log(nrow(Y))*(sum(A!=0)+sum(B!=0))/nrow(Y) 
}
AIC_value = function(Y,A,B){
  log(det(crossprod(diff(Y) - Y[-nrow(Y),]%*%B%*%t(A))/nrow(Y))) + 2*(sum(A!=0)+sum(B!=0))/nrow(Y) 
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
prox_nuclear = function(X,lambda){
  X_svd = svd(X)
  d_softt = softt_matrix(X_svd$d,lambda)
  X_new = X_svd$u%*%diag(d_softt)%*%t(X_svd$v)
  list(X_new=X_new,d=d_softt)
}

#DGP functions
draw_VECM = function(t,spec="HD_sparse",a=-0.4){
  
  #Specify coefficients
  if(spec=="HD_sparse"){
    r=4
    beta = matrix(0,12,r)
    for(k in 1:r){
      beta[1:3+(k-1)*3,k] = 1
    }
    beta = beta[-12,]
    alpha = a*beta
    n = nrow(beta)
  }
  
  #Generate data
  burnin = 50
  t_total = t+burnin
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
  
  list(Y=Y,A=alpha,B=beta,r=r)
}
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


#Functions estimating A and B

#Group penalty on columns + L1 penalty
AB_L2_L1 = function(Y,A_init,B_init,MA_init,MB_init,
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
                MA=MA_new,MB=MB_new,conv=conv,iterations=counter)
  output
}
                      
AB_L2_L1_tuned = function(Y,A_init,B_init,MA_init,MB_init,
                          lambda_grid,rho,omegas=NULL,
                          eps_abs=1e-3,eps_rel=1e-3,
                          max_iter=1000,adaptive_rho=T,
                          crit="AIC"){
  
  #Set omegas if needed
  if(is.null(omegas)){
    omegas = 1/eigen(crossprod(Y)/t^2)$values
  }
  
  n_lambdas = nrow(lambda_grid)
  lambda_L2_changes = which(diff(lambda_grid[,1])!=0) #Used for updating initialization
  if(ncol(lambda_grid)==1){lambda_L1 = 0}
  A_init_old = A_init; 
  B_init_old = B_init; 
  MA_init_old = MA_init;
  MB_init_old = MB_init;
  coefs = array(NA,dim=c(n,n,4,n_lambdas),
                dimnames = list(row=1:n,col=1:n,
                                mat=c("A","B","MA","MB"),
                                lambda = 1:n_lambdas))
  BICs = AICs = r = conv = iters = rep(0,n_lambdas);
  for(j in 1:n_lambdas){
    lambda_L2 = lambda_grid[j,1]
    if(ncol(lambda_grid)==2){ lambda_L1 = lambda_grid[j,2] }
    
    VECM_obj = AB_L2_L1(Y,A_init,B_init,MA_init,MB_init,
                        lambda_L2,lambda_L1,
                        rho,omegas=omegas,
                        eps_abs=eps_abs,eps_rel=eps_rel,max_iter=max_iter,
                        print_resid = FALSE,adaptive_rho=adaptive_rho)
    coefs[,,1,j] = VECM_obj$A
    coefs[,,2,j] = VECM_obj$B
    coefs[,,3,j] = VECM_obj$MA
    coefs[,,4,j] = VECM_obj$MB
    
    #Update initializers for warm start
    if(j %in% lambda_L2_changes){
      last_lambda_index = min(which(lambda_grid[,1]==lambda_grid[j,1]))
      A_init = coefs[,,1,last_lambda_index]
      B_init = coefs[,,2,last_lambda_index]
      MA_init = coefs[,,3,last_lambda_index]
      MB_init = coefs[,,4,last_lambda_index]
    }else{
      A_init = VECM_obj$A; B_init = VECM_obj$B; 
      MA_init = VECM_obj$MA; MB_init = VECM_obj$MB;
    }
    
    #Compute BICs
    r[j] = VECM_obj$rank
    BICs[j] = BIC_value(Y,VECM_obj$A,VECM_obj$B)
    AICs[j] = AIC_value(Y,VECM_obj$A,VECM_obj$B)
    conv[j] = VECM_obj$conv
    iters[j] = VECM_obj$iterations
    print(j)
  }
  
  #Choose optimal lambda
  if(crit=="AIC"){
    ind_opt = which.min(AICs)
  }else if(crit=="BIC"){
    ind_opt = which.min(BICs)
  }
  lambdas_opt = lambda_grid[ind_opt,]
  A_opt = coefs[,,"A",ind_opt]
  B_opt = coefs[,,"B",ind_opt]
  MA_opt = coefs[,,"MA",ind_opt]
  MB_opt = coefs[,,"MB",ind_opt]
  rank_opt = r[ind_opt]
    
  #return results
  list(lambda_opt = lambdas_opt,ind_opt=ind_opt,rank = rank_opt,
       A = A_opt, B = B_opt,MA = MA_opt,MB_opt = MB_opt,
       conv=conv,iterations=iters,
       As = coefs[,,"A",],Bs = coefs[,,"B",], 
       AICs = AICs,BICs = BICs)
  
}

#Nuclear norm penaly + hierarchical group penalty + L1 penalty

#Functions to update C
C_loss = function(C,U,DY,Y_lag,rho){
  sum((DY - Y_lag%*%t(C))^2)/(nrow(Y_lag)+1) + 
    rho*sum((C-U)^2)/2
}
C_step_size = function(C,nabla_C,U,DY,Y_lag,lambda,rho,
                       step_init=1,step_mult=0.5,step_max_iter = 1000){
  t = nrow(Y_lag)
  step_size = step_init/step_mult
  g_old = C_loss(C,U,DY,Y_lag,rho)
  count = 0
  while(TRUE){
    count = count+1
    step_size = step_mult*step_size
    step_lambda = step_size*lambda
    C_tmp = C - step_size*nabla_C
    C_prox_nuclear = prox_nuclear(C_tmp,step_lambda)
    C_new = C_prox_nuclear$X_new
    G_new = (C - C_new)/step_size
    g_new = C_loss(C_new,U,DY,Y_lag,rho)
    LB = g_old - step_size*sum(nabla_C*G_new) + step_size*sum(G_new^2)/2
    if(!(g_new > LB) | count> step_max_iter){
      break
    }
  }
  list(step_size=step_size,C_new=C_new,d=C_prox_nuclear$d,
       step_iterations=count)
}
C_prox = function(C,U,DY,Y_lag,YY,DYY,lambda,rho,
                  step_size="auto",step_init=1,step_mult=0.5,
                  step_max_iter=100,max_iter = 1000,thresh=1e-6){
  
  C_new = C
  n_par = nrow(C)^2
  iterations = 0
  while(TRUE){
    C_old = C_new
    nabla_C = 2*(C_old%*%YY - DYY) + rho*(C_old - U)
    
    if(step_size == "auto"){
      
      #Compute step size
      step_size_obj = C_step_size(C_old,nabla_C,U,DY,Y_lag,
                                  lambda,rho,step_init,
                                  step_mult,step_max_iter)
      C_new = step_size_obj$C_new
      d_new = step_size_obj$d
      
    }else if(is.numeric(step_size)){
      
      step_lambda = step_size*lambda
      C_tmp = C_old - step_size*nabla_C
      C_prox_obj = prox_nuclear(C_tmp,step_lambda)
      C_new = C_prox_obj$X_new
      d_new = C_prox_obj$d
    }
    
    #Check convergence
    iterations = iterations + 1
    if(norm(C_new-C_old,"F")/nrow(C_old) < thresh){
      convergence=0
      break
    }
    if(iterations > max_iter){
      convergence=1
      break
    }
  }
  list(C = C_new,iterations=iterations,convergence=convergence,d=d_new)
}

#Functions to update A and B
AB_loss = function(S,AB,rho){
  rho*sum((S-AB)^2)/2
}
AB_step_size = function(S,A,B,nabla_A,nabla_B,lambda_L2,lambda_L1=NULL,
                        omegas,rho,step_init=1,step_mult=0.5,
                        max_iter = 1000){
  
  n = nrow(A); n2 = 2*n
  AB = A%*%t(B)
  c = c(rbind(A,B))
  ind_A = c(sapply(1:n,function(j) 1:n + (j-1)*n2))
  ind_B = ind_A + n
  nabla_c = c(rbind(nabla_A,nabla_B))
  g_old = AB_loss(S,AB,rho)
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
    g_new = AB_loss(S,AB_new,rho)
    LB = g_old - step_size*t(nabla_c)%*%G_new + step_size*sum(G_new^2)/2
    if(!(g_new > LB) | count>max_iter){
      break
    }
  }
  list(step_size=step_size,A_new=A_new,B_new=B_new,
       AB_new = AB_new,step_iterations=count)
}

AB_prox = function(A,B,S,rho,lambda_L2,lambda_L1=NULL,omegas,
                   step_size="auto",step_init=1,step_mult=0.5,
                   step_max_iter=100,
                   max_iter=1000,thresh=1e-4){
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
    if(step_size == "auto"){
      
      AB_step_size_obj = AB_step_size(S,A_old,B_old,nabla_A,nabla_B,
                                      lambda_L2,lambda_L1,omegas,rho,
                                      step_init,step_mult,step_max_iter)
      A_new = AB_step_size_obj$A_new
      B_new = AB_step_size_obj$B_new
      AB_new = AB_step_size_obj$AB_new
      
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
    if(sqrt(sum((A_new-A_old)^2) + sum((B_new-B_old)^2))/(2*nrow(A_old)) < thresh){
      convergence=0
      break
    }
    if(iterations > max_iter){
      convergence=1
      break
    }
  }
  list(A=A_new,B=B_new,AB = AB_new,iterations=iterations,convergence=convergence)
}

#ADMM function to estimate the full VECM
AB_NN_L2_L1 = function(C_init,A_init,B_init,M_init,Y,
                         lambda_nuclear,lambda_L2,lambda_L1,
                         rho="auto",rho_init=0.1,rho_mult=0.5,mu=10,
                         step_size = "auto",step_init=1,step_mult=0.5,
                         step_max_iter=100,max_iter = 1000,
                         ADMM_thresh=1e-4,C_thresh=1e-4,AB_thresh=1e-4){
  
  #Transform data
  t=nrow(Y); n = ncol(Y)
  Y_lag = Y[-t,]
  Y_diff = diff(Y)
  Y_lag_cross = crossprod(Y_lag)/t
  DY_Y_lag = crossprod(Y_diff,Y_lag)/t
  omegas = 1/eigen(crossprod(Y)/t^2)$values
  if(lambda_nuclear==0){
    eye = diag(n)
  }
  
  #Initialize objects to update
  C_new = C_init
  A_new = A_init
  B_new = B_init
  AB_new = A_new%*%t(B_new)
  M_new = M_init
  if(rho=="auto"){
    rho_new = rho_init
  }else{
    rho_new = rho
  }
  
  
  
  #Run ADMM algorithm
  counter = 0
  while(TRUE){
    
    counter = counter+1
    C_old = C_new
    A_old = A_new
    B_old = B_new
    AB_old = AB_new
    M_old = M_new
    
    #Update C
    U = AB_old - M_old/rho_new
    if(lambda_nuclear>0){
      C_update = C_prox(C_old,U,Y_diff,Y_lag,
                        Y_lag_cross,DY_Y_lag,
                        lambda_nuclear,rho_new,
                        step_size,step_init,step_mult,
                        step_max_iter,max_iter,C_thresh)
      C_new = C_update$C
    }else{
      C_new = (2*DY_Y_lag + rho*U)%*%solve(2*Y_lag_cross + rho*eye)
    }
    
    
    #Update A,B
    S = C_new + M_old/rho_new
    AB_update = AB_prox(A_old,B_old,S,rho_new,lambda_L2,lambda_L1,omegas,
                        step_size,step_init,step_mult,step_max_iter,
                        max_iter,AB_thresh)
    A_new = AB_update$A
    B_new = AB_update$B
    AB_new = AB_update$AB        
    
    #Update M
    M_new = M_old + rho_new*(C_new - AB_new)
    
    #Check convergence
    if(norm(C_new-C_old,"F")/nrow(C_old) < ADMM_thresh){
      convergence=0
      break
    }
    print(norm(C_new-C_old,"F")/nrow(C_old))
    
    if(counter > max_iter){
      convergence=1
      break
    }
    
    #Update rho
    if(rho=="auto"){
      dual_resid = rho_new*(AB_new - AB_old)
      dual_resid_norm = sqrt(sum(dual_resid^2))
      primal_resid = C_new - AB_new
      primal_resid_norm = sqrt(sum(primal_resid^2))
      if(primal_resid_norm>mu*dual_resid_norm){
        rho_new = rho_new/rho_mult
      }else if(primal_resid_norm < dual_resid_norm/mu){
        rho_new = rho_mult*rho_new
      }
      #print(c(rho_new,dual_resid_norm,primal_resid_norm))
    }
  }
  if(lambda_nuclear>0){
    d = C_update$d
  }else{
    d = n
  }
  list(C=C_new,A=A_new,B=B_new,M=M_new,d = d,
       convergence=convergence,iterations=counter)
}
                       

AB_NN_L2_L1_tuned = function(C_init,A_init,B_init,M_init,Y,
                              lambda_nuclear,lambda_grid,crit="MSE",
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
  BICs = AICs = PI_MSEs = r_C = r_AB = rep(0,n_lambdas,2);
  for(j in 1:n_lambdas){
    lambda_L2 = lambda_grid[j,1]
    if(ncol(lambda_grid)==2){ lambda_L1 = lambda_grid[j,2] }
    
    VECM_obj = AB_NN_L2_L1(C_init = C_init,A_init = A_init,B_init = B_init,
                            M_init = M_init,Y = Y,lambda_nuclear = lambda_nuclear,
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
    BICs[j] = BIC_value(Y,VECM_obj$A,VECM_obj$B)
    AICs[j] = AIC_value(Y,VECM_obj$A,VECM_obj$B)
    PI_MSEs[j] = sqrt(sum((Pi-VECM_obj$A%*%t(VECM_obj$B))^2))
    print(j)
  }
  if(crit=="MSE"){
    ind_opt = which.min(PI_MSEs)
  }else if(crit=="AIC"){
    ind_opt = which.min(AICs)
  }else if(crit=="BIC"){
    ind_opt = which.min(BICs)
  }
  lambdas_opt = lambda_grid[ind_opt,]
  A_opt = coefs[,,"A",ind_opt]
  B_opt = coefs[,,"B",ind_opt]
  C_opt = coefs[,,"C",ind_opt]
  M_opt = coefs[,,"M",ind_opt]
  
  #return results
  list(lambda_opt = lambdas_opt,ind_opt=ind_opt,
       A = A_opt, B = B_opt,C = C_opt, 
       As = coefs[,,"A",],Bs = coefs[,,"B",],M = M_opt, 
       AICs = AICs,BICs = BICs, PI_MSEs = PI_MSEs)
}


#Functions estimating C
C_NN_L1 = function(Y,C_init,M_init,lambda_NN,lambda_L1,rho,
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
  C1_cst = DY_Y_lag%*%C1_den
  
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
    C1_new = (rho/2)*U1%*%C1_den + C1_cst
    
    #Update C2
    U2 = C_old - M2_old/rho
    U2_SVD = svd(U2)
    d_softt = softt_matrix(U2_SVD$d,lambda_NN/rho,omegas)
    C2_new = U2_SVD$u%*%diag(d_softt)%*%t(U2_SVD$v)
    
    #Update C3
    W_tilde = (C1_new + M1_old/rho + C2_new + M2_old/rho)/2
    C_new = softt_matrix(W_tilde,lambda_L1/(2*rho))
    
    #Update M
    M1_new = M1_old + rho*(C1_new - C_new)
    M2_new = M2_old + rho*(C2_new - C_new)
    
    #Compute primal and dual residuals
    dual_resid = rho*(C_new-C_old)
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
                      
