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

#VECM-ADMM
VECM_ADMM = function(Y,rho,lambda_L1,lambda_L2,lambda_F,A_init=NULL,B_init=NULL,
                     adaptive=T,max_iter = 1000,thresh=1e-4){
  
  t = nrow(Y)
  n = ncol(Y)
  
  #compute initial parameters
  if(is.null(A_init) | is.null(B_init)){
    colnames(Y) = paste("x",1:n,sep="")
    vecm_fit = ca.jo(Y,K=2)
    r = choose_r(vecm_fit,0.05)
    if(is.null(A_init)){
      A_init = vecm_fit@W
    }
    if(is.null(B_init)){
      B_init = vecm_fit@V
    }
  }
  
  #Declare storage object
  start_indices = c(rep(1,4),2:n)
  local_names = c("LS","L2",paste("G",rep(1:n),sep=""))
  As_k = MAs_k = Bs_k = MBs_k = list()
  Zeros = matrix(0,n,n)
  for(k in 1:length(start_indices)){
    j = start_indices[k]
    As_k[[k]] = A_init[,j:n]
    MAs_k[[k]] = Zeros[,j:n]
    Bs_k[[k]] = B_init[,j:n]
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
  
  #Compute statistics used in ADMM
  Y_lag = Y[-t,]
  Y_diff = diff(Y)
  t = t-1
  Y_lag_cross = crossprod(Y_lag)/t
  DY_Y_lag = crossprod(Y_diff,Y_lag)/t
  I_rhod2 = diag(rho/2,n)
  I_rhod2_kron = diag(rho/2,n^2)
  
  #Compute weights
  if(adaptive){
    eigen_values = eigen(crossprod(Y)/t^2)$values
    omegas = 1/eigen_values
    # eigen_sums = cumsum(eigen_values[n:1])[n:1]
    # omegas = 1/eigen_sums
  }else{
    omegas=rep(1,n)
  }
  
  #ADMM iterations
  counter = 0
  loss = loss_function(Y_diff,Y_lag,A_init,B_init,lambda_L1,lambda_L2,lambda_F,omegas)
  loss_vec = loss_diff_vec = rep(1,max_iter)
  while(TRUE){
    
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
    # As_kp1$A = softt_matrix(A_kp1_sums,lambda_L1/rho)%*%W_j
    # Bs_kp1$B = softt_matrix(B_kp1_sums,lambda_L1/rho)%*%W_j
    As_kp1$A = softt_matrix(A_kp1_sums%*%W_j,lambda_L1/rho)
    Bs_kp1$B = softt_matrix(B_kp1_sums%*%W_j,lambda_L1/rho)
    
    #Update MA and MB
    MAs_kp1$MA_LS = MAs_k$MA_LS + rho*(As_kp1$A_LS - As_kp1$A)
    MAs_kp1$MA_L2 = MAs_k$MA_L2 + rho*(As_kp1$A_L2 - As_kp1$A)
    MBs_kp1$MB_LS = MBs_k$MB_LS + rho*(Bs_kp1$B_LS - Bs_kp1$B)
    MBs_kp1$MB_L2 = MBs_k$MB_L2 + rho*(Bs_kp1$B_L2 - Bs_kp1$B)
    for(j in 1:n+3){
      MAs_kp1[[j]] = MAs_k[[j]] + rho*(As_kp1[[j]] - As_kp1$A[,(j-3):n])
      MBs_kp1[[j]] = MBs_k[[j]] + rho*(Bs_kp1[[j]] - Bs_kp1$B[,(j-3):n])
    }
    
    #Evaluate loss and stopping criterion
    loss = loss_function(Y_diff,Y_lag,As_kp1$A,Bs_kp1$B,lambda_L1,lambda_L2,lambda_F,omegas)
    loss_vec[counter] = loss
    loss_diff = (loss-loss_old)/loss_old
    loss_diff_vec[counter] = loss_diff
    if(counter > 3){
      if(all(abs(loss_diff_vec[(counter-3):counter]) < thresh) | counter > (max_iter-1)){
        break
      }
    }
    
    #Update old values
    As_k = As_kp1; Bs_k = Bs_kp1
    MAs_k = MAs_kp1; MBs_k = MBs_kp1
  }
  
  #Determine rank
  # zero_ind = which(sapply(1:n,function(j){
  #   all(c(As_kp1$A[,j:n],Bs_kp1$B[,j:n])==0)
  # }))
  zero_ind = which(sapply(1:n,function(j){
    all(c(As_kp1[[j+3]],Bs_kp1[[j+3]])==0)
  }))
  if(sum(zero_ind)>0){
    r_hat = min(zero_ind) - 1
  }else{
    r_hat = n
  }
  
  list(A = As_kp1$A,B = Bs_kp1$B, As = As_kp1, Bs = Bs_kp1,
       r = r_hat,omegas=omegas,
       losses = loss_vec,improvements = loss_diff_vec,iterations = counter)
}
