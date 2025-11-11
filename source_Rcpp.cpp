// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

void check(int i){
  Rcpp::Rcout << "Check " << i << std::endl;
}

int count_zero_columns(const arma::mat& B) {
  int count = 0;
  for (arma::uword j = 0; j < B.n_cols; ++j) {
    if (all(B.col(j) == 0)) {
      ++count;
    }
  }
  return count;
}

arma::uvec change_indices(const arma::vec& x) {
  
  // handle edge case
  if (x.n_elem < 2)
    return arma::uvec();  // empty
  
  // compute differences
  arma::vec diffx = arma::diff(x);
  
  // find where diff != 0
  arma::uvec change_pos = arma::find(diffx != 0);
  
  return change_pos;
}

// [[Rcpp::export]]
int rank_Rcpp(arma::mat A){
  
  int rank_A = A.n_cols;
  for(unsigned int j=0; j<A.n_cols; j++){
    if(arma::all(A.col(j)==0)){
      rank_A = rank_A - 1;
    }
  }
  return rank_A;
}

// [[Rcpp::export]]
double BIC_Rcpp(const arma::mat& Y,const arma::mat& A,const arma::mat& B){
  double t = Y.n_rows;
  arma::mat resid = Y.rows(1,t-1) - Y.rows(0,t-2) - Y.rows(0,t-2)*B*A.t();
  double SSR = arma::log_det_sympd(resid.t()*resid/t);
  double penalty = log(t)*(arma::accu(A!=0)+arma::accu(B!=0))/t;
  return SSR + penalty;
}

// [[Rcpp::export]]
double AIC_Rcpp(const arma::mat& Y,const arma::mat& A,const arma::mat& B){
  double t = Y.n_rows;
  arma::mat resid = Y.rows(1,t-1) - Y.rows(0,t-2) - Y.rows(0,t-2)*B*A.t();
  double SSR = arma::log_det_sympd(resid.t()*resid/t);
  double penalty = 2*(arma::accu(A!=0)+arma::accu(B!=0))/t;
  return SSR + penalty;
}

// [[Rcpp::export]]
arma::vec softt_vec(arma::vec x,const double lambda){
  unsigned int n = x.n_elem;
  for(unsigned int j = 0; j<n; j++){
    x(j) = arma::sign(x(j))*std::max(0.0,std::abs(x(j))-lambda);
  }
  return x;
}

// [[Rcpp::export]]
arma::vec prox_L2_Rcpp(arma::vec x,int n,double lambda,arma::vec omegas){
  int n2 = 2*n;
  for(int j=0; j<n; j++){
    int j_lb = j*n2;
    int j_ub = j_lb + n2 - 1;
    double c_mult = 1 - omegas(j)*lambda/std::sqrt(arma::sum(arma::pow(x.subvec(j_lb,j_ub),2)));
    if(c_mult>0){
      x.subvec(j_lb,j_ub) = c_mult*x.subvec(j_lb,j_ub);
    }else{
      x.subvec(j_lb,j_ub) = arma::vec(n2,arma::fill::zeros);
    }
  }
  return x;
}

// [[Rcpp::export]]
arma::vec prox_L2_L1_Rcpp(arma::vec x,int n,double lambda_L2,double lambda_L1,arma::vec omegas){
  int n2 = 2*n;
  for(int j=0; j<n; j++){
    int j_lb = j*n2;
    int j_ub = j_lb + n2 - 1;
    arma::vec x_tmp = softt_vec(x.subvec(j_lb,j_ub),lambda_L1);
    double c_mult = 1 - omegas(j)*lambda_L2/std::sqrt(arma::sum(arma::pow(x_tmp,2)));
    if(c_mult>0){
      x.subvec(j_lb,j_ub) = c_mult*x_tmp;
    }else{
      x.subvec(j_lb,j_ub) = arma::vec(n2,arma::fill::zeros);
    }
  }
  return x;
}

double AB_loss_Rcpp(const arma::mat& DY,const arma::mat& Y_lag,
                    const arma::mat& AB,double lambda_R,
                    const arma::mat& A = arma::mat(), 
                    const arma::mat& B = arma::mat()){
  double loss = arma::as_scalar(arma::accu(arma::pow(DY-Y_lag*AB.t(),2))/Y_lag.n_rows);
  if(lambda_R>0){
    loss = loss + lambda_R*(arma::accu(arma::pow(A,2))+arma::accu(arma::pow(B,2)));
  }
  return loss;
}

Rcpp::List AB_step_size_Rcpp(const arma::mat& DY,const arma::mat& Y_lag,
                             const arma::mat& A,const arma::mat& B,
                             const arma::mat& nabla_A,const arma::mat& nabla_B,
                             const arma::vec& omegas,double lambda_L2,
                             double lambda_L1=0,double lambda_R=0,
                             double step_init=1,double step_mult=0.5,
                             int max_iter = 1000){
  
  int n = A.n_rows; 
  int n2 = 2*n;
  arma::mat AB = A*B.t();
  arma::vec c = arma::vectorise(arma::join_cols(A,B));
  arma::vec nabla_c = arma::vectorise(arma::join_cols(nabla_A,nabla_B));
  
  //Construct indices to extract columns of A and B from c
  arma::uvec ind_A = arma::uvec(n*n,arma::fill::zeros);
  for(int j=0; j<n; j++){
    int a_lb = j*n2;
    ind_A.subvec(j*n,j*n + n-1) = arma::regspace<uvec>(a_lb,a_lb + n-1);
  }
  arma::uvec ind_B = ind_A + n;
  
  double g_old = AB_loss_Rcpp(DY,Y_lag,AB,lambda_R,A,B);
  double step_size = step_init/step_mult;
  int count = 0;
  arma::vec c_new;
  arma::mat A_new, B_new, AB_new;
  int step_convergence;
  while(TRUE){
    count++;
    step_size = step_mult*step_size;
    double step_lambda_L2 = step_size*lambda_L2;
    arma::vec c_tmp = c - step_size*nabla_c;
    
    //proximal update
    if(lambda_L1==0){
      c_new = prox_L2_Rcpp(c_tmp,n,step_lambda_L2,omegas);
    }else{
      double step_lambda_L1 = step_size*lambda_L1;
      c_new = prox_L2_L1_Rcpp(c_tmp,n,step_lambda_L2,step_lambda_L1,omegas);
    }
    
    //#Update parameters
    A_new = arma::reshape(c_new(ind_A),n,n);
    B_new = arma::reshape(c_new(ind_B),n,n);
    AB_new = A_new*B_new.t();
    
    //Backtracking line search
    arma::vec G_new = (c - c_new)/step_size;
    double g_new = AB_loss_Rcpp(DY,Y_lag,AB_new,lambda_R,A_new,B_new);
    double LB = arma::as_scalar(g_old - step_size*nabla_c.t()*G_new + step_size*arma::accu(arma::pow(G_new,2))/2);
    if(!(g_new > LB)){
      step_convergence = 0;
      break;
    }else if(count>max_iter){
      step_convergence = 1;
      break;
    }
  }
  
  Rcpp::List output = Rcpp::List::create(Rcpp::_["step_size"] = step_size,
                                         Rcpp::_["A_new"]=A_new, Rcpp::_["B_new"]=B_new,
                                         Rcpp::_["AB_new"]=AB_new,
                                         Rcpp::_["step_convergence"]=step_convergence);
  return output;
}

// [[Rcpp::export]]
Rcpp::List VECM_SG_Rcpp(const arma::mat& A,const arma::mat& B,
                        const arma::mat& Y,const double lambda_L2,
                        arma::vec omegas,const double lambda_L1=0,
                        const double lambda_R=0,
                        double step_init=1,double step_mult=0.5,
                        double step_max_iter=100,
                        int max_iter=1000,double thresh=1e-5,
                        bool print_dist=true,bool print_loss=false){
  
  //Transform data
  int n = Y.n_cols;
  int t = Y.n_rows;
  arma::mat Y_lag = Y.rows(0,t-2);
  arma::mat DY = Y.rows(1,t-1) - Y_lag;
  t = t-1;
  arma::mat DYYx2 = 2*DY.t()*Y_lag/t;
  arma::mat YYx2 = 2*Y_lag.t()*Y_lag/t;
  
  //Compute omegas if not user-supplied
  if (omegas(0) == -999) {
    arma::vec eigval = arma::eig_sym(Y.t()*Y);
    eigval = arma::reverse(eigval);
    omegas = t*t/eigval;
  }
  
  //Conduct PGD
  arma::vec step_convergence = arma::vec(max_iter);
  arma::mat A_new = A;
  arma::mat A_old = A;
  arma::mat B_old = B;
  arma::mat B_new = B;
  arma::mat AB_new = A*B.t();
  arma::mat AB_old = AB_new;
  arma::mat nabla_A, nabla_B;
  int iterations = 0;
  int convergence = 0;
  while(TRUE){
    Rcpp::checkUserInterrupt();
    
    //Update old coefficients
    A_old = A_new;
    B_old = B_new;
    AB_old = AB_new;
    
    //Update gradient
    if(lambda_R > 0){
      nabla_A = (AB_old*YYx2 - DYYx2)*B_old + 2*lambda_R*A_old;
      nabla_B = (YYx2*AB_old.t() - DYYx2.t())*A_old + 2*lambda_R*B_old;
    }else{
      nabla_A = (AB_old*YYx2 - DYYx2)*B_old;
      nabla_B = (YYx2*AB_old.t() - DYYx2.t())*A_old;
    }
    
    //Perform proximal update with backtracking
    Rcpp::List AB_step_size_obj = AB_step_size_Rcpp(DY,Y_lag,
                                                    A_old,B_old,nabla_A,nabla_B,omegas,
                                                    lambda_L2,lambda_L1,lambda_R,
                                                    step_init,step_mult,step_max_iter);
    A_new = Rcpp::as<arma::mat>(AB_step_size_obj["A_new"]);
    B_new = Rcpp::as<arma::mat>(AB_step_size_obj["B_new"]);
    AB_new = Rcpp::as<arma::mat>(AB_step_size_obj["AB_new"]);
    step_convergence(iterations) = AB_step_size_obj["step_convergence"];
    
    //Check convergence
    iterations++;
    double AB_dist = std::sqrt(arma::accu(arma::pow(AB_new-AB_old,2)))/n;
    if(AB_dist < thresh){
      break;
    }
    if(iterations > max_iter){
      convergence=1;
      break;
    }
    if(print_dist & (iterations % 100 == 0)){
      Rcpp::Rcout << AB_dist << std::endl;
    }
    if(print_loss & (iterations % 100 == 0)){
      double loss_new = AB_loss_Rcpp(DY,Y_lag,AB_new,lambda_R,A_new,B_new);
      Rcpp::Rcout << loss_new << std::endl;
    }
  }
  step_convergence = step_convergence.head(iterations-1);
  
  Rcpp::List output = Rcpp::List::create(Rcpp::_["A"] = A_new,Rcpp::_["B"]=B_new,
                                         Rcpp::_["AB"]=AB_new,
                                         Rcpp::_["iterations"]=iterations,
                                         Rcpp::_["convergence"]=convergence,
                                         Rcpp::_["step_convergence"]=step_convergence);
  return output;
}

// [[Rcpp::export]]
Rcpp::List VECM_SG_tuned_Rcpp(const arma::mat& A,const arma::mat& B,
                              const arma::mat& Y,const arma::mat lambda_grid,
                              arma::vec omegas,
                              const double lambda_R=0,
                              double step_init=1,double step_mult=0.5,
                              double step_max_iter=100,
                              int max_iter=1000,double thresh=1e-5,
                              std::string crit="AIC"){
  
  int t = Y.n_rows;
  int n = Y.n_cols;
  int n_lambdas = lambda_grid.n_rows;
  arma::uvec lambda_L2_changes = change_indices(lambda_grid.col(0));
  double lambda_L1, lambda_L2;
  if(lambda_grid.n_cols == 1 ){
    lambda_L1 = 0;
  }
  
  //Compute omegas if not user-supplied
  if (omegas(0) == -999) {
    arma::vec eigval = arma::eig_sym(Y.t()*Y);
    eigval = arma::reverse(eigval);
    omegas = t*t/eigval;
  }
  
  //Initialize storage objects
  arma::mat A_old = A;
  arma::mat B_old = B;
  arma::cube A_cube = arma::cube(n,n,n_lambdas);
  arma::cube B_cube = arma::cube(n,n,n_lambdas);
  arma::vec BICs = arma::vec(n_lambdas);
  arma::vec AICs = BICs; 
  arma::vec r_AB = BICs;
  
  //Run loop over lambdas
  for(int j=0; j<n_lambdas; j++){
    Rcpp::checkUserInterrupt();
    lambda_L2 = lambda_grid(j,0);
    if(lambda_grid.n_cols == 2){ 
      lambda_L1 = lambda_grid(j,1);
    }
    Rcpp::List VECM_obj = VECM_SG_Rcpp(A_old,B_old,Y,lambda_L2,omegas,lambda_L1,
                                       lambda_R,step_init,step_mult,
                                       step_max_iter,max_iter,thresh,
                                       false,false);
    A_cube.slice(j) = Rcpp::as<arma::mat>(VECM_obj["A"]);
    B_cube.slice(j) = Rcpp::as<arma::mat>(VECM_obj["B"]);
    
    //Update initializers for warm start
    if(arma::any(lambda_L2_changes == j)){
      int last_lambda_index = arma::min(arma::find(lambda_grid.col(0)==lambda_grid(j-1,0)));
      A_old = A_cube.slice(last_lambda_index);
      B_old = B_cube.slice(last_lambda_index);
      // return Rcpp::List::create(Rcpp::_["lambda_index"] = last_lambda_index,
      //                           Rcpp::_["lambda_changes"] = lambda_L2_changes,
      //                           Rcpp::_["j"]=j,
      //                           Rcpp::_["A_old"]=A_old,
      //                           Rcpp::_["B_old"]=B_old);
    }else{
      A_old = A_cube.slice(j);
      B_old = B_cube.slice(j);
    }
    
    //Compute BICs
    r_AB(j) = rank_Rcpp(A_cube.slice(j));
    BICs(j) = BIC_Rcpp(Y,A_cube.slice(j),B_cube.slice(j));
    AICs(j) = AIC_Rcpp(Y,A_cube.slice(j),B_cube.slice(j));
    Rcpp::Rcout << j << std::endl;
  }
  
  //Extract tuned coefficients
  int ind_opt=0;
  if(crit=="AIC"){
    ind_opt = arma::index_min(AICs);
  }else if(crit=="BIC"){
    ind_opt = arma::index_min(BICs);
  }
  arma::rowvec lambdas_opt = lambda_grid.row(ind_opt);
  
  //return results
  Rcpp::List output = Rcpp::List::create(Rcpp::_["lambda_opt"]=lambdas_opt,
                                         Rcpp::_["ind_opt"]=ind_opt,
                                         Rcpp::_["A"] = A_cube.slice(ind_opt),
                                         Rcpp::_["B"] = B_cube.slice(ind_opt),
                                         Rcpp::_["r"] = r_AB(ind_opt),
                                         Rcpp::_["AICs"] = AICs,
                                         Rcpp::_["BICs"] = BICs,
                                         Rcpp::_["ranks"] = r_AB,
                                         Rcpp::_["As"] = A_cube,
                                         Rcpp::_["Bs"] = B_cube);
  return output;
}

// [[Rcpp::export]]
Rcpp::List AB_L2_L1_Rcpp(const arma::mat& Y,arma::mat A_init,arma::mat B_init,
                         arma::mat MA_init,arma::mat MB_init,
                         const double lambda_L2,const double lambda_L1,
                         double rho,arma::vec omegas,
                         const double eps_abs=1e-3,const double eps_rel=1e-3,
                         const int max_iter=1000,const bool print_resid = false,
                         const bool adaptive_rho = true){
  
  //Transform data
  int t=Y.n_rows;
  int n = Y.n_cols;
  int n_sq = n*n;
  arma::mat Y_lag = Y.rows(0,t-2);
  arma::mat Y_diff = Y.rows(1,t-1) - Y_lag;
  arma::mat Y_lag_cross = Y_lag.t()*Y_lag/t;
  arma::mat DY_Y_lag = Y_diff.t()*Y_lag/t;

  //Compute omegas if not user-supplied
  if (omegas(0) == -999) {
    arma::vec eigval = arma::eig_sym(Y.t()*Y);
    eigval = arma::reverse(eigval);
    omegas = t*t/eigval;
  }

  //Transform desired accuracies
  double eps_abs_scaled = eps_abs*n;
  
  //Initialize objects to update
  arma::mat A1_new = A_init;
  arma::mat A_new = A_init;
  arma::mat B1_new = B_init;
  arma::mat B_new = B_init;
  arma::mat C_new = join_cols(A_init,B_init);
  arma::mat MA_new = MA_init;
  arma::mat MB_new = MB_init;
  
  arma::mat rho_div_2(n,n);
  arma::mat rho_div_2_kron(n_sq,n_sq);

  //Run ADMM algorithm
  int counter = 0;
  int conv = 0;
  while(true){
    counter++;
    
    //rho quantities
    double lambda_L1_div_rho = lambda_L1/rho;
    rho_div_2.diag() = arma::vec(n,arma::fill::value(rho/2));
    rho_div_2_kron.diag() = arma::vec(n_sq,arma::fill::value(rho/2));

    //Update old parameters
    arma::mat A1_old = A1_new;
    arma::mat B1_old = B1_new;
    arma::mat A_old = A_new;
    arma::mat B_old = B_new;
    arma::mat C_old = C_new;
    arma::mat MA_old = MA_new;
    arma::mat MB_old = MB_new;
      
    //Update A1
    arma::mat UA = A_old - MA_old/rho;
    A1_new = (DY_Y_lag*B1_old + (rho/2)*UA)*
        arma::inv(B1_old.t()*Y_lag_cross*B1_old +rho_div_2);
    A1_new.elem(find(abs(A1_new) < 1e-10)).zeros();

    //Update B1
    arma::mat UB = B_old - MB_old/rho;
    arma::vec B1_new_vec = arma::inv(arma::kron(A1_new.t()*A1_new,Y_lag_cross) + rho_div_2_kron)*
      arma::vectorise(DY_Y_lag.t()*A1_new + (rho/2)*UB);
    B1_new = arma::reshape(B1_new_vec,n,n);
    B1_new.elem(find(abs(B1_new) < 1e-10)).zeros();

    //Update A and B jointly
    arma::mat C1_new = arma::join_cols(A1_new,B1_new);
    arma::mat M_old = arma::join_cols(MA_old,MB_old);
    arma::mat WC = C1_new + M_old/rho;
    for(int j = 0; j<n; j++){
      arma::vec w_j = WC.col(j);
      arma::vec c_softt_j = softt_vec(w_j,lambda_L1_div_rho);
      double c_mult = 1-(lambda_L2*omegas(j)/rho)/std::sqrt(arma::sum(arma::pow(c_softt_j,2)));
      if(c_mult>0){
        C_new.col(j) = c_mult*c_softt_j;
      }else{
        C_new.col(j).zeros();
      }
    }
    A_new = C_new.rows(0,(n-1));
    B_new = C_new.rows(n,2*n-1);

    //Update M
    MA_new = MA_old + rho*(A1_new - A_new);
    MB_new = MB_old + rho*(B1_new - B_new);
    arma::mat M_new = arma::join_cols(MA_new,MB_new);
        
    //Compute primal and dual residuals
    arma::mat dual_resid = rho*(C_new-C_old);
    double dual_resid_norm = std::sqrt(arma::accu(arma::pow(dual_resid,2)));
    arma::mat primal_resid = C1_new - C_new;
    double primal_resid_norm = std::sqrt(arma::accu(arma::pow(primal_resid,2)));
    if(print_resid){
      Rcpp::Rcout << primal_resid_norm << " " <<  dual_resid_norm << std::endl;
    }
          
    //Check convergence
    arma::vec C_norms = {std::sqrt(arma::accu(arma::pow(C1_new,2))),
                         std::sqrt(arma::accu(arma::pow(C_new,2)))};
    double eps_primal = eps_abs_scaled + eps_rel*arma::max(C_norms);
    double M_norm = std::sqrt(arma::accu(arma::pow(M_new,2)));
    double eps_dual = eps_abs_scaled + eps_rel*M_norm;
    
    if((primal_resid_norm < eps_primal) & (dual_resid_norm < eps_dual)){
      break;
    }
    
    if(counter == max_iter){
      conv=1;
      break;
    }
          
    //Update rho
    if(adaptive_rho){
      if(primal_resid_norm > 10*dual_resid_norm){
        rho = 2*rho;
      }else if(dual_resid_norm > 10*primal_resid_norm){
        rho = rho/2;
      }
    }
  }
  
  int rank = count_zero_columns(C_new);
  Rcpp::List output = Rcpp::List::create(Rcpp::_["A"] = A_new,
                                         Rcpp::_["B"] = B_new,
                                         Rcpp::_["A1"] = A1_new,
                                         Rcpp::_["B1"] = B1_new,
                                         Rcpp::_["rank"] = rank,
                                         Rcpp::_["MA"] = MA_new,
                                         Rcpp::_["MB"] = MB_new,
                                         Rcpp::_["conv"] = conv,
                                         Rcpp::_["iterations"] = counter);
    return output;
}

// [[Rcpp::export]]
Rcpp::List AB_L2_L1_tuned_Rcpp(const arma::mat& Y,arma::mat A_init,arma::mat B_init,
                          arma::mat MA_init,arma::mat MB_init,
                          const arma::mat lambda_grid,
                          double rho,arma::vec omegas,
                          const double eps_abs=1e-3,const double eps_rel=1e-3,
                          const int max_iter=1000,const bool adaptive_rho = true,
                          const std::string crit = "AIC"){
  
  int t = Y.n_rows;
  int n = Y.n_cols;
  int n_lambdas = lambda_grid.n_rows;
  arma::uvec lambda_L2_changes = change_indices(lambda_grid.col(0));
  double lambda_L1, lambda_L2;
  if(lambda_grid.n_cols == 1 ){
    lambda_L1 = 0;
  }
  
  //Compute omegas if not user-supplied
  if (omegas(0) == -999) {
    arma::vec eigval = arma::eig_sym(Y.t()*Y);
    eigval = arma::reverse(eigval);
    omegas = t*t/eigval;
  }

  arma::field<arma::cube> coefs(n_lambdas);
  arma::vec BICs = arma::vec(n_lambdas);
  arma::vec AICs = BICs;
  arma::vec r = BICs;
  arma::vec conv = BICs;
  arma::vec iters = BICs;
  arma::cube coef_j = arma::cube(n,n,4);
  for(int j=0; j<n_lambdas; j++){
    
    lambda_L2 = lambda_grid(j,0);
    if(lambda_grid.n_cols==2){ 
      lambda_L1 = lambda_grid(j,1);
    }
    
    Rcpp::List VECM_obj = AB_L2_L1_Rcpp(Y,A_init,B_init,MA_init,MB_init,
                                        lambda_L2,lambda_L1,rho,omegas,
                                        eps_abs,eps_rel,max_iter,false,adaptive_rho);
                                   
    
    coef_j.slice(0) = Rcpp::as<arma::mat>(VECM_obj["A"]);
    coef_j.slice(1) = Rcpp::as<arma::mat>(VECM_obj["B"]);
    coef_j.slice(2) = Rcpp::as<arma::mat>(VECM_obj["MA"]);
    coef_j.slice(3) = Rcpp::as<arma::mat>(VECM_obj["MA"]);
    coefs(j) = coef_j;
    
    //Update initializers for warm start
    if(arma::any(lambda_L2_changes == j)){
      int last_lambda_index = arma::min(arma::find(lambda_grid.col(0)==lambda_grid(j-1,0)));
      A_init = coefs(last_lambda_index).slice(0);
      B_init = coefs(last_lambda_index).slice(1);
      MA_init = coefs(last_lambda_index).slice(2);
      MB_init = coefs(last_lambda_index).slice(3);
    }else{
      A_init = coef_j.slice(0);
      B_init = coef_j.slice(1);
      MA_init = coef_j.slice(2);
      MB_init = coef_j.slice(3);
    }
    
    //Compute BICs
    r(j) = VECM_obj["rank"];
    BICs(j) = BIC_Rcpp(Y,coef_j.slice(0),coef_j.slice(1));
    AICs(j) = AIC_Rcpp(Y,coef_j.slice(0),coef_j.slice(1));
    conv(j) = VECM_obj["conv"];
    iters(j) = VECM_obj["iterations"];
    Rcpp::Rcout << j << std::endl;
  }
  
  //Extract tuned coefficients
  unsigned int ind_opt=0;
  if(crit=="AIC"){
    ind_opt = arma::index_min(AICs);
  }else if(crit=="BIC"){
    ind_opt = arma::index_min(BICs);
  }
  arma::rowvec lambdas_opt = lambda_grid.row(ind_opt);
  
  //return results
  Rcpp::List output = Rcpp::List::create(Rcpp::_["lambda_opt"]=lambdas_opt,
                                         Rcpp::_["ind_opt"]=ind_opt,
                                         Rcpp::_["A"] = coefs(ind_opt).slice(0),
                                         Rcpp::_["B"] = coefs(ind_opt).slice(1),
                                         Rcpp::_["r"] = r(ind_opt),
                                         Rcpp::_["AICs"] = AICs,
                                         Rcpp::_["BICs"] = BICs,
                                         Rcpp::_["ranks"] = r,
                                         Rcpp::_["coefs"] = coefs);
  return output;
}