//#########################################################
//#                                                       #
//#            Bliss method : rcpp code                   #
//#                                                       #
//#########################################################
#define ARMA_DONT_PRINT_ERRORS
// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <string>
#include <iostream>
#include <vector>
#include <cstring>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//#############################################################################
//############################ basic functions  ###############################
//#############################################################################

// The R function : ginv (generalized matrix inversion using SVD decomposition)
// [[Rcpp::export]]
arma::mat ginv_cpp (arma::mat & x, double tol){
  int p;
  mat u;
  vec s;
  mat v;

  svd(u,s,v,x);
  p = s.size();

  tol   = tol * s(0);
  mat S = zeros<mat>(p,p);

  for( unsigned i=0 ; i<p; ++i){
    if( s(i) > tol ) S(i,i) = 1/s(i);
  }

  return( v * (S * trans(u)) );
}

// Sample in 0:n-1.
int sample_cpp (int n){
 double u = R::runif(0,1) * n;
 int res = trunc(u) ;
 if(res == n) res = n-1;
 return res;
}

// Vector sample in 0:n-1.
arma::vec sample_cpp (int nbre, int n){
 vec res = zeros<vec>(nbre);
 for (unsigned i = 0 ; i < nbre ; ++i) {
  res(i) = sample_cpp(n);
 }
 return res;
}

// Weighted sample in 0:n-1.
int sample_weight (arma::vec proba){
  if(sum(proba) == 0)   proba = ones<vec>(proba.n_rows) ;
  vec proba_cum = cumsum(proba)/sum(proba) ;

  int ret = 0;
  double u = R::runif(0,1);
  while (ret <= proba_cum.n_rows and u > proba_cum(ret)) {
    ret++;
  }
  return ret;
}

// Vector weighted sample in 0:n-1.
arma::vec sample_weight (int n, arma::vec proba){
 vec res = zeros<vec>(n);
 for (unsigned i = 0 ; i < n ; ++i) {
  res(i) = sample_weight(proba);
 }
 return res;
}

// Return the vector vec[-k].
arma::vec vec_drop_k(arma::vec vecteur, int k){
 vecteur.shed_row(k);
 return vecteur;
}

// Return the matrix mat[,-k].
arma::mat mat_drop_col_k(arma::mat matrix, int k){
 matrix.shed_col(k);
 return matrix;
}

// seq.
arma::vec sequence(int a,int b,double by){
  int range = floor((b-a)/by + 1) ;
  vec res = zeros<vec>(range);
  for(int i=0 ; i<range ; i++){
    res(i) = a + i*by;
  }
  return res;
}

// Extract an element from a cube.
double cube_extract(NumericVector & cube, int x , int y, int z, arma::vec & dims){
  double res;
  res = cube[x + y*dims(0) + z*dims(1)*dims(0)];
  return res;
}

// Compute the square root matrix using the SVD decomposition
// [[Rcpp::export]]
arma::mat sqrt_mat (arma::mat & x){
  int p;
  mat u;
  vec s;
  mat v;

  svd(u,s,v,x);
  p = s.size();

  mat S = zeros<mat>(p,p);

  for( unsigned i=0 ; i<p; ++i){
    S(i,i) = sqrt(s(i));
  }

  return( u * (S * trans(v)) );
}

// Simulate from a multidimensional gaussian.
// [[Rcpp::export]]
arma::vec mvrnormArma(arma::vec mu, arma::mat VarCovar, double sigma_sq) {
  int ncols = VarCovar.n_cols;
  vec y = randn<vec>(ncols);
  VarCovar = chol(VarCovar); // xxxxxxxxxxxx

  return  mu + sqrt(sigma_sq) * trans(trans(y) * VarCovar);
}

// Compute a trapezoidal approximation of area under curve.
// [[Rcpp::export]]
double integrate_trapeze_cpp (arma::vec & x, arma::vec & y){
 vec diff_x = vec_drop_k(x,0) - vec_drop_k(x,x.size()-1);
 vec cumu_y = vec_drop_k(y,0) + vec_drop_k(y,y.size()-1);
  return sum( diff_x % cumu_y  )/2 ;
}

//  Compute the norm of a vector.
double norm_fct(arma::vec & x,arma::vec & y){
  vec tmp = zeros<vec>(x.size());
  for(int i=0 ; i<tmp.size() ; ++i){
    tmp(i) = pow( y(i),2 );
  }
  double res;
  res = sqrt(integrate_trapeze_cpp(x,tmp));

  return res;
}

// Use to compute an uniform function, see function compute_beta.
// [[Rcpp::export]]
arma::vec uniform_cpp (int m, int l, arma::vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  vec index = sequence(m-l,m+l,1);
  int tmp;
  for(int i=0 ; i<index.size() ; ++i){
    tmp = index(i);
    if( (tmp <= p) && (tmp >= 1) ){
      res(index(i)-1) = 1;
    }
  }
  double res_norm = norm_fct(grid,res);
  res = res / res_norm ;
  return res;
}

// Use to compute a triangular function, see function compute_beta.
// [[Rcpp::export]]
arma::vec triangular_cpp (int m, int l, arma::vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  1 - i/l_double ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
  double res_norm = norm_fct(grid,res);
  res = res / res_norm ;
  return res;
}

// Use to compute a gaussian function, see function compute_beta.
// [[Rcpp::export]]
arma::vec gaussian_cpp (int m, int l, arma::vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  exp( - 9*pow(i/l_double,2)/2) ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
  double res_norm = norm_fct(grid,res);
  res = res / res_norm ;
  return res;
}

// Use to compute a gaussian function, see function compute_beta.
// [[Rcpp::export]]
arma::vec Epanechnikov_cpp (int m, int l, arma::vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  1-pow(i/l_double,2) ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
  double res_norm = norm_fct(grid,res);
  res = res / res_norm ;
  return res;
}

// compute_beta in cpp.
// [[Rcpp::export]]
arma::vec compute_beta_cpp (arma::vec & b, arma::vec & m, arma::vec & l,
                          arma::vec & grid, int p, int K, std::string basis, 
                          arma::mat & normalization_values ){
  vec res = zeros<vec>(p) ;

  if(basis == "Uniform"){
    for(int i=0 ; i<K ; ++i){
      res = res + b(i)/normalization_values( m(i)-1 , l(i)-1 ) *
        uniform_cpp(m(i),l(i),grid);
    }
  }
  if(basis == "Triangular"){
    for(int i=0 ; i<K ; ++i){
      res = res + b(i)/normalization_values( m(i)-1 , l(i)-1 )  *
        triangular_cpp(m(i),l(i),grid);
    }
  }
  if(basis == "Gaussian"){
    for(int i=0 ; i<K ; ++i){
      res = res + b(i)/normalization_values( m(i)-1 , l(i)-1 )  *
        gaussian_cpp(m(i),l(i),grid);
    }
  }
  if(basis == "Epanechnikov"){
    for(int i=0 ; i<K ; ++i){
      res = res + b(i)/normalization_values( m(i)-1 , l(i)-1 )  *
        Epanechnikov_cpp(m(i),l(i),grid);
    }
  }
  return res;
}

// Compute the functions beta_i for each iteration i.
// [[Rcpp::export]]
arma::mat compute_beta_sample_cpp (arma::mat &  trace, int p, int K, arma::vec & grid,
                                std::string basis, arma::mat & normalization_values){
  mat res = zeros<mat>(trace.n_rows,p);
  vec b;
  vec m   ;
  vec l   ;
  vec tmp ;
  
  for(int i=0 ; i<res.n_rows ; ++i){
    tmp  = trans(trace.row(i))     ;
    b = tmp.subvec(0,K-1)    ;
    m    = tmp.subvec(K,2*K-1)   ;
    l    = tmp.subvec(2*K,3*K-1) ;

    res.row(i) = trans(compute_beta_cpp(b,m,l,grid,p,K,basis,normalization_values)) ;
  }
  return res ;
}

// Compute all the alternative for the value of the intergral for all m and l.
// [[Rcpp::export]]
arma::cube potential_intervals_List(List & x_list, List & grids,arma::vec & l_max_vec,
                                    CharacterVector & basis_vec, int q){
  mat x = as<mat>(x_list[q]);
  vec grid = as<vec>(grids[q]);
  int l_max = l_max_vec(q);
  std::string basis = as<std::string>(basis_vec(q));

  int n = x.n_rows ;
  int p = x.n_cols ;
  vec tub;

  arma::cube res(p,l_max,n+1);
  vec tmp;
  vec x_tmp;
  vec tmp2;
  if(basis == "Uniform"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = uniform_cpp(i+1,j+1,grid);
          x_tmp = trans(x.row(k)) ;
          tmp2  =  x_tmp % tmp ;

          res(i,j,k) = integrate_trapeze_cpp(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "Triangular"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = triangular_cpp(i+1,j+1,grid);
          x_tmp = trans(x.row(k)) ;
          tmp2  =  x_tmp % tmp ;

          res(i,j,k) = integrate_trapeze_cpp(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "Gaussian"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = gaussian_cpp(i+1,j+1,grid);
          x_tmp = trans(x.row(k)) ;
          tmp2  =  x_tmp % tmp ;

          res(i,j,k) = integrate_trapeze_cpp(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "Epanechnikov"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = Epanechnikov_cpp(i+1,j+1,grid);
          x_tmp = trans(x.row(k)) ;
          tmp2  =  x_tmp % tmp ;

          res(i,j,k) = integrate_trapeze_cpp(grid, tmp2 );
        }
      }
    }
  }

  for(int i=0 ; i<p ; ++i){
    for(int j=0 ; j<l_max ; ++j){
      // normalize by the scale \hat{s}_k
      tub = res.tube(i,j);
      tub = tub.subvec(0,n-1) ;
      res(i,j,n) = stddev( tub );
      for( int k=0 ; k<n ; ++k){
        res(i,j,k) = res(i,j,k) / res(i,j,n);
      }
    }
  }
  return res;
}

// Find the first element of a vector which equals to n.
int which_first (NumericVector & v, int n){
  for(int i=0 ; i<v.size() ; ++i){
    if( v(i)==n ) return i;
  }
  return -1;
}

// Compute a moving average on the vector v.
// [[Rcpp::export]]
arma::vec moving_average_cpp (arma::vec & v, int range){
  int n = v.size();
  vec res = zeros<vec>(n) ;
  int b_inf;
  int b_sup;

  for( int i=0; i<n; ++i){
    if(i-range < 0  ) b_inf = 0   ; else b_inf = i-range;
    if(i+range > n-1) b_sup = n-1 ; else b_sup = i+range;
    res(i) = mean( v.subvec(b_inf,b_sup) )  ;
  }

  return res;
}

//#############################################################################
//############################ Auxiliary functions ############################
//#############################################################################

// Compute the matrix V (for a Ridge Zellner prior)
// (for Q functional covaribles)
arma::mat compute_W_inv_List (int Q, arma::vec & K, double g, arma::mat & x_tilde, 
                              int sum_K, arma::mat & lambda_id0){
 mat W_inv = zeros<mat>(sum_K+1,sum_K+1);
 mat lambda_id = lambda_id0 ;
 
 mat x_tilde_temp = mat_drop_col_k(x_tilde,0);
 mat u;
 vec s;
 mat v;
 svd(u,s,v,x_tilde_temp);
 
 W_inv(0,0) = 1/lambda_id(0,0);
 for( unsigned i=1 ; i<sum_K+1; ++i){
  lambda_id(i,i) = lambda_id(i,i) * max(s); // try with min 
 }
 
 int count = 0;
 for( unsigned q=0 ; q<Q ; ++q){
  W_inv.submat(1+count,1+count,K(q)+count,K(q)+count) =
   ( trans(x_tilde.cols(1+count,K(q)+count)) *
   x_tilde.cols(1+count,K(q)+count) +
   lambda_id.submat(1+count,1+count,K(q)+count,K(q)+count) )  /g;
  count = count + K(q);
 }
 
 return W_inv;
}

// Extract a subvector from the cube potential_intervals with a m_k and a l_k.
// [[Rcpp::export]]
arma::vec potential_intervals_extract (NumericVector & potential_intervals, int mk , 
                                       int lk, arma::vec & dims) {
  vec res = zeros<vec>(dims(2));
  for (int i = 0; i < dims(2); i++) {
    res(i) = cube_extract(potential_intervals, mk - 1, lk - 1, i, dims);
  }
  return res;
}

// Extract a submatrix from the cube potential_intervals with the vectors m and l.
arma::mat extraire_x_tilde(arma::vec & m, arma::vec & l, NumericVector & potential_intervals,
                           arma::vec & dims){
  int K = m.size();
  mat res = ones<mat>(dims(2), K + 1);
  for(int i=0; i<K ; i++){
    res.col(i+1) = potential_intervals_extract(potential_intervals,m(i),l(i),dims);
  }
  return res;
}


// Update the parameter m_k
// [[Rcpp::export]]
void update_mqk (int count, int k, arma::vec & y, arma::vec & b_tilde, double sigma_sq,
                 arma::vec & m_q, arma::vec & l_q, arma::mat x_tilde,
                 NumericVector & potential_intervals_q, arma::vec & potential_intervals_dims_q,
                 arma::vec & m_alternative_q, int p_q, int Q,
                 arma::vec K, double g, int sum_K, arma::mat & lambda_id0) {
  vec aux = zeros<vec>(p_q);
  vec aux2 = zeros<vec>(p_q);
  mat W_inv_temp;

  // Compute the probabilities
  vec probs = ones<vec>(p_q);
  vec x_tilde_qki = zeros<vec>(potential_intervals_dims_q(2)) ;
  for(int  i=0 ; i<p_q ; ++i){
    x_tilde_qki = potential_intervals_extract(potential_intervals_q,m_alternative_q(i),
                                        l_q(k),potential_intervals_dims_q);

   x_tilde.col(count + k + 1) = x_tilde_qki;
   aux(i) = dot( y - x_tilde * b_tilde ,
       y - x_tilde * b_tilde ) /(2*sigma_sq) ;

    W_inv_temp = compute_W_inv_List(Q,K,g,x_tilde,sum_K,lambda_id0);
    
    aux(i) = aux(i) +
     dot( b_tilde , W_inv_temp * b_tilde ) / (2*sigma_sq);
    aux2(i) = sqrt( det(W_inv_temp) );
  }

  double min_aux = min(aux);
  for(int  i=0 ; i<p_q ; ++i){
    aux(i)  = aux(i) - min_aux;
    probs(i) = aux2(i) * exp( - aux(i) ) ;
  }
  // Simulate a mk
  m_q(k) = sample_weight(probs) + 1 ;
}


// Update the parameter l_k
// [[Rcpp::export]]
void update_lqk (int count, int k, arma::vec & y, arma::vec & b_tilde, double sigma_sq,
                 arma::vec & m_q, arma::vec & l_q, arma::mat x_tilde,
                 NumericVector & potential_intervals_q, arma::vec & potential_intervals_dims_q,
                 arma::vec & l_alternative_q, arma::vec & phi_l_q, int l_values_length_q, 
                 int Q, arma::vec K, double g, int sum_K, arma::mat & lambda_id0) {
  vec aux = zeros<vec>(l_values_length_q);
  vec aux2 = zeros<vec>(l_values_length_q);
  mat W_inv_temp;

  // Compute the probabilities
  vec probs = ones<vec>(l_values_length_q);
  vec x_tilde_qki = zeros<vec>(potential_intervals_dims_q(2)) ;
  for(int  i=0 ; i<l_values_length_q ; ++i){
    x_tilde_qki = potential_intervals_extract(potential_intervals_q,m_q(k),
                                        l_alternative_q(i),
                                        potential_intervals_dims_q);

    x_tilde.col(count + k + 1) = x_tilde_qki;
    aux(i) = dot( y - x_tilde * b_tilde ,
         y - x_tilde * b_tilde ) /(2*sigma_sq) ;
    
    W_inv_temp = compute_W_inv_List(Q,K,g,x_tilde,sum_K,lambda_id0);
    
    aux(i) = aux(i) +
     dot( b_tilde , W_inv_temp * b_tilde ) / (2*sigma_sq);
    aux2(i) = sqrt( det(W_inv_temp) );
  }

  double min_aux = min(aux);
  for(int  i=0 ; i<l_values_length_q ; ++i){
    aux(i)  = aux(i) - min_aux;
    probs(i) = aux2(i) * exp( - aux(i) ) * phi_l_q(i);
  }
  // Simulate a lk
  m_q(k) = sample_weight(probs) + 1 ;
}

// update the parameter sigma_sq
void update_sigma_sq (arma::vec & y, arma::vec & b_tilde, arma::mat & W_inv, 
                      arma::mat & x_tilde, int n, int sum_K, double & sigma_sq) {
  vec y_tmp     = y - x_tilde * b_tilde ;
  double y_tmp2 = dot(y_tmp,y_tmp) ;
  double b_tilde_tmp = dot(b_tilde, W_inv * b_tilde) ;

  double a_star = (sum_K+n+1)/2 ;
  double b_star = 0.5*( y_tmp2 + b_tilde_tmp);
  
  sigma_sq = 1. / (R::rgamma(a_star, 1/b_star) );
}

// update the parameter b
// [[Rcpp::export]]
void update_b_tilde (arma::vec & y, double sigma_sq, arma::mat & x_tilde, 
                          arma::mat & Sigma_b_tilde_inv, double tol,
                          arma::vec & b_tilde) {
  vec mu_b_tilde = trans(x_tilde) * y;
  b_tilde = mvrnormArma( ginv_cpp(Sigma_b_tilde_inv,tol) * mu_b_tilde ,
                         ginv_cpp(Sigma_b_tilde_inv,tol), sigma_sq);
}

// Compute the loss function for a proposal d
// [[Rcpp::export]]
double loss_cpp (arma::vec & d, arma::vec & grid, arma::vec & posterior_expe){
  vec tmp  = d-posterior_expe ;

  return pow(norm_fct(grid, tmp),2 );
}

// Compute the decrease of the Temperature
double cooling_cpp (int i, double Temp){
  double res = Temp / log( ( i / 10)*10 + exp(1));
  return res;
}

void update_x_tilde (int Q, arma::vec & K, List & potential_intervals,
                    List & potential_intervals_dims, List & m, List & l,
                    arma::mat & x_tilde){
  int count = 0;
  for( unsigned q=0 ; q<Q ; ++q){
    for(unsigned k=0 ; k<K(q) ; ++k) {
      vec m_temp = m[q] ;
      vec l_temp = l[q] ;
      NumericVector potential_intervals_temp = potential_intervals[q];
      vec potential_intervals_dims_temp = potential_intervals_dims[q];

      x_tilde.col(k+1+count) = potential_intervals_extract(potential_intervals_temp,
                  m_temp(k),l_temp(k),potential_intervals_dims_temp);
    }
    count = count + K(q);
  }
}


//##############################################################################
//#################### Gibbs Sampler and Simulated Annealing ###################
//##############################################################################
//

// Perform the Gibbs Sampler algorithm for the Bliss model
// returned values; trace and param.
// trace : a matrix, with the different parameters in columns and the
// iterations in rows.
// The different parameters are : b_1, m1, l1, b_2, m2,
//l2, ..., b_Q, mQ, lQ, mu, sigma2.
// Hence the trace has (iter + 1) rows (+1 for the initailisation), and
// 3*(K1+K2+...KQ)+2 columns. Indeed,
// b_1, m1 and l1 are of length K, ..., b_Q, mQ and lQ are of
// length KQ.
// [[Rcpp::export]]
List Bliss_Gibbs_Sampler_cpp (int Q, arma::vec & y, List & x, List & grids,
                              int iter, arma::vec & K, CharacterVector & basis, 
                              double g, double lambda ,arma::mat & V_tilde, 
                              List & l_values, arma::vec & l_values_length,List & probs_l, 
                              bool progress, double tol) {
  if(progress) Rcpp::Rcout << "Gibbs Sampler: " <<  std::endl;
  if(progress) Rcpp::Rcout << "\t Initialization." <<  std::endl;

  // Compute the value of n and the p's
  int n = as<mat>(x[0]).n_rows ;

  vec p = zeros<vec>(Q)        ;
  for(int i=0 ; i<Q ; ++i){
    p(i) = as<mat>(x[i]).n_cols;
  }

  // Compute projection of the x_i on all the intervals
  List normalization_values(Q) ;    // normalization_values is used to normalize the predictors
  List potential_intervals(Q)        ;    // will be contain all the projections
  List potential_intervals_dims(Q)   ;    // will be contain the dim of the potential_intervals's
  for( int q=0 ; q<Q ; ++q){
    arma::cube temp = potential_intervals_List (x, grids, l_values_length, basis,q) ;
    normalization_values[q] = temp.slice(n);

    temp = temp.subcube(0,0,0,p(q)-1,l_values_length(q)-1,n-1);
    potential_intervals[q] = temp;

    vec temp2 = zeros<vec>(3);
    temp2(0) = p(q)    ;
    temp2(1) = l_values_length(q) ;
    temp2(2) = n    ;
    potential_intervals_dims[q] = temp2;
  }

  // Compute the matrix of lambda for the Ridge penalty of the Rigde
  // Zellner prior...
  int sum_K = sum(K);
  mat lambda_id0  = zeros<mat>(sum_K+1,sum_K+1) ;
  lambda_id0(0,0) = 100*var(y);  // XXXXXXXXXXx
  for( unsigned i=1 ; i<sum_K+1; ++i){
   lambda_id0(i,i) = lambda ;
  }
  // ... or the constant matrix if V does not depend on the intervals

  // Determine the start point
  if(progress) Rcpp::Rcout << "\t Determine the starting point." <<  std::endl;
  double sigma_sq       ;
  vec b_tilde           ;
  mat Sigma_b_tilde_inv ;
  mat W_inv             ; 

  bool success = false ;
  mat R                ;
  mat test             ;

  List m(Q) ;
  List l(Q) ;
  mat x_tilde = ones<mat>(n,sum_K+1) ;
  
  // Try to determine a starting point which not leads to a non-invertible
  // matrix problem
  while(success == false){
    // Initialization of sigma_sq
    sigma_sq  = var(y) ;

    // Initialization of the middle and length of the intervals
    for( unsigned q=0 ; q<Q ; ++q){
      vec probs_l_temp = probs_l[q];
      m[q] = sample_cpp(K(q),p(q)) + 1 ;
      l[q] = sample_weight(K(q),probs_l_temp) + 1 ;
    }

    // Initialize the current x_tilde matrix (which depend on the intervals)
    int count = 0;
    for( unsigned q=0 ; q<Q ; ++q){
      for(unsigned k=0 ; k<K(q) ; ++k) {
        vec m_temp = m[q] ;
        vec l_temp = l[q] ;
        NumericVector potential_intervals_temp = potential_intervals[q];
        vec potential_intervals_dims_temp = potential_intervals_dims[q];
        x_tilde.col(k+1+count) = potential_intervals_extract(potential_intervals_temp,
                    m_temp(k),l_temp(k),potential_intervals_dims_temp);
      }
      count = count + K(q);
    }

    // Initialize the current W_inv matrix (which depend on the intervals)
    // (W_inv is the covariance matrix of the Ridge Zellner prior) (without sig)
    W_inv = compute_W_inv_List (Q,K, g, x_tilde,sum_K,lambda_id0) ;

    // Check if there is a non-invertible matrix problem
    Sigma_b_tilde_inv = W_inv + trans(x_tilde) * x_tilde ;
    test            = ginv_cpp(Sigma_b_tilde_inv,tol)    ;
    success         = accu(abs(test)) != 0                  ;
  }
  
  // Initialization of b_tilde 
  // peut etre definir le vecteur de 0 avant XXXXXXXXXXXXXXXXXXXXX
  b_tilde = mvrnormArma( zeros<vec>(sum_K+1) , ginv_cpp(W_inv,tol) , sigma_sq) ;
  
  // Initialize the matrix trace
  mat trace = zeros<mat>(iter+1,3*sum_K+2);
  int count = 0;
  for( unsigned q=0 ; q<Q ; ++q){
    vec m_temp = m[q] ;
    vec l_temp = l[q] ;
    trace.row(0).subvec( 3*count        , 3*count+  K(q)-1) =
      trans(b_tilde.subvec( 1+count , K(q)+count )) ;
    trace.row(0).subvec( 3*count+  K(q) , 3*count+2*K(q)-1)   = trans(m_temp) ;
    trace.row(0).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1)   = trans(l_temp) ;

    trace(0,3*sum_K  ) = b_tilde(0) ;
    trace(0,3*sum_K+1) = sigma_sq      ;
    count = count + K(q) ;
  }

  // Initialize some variable used in the Gibbs loop
  List l_alternative(Q);
  for( int q=0 ; q<Q ; ++q){
    l_alternative[q] = sequence(1,l_values_length(q),1);
  }
  List m_alternative(Q);
  for( int q=0 ; q<Q ; ++q){
    m_alternative[q] = sequence(1,p(q),1);
  }

  // The Gibbs loop
  if(progress) Rcpp::Rcout << "\t Start the Gibbs Sampler." <<  std::endl;
  for(unsigned i=1  ; i < iter+1 ; ++i ) {
    if( i % (iter / 10)  == 0)
      if(progress) Rcpp::Rcout << "\t " << i / (iter / 100) << "%" << std::endl;

    // update sigma_sq
    update_sigma_sq(y,b_tilde,W_inv,x_tilde,n,sum_K,sigma_sq) ;

    // update m
    count = 0 ;
    // count is used to browse some vec/mat when p(q) is not constant wrt q.
    for( unsigned q=0 ; q<Q ; ++q ){
      // Compute some quantities which do not vary with k
      vec m_q = m[q];
      vec l_q = l[q];
      int p_q = p(q);
      NumericVector potential_intervals_q = potential_intervals[q];
      vec potential_intervals_dims_q      = potential_intervals_dims[q];
      vec m_alternative_q = sequence(1,p_q,1) ;

      for(int k=0 ; k<K(q) ; ++k){
        // update m_k
        update_mqk(count,k,y,b_tilde,sigma_sq,m_q,l_q,x_tilde,
            potential_intervals_q,potential_intervals_dims_q,m_alternative_q,p_q,Q,K,g,
            sum_K,lambda_id0);
       
        // update the value "x_tilde"
        update_x_tilde(Q,K,potential_intervals,potential_intervals_dims,m,l,
                       x_tilde);
      }

      // Update the m_q value
      m[q] = m_q;
      // Update count
      count = count + K(q);
    }

    // update l
    count = 0 ;
    // count is used to browse some vec/mat when p(q) is not constant wrt q.
    for( unsigned q=0 ; q<Q ; ++q ){
      // Compute some quantities which do not vary with k
      vec m_q = m[q];
      vec l_q = l[q];
      int l_values_length_q = l_values_length(q);
      NumericVector potential_intervals_q = potential_intervals[q];
      vec potential_intervals_dims_q      = potential_intervals_dims[q];
      vec l_alternative_q = sequence(1,l_values_length_q,1) ;
      vec probs_l_q         = probs_l[q];

      for(int k=0 ; k<K(q) ; ++k){
        // update l_k
        update_lqk(count,k,y,b_tilde,sigma_sq,m_q,l_q,x_tilde,
            potential_intervals_q,potential_intervals_dims_q,l_alternative_q,probs_l_q,
            l_values_length_q, Q,K,g,sum_K,lambda_id0);
       

        // update the value "x_tilde"
        update_x_tilde(Q,K,potential_intervals,potential_intervals_dims,m,l,
                                 x_tilde);
      }

      // Update the m_q value
      l[q] = l_q;
      // Update count
      count = count + K(q);
    }

    // update the value "W_inv" (only after the updating of the m's and l's)
    W_inv = compute_W_inv_List (Q,K, g, x_tilde,sum_K,lambda_id0) ;

    // update the matrix Sigma_b_tilde (only after the updating of
    // the m's and l's)
    Sigma_b_tilde_inv = W_inv + trans(x_tilde) * x_tilde   ;
    test                 = ginv_cpp(Sigma_b_tilde_inv,tol) ;
    success              = accu(abs(test)) != 0               ;

    // Try to determine an update which not leads to a non-invertible
    // matrix problem. If there is a problem, go back to the beginning of the
    // updating process.
    if(success){
      // update the b_tilde
      update_b_tilde(y,sigma_sq,x_tilde,Sigma_b_tilde_inv,tol,b_tilde) ;

      // update the matrix trace
      count = 0;
      for( unsigned q=0 ; q<Q ; ++q){
        vec m_temp = m[q] ;
        vec l_temp = l[q] ;
        trace.row(i).subvec( 3*count        , 3*count+  K(q)-1) =
          trans(b_tilde.subvec( 1+count , K(q)+count )) ;
        trace.row(i).subvec( 3*count+  K(q) , 3*count+2*K(q)-1) = trans(m_temp);
        trace.row(i).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1) = trans(l_temp);

        trace(i,3*sum_K  ) = b_tilde(0) ;
        trace(i,3*sum_K+1) = sigma_sq      ;
        count = count + K(q) ;
      }
    }else{ //... go back to the beginning of the updating process.
      i     = i - 1 ;
      count = 0;
      for( unsigned q=0 ; q<Q ; ++q){
        b_tilde.subvec( 1+count , K(q)+count ) =
          trans(trace.row(i).subvec( 3*count , 3*count+  K(q)-1))  ;
        m[q] = trans(trace.row(i).subvec( 3*count+  K(q) , 3*count+2*K(q)-1)) ;
        l[q] = trans(trace.row(i).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1)) ;

        b_tilde(0) = trace(i,3*sum_K  ) ;
        sigma_sq      = trace(i,3*sum_K+1) ;
        count = count + K(q) ;
      }

      // update the value "x_tilde"
      update_x_tilde(Q,K,potential_intervals,potential_intervals_dims,m,l,
                     x_tilde);

      // update the value "W_inv"
      W_inv = compute_W_inv_List (Q,K, g, x_tilde,sum_K,lambda_id0) ;
  }
  }

  // return the trace and the parameters
  if(progress) Rcpp::Rcout << "\t Return the result." <<  std::endl;
  return  List::create(_["trace"]=trace,
                       _["param"]=List::create(_["phi_l"]=probs_l,
                                          _["K"]=K,
                                          _["l_values_length"]=l_values_length,
                                          _["potential_intervals"]=potential_intervals,
                                          _["grids"]=grids,
                                          _["normalization_values"]=normalization_values
                       ));
  }

// Perform the Simulated Annealing algorithm to minimize the loss function
// [[Rcpp::export]]
List Bliss_Simulated_Annealing_cpp (int iter, arma::mat & beta_sample, arma::vec & grid,
                                    int burnin, double Temp,int k_max,
                                    int l_max, int dm, int dl,
                                    int p,std::string basis, arma::mat & normalization_values,
                                    bool progress){
  if(progress) Rcpp::Rcout << "Simulated Annealing:" <<  std::endl;
  // Initialization
  if(progress) Rcpp::Rcout << "\t Initialization." <<  std::endl;
  int N = beta_sample.n_rows;
  vec posterior_expe = zeros<vec>(p);
  vec posterior_var  = zeros<vec>(p); // 
  for(int i=0 ; i<p ; ++i){
    posterior_expe(i) = mean(beta_sample.col(i));
    posterior_var(i)  =  var(beta_sample.col(i));
  }

  vec probs;
  int k;
  vec m;
  vec l;
  vec b;
  vec d;
  double d_loss;
  int boundary_min;
  int boundary_max;
  vec difference;

  vec d_tmp;
  double d_loss_tmp;
  double proba_acceptance;
  double u;
  int j;
  double Temperature;
  vec b_tmp;
  vec m_tmp;
  vec l_tmp;
  int k_tmp;
  int accepted;
  vec choice_prob_interval;
  vec choice_prob;
  int choice;
  int interval_min;
  int interval_max;
  double var_b;
  vec boundaries_min;
  vec boundaries_max;
  vec boundaries;
  int new_m;
  int new_l;
  double new_b;

  // Initialize the matrix trace
  mat trace = zeros<mat>(iter+1,3*k_max+3);

  // Determine the start point
  if(progress) Rcpp::Rcout << "\t Determine the starting point." <<  std::endl;
  probs = ones<vec>(k_max);
  k     = sample_weight( probs )+1;
  m     = zeros<vec>(k);
  l     = zeros<vec>(k);
  b     = zeros<vec>(k);

  probs = ones<vec>(p);
  m(0)  = sample_weight( probs )+1;
  probs = ones<vec>(l_max); // besoin de lmax: a remplacer par du calcul de l_values
  l(0)  = sample_weight( probs )+1;

  boundary_min = m(0)-l(0)-1;
  boundary_max = m(0)+l(0)-1;
  if(boundary_min < 0   ) boundary_min = 0   ;
  if(boundary_max > p-1 ) boundary_max = p-1 ;
  b(0) = mean(posterior_expe.subvec( boundary_min , boundary_max ));
  d = compute_beta_cpp(b,m,l,grid,p,1,basis,normalization_values);

  if(k > 1){
    for(int i=1 ; i<k ; ++i ){
      // Compute the difference ...
      difference = abs(posterior_expe - d);
      // ... and its smoothed version.
      difference = moving_average_cpp(difference,floor(p/10));

      // Which intervals are possible ?
      for(int o=0 ; o<i ; ++o){
        if( m(o) - l(o) -1 > 0  ) boundary_min = m(o) - l(o) -1; else boundary_min = 1;
        if( m(o) + l(o) +1 < p+1) boundary_max = m(o) + l(o) +1; else boundary_max = p;
        for(int z=boundary_min; z < boundary_max+1 ; ++z){
          difference(z-1) = 0;
        }
      }

      // Simulate an interval
      if(sum(difference) > 0){
        vec boundaries;
        vec boundaries_min;
        vec boundaries_max;
        int boundary_min;
        int boundary_max;

        // Simulate a m
        m(i) = sample_weight(difference) +1;

        // Simulate a l
        boundaries_max = zeros<vec>(i);
        boundaries_min = zeros<vec>(i);

        boundaries_min = abs(m.subvec(0,i-1) - l.subvec(0,i-1) - m(i))-1 ;
        boundaries_max = abs(m.subvec(0,i-1) + l.subvec(0,i-1) - m(i))-1 ;

        boundaries = zeros<vec>(2*i+1);
        boundaries(0) = l_max;
        boundaries.subvec(1  ,i  ) = boundaries_min;
        boundaries.subvec(1+i,i+i) = boundaries_max;
        boundaries = sort(boundaries);
        boundary_max = boundaries(0);

        if(boundary_max < 1) boundary_max = 1;

        l(i) = sample_weight( ones<vec>(boundary_max) ) + 1 ;
        // Simulate a b (from the smoothed difference)
        if( m(i) - l(i) -1 > 0  ) boundary_min = m(i) - l(i) -1; else boundary_min = 1;
        if( m(i) + l(i) +1 < p+1) boundary_max = m(i) + l(i) +1; else boundary_max = p;
        b(i) = mean( difference.subvec(boundary_min-1 , boundary_max-1) );
        // Compute the function with these intervals
        d = compute_beta_cpp(b,m,l,grid,p,i+1,basis,normalization_values);
      }else{
        // sortir de la boucle avec une bonne valeur de k
        unsigned i_tmp = i;
        i = k ;
        k = i_tmp;
      }
    }
  }

  // Compute the first function with K intervals (and its loss)
  d      = compute_beta_cpp(b,m,l,grid,p,k,basis,normalization_values);
  d_loss = loss_cpp(d,grid,posterior_expe);

  // Update the trace with the start point
  trace.row(0).subvec( 0       ,           k-1) = trans(b.subvec(0,k-1)) ;
  trace.row(0).subvec( k_max   , k_max   + k-1) = trans(m.subvec(0,k-1))         ;
  trace.row(0).subvec( k_max*2 , k_max*2 + k-1) = trans(l.subvec(0,k-1))         ;
  trace(0,3*k_max)   = 1      ;
  trace(0,3*k_max+1) = k      ;
  trace(0,3*k_max+2) = d_loss ;

  // Start the loop
  if(progress) Rcpp::Rcout << "\t Start the loop." <<  std::endl;
  for(int i=0 ; i<iter ; ++i){
    Temperature = cooling_cpp(i,Temp);
    // Progress
    if( (i+1) % (iter / 10)  == 0)
      if(progress) Rcpp::Rcout << "\t " << (i+1) / (iter / 100) << "%" << std::endl;
    // Initialize the proposal
    b_tmp     = b ;
    m_tmp     = m ;
    l_tmp     = l ;
    k_tmp     = k ;
    accepted  = 0 ;
    choice_prob_interval = ones<vec>(k);

    // Choose a move
    choice_prob = ones<vec>(5);
    if(k == k_max) choice_prob(3) = choice_prob(3) -1;
    if(k == 1    ) choice_prob(4) = choice_prob(4) -1;

    choice = sample_weight( choice_prob ) + 1;
    // change a b
    if(choice == 1){
      // choose an interval
      j = sample_weight(choice_prob_interval);

      // Simulate a new b_k
      interval_min = m(j)-l(j) -1;
      interval_max = m(j)+l(j) -1;
      if(interval_min < 0   ) interval_min = 0   ;
      if(interval_max > p-1 ) interval_max = p-1 ;

      var_b = mean(posterior_var.subvec( interval_min , interval_max));
      b_tmp(j) = R::rnorm( b(j), sqrt(var_b) );
    }
    // change a m_k
    if(choice == 2){
      // choose an interval
      j = sample_weight(choice_prob_interval);

      // Simulate a new m_k
      if(k > 1){
        boundaries_max = zeros<vec>(k);
        boundaries_max.subvec(0,k-2) = vec_drop_k(m,j).subvec(0,k-2) -
          vec_drop_k(l,j).subvec(0,k-2) - l(j)-1 ;
        boundaries_max(k-1)          = p ;
        boundaries_max               = sort(boundaries_max);
        boundary_max                 = boundaries_max(boundaries_max.size()-1);

        for(int o=0 ; o<boundaries_max.size() ; ++o){
          if(m(j) < boundaries_max(o) ){
            boundary_max = boundaries_max(o) ;
            break;
          }
        }

        boundaries_min = zeros<vec>(k);
        boundaries_min.subvec(0,k-2) = vec_drop_k(m,j).subvec(0,k-2) +
          vec_drop_k(l,j).subvec(0,k-2) + l(j)+1 ;
        boundaries_min(k-1)          = 1 ;
        boundaries_min               = sort(boundaries_min);
        boundary_min                 = boundaries_min(0);

        for(int o=boundaries_max.size()-1 ; o>=0 ; --o){
          if(m(j) > boundaries_min(o) ){
            boundary_min = boundaries_min(o) ;
            break;
          }
        }
      }else{
        boundary_max = m(j) + dm ;
        boundary_min = m(j) - dm ;
        if(boundary_max > p) boundary_max = p ;
        if(boundary_min < 1) boundary_min = 1 ;
      }

      if(boundary_max - boundary_min + 1 > 0){
        probs = ones<vec>( boundary_max - boundary_min + 1 ) ;
        m_tmp(j)  = sample_weight( probs ) + boundary_min ;
      }
    }
    // change a l_k
    if(choice == 3){
      // choose an interval
      j = sample_weight(choice_prob_interval);

      // Simulate a new l_k
      if(k > 1){
        boundaries_max = zeros<vec>(k-1);
        boundaries_min = zeros<vec>(k-1);
        boundaries_max = abs(vec_drop_k(m,j) + vec_drop_k(l,j) - m(j))-1 ;
        boundaries_min = abs(vec_drop_k(m,j) - vec_drop_k(l,j) - m(j))-1 ;

        boundaries = zeros<vec>(2*k-1);
        boundaries(0) = dl;
        boundaries.subvec(1,k-1)   = boundaries_min.subvec(0,k-2);
        boundaries.subvec(k,2*(k-1)) = boundaries_max.subvec(0,k-2);
        boundaries = sort(boundaries);
        boundary_max = boundaries(0);

        if(boundary_max > 1){
          l_tmp(j) = sample_weight( ones<vec>(boundary_max) ) + 1 ;
        }
      }
    }
    // birth
    if(choice == 4){
      // compute the difference ...
      difference = posterior_expe - d;
      // ... and its smoothed version
      difference = moving_average_cpp(difference,floor(p/10));

      // Which intervals are possible ?
      for(int o=0 ; o<k ; ++o){
        if( m(o) - l(o) -1 > 0  ) boundary_min = m(o) - l(o) -1; else
          boundary_min = 1;
        if( m(o) + l(o) +1 < p+1) boundary_max = m(o) + l(o) +1; else
          boundary_max = p;
        for(int z=boundary_min; z < boundary_max+1 ; ++z){
          difference(z-1) = 0;
        }
      }
      if(sum(abs(difference)) > 0){
        // update k
        k_tmp = k+1;
        // Simulate a new m
        new_m = sample_weight(abs(difference)) +1;
        m_tmp = zeros<vec>(k_tmp);
        m_tmp.subvec(0,k-1) = m;
        m_tmp(k_tmp-1)      = new_m;
        // Simulate a new l
        boundaries_max = zeros<vec>(k_tmp-1);
        boundaries_min = zeros<vec>(k_tmp-1);

        boundaries_min = abs(m - l - new_m)-1 ;
        boundaries_max = abs(m + l - new_m)+1 ;

        boundaries = zeros<vec>(2*k+1);
        boundaries(0) = l_max;
        boundaries.subvec(1  ,k  ) = boundaries_min;
        boundaries.subvec(1+k,k+k) = boundaries_max;
        boundaries = sort(boundaries);
        boundary_max = boundaries(0);

        new_l = sample_weight( ones<vec>(boundary_max) ) + 1 ;
        l_tmp = zeros<vec>(k_tmp);
        l_tmp.subvec(0,k-1) = l;
        l_tmp(k_tmp-1)      = new_l;

        // Simulate a new b (from the smoothed difference)
        if( new_m - new_l -1 > 0  ) boundary_min = new_m - new_l -1; else boundary_min = 1;
        if( new_m + new_l +1 < p+1) boundary_max = new_m + new_l +1; else boundary_max = p;
        new_b = mean( difference.subvec(boundary_min-1 , boundary_max-1) );
        b_tmp = zeros<vec>(k_tmp);
        b_tmp.subvec(0,k_tmp-2) = b;
        b_tmp(k_tmp-1)          = new_b;
      }
    }
    // death
    if(choice == 5){
      // Choose an interval to drop
      j = sample_weight(choice_prob_interval);

      // Drop the interval
      k_tmp = k-1;
      b_tmp = zeros<vec>(k_tmp);
      m_tmp    = zeros<vec>(k_tmp);
      l_tmp    = zeros<vec>(k_tmp);

      b_tmp = vec_drop_k(b,j);
      m_tmp    = vec_drop_k(m   ,j);
      l_tmp    = vec_drop_k(l   ,j);
    }

    // Compute the acceptance probability
    d_tmp      = compute_beta_cpp(b_tmp,m_tmp,l_tmp,grid,p,k_tmp,basis,
                                normalization_values);
    d_loss_tmp = loss_cpp(d_tmp,grid,posterior_expe);

    proba_acceptance = exp( -( d_loss_tmp-d_loss )/ Temperature );

    // Accept/reject
    u = R::runif(0,1) ;
    if(u < proba_acceptance){
      b = zeros<vec>(k_tmp)     ;
      l = zeros<vec>(k_tmp)     ;
      m = zeros<vec>(k_tmp)     ;
      accepted  = 1             ;
      b = b_tmp.subvec(0,k_tmp-1) ;
      m         = m_tmp.subvec(0,k_tmp-1)         ;
      l         = l_tmp.subvec(0,k_tmp-1)         ;
      k         = k_tmp         ;
      d         = d_tmp         ;
      d_loss    = d_loss_tmp    ;
    }

    // Update the trace
    trace.row(i+1).subvec( 0       ,           k-1) = trans(b.subvec( 0,k-1)) ;
    trace.row(i+1).subvec( k_max   , k_max   + k-1) = trans(m.subvec( 0,k-1))         ;
    trace.row(i+1).subvec( k_max*2 , k_max*2 + k-1) = trans(l.subvec( 0,k-1))         ;
    trace(i+1,3*k_max  ) = accepted ;
    trace(i+1,3*k_max+1) = k        ;
    trace(i+1,3*k_max+2) = d_loss   ;
  }

  // Return the result
  if(progress) Rcpp::Rcout << "\t Return the result." <<  std::endl;
  return  List::create(_["trace"]         =trace,
                       _["posterior_expe"]=posterior_expe);
}





// Perform the Simulated Annealing algorithm to minimize the loss function
// [[Rcpp::export]]
arma::mat dposterior_cpp (arma::mat & rposterior, arma::vec & y, unsigned N, unsigned K,
                    NumericVector & potential_intervals, arma::vec & potential_intervals_dims,
                    double lambda, double l_max){
  mat      res = zeros<mat>(N,6);
  unsigned n   = y.size()       ;
  double   RSS = 0              ;

  vec    b_tilde = zeros<vec>(K+1) ;
  double sigma_sq ;
  vec    m        ;
  vec    l        ;

  double posterior_d      ;
  double log_posterior_d  ;
  double prior_d      ;
  double log_prior_d  ;
  double lkh     ;
  double log_lkh ;


  mat    W_inv  ;
  double a_l = 5*K ;

  mat    x_tilde  = ones<mat>(n,K+1) ;
  mat lambda_id0  = zeros<mat>(K+1,K+1) ;
  lambda_id0(0,0) = 100*var(y);
  for( unsigned i=1 ; i<K+1; ++i){
    lambda_id0(i,i) = lambda ;
  }
  vec K_vec = ones<vec>(1) ;
  K_vec(0)  = K ;

  for(unsigned j=0 ; j<N ; ++j ){
    // load the parameter value
    b_tilde(0)          = rposterior(j,3*K  ) ;
    b_tilde.subvec(1,K) = trans(rposterior.row(j).subvec(  0,  K-1)) ;
    sigma_sq = rposterior(j,3*K+1) ;
    m        = trans(rposterior.row(j).subvec(  K,2*K-1)) ;
    l        = trans(rposterior.row(j).subvec(2*K,3*K-1)) ;

    // compute x_tilde
    for(unsigned k=0 ; k<K ; ++k) {
      x_tilde.col(k+1) = potential_intervals_extract(potential_intervals,m(k),l(k),
                  potential_intervals_dims);
    }
    // compute Sigma
    W_inv = compute_W_inv_List (1,K_vec,n, x_tilde,K,lambda_id0) ;

    RSS   = sum(square(y - x_tilde * b_tilde)) ;
    // compute the (log) likelihood
    log_lkh = -1./(2.*sigma_sq) * RSS - n/2. * log(sigma_sq);
    lkh     = exp(log_lkh);

    // compute the (log) prior density
    log_prior_d = -1./(2.*sigma_sq) * dot(b_tilde, W_inv * b_tilde) -
      (K+3.)/2. * log(sigma_sq) - a_l * sum(l) / l_max + log(det(W_inv)) / 2.;
    prior_d = exp(log_prior_d);
    // compute the (log) posterior density

    log_posterior_d = log_lkh + log_prior_d;
    posterior_d     = exp(log_posterior_d);

    // update the result
    res(j,0) = posterior_d     ;
    res(j,1) = log_posterior_d ;
    res(j,2) = lkh     ;
    res(j,3) = log_lkh ;
    res(j,4) = prior_d     ;
    res(j,5) = log_prior_d ;
  }

  // return the (log) densities
  return(res);
}
