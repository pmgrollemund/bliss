//#########################################################
//#                                                       #
//#            Bliss method : rcpp code                   #
//#                                                       #
//#########################################################
#define ARMA_DONT_PRINT_ERRORS
// #include <Rcpp.h>

#define ARMA_WARN_LEVEL 0
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <string>
#include <iostream>
#include <vector>
#include <cstring>
// #include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// //[[Rcpp::depends(RcppClock)]]
// #include <RcppClock.h>
// #include <thread>

using namespace Rcpp;
using namespace arma;

//#############################################################################
//############################ basic functions  ###############################
//#############################################################################

// The R function : ginv (generalized matrix inversion using SVD decomposition)
// [[Rcpp::export]]
arma::mat ginv_cpp (arma::mat & x, double tol){
  unsigned p;
  arma::mat u;
  arma::vec s;
  arma::mat v;

  svd(u,s,v,x);
  p = s.size();

  // tol   = tol * s(0);
  arma::mat S = arma::zeros<arma::mat>(p,p);

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
  arma::vec res = arma::zeros<arma::vec>(nbre);
  for (int i = 0 ; i < nbre ; ++i) {
    res(i) = sample_cpp(n);
  }
  return res;
}

// Weighted sample in 0:n-1.
int sample_weight (arma::vec proba){
  if(sum(proba) == 0)   proba = ones<arma::vec>(proba.n_rows) ;
  arma::vec proba_cum = cumsum(proba)/sum(proba) ;

  unsigned ret = 0;
  double u = R::runif(0,1);
  while (ret <= proba_cum.n_rows and u > proba_cum(ret)) {
    ret++;
  }
  return ret;
}

// Vector weighted sample in 0:n-1.
arma::vec sample_weight (int n, arma::vec proba){
  arma::vec res = arma::zeros<arma::vec>(n);
  for (int i = 0 ; i < n ; ++i) {
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

// Function seq
arma::vec sequence(int a,int b,double by){
  int range = floor((b-a)/by + 1) ;
  arma::vec res = arma::zeros<arma::vec>(range);
  for(int i=0 ; i<range ; i++){
    res(i) = a + i*by;
  }
  return res;
}

// Extract an element from a cube
double cube_extract(NumericVector & cube, int x , int y, int z, arma::vec & dims){
  double res;
  res = cube[x + y*dims(0) + z*dims(1)*dims(0)];
  return res;
}

// Simulate from a multidimensional gaussian.
arma::vec mvrnormArma(arma::vec mu, arma::mat VarCovar, double sigma_sq, double tol) {
  int ncols = VarCovar.n_cols;
  arma::vec y = arma::randn<arma::vec>(ncols);

  if(VarCovar(0,0) < tol) VarCovar(0,0) = 1.1*tol;
  VarCovar = chol(VarCovar);

  return  mu + std::sqrt(sigma_sq) * trans(trans(y) * VarCovar);
}

// Compute a trapezoidal approximation of area under curve.
// [[Rcpp::export]]
double integrate_trapeze_cpp (arma::vec & x, arma::vec & y){
  arma::vec diff_x = vec_drop_k(x,0) - vec_drop_k(x,x.size()-1);
  arma::vec cumu_y = vec_drop_k(y,0) + vec_drop_k(y,y.size()-1);
  return sum( diff_x % cumu_y  )/2 ;
}

//  Compute the norm of a vector.
double L2_norm(arma::vec & x,arma::vec & y){
  arma::vec tmp = arma::zeros<arma::vec>(x.size());
  for(unsigned i=0 ; i<tmp.size() ; ++i){
    tmp(i) = std::pow( y(i),2 );
  }
  double res;
  res = sqrt(integrate_trapeze_cpp(x,tmp));

  return res;
}

// Use to compute an uniform function, see function compute_beta.
// [[Rcpp::export]]
arma::vec uniform_cpp (int m, int l, arma::vec & grid){
  int p = grid.size();
  arma::vec res = arma::zeros<arma::vec>(p);
  arma::vec index = sequence(m-l,m+l,1);
  int tmp;
  for(unsigned i=0 ; i<index.size() ; ++i){
    tmp = index(i);
    if( (tmp <= p) && (tmp >= 1) ){
      res(index(i)-1) = 1;
    }
  }
  double res_norm = L2_norm(grid,res);
  res = res / res_norm ;
  return res;
}

// Use to compute a triangular function, see function compute_beta.
// [[Rcpp::export]]
arma::vec triangular_cpp (int m, int l, arma::vec & grid){
  int p = grid.size();
  arma::vec res = arma::zeros<arma::vec>(p);
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
  double res_norm = L2_norm(grid,res);
  res = res / res_norm ;
  return res;
}

// Use to compute a gaussian function, see function compute_beta.
// [[Rcpp::export]]
arma::vec gaussian_cpp (int m, int l, arma::vec & grid){
  int p = grid.size();
  arma::vec res = arma::zeros<arma::vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  exp( - 9*std::pow(i/l_double,2)/2) ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
  double res_norm = L2_norm(grid,res);
  res = res / res_norm ;
  return res;
}

// Use to compute a gaussian function, see function compute_beta.
// [[Rcpp::export]]
arma::vec Epanechnikov_cpp (int m, int l, arma::vec & grid){
  int p = grid.size();
  arma::vec res = arma::zeros<arma::vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  1-std::pow(i/l_double,2) ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
  double res_norm = L2_norm(grid,res);
  res = res / res_norm ;
  return res;
}

// compute_beta in cpp.
// [[Rcpp::export]]
arma::vec compute_beta_cpp (arma::vec & b, arma::vec & m, arma::vec & l,
                            arma::vec & grid, int p, int K, std::string basis,
                            arma::mat & normalization_values ){
  arma::vec res = arma::zeros<arma::vec>(p) ;

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
arma::mat compute_beta_sample_cpp (arma::mat & posterior_sample,
                                   int K, arma::vec & grid, int p, std::string & basis,
                                   arma::mat & normalization_values){
  arma::mat res = zeros<mat>(posterior_sample.n_rows,p) ;
  arma::vec b   ;
  arma::vec m   ;
  arma::vec l   ;
  arma::vec tmp ;

  for(unsigned i=0 ; i<res.n_rows ; ++i){
    tmp = trans(posterior_sample.row(i))     ;
    b   = tmp.subvec(0,K-1)    ;
    m   = tmp.subvec(K,2*K-1)   ;
    l   = tmp.subvec(2*K,3*K-1) ;

    res.row(i) = trans(compute_beta_cpp(b,m,l,grid,p,K,basis,normalization_values)) ;
  }
  return res ;
}

// Compute all the alternative for the value of the intergral for all m and l.
// [[Rcpp::export]]
arma::cube potential_intervals_List(List & x_list, List & grids,arma::vec & p_l_vec,
                                    CharacterVector & basis_vec, int q){
  arma::mat x = as<mat>(x_list[q]);
  arma::vec grid = as<arma::vec>(grids[q]);
  int p_l = p_l_vec(q);
  std::string basis = as<std::string>(basis_vec(q));

  int n = x.n_rows ;
  int p = x.n_cols ;
  arma::vec tub;

  arma::cube res(p,p_l,n+1);
  arma::vec tmp;
  arma::vec x_tmp;
  arma::vec tmp2;
  if(basis == "Uniform"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<p_l ; ++j){
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
      for(int j=0 ; j<p_l ; ++j){
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
      for(int j=0 ; j<p_l ; ++j){
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
      for(int j=0 ; j<p_l ; ++j){
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
    for(int j=0 ; j<p_l ; ++j){
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

// Compute a moving average on the vector v.
// [[Rcpp::export]]
arma::vec moving_average_cpp (arma::vec & v, int range){
  int n = v.size();
  arma::vec res = arma::zeros<arma::vec>(n) ;
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

// Compute the largest singular value
double compute_lambda_max (arma::mat G){
  arma::mat u;
  arma::vec s;
  arma::mat v;

  svd(u,s,v,G);

  double lambda_max = max(s);

  return(lambda_max) ;
}

// Compute the matrix V (for a Ridge Zellner prior)
// (for Q functional covaribles)
arma::mat compute_W_inv_List_old (int Q, arma::vec & K, double g, arma::mat & x_tilde,
                                  int sum_K, arma::mat & lambda_id0){
  arma::mat W_inv = arma::zeros<arma::mat>(sum_K+1,sum_K+1);
  arma::mat lambda_id = lambda_id0 ;

  arma::mat x_tilde_temp = mat_drop_col_k(x_tilde,0);
  arma::mat u;
  arma::vec s;
  arma::mat v;
  arma::mat G_0;
  G_0 = trans(x_tilde_temp) * x_tilde_temp; //2024-07-10
  svd(u,s,v,G_0);

  W_inv(0,0) = 1/lambda_id(0,0);
  for( int i=1 ; i<sum_K+1; ++i){
    lambda_id(i,i) = lambda_id(i,i) * max(s); // try with min
  }

  Rcpp::Rcout << s <<  std::endl;
  Rcpp::Rcout << lambda_id <<  std::endl;

  int count = 0;
  for( int q=0 ; q<Q ; ++q){
    W_inv.submat(1+count,1+count,K(q)+count,K(q)+count) =
      ( trans(x_tilde.cols(1+count,K(q)+count)) *
      x_tilde.cols(1+count,K(q)+count) +
      lambda_id.submat(1+count,1+count,K(q)+count,K(q)+count) )  /g;
    count = count + K(q);
  }

  return W_inv;
}
arma::mat compute_Sigma_b_tilde_inv_prior (int Q, arma::vec & K, double g,
                                           arma::mat & x_tilde,int sum_K,
                                           double v0, double v){
  arma::mat Sigma_b_tilde_inv_prior = arma::zeros<arma::mat>(sum_K+1,sum_K+1);

  arma::mat x_tilde_temp = mat_drop_col_k(x_tilde,0);
  arma::mat G = trans(x_tilde_temp) * x_tilde_temp; //2024-07-10
  double lambda = v * compute_lambda_max(G);

  Sigma_b_tilde_inv_prior(0,0) = 1/v0;

  Sigma_b_tilde_inv_prior.submat(1,1,sum_K,sum_K) = G / g ;
  for( int k=1 ; k<sum_K+1 ; ++k){
    Sigma_b_tilde_inv_prior(k,k) = Sigma_b_tilde_inv_prior(k,k) + lambda/g ;
  }

  return Sigma_b_tilde_inv_prior;
}
// Sigma_b_tilde_inv_prior = n^{-1} ( G + v lambda_max(G) I)

// Extract a subvector from the cube potential_intervals with a m_k and a l_k.
// [[Rcpp::export]]
arma::vec potential_intervals_extract (NumericVector & potential_intervals, int mk ,
                                       int lk, arma::vec & dims) {
  arma::vec res = arma::zeros<arma::vec>(dims(2));
  for (int i = 0; i < dims(2); i++) {
    res(i) = cube_extract(potential_intervals, mk - 1, lk - 1, i, dims);
  }
  return res;
}

// Update the parameter m_k
// [[Rcpp::export]]
void update_mqk (int count, int k, arma::vec & y, arma::vec & b_tilde, double sigma_sq,
                 arma::vec & m_q, arma::vec & l_q, arma::mat x_tilde,
                 NumericVector & potential_intervals_q, arma::vec & potential_intervals_dims_q,
                 arma::vec & m_possible_q, int p_q, int Q,
                 arma::vec K, double g, int sum_K, double v0, double v) {
  arma::vec aux = arma::zeros<arma::vec>(p_q);
  arma::vec aux2 = arma::zeros<arma::vec>(p_q);
  arma::mat Sigma_b_tilde_inv_prior_temp;

  // Compute the probabilities
  arma::vec probs = arma::ones<arma::vec>(p_q);
  arma::vec x_tilde_qki = arma::zeros<arma::vec>(potential_intervals_dims_q(2)) ;
  for(int  i=0 ; i<p_q ; ++i){
    x_tilde_qki = potential_intervals_extract(potential_intervals_q,m_possible_q(i),
                                              l_q(k),potential_intervals_dims_q);

    x_tilde.col(count + k + 1) = x_tilde_qki;
    aux(i) = dot( y - x_tilde * b_tilde ,
        y - x_tilde * b_tilde ) /(2*sigma_sq) ;

    Sigma_b_tilde_inv_prior_temp =
      compute_Sigma_b_tilde_inv_prior(Q,K,g,x_tilde,sum_K,v0,v);

    aux(i) = aux(i) +
      dot( b_tilde , Sigma_b_tilde_inv_prior_temp * b_tilde ) / (2*sigma_sq);
    aux2(i) = sqrt( det(Sigma_b_tilde_inv_prior_temp) );
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
                 arma::vec & l_possible_q, arma::vec & phi_l_q, int l_values_length_q,
                 int Q, arma::vec K, double g, int sum_K, double v0, double v) {
  arma::vec aux = arma::zeros<arma::vec>(l_values_length_q);
  arma::vec aux2 = arma::zeros<arma::vec>(l_values_length_q);
  arma::mat Sigma_b_tilde_inv_prior_temp;

  // Compute the probabilities
  arma::vec probs = arma::ones<arma::vec>(l_values_length_q);
  arma::vec x_tilde_qki = arma::zeros<arma::vec>(potential_intervals_dims_q(2)) ;
  for(int  i=0 ; i<l_values_length_q ; ++i){
    x_tilde_qki = potential_intervals_extract(potential_intervals_q,m_q(k),
                                              l_possible_q(i),
                                              potential_intervals_dims_q);

    x_tilde.col(count + k + 1) = x_tilde_qki;
    aux(i) = dot( y - x_tilde * b_tilde ,
        y - x_tilde * b_tilde ) /(2*sigma_sq) ;

    Sigma_b_tilde_inv_prior_temp =
      compute_Sigma_b_tilde_inv_prior(Q,K,g,x_tilde,sum_K,v0,v);

    aux(i) = aux(i) +
      dot( b_tilde , Sigma_b_tilde_inv_prior_temp * b_tilde ) / (2*sigma_sq);
    aux2(i) = sqrt( det(Sigma_b_tilde_inv_prior_temp) );
  }

  double min_aux = min(aux);
  for(int  i=0 ; i<l_values_length_q ; ++i){
    aux(i)  = aux(i) - min_aux;
    probs(i) = aux2(i) * exp( - aux(i) ) * phi_l_q(i);
  }
  // Simulate a lk
  l_q(k) = sample_weight(probs) + 1 ;
}

// update the parameter sigma_sq
void update_sigma_sq (arma::vec & y, arma::vec & b_tilde, arma::mat & Sigma_b_tilde_inv_prior,
                      arma::mat & x_tilde, double & a_star, double & sigma_sq) {
  arma::vec y_tmp= y - x_tilde * b_tilde ;
  double y_tmp2 = dot(y_tmp,y_tmp) ;
  double b_tilde_tmp = dot(b_tilde, Sigma_b_tilde_inv_prior * b_tilde) ;
  double b_star = 0.5*( y_tmp2 + b_tilde_tmp);

  sigma_sq = 1. / (R::rgamma(a_star, 1/b_star) );
}

// update the parameter sigma_sq
void update_sigma_sq_prior (arma::vec & y, arma::vec & b_tilde, arma::mat & W_inv,
                            arma::mat & x_tilde, int n, int sum_K, double & sigma_sq) {
  double b_tilde_tmp = dot(b_tilde, W_inv * b_tilde) ;

  double a_star = (sum_K+1)/2 ;
  double b_star = 0.5*( b_tilde_tmp);

  sigma_sq = 1. / (R::rgamma(a_star, 1/b_star) );
}


// update the parameter b
// [[Rcpp::export]]
void update_b_tilde (arma::vec & y, double sigma_sq, arma::mat & x_tilde,
                     arma::mat & Sigma_b_tilde_inv, double tol,
                     arma::vec & b_tilde) {
  arma::vec mu_b_tilde = trans(x_tilde) * y;
  b_tilde = mvrnormArma( ginv_cpp(Sigma_b_tilde_inv,tol) * mu_b_tilde ,
                         ginv_cpp(Sigma_b_tilde_inv,tol), sigma_sq,tol);
}

// update the parameter b
void update_b_tilde_prior (arma::vec & y, double sigma_sq, arma::mat & x_tilde,
                     arma::mat & Sigma_b_tilde_inv, double tol,
                     arma::vec & b_tilde) {
  arma::vec mu_b_tilde = arma::zeros<arma::vec>(Sigma_b_tilde_inv.n_cols);
  b_tilde = mvrnormArma( mu_b_tilde ,
                         ginv_cpp(Sigma_b_tilde_inv,tol), sigma_sq,tol);
}

// Compute the loss function for a proposal d
// [[Rcpp::export]]
double loss_cpp (arma::vec & d, arma::vec & grid, arma::vec & posterior_expe){
  arma::vec tmp  = d-posterior_expe ;

  return std::pow(L2_norm(grid, tmp),2 );
}

// Compute the decrease of the Temperature
double cooling_cpp (int i, double Temp){
  double tmp = 1.0;
  double res = Temp / log( ( i / 10)*10 + exp(tmp));
  return res;
}

// Update the matrix x_tilde with respect to current intervals
void update_x_tilde (int Q, arma::vec & K, List & potential_intervals,
                     List & potential_intervals_dims, List & m, List & l,
                     arma::mat & x_tilde){
  int count = 0;
  for( int q=0 ; q<Q ; ++q){
    for(int k=0 ; k<K(q) ; ++k) {
      arma::vec m_temp = m[q] ;
      arma::vec l_temp = l[q] ;
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
List Bliss_Gibbs_Sampler_cpp (
    int Q, arma::vec & y, List & x, List & grids,
    int iter, arma::vec & K, CharacterVector & basis,
    double g, double v ,
    arma::vec & l_values_length,List & probs_l,
    bool verbose, double tol) {
  ////////////////////////////////////////////////////////////////
  //////////////// Initialization
  ////////////////////////////////////////////////////////////////
  if(verbose) Rcpp::Rcout << "\t Initialization." <<  std::endl;

  //////////////// Initialization - Primary objects
  // Compute the value of n and the p
  int n = as<mat>(x[0]).n_rows ;
  arma::vec p = arma::zeros<arma::vec>(Q) ;
  for(int i=0 ; i<Q ; ++i){
    p(i) = as<mat>(x[i]).n_cols;
  }
  int sum_K = sum(K);

  //////////////// Initialization - Declare types
  double sigma_sq             ;
  arma::vec b_tilde           ;
  arma::mat Sigma_b_tilde_inv       ;
  arma::mat Sigma_b_tilde_inv_prior ;
  double v0 = 100*var(y); // Weakly informative prior
  if(v0 == 0){
    v0 = sqrt(1/tol);
  }
  if(v0 >= 1/tol){
    v0 = sqrt(1/tol);
  }

  List m(Q) ;
  List l(Q) ;
  List l_possible(Q);
  List m_possible(Q);

  List normalization_values(Q)     ;                                            // normalization_values is used to normalize the predictors
  List potential_intervals(Q)      ;                                            // will be contain all the projections
  List potential_intervals_dims(Q) ;                                            // will be contain the dim of the potential_intervals's

  bool success = false ;
  arma::mat R          ;
  arma::mat test       ;

  arma::mat x_tilde = arma::ones<arma::mat>(n,sum_K+1) ;
  mat trace = zeros<mat>(iter+1,3*sum_K+2);

  //////////////// Initialization - Related to prior
  double a_star = (n+sum_K+1)/2 ;

  //////////////// Initialization - Compute projection of the x_i on all the intervals
  for( int q=0 ; q<Q ; ++q){
    arma::cube temp = potential_intervals_List (x, grids, l_values_length, basis,q) ;
    normalization_values[q] = temp.slice(n);

    temp = temp.subcube(0,0,0,p(q)-1,l_values_length(q)-1,n-1);
    potential_intervals[q] = temp;

    arma::vec temp2 = arma::zeros<arma::vec>(3);
    temp2(0) = p(q)    ;
    temp2(1) = l_values_length(q) ;
    temp2(2) = n    ;
    potential_intervals_dims[q] = temp2;
  }

  ////////////////////////////////////////////////////////////////
  //////////////// Determine the starting point
  ////////////////////////////////////////////////////////////////
  if(verbose) Rcpp::Rcout << "\t Determine the starting point." <<  std::endl;

  //////////////// Try to determine a starting point
  // which not leads to a non-invertible matrix problem
  while(success == false){
    // Starting point for sigma_sq
    sigma_sq  = var(y) ;

    // Starting point for the middle and length of the intervals
    for( int q=0 ; q<Q ; ++q){
      arma::vec probs_l_temp = probs_l[q];
      m[q] = sample_cpp(K(q),p(q)) + 1 ;
      l[q] = sample_weight(K(q),probs_l_temp) + 1 ;
    }

    // Initialize the current x_tilde matrix
    // (which depend on the intervals)
    int count = 0;
    for( int q=0 ; q<Q ; ++q){
      for(int k=0 ; k<K(q) ; ++k) {
        arma::vec m_temp = m[q] ;
        arma::vec l_temp = l[q] ;
        NumericVector potential_intervals_temp = potential_intervals[q];
        arma::vec potential_intervals_dims_temp = potential_intervals_dims[q];
        x_tilde.col(k+1+count) = potential_intervals_extract(potential_intervals_temp,
                    m_temp(k),l_temp(k),potential_intervals_dims_temp);
      }
      count = count + K(q);
    }

    // Initialize the current Sigma_b_tilde_inv_prior matrix
    // (which depend on the intervals)
    // (Sigma_b_tilde_inv_prior is the covariance matrix of the Ridge Zellner prior) (without sig)
    Sigma_b_tilde_inv_prior =
      compute_Sigma_b_tilde_inv_prior(Q,K, g, x_tilde,sum_K,v0,v) ;

    // Check if there is a non-invertible matrix problem
    Sigma_b_tilde_inv = Sigma_b_tilde_inv_prior + trans(x_tilde) * x_tilde ;
    test            = ginv_cpp(Sigma_b_tilde_inv,tol)    ;
    success         = accu(abs(test)) != 0               ;
  }

  // Starting point for b_tilde
  b_tilde = mvrnormArma( zeros<vec>(sum_K+1) ,
                         ginv_cpp(Sigma_b_tilde_inv_prior,tol) ,
                         sigma_sq,tol) ;

  //////////////// Initialize the matrix trace
  int count = 0;
  for( int q=0 ; q<Q ; ++q){
    arma::vec m_temp = m[q] ;
    arma::vec l_temp = l[q] ;
    trace.row(0).subvec( 3*count        , 3*count+  K(q)-1) =
      trans(b_tilde.subvec( 1+count , K(q)+count )) ;
    trace.row(0).subvec( 3*count+  K(q) , 3*count+2*K(q)-1)   = trans(m_temp) ;
    trace.row(0).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1)   = trans(l_temp) ;

    trace(0,3*sum_K  ) = b_tilde(0) ;
    trace(0,3*sum_K+1) = sigma_sq      ;
    count = count + K(q) ;
  }

  //////////////// Initialize some variable used in the Gibbs loop
  for( int q=0 ; q<Q ; ++q){
    l_possible[q] = sequence(1,l_values_length(q),1);
  }
  for( int q=0 ; q<Q ; ++q){
    m_possible[q] = sequence(1,p(q),1);
  }

  ////////////////////////////////////////////////////////////////
  //////////////// The Gibbs loop
  ////////////////////////////////////////////////////////////////
  if(verbose){
    Rcpp::Rcout << "\t Start the Gibbs Sampler loop." <<  std::endl;
  }
  Progress progress_bar(iter, verbose);

  for(int i=1  ; i < iter+1 ; ++i ) {
    // Progress bar
    if(verbose){
      progress_bar.increment();
    }

    // update sigma_sq
    update_sigma_sq(y,b_tilde,Sigma_b_tilde_inv_prior,x_tilde,a_star,sigma_sq) ;

    // update m
    count = 0 ;                                                                 // count is used to browse some vec/mat when p(q) is not constant wrt q.
    for( int q=0 ; q<Q ; ++q ){
      // Compute some quantities which do not vary with k
      arma::vec m_q = m[q];
      arma::vec l_q = l[q];
      int p_q = p(q);
      NumericVector potential_intervals_q = potential_intervals[q];
      arma::vec potential_intervals_dims_q      = potential_intervals_dims[q];
      arma::vec m_possible_q = sequence(1,p_q,1) ;

      for(int k=0 ; k<K(q) ; ++k){
        // update m_k
        update_mqk(count,k,y,b_tilde,sigma_sq,m_q,l_q,x_tilde,
                   potential_intervals_q,potential_intervals_dims_q,m_possible_q,p_q,Q,K,g,
                   sum_K,v0,v);

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
    count = 0 ;                                                                 // count is used to browse some vec/mat when p(q) is not constant wrt q.
    for( int q=0 ; q<Q ; ++q ){
      // Compute some quantities which do not vary with k
      arma::vec m_q = m[q];
      arma::vec l_q = l[q];
      int l_values_length_q = l_values_length(q);
      NumericVector potential_intervals_q = potential_intervals[q];
      arma::vec potential_intervals_dims_q      = potential_intervals_dims[q];
      arma::vec l_possible_q = sequence(1,l_values_length_q,1) ;
      arma::vec probs_l_q         = probs_l[q];

      for(int k=0 ; k<K(q) ; ++k){
        // update l_k
        update_lqk(count,k,y,b_tilde,sigma_sq,m_q,l_q,x_tilde,
                   potential_intervals_q,potential_intervals_dims_q,l_possible_q,probs_l_q,
                   l_values_length_q, Q,K,g,sum_K,v0,v);

        // update the value "x_tilde"
        update_x_tilde(Q,K,potential_intervals,potential_intervals_dims,m,l,
                       x_tilde);
      }

      // Update the m_q value
      l[q] = l_q;
      // Update count
      count = count + K(q);
    }
    // update the value "x_tilde"
    update_x_tilde(Q,K,potential_intervals,potential_intervals_dims,m,l,
                   x_tilde);

    // Update the value "Sigma_b_tilde_inv_prior"
    // (only after the updating of the m's and l's)
    Sigma_b_tilde_inv_prior =
      compute_Sigma_b_tilde_inv_prior (Q,K, g, x_tilde,sum_K,v0,v) ;

    // Update the matrix Sigma_b_tilde
    // (only after the updating of the m's and l's)
    Sigma_b_tilde_inv = Sigma_b_tilde_inv_prior + trans(x_tilde) * x_tilde   ;
    test                 = ginv_cpp(Sigma_b_tilde_inv,tol) ;
    success              = accu(abs(test)) != 0               ;

    // Try to determine an update
    // ()which not leads to a non-invertible matrix problem)
    // If there is a problem, go back to the beginning of the updating process.
    if(success){
      // Update the b_tilde
      update_b_tilde(y,sigma_sq,x_tilde,Sigma_b_tilde_inv,tol,b_tilde) ;

      // Update the matrix trace
      count = 0;
      for( int q=0 ; q<Q ; ++q){
        arma::vec m_temp = m[q] ;
        arma::vec l_temp = l[q] ;
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
      for( int q=0 ; q<Q ; ++q){
        b_tilde.subvec( 1+count , K(q)+count ) =
          trans(trace.row(i).subvec( 3*count , 3*count+  K(q)-1))  ;
        m[q] = trans(trace.row(i).subvec( 3*count+  K(q) , 3*count+2*K(q)-1)) ;
        l[q] = trans(trace.row(i).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1)) ;

        b_tilde(0) = trace(i,3*sum_K  ) ;
        sigma_sq      = trace(i,3*sum_K+1) ;
        count = count + K(q) ;
      }

      // Update the value "x_tilde"
      update_x_tilde(Q,K,potential_intervals,potential_intervals_dims,m,l,
                     x_tilde);

      // Update the value "Sigma_b_tilde_inv_prior"
      Sigma_b_tilde_inv_prior =
        compute_Sigma_b_tilde_inv_prior (Q,K, g, x_tilde,sum_K,v0,v) ;
    }
  }

  ////////////////////////////////////////////////////////////////
  //////////////// Return the trace and the parameters
  ////////////////////////////////////////////////////////////////
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
                                    double Temp,int k_max,
                                    int p_l, int dm, int dl, int p,
                                    std::string basis, arma::mat & normalization_values,
                                    bool verbose, arma::mat & starting_point){
  ////////////////////////////////////////////////////////////////
  //////////////// Initialization
  ////////////////////////////////////////////////////////////////
  if(verbose) Rcpp::Rcout << "\t Initialization." <<  std::endl;

  //////////////// Initialization - Primary objects
  // Approximate the posterior expectation
  vec posterior_expe = zeros<vec>(p);
  for(int i=0 ; i<p ; ++i){
    posterior_expe(i) = mean(beta_sample.col(i));
  }

  //////////////// Initialization - Declare types
  int j;

  arma::vec b;
  arma::vec m;
  arma::vec l;
  int k;
  arma::vec b_tmp;
  arma::vec m_tmp;
  arma::vec l_tmp;
  int k_tmp;

  double mean_b;
  double var_b;
  double new_b;
  int new_m;
  int new_l;

  int value_min;
  int value_max;

  double Temperature;
  arma::vec d;
  arma::vec d_tmp;
  double d_loss;
  double d_loss_tmp;

  int choice;
  arma::vec choice_prob;
  arma::vec choice_prob_interval;
  arma::vec probs;
  double proba_acceptance;
  double u;
  unsigned accepted;

  arma::vec difference;
  arma::vec difference_l;

  unsigned smooth_param = 1;
  if(smooth_param == 0){
    smooth_param = 1;
  }
  arma::vec tmp_basis;

  //////////////// Initialize the matrix trace
  arma::mat trace = arma::zeros<arma::mat>(iter+1,3*k_max+4);
  Temp = std::pow(L2_norm(grid, posterior_expe),2 )*(std::floor(iter)/1000);
  difference = abs(posterior_expe);
  difference = moving_average_cpp(difference,smooth_param);

  ////////////////////////////////////////////////////////////////
  //////////////// Determine the starting point
  ////////////////////////////////////////////////////////////////
  if(verbose) Rcpp::Rcout << "\t Determine the starting point." <<  std::endl;
  m = starting_point.col(0);
  l = starting_point.col(1);
  b = starting_point.col(2);
  k = m.size()   ;
  for(int i=0 ; i<k ; ++i){
    tmp_basis = uniform_cpp(m(i),l(i),grid);
    b(i) = b(i) * normalization_values( m(i)-1 , l(i)-1 ) / tmp_basis(m(i)-1);
  }

  // Compute the first function with K intervals (and its loss)
  d      = compute_beta_cpp(b,m,l,grid,p,k,basis,normalization_values);
  d_loss = loss_cpp(d,grid,posterior_expe);

  // Update the trace with the starting point
  trace.row(0).subvec( 0       ,           k-1) = trans(b.subvec(0,k-1)) ;
  trace.row(0).subvec( k_max   , k_max   + k-1) = trans(m.subvec(0,k-1))         ;
  trace.row(0).subvec( k_max*2 , k_max*2 + k-1) = trans(l.subvec(0,k-1))         ;
  trace(0,3*k_max)   = 1      ;
  trace(0,3*k_max+1) = 0      ;
  trace(0,3*k_max+2) = k      ;
  trace(0,3*k_max+3) = d_loss ;

  ////////////////////////////////////////////////////////////////
  //////////////// The SANN loop
  ////////////////////////////////////////////////////////////////
  if(verbose){
    Rcpp::Rcout << "\t Start the SANN loop." <<  std::endl;
  }
  Progress progress_bar(iter, verbose);

  for(int i=0 ; i<iter ; ++i){
    // Progress bar
    if(verbose){
      progress_bar.increment();
    }

    //Get the current temperature
    Temperature = cooling_cpp(i,Temp);

    // Initialize the proposal
    b_tmp     = b ;
    m_tmp     = m ;
    l_tmp     = l ;
    k_tmp     = k ;
    accepted  = 0 ;
    choice_prob_interval = ones<vec>(k);

    // Choose a move
    choice_prob = ones<vec>(5);
    if(k == k_max){
      choice_prob(3) = choice_prob(3) -1;
    }
    if(k == 1    ){
      choice_prob(4) = choice_prob(4) -1;
    }
    choice = sample_weight( choice_prob ) + 1;

    // change a b_j
    if(choice == 1){
      // choose an interval
      j = sample_weight(choice_prob_interval);

      // Simulate a new b_j
      value_min = m(j)-l(j) -1;
      value_max = m(j)+l(j) -1;
      if(value_min < 0   ){
        value_min = 0   ;
      }
      if(value_max > p-1 ){
        value_max = p-1 ;
      }

      // Compute the difference ...
      difference = (posterior_expe - d);
      difference = moving_average_cpp(difference,smooth_param);
      tmp_basis  = uniform_cpp(m(j),l(j),grid);
      mean_b = mean(difference.subvec( value_min , value_max)*
        normalization_values( m(j)-1 , l(j)-1 )) / tmp_basis(m(j)-1);

      // ... and its smoothed version.
      difference = abs(posterior_expe - d);
      var_b = var(difference.subvec( value_min , value_max)*
        normalization_values( m(j)-1 , l(j)-1 )) / tmp_basis(m(j)-1);

      b_tmp(j) = R::rnorm( mean_b, sqrt(var_b) );
    }
    // change a m_j and b_j
    if(choice == 2){
      // choose an interval
      j = sample_weight(choice_prob_interval);

      // Simulate a new m_j
      difference = (posterior_expe - d);
      difference = moving_average_cpp(difference,smooth_param);
      m_tmp(j) = sample_weight(abs(difference)) +1;

      // Simulate a new b_j
      value_min = m(j)-l(j) -1;
      value_max = m(j)+l(j) -1;
      if(value_min < 0   ){
        value_min = 0   ;
      }
      if(value_max > p-1 ){
        value_max = p-1 ;
      }
  tmp_basis  = uniform_cpp(m(j),l(j),grid);
      b_tmp(j) = mean(difference.subvec( value_min , value_max)*
        normalization_values( m(j)-1 , l(j)-1 )) / tmp_basis(m(j)-1);
    }
    // change a l_j
    if(choice == 3){
      // choose an interval
      j = sample_weight(choice_prob_interval);

      // Simulate a new l_j
      value_max = l(j) + dl ;
      value_min = l(j) - dl ;
      if(value_max > p_l){
        value_max = p_l ;
        value_min = p_l - 2*dl ; // PMG - 2024-07-15
      }
      if(value_min < 1){
        value_min = 1 ;
      }
      probs = ones<vec>( value_max - value_min + 1 ) ;
      l_tmp(j)  = sample_weight( probs ) + value_min ;
    }
    // birth
    if(choice == 4){
      // compute the difference ...
      difference = abs(posterior_expe - d);
      // ... and its smoothed version
      difference = moving_average_cpp(difference,smooth_param);

      if(sum(abs(difference)) > 0){
        // update k
        k_tmp = k+1;
        // Simulate a new m
        new_m = sample_weight(abs(difference)) +1;
        m_tmp = zeros<vec>(k_tmp);
        m_tmp.subvec(0,k-1) = m;
        m_tmp(k_tmp-1)      = new_m;

        // Simulate a new l
        difference_l = zeros<vec>(dl);
        for(int o=0 ; o<dl ; ++o){
          value_min = new_m-o-1; if(value_min < 1) value_min = 1;
          value_max = new_m+o+1; if(value_max > p) value_max = p;
          difference_l(o) = exp(-std::pow(mean(abs(
            difference.subvec(value_min-1,value_max-1) - difference(new_m-1)
          )),2));
        }

        new_l = sample_weight(difference_l) + 1 ;

        l_tmp = zeros<vec>(k_tmp);
        l_tmp.subvec(0,k-1) = l;
        l_tmp(k_tmp-1)      = new_l;

        // Simulate a new b (from the smoothed difference)
        if( new_m - new_l -1 > 0  ) value_min = new_m - new_l -1; else value_min = 1;
        if( new_m + new_l +1 < p+1) value_max = new_m + new_l +1; else value_max = p;

        tmp_basis  = uniform_cpp(new_m,new_l,grid);
        difference = (posterior_expe - d);
        difference = moving_average_cpp(difference,smooth_param);
        new_b = mean( difference.subvec(value_min-1 , value_max-1) )*
          normalization_values( new_m-1 , new_l-1 )/ tmp_basis(new_m-1);

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
      m_tmp = zeros<vec>(k_tmp);
      l_tmp = zeros<vec>(k_tmp);

      b_tmp = vec_drop_k(b,j);
      m_tmp = vec_drop_k(m   ,j);
      l_tmp = vec_drop_k(l   ,j);
    }

    // Compute the acceptance probability
    d_tmp      = compute_beta_cpp(b_tmp,m_tmp,l_tmp,grid,p,k_tmp,basis,
                                  normalization_values);
    d_loss_tmp = loss_cpp(d_tmp,grid,posterior_expe);

    proba_acceptance = exp( -( d_loss_tmp-d_loss )/ Temperature );

    // Accept/reject
    u = R::runif(0,1) ;
    if(u < proba_acceptance){
      b = zeros<vec>(k_tmp) ;
      l = zeros<vec>(k_tmp) ;
      m = zeros<vec>(k_tmp) ;
      accepted  = 1         ;
      b = b_tmp.subvec(0,k_tmp-1) ;
      m = m_tmp.subvec(0,k_tmp-1) ;
      l = l_tmp.subvec(0,k_tmp-1) ;
      k = k_tmp ;
      d = d_tmp ;
      d_loss    = d_loss_tmp ;
    }

    // Update the trace
    trace.row(i+1).subvec( 0       ,           k-1) = trans(b.subvec( 0,k-1)) ;
    trace.row(i+1).subvec( k_max   , k_max   + k-1) = trans(m.subvec( 0,k-1))         ;
    trace.row(i+1).subvec( k_max*2 , k_max*2 + k-1) = trans(l.subvec( 0,k-1))         ;
    trace(i+1,3*k_max  ) = accepted ;
    trace(i+1,3*k_max+1) = choice   ;
    trace(i+1,3*k_max+2) = k        ;
    trace(i+1,3*k_max+3) = d_loss   ;
  }

  ////////////////////////////////////////////////////////////////
  //////////////// Return the trace and the parameters
  ////////////////////////////////////////////////////////////////
  if(verbose) Rcpp::Rcout << "\t Return the result." <<  std::endl;
  return  List::create(_["trace"]         =trace,
                       _["posterior_expe"]=posterior_expe);
}

// Perform the Simulated Annealing algorithm to minimize the loss function
// [[Rcpp::export]]
arma::mat dposterior_cpp (arma::mat & rposterior, arma::vec & y, unsigned N,
                          arma::vec & K, List & potential_intervals,
                          List & potential_intervals_dims, arma::vec & p_l,
                          unsigned Q, double tol){
  ////////////////////////////////////////////////////////////////
  //////////////// Initialization
  ////////////////////////////////////////////////////////////////
  int      sum_K = sum(K)       ;
  arma::mat      res = zeros<mat>(N,6);
  unsigned n   = y.size()       ;
  double   RSS = 0              ;

  arma::vec    b_tilde = zeros<vec>(sum_K+1) ;
  double sigma_sq ;
  arma::vec    m = zeros<vec>(sum_K) ;
  arma::vec    l = zeros<vec>(sum_K) ;
  arma::vec    m_temp   ;
  arma::vec    l_temp   ;

  double posterior_d     ;
  double log_posterior_d ;
  double prior_d         ;
  double log_prior_d     ;
  double lkh             ;
  double log_lkh         ;
  double v = 5      ;
  unsigned count         ;
  unsigned count2        ;

  arma::mat Sigma_b_tilde_inv_prior     ;
  arma::vec a_l = 5*K ; // fixed hyperparameter

  arma::mat    x_tilde  = arma::ones<arma::mat>(n,sum_K+1) ;
  // arma::mat lambda_id0  = arma::zeros<arma::mat>(sum_K+1,sum_K+1) ;
  // lambda_id0(0,0) = 100*var(y); // Weakly information prior: should be the same that in the Bliss_Gibbs_Sampler_cpp function
  // for( int i=1 ; i<sum_K+1; ++i){
  //   lambda_id0(i,i) = lambda ;
  // }

  double v0 = 100*var(y); // Weakly informative prior
  if(v0 > 1/tol){
    v0 = 0.5/tol;
  }
  ////////////////////////////////////////////////////////////////
  //////////////// For all data
  ////////////////////////////////////////////////////////////////
  for(unsigned j=0 ; j<N ; ++j ){
    // load the parameter value
    b_tilde(0) = rposterior(j,3*sum_K  ) ;
    sigma_sq   = rposterior(j,3*sum_K+1) ;
    count = 0 ;
    count2 = 0;
    for(unsigned q=0 ; q<Q ; ++q){
      b_tilde.subvec(count+1,count+K(q)) =
        trans(rposterior.row(j).subvec(  count2,  count2 + K(q)-1)) ;
      count2 += K(q);
      m.subvec(count,count+K(q)-1) =
        trans(rposterior.row(j).subvec(  count2,  count2 + K(q)-1)) ;
      count2 += K(q);
      l.subvec(count,count+K(q)-1) =
        trans(rposterior.row(j).subvec(  count2,  count2 + K(q)-1)) ;
      count2 += K(q);
      count += K(q);
    }

    // compute x_tilde
    int count = 0;
    for( unsigned q=0 ; q<Q ; ++q){
      for(unsigned k=0 ; k<K(q) ; ++k) {
        m_temp = m.subvec(count,count+K(q)-1) ;
        l_temp = l.subvec(count,count+K(q)-1) ;
        NumericVector potential_intervals_temp = potential_intervals[q];
        vec potential_intervals_dims_temp = potential_intervals_dims[q];
        x_tilde.col(k+1+count) = potential_intervals_extract(potential_intervals_temp,
                    m_temp(k),l_temp(k),potential_intervals_dims_temp);
      }
      count += K(q);
    }

    // compute Sigma
    Sigma_b_tilde_inv_prior =
      compute_Sigma_b_tilde_inv_prior (Q,K,n, x_tilde,sum_K,v0,v) ;

    RSS   = sum(square(y - x_tilde * b_tilde)) ;
    // compute the (log) likelihood
    log_lkh = -1./(2.*sigma_sq) * RSS - n/2. * log(sigma_sq);
    lkh     = exp(log_lkh);

    // compute the (log) prior density
    log_prior_d = -1./(2.*sigma_sq) * dot(b_tilde, Sigma_b_tilde_inv_prior * b_tilde) -
      (sum_K+3.)/2. * log(sigma_sq)  + log(det(Sigma_b_tilde_inv_prior)) / 2.;

    count =0;
    for(unsigned q=0 ; q<Q ; ++q){
      l_temp = l.subvec(count,count+K(q)-1);
      log_prior_d -=  a_l(q) * sum(l_temp) / p_l(q);
      count += K(q) ;
    }
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

  ////////////////////////////////////////////////////////////////
  //////////////// Output
  ////////////////////////////////////////////////////////////////
  return(res);
}


