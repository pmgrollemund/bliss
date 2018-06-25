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

// square
inline double sq(double x){
  return x * x;
}

// The R function : ginv (generalized matrix inversion using SVD decomposition)
// [[Rcpp::export]]
arma::mat ginv_cpp (arma::mat & X, double tol){
  int p;
  mat u;
  vec s;
  mat v;

  svd(u,s,v,X);
  p = s.size();

  tol   = tol * s(0);
  mat S = zeros<mat>(p,p);

  for( unsigned i=0 ; i<p; ++i){
    if( s(i) > tol ) S(i,i) = 1/s(i);
  }

  return( v * (S * trans(u)) );
}

// Termwise product
arma::vec termwise_product (arma::vec & v, arma::vec & u){
  vec res = v;
  for(int i=0 ; i<v.size() ; ++i){
    res(i) = res(i) * u(i);
  }
  return(res);
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
arma::vec vec_drop_k(arma::vec & vecteur, int k){
  unsigned n = vecteur.n_rows;
  if (k==0   && n > 1) return vecteur.subvec(1,n-1);
  if (k==n-1 && n > 1) return vecteur.subvec(0,n-2);
  vec res = zeros<vec>(n-1);

  res.subvec( 0  , k-1 ) = vecteur.subvec(0,k-1);
  res.subvec( k  , n-2 ) = vecteur.subvec(k+1,n-1);
  return res;
}

// Return the matrix mat[,-k].
arma::mat mat_drop_col_k(arma::mat & matrix, int k){
  unsigned n = matrix.n_rows;
  unsigned p = matrix.n_cols;
  if(k < 0  || k > (p-1)) return matrix;
  if(k==0   && p>1) return matrix.cols(1,p-1);
  if(k==p-1 && p>1) return matrix.cols(0,p-2);

  mat res = zeros<mat>(n,p-1);
  res.cols( 0 , k-1 ) = matrix.cols(0   , k-1);
  res.cols( k , p-2 ) = matrix.cols(k+1 , p-1);

  return res;
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
arma::mat sqrt_mat (arma::mat & X){
  int p;
  mat u;
  vec s;
  mat v;

  svd(u,s,v,X);
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
  vec Y = randn<vec>(ncols);
  VarCovar = chol(VarCovar); // XXXXXXXXXXXX

  return  mu + sqrt(sigma_sq) * trans(trans(Y) * VarCovar);
}

// Compute a trapezoidal approximation of area under curve.
// [[Rcpp::export]]
double integrate_trapeze (arma::vec & x, arma::vec & y){
  vec diff_x = vec_drop_k(x,0) - vec_drop_k(x,x.size()-1);
  vec cumu_y = vec_drop_k(y,0) + vec_drop_k(y,y.size()-1);
  return sum( termwise_product (diff_x  , cumu_y ) )/2 ;
}

//  Compute the norm of a vector.
double norm_fct(arma::vec & x,arma::vec & y){
  vec tmp = zeros<vec>(x.size());
  for(int i=0 ; i<tmp.size() ; ++i){
    tmp(i) = sq( y(i) );
  }
  double res;
  res = sqrt(integrate_trapeze(x,tmp));

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

// Use to compute an uniform function, see function compute_beta.
// [[Rcpp::export]]
arma::vec uniform_cpp_unnormalized (int m, int l, arma::vec & grid){
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


// Use to compute a triangular function, see function compute_beta.
// [[Rcpp::export]]
arma::vec triangular_cpp_unnormalized (int m, int l, arma::vec & grid){
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
    tmp2 =  exp( - 9*sq(i/l_double)/2) ;
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
arma::vec gaussian_cpp_unnormalized (int m, int l, arma::vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  exp( - 9*sq(i/l_double)/2) ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
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
    tmp2 =  1-sq(i/l_double) ;
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
arma::vec Epanechnikov_cpp_unnormalized (int m, int l, arma::vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  1-sq(i/l_double) ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
  return res;
}

// compute_beta in cpp.
// [[Rcpp::export]]
arma::vec compute_beta_cpp (arma::vec & b, arma::vec & m, arma::vec & l,
                          arma::vec & grid, int p, int K, std::string basis, arma::mat & normalization_values ){
  vec res = zeros<vec>(p) ;

  if(basis == "uniform"){
    for(int i=0 ; i<K ; ++i){
      res = res + b(i)/normalization_values( m(i)-1 , l(i)-1 ) *
        uniform_cpp(m(i),l(i),grid);
    }
  }
  if(basis == "uniform_unnormalized"){
    for(int i=0 ; i<K ; ++i){
      res = res + b(i) * uniform_cpp_unnormalized(m(i),l(i),grid);
    }
  }
  if(basis == "triangular"){
    for(int i=0 ; i<K ; ++i){
      res = res + b(i)/normalization_values( m(i)-1 , l(i)-1 )  *
        triangular_cpp(m(i),l(i),grid);
    }
  }
  if(basis == "triangular_unnormalized"){
    for(int i=0 ; i<K ; ++i){
      res = res + b(i) * triangular_cpp_unnormalized(m(i),l(i),grid);
    }
  }
  if(basis == "gaussian"){
    for(int i=0 ; i<K ; ++i){
      res = res + b(i)/normalization_values( m(i)-1 , l(i)-1 )  *
        gaussian_cpp(m(i),l(i),grid);
    }
  }
  if(basis == "gaussian_unnormalized"){
    for(int i=0 ; i<K ; ++i){
      res = res + b(i) * gaussian_cpp_unnormalized(m(i),l(i),grid);
    }
  }
  if(basis == "Epanechnikov"){
    for(int i=0 ; i<K ; ++i){
      res = res + b(i)/normalization_values( m(i)-1 , l(i)-1 )  *
        Epanechnikov_cpp(m(i),l(i),grid);
    }
  }
  if(basis == "Epanechnikov_unnormalized"){
    for(int i=0 ; i<K ; ++i){
      res = res + b(i) * Epanechnikov_cpp_unnormalized(m(i),l(i),grid);
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
arma::cube potential_intervals (arma::mat & X, arma::vec & grid, int l_max,
                                std::string basis){
  int n = X.n_rows ;
  int p = X.n_cols ;
  vec tub;

  arma::cube res(p,l_max,n+1);
  vec tmp;
  vec X_tmp;
  vec tmp2;
  if(basis == "uniform"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = uniform_cpp(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "uniform_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = uniform_cpp_unnormalized(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "triangular"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = triangular_cpp(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "triangular_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = triangular_cpp_unnormalized(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "gaussian"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = gaussian_cpp(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "gaussian_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = gaussian_cpp_unnormalized(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "Epanechnikov"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = Epanechnikov_cpp(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "Epanechnikov_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = Epanechnikov_cpp_unnormalized(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }

  // pour plus de clarete, je recommence une nouvelle double boucle plutot
  // que de copier ce bout de code dans chacun des if. A changer ?
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

// Compute all the alternative for the value of the intergral for all m and l.
// [[Rcpp::export]]
arma::cube potential_intervals_List(List & X_list, List & grids,arma::vec & l_max_vec,
                                    CharacterVector & basis_vec, int q){
  mat X = as<mat>(X_list[q]);
  vec grid = as<vec>(grids[q]);
  int l_max = l_max_vec(q);
  std::string basis = as<std::string>(basis_vec(q));

  int n = X.n_rows ;
  int p = X.n_cols ;
  vec tub;

  arma::cube res(p,l_max,n+1);
  vec tmp;
  vec X_tmp;
  vec tmp2;
  if(basis == "uniform"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = uniform_cpp(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "uniform_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = uniform_cpp_unnormalized(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "triangular"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = triangular_cpp(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "triangular_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = triangular_cpp_unnormalized(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "gaussian"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = gaussian_cpp(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "gaussian_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = gaussian_cpp_unnormalized(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "Epanechnikov"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = Epanechnikov_cpp(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }
  if(basis == "Epanechnikov_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = Epanechnikov_cpp_unnormalized(i+1,j+1,grid);
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;

          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }
  }

  // pour plus de clarete, je recommence une nouvelle double boucle plutot
  // que de copier ce bout de code dans chacun des if. A changer ?
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

// Concatenate two vectors.
arma::vec concatenate (arma::vec & v, arma::vec & u){
  vec res = zeros<vec>( v.size()+u.size() ) ;

  res.subvec( 0        ,   v.size()-1 ) = v;
  res.subvec( v.size() , res.size()-1 ) = u;
  return res;
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

// Compute the matrix V (for a Ridge Zellner prior)
arma::mat compute_W_inv_RZ (int K, double v0, double g, arma::mat & X_tilde, arma::mat & V){
  mat V_inv = zeros<mat>(K+1,K+1);

  V_inv(0,0) = 1/v0;
  V_inv.submat(1,1,K,K) = ( trans(X_tilde.cols(1,K)) *
    X_tilde.cols(1,K) + V.submat(1,1,K,K) )  /g;

  return V_inv;
}

// Old wrong version of the compute_W_inv_RZ_List
// Compute the matrix V (for a Ridge Zellner prior)
// (for Q functional covaribles)
arma::mat compute_W_inv_RZ2 (int Q, arma::vec K, double g, arma::mat & X_tilde, arma::mat & lambda_id,
                       int sum_K){
  mat W_inv = zeros<mat>(sum_K+1,sum_K+1);

  W_inv(0,0) = 1/lambda_id(0,0);
  int count = 0;
  for( unsigned q=0 ; q<Q ; ++q){
    W_inv.submat(1+count,1+count,K(q)+count,K(q)+count) =
      ( trans(X_tilde.cols(1+count,K(q)+count)) *
      X_tilde.cols(1+count,K(q)+count) +
      lambda_id.submat(1+count,1+count,K(q)+count,K(q)+count) )  /g;
    count = count + K(q);
  }

  return W_inv;
}

// Compute the matrix V (for a Ridge Zellner prior)
// (for Q functional covaribles)
arma::mat compute_W_inv_RZ_List (int Q, arma::vec & K, double g, arma::mat & X_tilde, int sum_K,
                                 arma::mat & lambda_id0){
  mat W_inv = zeros<mat>(sum_K+1,sum_K+1);
  mat lambda_id = lambda_id0 ;

  mat X_tilde_temp = mat_drop_col_k(X_tilde,0);
  mat u;
  vec s;
  mat v;
  svd(u,s,v,X_tilde_temp);

  W_inv(0,0) = 1/lambda_id(0,0);
  for( unsigned i=1 ; i<sum_K+1; ++i){
    lambda_id(i,i) = lambda_id(i,i) * max(s);
  }

  int count = 0;
  for( unsigned q=0 ; q<Q ; ++q){
    W_inv.submat(1+count,1+count,K(q)+count,K(q)+count) =
      ( trans(X_tilde.cols(1+count,K(q)+count)) *
      X_tilde.cols(1+count,K(q)+count) +
      lambda_id.submat(1+count,1+count,K(q)+count,K(q)+count) )  /g;
    count = count + K(q);
  }

  return W_inv;
}

//#############################################################################
//############################ Auxiliary functions ############################
//#############################################################################

// Extract a subvector from the cube all_intervals with a m_k and a l_k.
// [[Rcpp::export]]
arma::vec all_intervals_extract (NumericVector & all_intervals, int mk , int lk,
                                 arma::vec & dims) {
  vec res = zeros<vec>(dims(2));
  for (int i = 0; i < dims(2); i++) {
    res(i) = cube_extract(all_intervals, mk - 1, lk - 1, i, dims);
  }
  return res;
}

// Extract a submatrix from the cube all_intervals with the vectors m and l.
arma::mat extraire_X_tilde(arma::vec & m, arma::vec & l, NumericVector & all_intervals,
                           arma::vec & dims){
  int K = m.size();
  mat res = ones<mat>(dims(2), K + 1);
  for(int i=0; i<K ; i++){
    res.col(i+1) = all_intervals_extract(all_intervals,m(i),l(i),dims);
  }
  return res;
}

// Update the parameter l_k                                                     // single covariate version : the matrix W_inv does not change here when all the possibility of l_k are evaluated
int lk_update (int k, arma::vec & Y, arma::vec & b, double  sigma_sq, arma::vec & m,
               arma::vec & l, arma::mat & X_tilde, NumericVector & all_intervals,
               arma::vec & all_intervals_dims, arma::vec & l_alternative, arma::mat & V_inv,
               arma::vec & eta, arma::vec & phi_l) {
  double aux;
  const int lmax = l_alternative.size();
  // Compute mu.mlk
  vec mu_mlk = (Y - mat_drop_col_k(X_tilde,k+1) *
    vec_drop_k(b,k+1) ) / b(k+1);
  // Compute the probabilities
  vec probs = zeros<vec>(lmax);
  for(int  i=0 ; i<lmax ; i++){
    vec intervals_lki = all_intervals_extract(all_intervals,m(k),
                                              l_alternative(i),all_intervals_dims);
    aux = sq(b(k+1)) * dot(intervals_lki - mu_mlk ,
             intervals_lki - mu_mlk) / (2*sigma_sq) +
               1/(2*sigma_sq) * dot(b - eta, V_inv *
               (b - eta));
    probs(i) = exp( - aux ) * phi_l(i);
  }
  // Simulate a lk
  int lk = sample_weight(probs) + 1 ;

  return lk;
}

//Update the parameter m_k                                                      // single covariate version : the matrix W_inv does not change here when all the possibility of m_k are evaluated
int mk_update (int k, arma::vec & Y, arma::vec & b, double  sigma_sq, arma::vec & m, arma::vec & l,
               arma::mat & X_tilde, NumericVector & all_intervals, arma::vec & all_intervals_dims,
               arma::vec & m_alternative, arma::mat & V_inv, arma::vec & eta, arma::vec & phi_m) {
  double aux;
  const int p = m_alternative.size();
  // Compute mu.mmk
  vec mu_mmk = (Y - mat_drop_col_k(X_tilde,k+1) * vec_drop_k(b,k+1) )
    / b(k+1);
  // Compute the probabilities
  vec probs = ones<vec>(p);
  for(int  i=0 ; i<p ; ++i){
    vec intervals_mki = all_intervals_extract(all_intervals,m_alternative(i),
                                              l(k),all_intervals_dims);
    aux = sq(b(k+1)) * dot(intervals_mki - mu_mmk ,
             intervals_mki - mu_mmk) / (2*sigma_sq)  +
               1/(2*sigma_sq) * dot(b - eta, V_inv *
               (b - eta));
    probs(i) = exp( - aux ) * phi_m(i);
  }
  // Simulate a mk
  int mk = sample_weight(probs) + 1 ;

  return mk;
}

// Old wrong version of the mk_update_List
// Update the parameter m_k
// [[Rcpp::export]]
int mk_update2 (int count, int k, arma::vec & Y, arma::vec & beta_tilde, double sigma_sq,
                arma::vec & m_q, arma::vec & l_q, arma::mat & X_tilde,
                NumericVector & all_intervals_q, arma::vec & all_intervals_dims_q,
                arma::vec & m_alternative_q, arma::mat & W_inv, arma::vec & eta_tilde,
                arma::vec & phi_m_q, int p_q, std::string prior_beta) {
  double aux;
  vec aux2 = zeros<vec>(p_q);
  double beta_qk  = beta_tilde(count + k + 1); // mis des + 1
  colvec beta_mqk = vec_drop_k(beta_tilde,count + k + 1); // mis des + 1
  mat X_tilde_mqk = mat_drop_col_k(X_tilde,count + k + 1);  // mis des + 1

  // Compute the probabilities
  vec probs = ones<vec>(p_q);
  vec X_tilde_qki = zeros<vec>(all_intervals_dims_q(2)) ;
  for(int  i=0 ; i<p_q ; ++i){
    X_tilde_qki = all_intervals_extract(all_intervals_q,m_alternative_q(i),
                                        l_q(k),all_intervals_dims_q);
    aux = (dot(X_tilde_qki,X_tilde_qki) -2/beta_qk *
      dot(Y - X_tilde_mqk*beta_mqk, X_tilde_qki)) / (2 * sigma_sq);
    aux = aux * sq(beta_qk) ;
    if(prior_beta == "Ridge_Zellner")
      aux = aux + 1/(2*sigma_sq) * dot(beta_tilde - eta_tilde, W_inv *
        (beta_tilde - eta_tilde));
    aux2(i) = aux;
  }
  double min_aux = min(aux2);
  for(int  i=0 ; i<p_q ; ++i){
    aux2(i)  = aux2(i) - min_aux;
    probs(i) = exp( - aux2(i) ) * phi_m_q(i);
  }
  // Simulate a mk
  int mk = sample_weight(probs) + 1 ;

  return mk;
}


// Update the parameter m_k
// [[Rcpp::export]]
int mk_update_List (int count, int k, arma::vec & Y, arma::vec & beta_tilde, double sigma_sq,
                    arma::vec & m_q, arma::vec & l_q, arma::mat & X_tilde,
                    NumericVector & all_intervals_q, arma::vec & all_intervals_dims_q,
                    arma::vec & m_alternative_q, arma::vec & eta_tilde,
                    arma::vec & phi_m_q, int p_q, std::string prior_beta, int Q,
                    arma::vec K, double g, int sum_K, arma::mat & lambda_id0) {
  double aux;
  vec aux2 = zeros<vec>(p_q);
  vec aux3 = zeros<vec>(p_q);
  mat X_tilde_mqk = mat_drop_col_k(X_tilde,count + k + 1);  // mis des + 1
  mat X_tilde_temp = X_tilde;
  mat W_inv_temp;

  // Compute the probabilities
  vec probs = ones<vec>(p_q);
  vec X_tilde_qki = zeros<vec>(all_intervals_dims_q(2)) ;
  for(int  i=0 ; i<p_q ; ++i){
    X_tilde_qki = all_intervals_extract(all_intervals_q,m_alternative_q(i),
                                        l_q(k),all_intervals_dims_q);

    X_tilde_temp.col(count + k + 1) = X_tilde_qki;
    aux2(i) = dot( Y - X_tilde_temp * beta_tilde ,
         Y - X_tilde_temp * beta_tilde ) /(2*sigma_sq) ;

    aux3(i) = 1 ; // if prior is not the Rigde Zellner prior
    if(prior_beta == "Ridge_Zellner"){
      W_inv_temp = compute_W_inv_RZ_List(Q,K,g,X_tilde_temp,sum_K,lambda_id0);

      aux2(i) = aux2(i) +
        dot( beta_tilde , W_inv_temp * beta_tilde ) / (2*sigma_sq);
      aux3(i) = sqrt( det(W_inv_temp) );
    }
  }

  double min_aux = min(aux2);
  for(int  i=0 ; i<p_q ; ++i){
    aux2(i)  = aux2(i) - min_aux;
    // aux3(i)  = aux3(i) / max_aux;
    probs(i) = aux3(i) * exp( - aux2(i) ) * phi_m_q(i);
  }
  // Simulate a mk
  int mk = sample_weight(probs) + 1 ;

  return mk;
}


// Old wrong version of the lk_update_List
// Update the parameter l_k
int lk_update2 (int count, int k, arma::vec & Y, arma::vec & beta_tilde, double sigma_sq,
                arma::vec & m_q, arma::vec & l_q, arma::mat & X_tilde,
                NumericVector & all_intervals_q, arma::vec & all_intervals_dims_q,
                arma::vec & l_alternative_q, arma::mat & W_inv, arma::vec & eta_tilde,
                arma::vec & phi_l_q, int lmax_q, std::string prior_beta) {
  double aux;
  vec aux2 = zeros<vec>(lmax_q);
  double beta_qk  = beta_tilde(count + k + 1);
  colvec beta_mqk = vec_drop_k(beta_tilde,count + k + 1);
  mat X_tilde_mqk = mat_drop_col_k(X_tilde,count + k + 1);

  // Compute the probabilities
  vec probs = ones<vec>(lmax_q);
  vec X_tilde_qki = zeros<vec>(all_intervals_dims_q(2)) ;
  for(int  i=0 ; i<lmax_q ; ++i){
    X_tilde_qki = all_intervals_extract(all_intervals_q,m_q(k),
                                        l_alternative_q(i),all_intervals_dims_q);

    aux = sq(beta_qk)* (dot(X_tilde_qki,X_tilde_qki) -2/beta_qk *
      dot(Y - X_tilde_mqk*beta_mqk, X_tilde_qki)) / (2 * sigma_sq);
    if(prior_beta == "Ridge_Zellner")
      aux = aux + 1/(2*sigma_sq) * dot(beta_tilde - eta_tilde, W_inv *
        (beta_tilde - eta_tilde));
    aux2(i) = aux;
  }
  double min_aux = min(aux2);
  for(int  i=0 ; i<lmax_q ; ++i){
    aux2(i)  = aux2(i) - min_aux;
    probs(i) = exp( - aux2(i) ) * phi_l_q(i);
  }
  // Simulate a lk
  int lk = sample_weight(probs) + 1 ;

  return lk;
}


// Update the parameter l_k
// [[Rcpp::export]]
int lk_update_List (int count, int k, arma::vec & Y, arma::vec & beta_tilde, double sigma_sq,
                    arma::vec & m_q, arma::vec & l_q, arma::mat & X_tilde,
                    NumericVector & all_intervals_q, arma::vec & all_intervals_dims_q,
                    arma::vec & l_alternative_q, arma::vec & eta_tilde,
                    arma::vec & phi_l_q, int lmax_q, std::string prior_beta, int Q,
                    arma::vec K, double g, int sum_K, arma::mat & lambda_id0) {
  double aux;
  vec aux2 = zeros<vec>(lmax_q);
  vec aux3 = zeros<vec>(lmax_q);
  mat X_tilde_mqk = mat_drop_col_k(X_tilde,count + k + 1);  // mis des + 1
  mat X_tilde_temp = X_tilde;
  mat W_inv_temp;

  // Compute the probabilities
  vec probs = ones<vec>(lmax_q);
  vec X_tilde_qki = zeros<vec>(all_intervals_dims_q(2)) ;
  for(int  i=0 ; i<lmax_q ; ++i){
    X_tilde_qki = all_intervals_extract(all_intervals_q,m_q(k),
                                        l_alternative_q(i),
                                        all_intervals_dims_q);

    X_tilde_temp.col(count + k + 1) = X_tilde_qki;
    aux2(i) = dot( Y - X_tilde_temp * beta_tilde ,
         Y - X_tilde_temp * beta_tilde ) /(2*sigma_sq) ;

    aux3(i) = 1 ; // if prior is not the Rigde Zellner prior
    if(prior_beta == "Ridge_Zellner"){
      W_inv_temp = compute_W_inv_RZ_List(Q,K,g,X_tilde_temp,sum_K,lambda_id0);

      aux2(i) = aux2(i) +
        dot( beta_tilde , W_inv_temp * beta_tilde ) / (2*sigma_sq);
      aux3(i) = sqrt( det(W_inv_temp) );
    }
  }

  double min_aux = min(aux2);
  for(int  i=0 ; i<lmax_q ; ++i){
    aux2(i)  = aux2(i) - min_aux;
    // aux3(i)  = aux3(i) / max_aux;
    probs(i) = aux3(i) * exp( - aux2(i) ) * phi_l_q(i);
  }
  // Simulate a lk
  int lk = sample_weight(probs) + 1 ;

  return lk;
}

// update the parameter sigma_sq
double sigma_sq_update (arma::vec & Y, arma::vec & b, arma::vec & eta, arma::mat & V_inv,
                        arma::mat & X_tilde, double a, double b, int n, int K) {
  double a_tmp     = K+n+1 ;
  double a_star    = a + a_tmp/2 ;

  vec Y_tmp        = Y - X_tilde * b ;
  double Y_tmp2    = dot(Y_tmp,Y_tmp) ;
  vec b_tmp     = b - eta ;
  double b_tmp2 = dot(b_tmp, V_inv * b_tmp) ;

  double b_star    = b + 0.5*( Y_tmp2 + b_tmp2 ) ;

  double res = 1. / (R::rgamma(a_star, 1/b_star) );

  return res ;
}
// update the parameter sigma_sq
double sigma_sq_update_List (arma::vec & Y, arma::vec & beta_tilde, arma::vec & eta_tilde,
                             arma::mat & W_inv, arma::mat & X_tilde, double a, double b,
                             int n, int sum_K) {
  double a_tmp     = sum_K+n+1 ;
  double a_star    = a + a_tmp/2 ;

  vec Y_tmp        = Y - X_tilde * beta_tilde ;
  double Y_tmp2    = dot(Y_tmp,Y_tmp) ;
  vec beta_tilde_tmp     = beta_tilde - eta_tilde ;
  double beta_tilde_tmp2 = dot(beta_tilde_tmp, W_inv * beta_tilde_tmp) ;

  double b_star    = b + 0.5*( Y_tmp2 + beta_tilde_tmp2);
  double res = 1. / (R::rgamma(a_star, 1/b_star) );

  return res ;
}

// update the parameter b
// [[Rcpp::export]]
arma::vec beta_tilde_update (arma::vec & Y, double sigma_sq, arma::vec & eta_tilde, arma::mat & W_inv,
                             arma::mat & X_tilde, arma::mat & Sigma_beta_tilde_inv,
                       double tol) {
  vec mu_beta_tilde = W_inv * eta_tilde + trans(X_tilde) * Y;
  return mvrnormArma( ginv_cpp(Sigma_beta_tilde_inv,tol) * mu_beta_tilde ,
                      ginv_cpp(Sigma_beta_tilde_inv,tol), sigma_sq);
}

// Compute the loss function for a proposal d
// [[Rcpp::export]]
double loss_cpp (arma::vec & d, arma::vec & grid, arma::vec & posterior_expe){
  vec tmp  = d-posterior_expe ;

  return sq(norm_fct(grid, tmp) );
}

// Compute the decrease of the Temperature
double cooling_cpp (int i, double Temp){
  double res;
  res = Temp / log( ( i / 10)*10 + exp(1));
  return res;
}

arma::mat update_X_tilde (int Q, arma::vec & K, List & all_intervals,
                    List & all_intervals_dims, List & m, List & l,
                    arma::mat & X_tilde){
  int count = 0;
  for( unsigned q=0 ; q<Q ; ++q){
    for(unsigned k=0 ; k<K(q) ; ++k) {
      vec m_temp = m[q] ;
      vec l_temp = l[q] ;
      NumericVector all_intervals_temp = all_intervals[q];
      vec all_intervals_dims_temp = all_intervals_dims[q];

      X_tilde.col(k+1+count) = all_intervals_extract(all_intervals_temp,
                  m_temp(k),l_temp(k),all_intervals_dims_temp);
    }
    count = count + K(q);
  }

  return X_tilde;
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
List Bliss_Gibbs_Sampler_multiple_cpp (int Q, arma::vec & y, List & x, List & grids,
                                       int iter, arma::vec & K, CharacterVector & basis, // rajouter du grid_l ?
                                       double g, double lambda,arma::mat & V_tilde,List & probs_l,
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
  List normalization_values(Q);    // normalization_values is used to normalize the predictors
  List all_intervals(Q);      // will be contain all the projections
  List all_intervals_dims(Q); // will be contain the dim of the all_intervals's
  for( int q=0 ; q<Q ; ++q){
    arma::cube temp = potential_intervals_List (X, grids, lmax, basis,q) ;
    normalization_values[q]      = temp.slice(n);

    temp = temp.subcube(0,0,0,p(q)-1,lmax(q)-1,n-1);
    all_intervals[q] = temp;

    vec temp2 = zeros<vec>(3);
    temp2(0) = p(q)    ;
    temp2(1) = lmax(q) ;
    temp2(2) = n    ;
    all_intervals_dims[q] = temp2;
  }

  // Compute the matrix of lambda for the Ridge penalty of the Rigde
  // Zellner prior...
  int sum_K = sum(K);
  mat lambda_id0  = zeros<mat>(sum_K+1,sum_K+1) ;
  if(prior_beta == "Ridge_Zellner"){
    lambda_id0(0,0) = 100*var(Y);
    for( unsigned i=1 ; i<sum_K+1; ++i){
      lambda_id0(i,i) = lambda ;
    }
  }
  // ... or the constant matrix if V does not depend on the intervals
  mat W_inv ;
  if(prior_beta == "diag") W_inv = ginv_cpp(V_tilde,tol);

  // Determine the start point
  Rcpp::Rcout << "\t Determine the start point." <<  std::endl;
  double sigma_sq           ;
  vec beta_tilde            ;
  mat Sigma_beta_tilde_inv  ;

  bool success = false      ;
  mat R                     ;
  mat test                  ;

  List m(Q)                 ;
  List l(Q)                 ;
  mat X_tilde = ones<mat>(n,sum_K+1) ;

  // Try to determine a starting point which not leads to a non-invertible
  // matrix problem
  while(success == false){
    // Initialization of sigma_sq
    sigma_sq  = var(Y) ;

    // Initialization of the middle and length of the intervals
    for( unsigned q=0 ; q<Q ; ++q){
      vec probs_m_temp = probs_m[q];
      m[q]         = sample_weight(K(q),probs_m_temp) + 1 ;
      vec probs_l_temp = probs_l[q];
      l[q]         = sample_weight(K(q),probs_l_temp) + 1 ;
    }

    // Initialize the current X_tilde matrix (which depend on the intervals)
    int count = 0;
    for( unsigned q=0 ; q<Q ; ++q){
      for(unsigned k=0 ; k<K(q) ; ++k) {
        vec m_temp = m[q] ;
        vec l_temp = l[q] ;
        NumericVector all_intervals_temp = all_intervals[q];
        vec all_intervals_dims_temp = all_intervals_dims[q];

        X_tilde.col(k+1+count) = all_intervals_extract(all_intervals_temp,
                    m_temp(k),l_temp(k),all_intervals_dims_temp);
      }
      count = count + K(q);
    }

    // Initialize the current W_inv matrix (which depend on the intervals)
    // (W_inv is the covariance matrix of the Ridge Zellner prior) (without sig)
    if(prior_beta == "Ridge_Zellner")
      W_inv = compute_W_inv_RZ_List (Q,K, g, X_tilde,sum_K,lambda_id0) ;

    // Check if there is a non-invertible matrix problem
    Sigma_beta_tilde_inv = W_inv + trans(X_tilde) * X_tilde ;
    test            = ginv_cpp(Sigma_beta_tilde_inv,tol)    ;
    success         = accu(abs(test)) != 0                  ;
  }

  // Initialization of beta_tilde
  beta_tilde = mvrnormArma( eta_tilde , ginv_cpp(W_inv,tol) , sigma_sq) ;

  // Initialize the matrix trace
  mat trace = zeros<mat>(iter+1,3*sum_K+2);
  int count = 0;
  for( unsigned q=0 ; q<Q ; ++q){
    vec m_temp = m[q] ;
    vec l_temp = l[q] ;
    trace.row(0).subvec( 3*count        , 3*count+  K(q)-1) =
      trans(beta_tilde.subvec( 1+count , K(q)+count )) ;
    trace.row(0).subvec( 3*count+  K(q) , 3*count+2*K(q)-1)   = trans(m_temp) ;
    trace.row(0).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1)   = trans(l_temp) ;

    trace(0,3*sum_K  ) = beta_tilde(0) ;
    trace(0,3*sum_K+1) = sigma_sq      ;
    count = count + K(q) ;
  }

  // Initialize some variable used in the Gibbs loop
  List l_alternative(Q);
  for( int q=0 ; q<Q ; ++q){
    l_alternative[q] = sequence(1,lmax(q),1);
  }
  List m_alternative(Q);
  for( int q=0 ; q<Q ; ++q){
    m_alternative[q] = sequence(1,p(q),1);
  }

  // The Gibbs loop
  Rcpp::Rcout << "\t Start the Gibbs loop." <<  std::endl;
  for(unsigned i=1  ; i < iter+1 ; ++i ) {
    // Progress
    // std::cout << "\t " << i << std::endl;
    if( i % (iter / 10)  == 0)
      Rcpp::Rcout << "\t " << i / (iter / 100) << "%" << std::endl;

    // update sigma_sq
    sigma_sq = sigma_sq_update_List(Y,beta_tilde,eta_tilde,W_inv,X_tilde,a,b,
                                    n,sum_K) ;

    // update m
    count = 0 ;
    // count is used to browse some vec/mat when p(q) is not constant wrt q.
    for( unsigned q=0 ; q<Q ; ++q ){
      // Compute some quantities which do not vary with k
      vec m_q = m[q];
      vec l_q = l[q];
      int p_q = p(q);
      NumericVector all_intervals_q = all_intervals[q];
      vec all_intervals_dims_q      = all_intervals_dims[q];
      vec m_alternative_q = sequence(1,p_q,1) ;
      vec probs_m_q       = probs_m[q];

      for(int k=0 ; k<K(q) ; ++k){
        // update m_k
        if(prior_beta == "Ridge_Zellner"){
          m_q(k) = mk_update_List(count,k,Y,beta_tilde,sigma_sq,m_q,l_q,X_tilde,
              all_intervals_q,all_intervals_dims_q,m_alternative_q,
              eta_tilde,probs_m_q,p_q, prior_beta,Q,K,g,sum_K,lambda_id0);
        }
        if(prior_beta == "diag"){
          m_q(k) = mk_update2(count,k,Y,beta_tilde,sigma_sq,m_q,l_q,X_tilde,
              all_intervals_q,all_intervals_dims_q,m_alternative_q,W_inv,
              eta_tilde,probs_m_q,p_q, prior_beta);
        }

        // update the value "X_tilde"
        X_tilde = update_X_tilde(Q,K,all_intervals,all_intervals_dims,m,l,
                                 X_tilde);

        // Do not need to update W_inv for now (because the "mk_update_List"
        // function does not require this matrix as an input. It have to be
        // internal computed for each possibility of m_k)
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
      int lmax_q = lmax(q);
      NumericVector all_intervals_q = all_intervals[q];
      vec all_intervals_dims_q      = all_intervals_dims[q];
      vec l_alternative_q = sequence(1,lmax_q,1) ;
      vec probs_l_q         = probs_l[q];

      for(int k=0 ; k<K(q) ; ++k){
        // update l_k
        if(prior_beta == "Ridge_Zellner"){
          l_q(k) = lk_update_List(count,k,Y,beta_tilde,sigma_sq,m_q,l_q,X_tilde,
              all_intervals_q,all_intervals_dims_q,l_alternative_q,
              eta_tilde,probs_l_q,lmax_q, prior_beta,Q,K,g,sum_K,lambda_id0);
        }
        if(prior_beta == "diag"){
          l_q(k) = lk_update2(count,k,Y,beta_tilde,sigma_sq,m_q,l_q,X_tilde,
              all_intervals_q,all_intervals_dims_q,l_alternative_q,W_inv,
              eta_tilde,probs_l_q,lmax_q, prior_beta);
        }

        // update the value "X_tilde"
        X_tilde = update_X_tilde(Q,K,all_intervals,all_intervals_dims,m,l,
                                 X_tilde);

        // Do not need to update W_inv for now (because the "mk_update_List"
        // function does not require this matrix as an input. It have to be
        // internal computed for each possibility of l_k)
      }

      // Update the m_q value
      l[q] = l_q;
      // Update count
      count = count + K(q);
    }

    // update the value "W_inv" (only after the updating of the m's and l's)
    if(prior_beta == "Ridge_Zellner")
      W_inv = compute_W_inv_RZ_List (Q,K, g, X_tilde,sum_K,lambda_id0) ;

    // update the matrix Sigma_beta_tilde (only after the updating of
    // the m's and l's)
    Sigma_beta_tilde_inv = W_inv + trans(X_tilde) * X_tilde   ;
    test                 = ginv_cpp(Sigma_beta_tilde_inv,tol) ;
    success              = accu(abs(test)) != 0               ;

    // Try to determine an update which not leads to a non-invertible
    // matrix problem. If there is a problem, go back to the beginning of the
    // updating process.
    if(success){
      // update the beta_tilde
      beta_tilde = beta_tilde_update(Y,sigma_sq,eta_tilde,W_inv,X_tilde,
                                     Sigma_beta_tilde_inv,tol) ;

      // update the matrix trace
      count = 0;
      for( unsigned q=0 ; q<Q ; ++q){
        vec m_temp = m[q] ;
        vec l_temp = l[q] ;
        trace.row(i).subvec( 3*count        , 3*count+  K(q)-1) =
          trans(beta_tilde.subvec( 1+count , K(q)+count )) ;
        trace.row(i).subvec( 3*count+  K(q) , 3*count+2*K(q)-1) = trans(m_temp);
        trace.row(i).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1) = trans(l_temp);

        trace(i,3*sum_K  ) = beta_tilde(0) ;
        trace(i,3*sum_K+1) = sigma_sq      ;
        count = count + K(q) ;
      }
    }else{ //... go back to the beginning of the updating process.
      i     = i - 1 ;
      count = 0;
      for( unsigned q=0 ; q<Q ; ++q){
        beta_tilde.subvec( 1+count , K(q)+count ) =
          trans(trace.row(i).subvec( 3*count , 3*count+  K(q)-1))  ;
        m[q] = trans(trace.row(i).subvec( 3*count+  K(q) , 3*count+2*K(q)-1)) ;
        l[q] = trans(trace.row(i).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1)) ;

        beta_tilde(0) = trace(i,3*sum_K  ) ;
        sigma_sq      = trace(i,3*sum_K+1) ;
        count = count + K(q) ;
      }

      // update the value "X_tilde"
      X_tilde = update_X_tilde(Q,K,all_intervals,all_intervals_dims,m,l,
                               X_tilde);

      // update the value "W_inv"
      if(prior_beta == "Ridge_Zellner")
        W_inv = compute_W_inv_RZ_List (Q,K, g, X_tilde,sum_K,lambda_id0) ;
    }
  }

  // return the trace and the parameters
  Rcpp::Rcout << "\t Return the result." <<  std::endl;
  return  List::create(_["trace"]=trace,
                       _["param"]=List::create(_["a"]=a,
                                          _["b"]=b,
                                          _["phi_m"]=probs_m,
                                          _["phi_l"]=probs_l,
                                          _["K"]=K,
                                          _["eta_tilde"]=eta_tilde,
                                          _["l_max"]=lmax,
                                          _["all_intervals"]=all_intervals,
                                          _["grids"]=grids,
                                          _["normalization_values"]=normalization_values
                       ));
}

// Perform the Simulated Annealing algorithm to minimize the loss function
// [[Rcpp::export]]
List Bliss_Simulated_Annealing_cpp (int iter, arma::mat beta_functions, arma::vec & grid,
                                    int burnin, double Temp,int k_max,
                                    int l_max, int dm, int dl,
                                    int p,std::string basis, arma::mat & normalization_values){
  Rcpp::Rcout << "Simulated Annealing:" <<  std::endl;
  // Initialization
  Rcpp::Rcout << "\t Initialization." <<  std::endl;
  int N = beta_functions.n_rows;
  vec posterior_expe = zeros<vec>(p);
  vec posterior_var  = zeros<vec>(p);
  for(int i=0 ; i<p ; ++i){
    posterior_expe(i) = mean(beta_functions.col(i));
    posterior_var(i)  =  var(beta_functions.col(i));
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
  Rcpp::Rcout << "\t Determine the start point." <<  std::endl;
  probs = ones<vec>(k_max);
  k      = sample_weight( probs )+1;
  m      = zeros<vec>(k);
  l      = zeros<vec>(k);
  b   = zeros<vec>(k);

  probs = ones<vec>(p);
  m(0)   = sample_weight( probs )+1;
  probs = ones<vec>(l_max);
  l(0)   = sample_weight( probs )+1;

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
      difference = moving_average_cpp(difference,4);

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
  Rcpp::Rcout << "\t Start the loop." <<  std::endl;
  for(int i=0 ; i<iter ; ++i){
    Temperature = cooling_cpp(i,Temp);
    // Progress
    if( (i+1) % (iter / 10)  == 0)
      Rcpp::Rcout << "\t " << (i+1) / (iter / 100) << "%" << std::endl;
    // Initialize the proposal
    b_tmp = b;
    m_tmp    = m   ;
    l_tmp    = l   ;
    k_tmp    = k   ;
    accepted  = 0   ;
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
      difference = moving_average_cpp(difference,4);

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
  Rcpp::Rcout << "\t Return the result." <<  std::endl;
  return  List::create(_["trace"]         =trace,
                       _["posterior_expe"]=posterior_expe,
                       _["posterior_var"] =posterior_var);
}

// Perform the Simulated Annealing algorithm to minimize the loss function
// [[Rcpp::export]]
arma::mat dposterior_cpp (arma::mat & rposterior, arma::vec & y, unsigned N, unsigned K,
                    NumericVector & all_intervals, arma::vec & all_intervals_dims,
                    double lambda, double l_max){
  mat      res = zeros<mat>(N,6);
  unsigned n   = y.size()       ;
  double   RSS = 0              ;

  vec    beta_tilde = zeros<vec>(K+1) ;
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

  mat    X_tilde  = ones<mat>(n,K+1) ;
  mat lambda_id0  = zeros<mat>(K+1,K+1) ;
  lambda_id0(0,0) = 100*var(y);
  for( unsigned i=1 ; i<K+1; ++i){
    lambda_id0(i,i) = lambda ;
  }
  vec K_vec = ones<vec>(1) ;
  K_vec(0)  = K ;

  for(unsigned j=0 ; j<N ; ++j ){
    // load the parameter value
    beta_tilde(0)          = rposterior(j,3*K  ) ;
    beta_tilde.subvec(1,K) = trans(rposterior.row(j).subvec(  0,  K-1)) ;
    sigma_sq = rposterior(j,3*K+1) ;
    m        = trans(rposterior.row(j).subvec(  K,2*K-1)) ;
    l        = trans(rposterior.row(j).subvec(2*K,3*K-1)) ;

    // compute X_tilde
    for(unsigned k=0 ; k<K ; ++k) {
      X_tilde.col(k+1) = all_intervals_extract(all_intervals,m(k),l(k),
                  all_intervals_dims);
    }
    // compute Sigma
    W_inv = compute_W_inv_RZ_List (1,K_vec,n, X_tilde,K,lambda_id0) ;

    RSS   = sum(square(y - X_tilde * beta_tilde)) ;
    // compute the (log) likelihood
    log_lkh = -1./(2.*sigma_sq) * RSS - n/2. * log(sigma_sq);
    lkh     = exp(log_lkh);

    // compute the (log) prior density
    log_prior_d = -1./(2.*sigma_sq) * dot(beta_tilde, W_inv * beta_tilde) -
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
