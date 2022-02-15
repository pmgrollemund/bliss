// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ginv_cpp
arma::mat ginv_cpp(arma::mat& x, double tol);
RcppExport SEXP _bliss_ginv_cpp(SEXP xSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(ginv_cpp(x, tol));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma
arma::vec mvrnormArma(arma::vec mu, arma::mat VarCovar, double sigma_sq);
RcppExport SEXP _bliss_mvrnormArma(SEXP muSEXP, SEXP VarCovarSEXP, SEXP sigma_sqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type VarCovar(VarCovarSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_sq(sigma_sqSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(mu, VarCovar, sigma_sq));
    return rcpp_result_gen;
END_RCPP
}
// integrate_trapeze_cpp
double integrate_trapeze_cpp(arma::vec& x, arma::vec& y);
RcppExport SEXP _bliss_integrate_trapeze_cpp(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(integrate_trapeze_cpp(x, y));
    return rcpp_result_gen;
END_RCPP
}
// uniform_cpp
arma::vec uniform_cpp(int m, int l, arma::vec& grid);
RcppExport SEXP _bliss_uniform_cpp(SEXP mSEXP, SEXP lSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(uniform_cpp(m, l, grid));
    return rcpp_result_gen;
END_RCPP
}
// triangular_cpp
arma::vec triangular_cpp(int m, int l, arma::vec& grid);
RcppExport SEXP _bliss_triangular_cpp(SEXP mSEXP, SEXP lSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(triangular_cpp(m, l, grid));
    return rcpp_result_gen;
END_RCPP
}
// gaussian_cpp
arma::vec gaussian_cpp(int m, int l, arma::vec& grid);
RcppExport SEXP _bliss_gaussian_cpp(SEXP mSEXP, SEXP lSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussian_cpp(m, l, grid));
    return rcpp_result_gen;
END_RCPP
}
// Epanechnikov_cpp
arma::vec Epanechnikov_cpp(int m, int l, arma::vec& grid);
RcppExport SEXP _bliss_Epanechnikov_cpp(SEXP mSEXP, SEXP lSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(Epanechnikov_cpp(m, l, grid));
    return rcpp_result_gen;
END_RCPP
}
// compute_beta_cpp
arma::vec compute_beta_cpp(arma::vec& b, arma::vec& m, arma::vec& l, arma::vec& grid, int p, int K, std::string basis, arma::mat& normalization_values);
RcppExport SEXP _bliss_compute_beta_cpp(SEXP bSEXP, SEXP mSEXP, SEXP lSEXP, SEXP gridSEXP, SEXP pSEXP, SEXP KSEXP, SEXP basisSEXP, SEXP normalization_valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type l(lSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< std::string >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type normalization_values(normalization_valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_beta_cpp(b, m, l, grid, p, K, basis, normalization_values));
    return rcpp_result_gen;
END_RCPP
}
// compute_beta_sample_cpp
arma::mat compute_beta_sample_cpp(arma::mat& posterior_sample, int K, arma::vec& grid, int p, std::string& basis, arma::mat& normalization_values);
RcppExport SEXP _bliss_compute_beta_sample_cpp(SEXP posterior_sampleSEXP, SEXP KSEXP, SEXP gridSEXP, SEXP pSEXP, SEXP basisSEXP, SEXP normalization_valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type posterior_sample(posterior_sampleSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< std::string& >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type normalization_values(normalization_valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_beta_sample_cpp(posterior_sample, K, grid, p, basis, normalization_values));
    return rcpp_result_gen;
END_RCPP
}
// potential_intervals_List
arma::cube potential_intervals_List(List& x_list, List& grids, arma::vec& p_l_vec, CharacterVector& basis_vec, int q);
RcppExport SEXP _bliss_potential_intervals_List(SEXP x_listSEXP, SEXP gridsSEXP, SEXP p_l_vecSEXP, SEXP basis_vecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List& >::type x_list(x_listSEXP);
    Rcpp::traits::input_parameter< List& >::type grids(gridsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type p_l_vec(p_l_vecSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type basis_vec(basis_vecSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(potential_intervals_List(x_list, grids, p_l_vec, basis_vec, q));
    return rcpp_result_gen;
END_RCPP
}
// moving_average_cpp
arma::vec moving_average_cpp(arma::vec& v, int range);
RcppExport SEXP _bliss_moving_average_cpp(SEXP vSEXP, SEXP rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type range(rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(moving_average_cpp(v, range));
    return rcpp_result_gen;
END_RCPP
}
// potential_intervals_extract
arma::vec potential_intervals_extract(NumericVector& potential_intervals, int mk, int lk, arma::vec& dims);
RcppExport SEXP _bliss_potential_intervals_extract(SEXP potential_intervalsSEXP, SEXP mkSEXP, SEXP lkSEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type potential_intervals(potential_intervalsSEXP);
    Rcpp::traits::input_parameter< int >::type mk(mkSEXP);
    Rcpp::traits::input_parameter< int >::type lk(lkSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type dims(dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(potential_intervals_extract(potential_intervals, mk, lk, dims));
    return rcpp_result_gen;
END_RCPP
}
// update_mqk
void update_mqk(int count, int k, arma::vec& y, arma::vec& b_tilde, double sigma_sq, arma::vec& m_q, arma::vec& l_q, arma::mat x_tilde, NumericVector& potential_intervals_q, arma::vec& potential_intervals_dims_q, arma::vec& m_possible_q, int p_q, int Q, arma::vec K, double g, int sum_K, arma::mat& lambda_id0);
RcppExport SEXP _bliss_update_mqk(SEXP countSEXP, SEXP kSEXP, SEXP ySEXP, SEXP b_tildeSEXP, SEXP sigma_sqSEXP, SEXP m_qSEXP, SEXP l_qSEXP, SEXP x_tildeSEXP, SEXP potential_intervals_qSEXP, SEXP potential_intervals_dims_qSEXP, SEXP m_possible_qSEXP, SEXP p_qSEXP, SEXP QSEXP, SEXP KSEXP, SEXP gSEXP, SEXP sum_KSEXP, SEXP lambda_id0SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type count(countSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_tilde(b_tildeSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_sq(sigma_sqSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type m_q(m_qSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type l_q(l_qSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_tilde(x_tildeSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type potential_intervals_q(potential_intervals_qSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type potential_intervals_dims_q(potential_intervals_dims_qSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type m_possible_q(m_possible_qSEXP);
    Rcpp::traits::input_parameter< int >::type p_q(p_qSEXP);
    Rcpp::traits::input_parameter< int >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type sum_K(sum_KSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type lambda_id0(lambda_id0SEXP);
    update_mqk(count, k, y, b_tilde, sigma_sq, m_q, l_q, x_tilde, potential_intervals_q, potential_intervals_dims_q, m_possible_q, p_q, Q, K, g, sum_K, lambda_id0);
    return R_NilValue;
END_RCPP
}
// update_lqk
void update_lqk(int count, int k, arma::vec& y, arma::vec& b_tilde, double sigma_sq, arma::vec& m_q, arma::vec& l_q, arma::mat x_tilde, NumericVector& potential_intervals_q, arma::vec& potential_intervals_dims_q, arma::vec& l_possible_q, arma::vec& phi_l_q, int l_values_length_q, int Q, arma::vec K, double g, int sum_K, arma::mat& lambda_id0);
RcppExport SEXP _bliss_update_lqk(SEXP countSEXP, SEXP kSEXP, SEXP ySEXP, SEXP b_tildeSEXP, SEXP sigma_sqSEXP, SEXP m_qSEXP, SEXP l_qSEXP, SEXP x_tildeSEXP, SEXP potential_intervals_qSEXP, SEXP potential_intervals_dims_qSEXP, SEXP l_possible_qSEXP, SEXP phi_l_qSEXP, SEXP l_values_length_qSEXP, SEXP QSEXP, SEXP KSEXP, SEXP gSEXP, SEXP sum_KSEXP, SEXP lambda_id0SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type count(countSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_tilde(b_tildeSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_sq(sigma_sqSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type m_q(m_qSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type l_q(l_qSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_tilde(x_tildeSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type potential_intervals_q(potential_intervals_qSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type potential_intervals_dims_q(potential_intervals_dims_qSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type l_possible_q(l_possible_qSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type phi_l_q(phi_l_qSEXP);
    Rcpp::traits::input_parameter< int >::type l_values_length_q(l_values_length_qSEXP);
    Rcpp::traits::input_parameter< int >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type sum_K(sum_KSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type lambda_id0(lambda_id0SEXP);
    update_lqk(count, k, y, b_tilde, sigma_sq, m_q, l_q, x_tilde, potential_intervals_q, potential_intervals_dims_q, l_possible_q, phi_l_q, l_values_length_q, Q, K, g, sum_K, lambda_id0);
    return R_NilValue;
END_RCPP
}
// update_b_tilde
void update_b_tilde(arma::vec& y, double sigma_sq, arma::mat& x_tilde, arma::mat& Sigma_b_tilde_inv, double tol, arma::vec& b_tilde);
RcppExport SEXP _bliss_update_b_tilde(SEXP ySEXP, SEXP sigma_sqSEXP, SEXP x_tildeSEXP, SEXP Sigma_b_tilde_invSEXP, SEXP tolSEXP, SEXP b_tildeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type sigma_sq(sigma_sqSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x_tilde(x_tildeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Sigma_b_tilde_inv(Sigma_b_tilde_invSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_tilde(b_tildeSEXP);
    update_b_tilde(y, sigma_sq, x_tilde, Sigma_b_tilde_inv, tol, b_tilde);
    return R_NilValue;
END_RCPP
}
// loss_cpp
double loss_cpp(arma::vec& d, arma::vec& grid, arma::vec& posterior_expe);
RcppExport SEXP _bliss_loss_cpp(SEXP dSEXP, SEXP gridSEXP, SEXP posterior_expeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type posterior_expe(posterior_expeSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_cpp(d, grid, posterior_expe));
    return rcpp_result_gen;
END_RCPP
}
// Bliss_Gibbs_Sampler_cpp
List Bliss_Gibbs_Sampler_cpp(int Q, arma::vec& y, List& x, List& grids, int iter, arma::vec& K, CharacterVector& basis, double g, double lambda, arma::mat& V_tilde, arma::vec& l_values_length, List& probs_l, bool progress, double tol);
RcppExport SEXP _bliss_Bliss_Gibbs_Sampler_cpp(SEXP QSEXP, SEXP ySEXP, SEXP xSEXP, SEXP gridsSEXP, SEXP iterSEXP, SEXP KSEXP, SEXP basisSEXP, SEXP gSEXP, SEXP lambdaSEXP, SEXP V_tildeSEXP, SEXP l_values_lengthSEXP, SEXP probs_lSEXP, SEXP progressSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< List& >::type x(xSEXP);
    Rcpp::traits::input_parameter< List& >::type grids(gridsSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type K(KSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< double >::type g(gSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type V_tilde(V_tildeSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type l_values_length(l_values_lengthSEXP);
    Rcpp::traits::input_parameter< List& >::type probs_l(probs_lSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(Bliss_Gibbs_Sampler_cpp(Q, y, x, grids, iter, K, basis, g, lambda, V_tilde, l_values_length, probs_l, progress, tol));
    return rcpp_result_gen;
END_RCPP
}
// Bliss_Simulated_Annealing_cpp
List Bliss_Simulated_Annealing_cpp(int iter, arma::mat& beta_sample, arma::vec& grid, int burnin, double Temp, int k_max, int p_l, int dm, int dl, int p, std::string basis, arma::mat& normalization_values, bool progress, arma::mat& starting_point);
RcppExport SEXP _bliss_Bliss_Simulated_Annealing_cpp(SEXP iterSEXP, SEXP beta_sampleSEXP, SEXP gridSEXP, SEXP burninSEXP, SEXP TempSEXP, SEXP k_maxSEXP, SEXP p_lSEXP, SEXP dmSEXP, SEXP dlSEXP, SEXP pSEXP, SEXP basisSEXP, SEXP normalization_valuesSEXP, SEXP progressSEXP, SEXP starting_pointSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type beta_sample(beta_sampleSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< double >::type Temp(TempSEXP);
    Rcpp::traits::input_parameter< int >::type k_max(k_maxSEXP);
    Rcpp::traits::input_parameter< int >::type p_l(p_lSEXP);
    Rcpp::traits::input_parameter< int >::type dm(dmSEXP);
    Rcpp::traits::input_parameter< int >::type dl(dlSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< std::string >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type normalization_values(normalization_valuesSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type starting_point(starting_pointSEXP);
    rcpp_result_gen = Rcpp::wrap(Bliss_Simulated_Annealing_cpp(iter, beta_sample, grid, burnin, Temp, k_max, p_l, dm, dl, p, basis, normalization_values, progress, starting_point));
    return rcpp_result_gen;
END_RCPP
}
// dposterior_cpp
arma::mat dposterior_cpp(arma::mat& rposterior, arma::vec& y, unsigned N, arma::vec& K, List& potential_intervals, List& potential_intervals_dims, arma::vec& p_l, unsigned Q);
RcppExport SEXP _bliss_dposterior_cpp(SEXP rposteriorSEXP, SEXP ySEXP, SEXP NSEXP, SEXP KSEXP, SEXP potential_intervalsSEXP, SEXP potential_intervals_dimsSEXP, SEXP p_lSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type rposterior(rposteriorSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< unsigned >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type K(KSEXP);
    Rcpp::traits::input_parameter< List& >::type potential_intervals(potential_intervalsSEXP);
    Rcpp::traits::input_parameter< List& >::type potential_intervals_dims(potential_intervals_dimsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type p_l(p_lSEXP);
    Rcpp::traits::input_parameter< unsigned >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(dposterior_cpp(rposterior, y, N, K, potential_intervals, potential_intervals_dims, p_l, Q));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bliss_ginv_cpp", (DL_FUNC) &_bliss_ginv_cpp, 2},
    {"_bliss_mvrnormArma", (DL_FUNC) &_bliss_mvrnormArma, 3},
    {"_bliss_integrate_trapeze_cpp", (DL_FUNC) &_bliss_integrate_trapeze_cpp, 2},
    {"_bliss_uniform_cpp", (DL_FUNC) &_bliss_uniform_cpp, 3},
    {"_bliss_triangular_cpp", (DL_FUNC) &_bliss_triangular_cpp, 3},
    {"_bliss_gaussian_cpp", (DL_FUNC) &_bliss_gaussian_cpp, 3},
    {"_bliss_Epanechnikov_cpp", (DL_FUNC) &_bliss_Epanechnikov_cpp, 3},
    {"_bliss_compute_beta_cpp", (DL_FUNC) &_bliss_compute_beta_cpp, 8},
    {"_bliss_compute_beta_sample_cpp", (DL_FUNC) &_bliss_compute_beta_sample_cpp, 6},
    {"_bliss_potential_intervals_List", (DL_FUNC) &_bliss_potential_intervals_List, 5},
    {"_bliss_moving_average_cpp", (DL_FUNC) &_bliss_moving_average_cpp, 2},
    {"_bliss_potential_intervals_extract", (DL_FUNC) &_bliss_potential_intervals_extract, 4},
    {"_bliss_update_mqk", (DL_FUNC) &_bliss_update_mqk, 17},
    {"_bliss_update_lqk", (DL_FUNC) &_bliss_update_lqk, 18},
    {"_bliss_update_b_tilde", (DL_FUNC) &_bliss_update_b_tilde, 6},
    {"_bliss_loss_cpp", (DL_FUNC) &_bliss_loss_cpp, 3},
    {"_bliss_Bliss_Gibbs_Sampler_cpp", (DL_FUNC) &_bliss_Bliss_Gibbs_Sampler_cpp, 14},
    {"_bliss_Bliss_Simulated_Annealing_cpp", (DL_FUNC) &_bliss_Bliss_Simulated_Annealing_cpp, 14},
    {"_bliss_dposterior_cpp", (DL_FUNC) &_bliss_dposterior_cpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_bliss(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
