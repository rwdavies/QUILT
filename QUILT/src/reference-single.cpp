// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

#include <sys/time.h>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <bitset>


Rcpp::IntegerVector rcpp_int_expand(arma::ivec& hapc, const int nSNPs);

arma::imat inflate_fhb(
    arma::imat& rhb,
    Rcpp::IntegerVector& haps_to_get,
    const int nSNPs
);


double print_times(
    double prev,
    int suppressOutput,
    std::string past_text,
    std::string next_text
);


int rcpp_simple_binary_matrix_search(
    int val,
    Rcpp::IntegerMatrix mat,
    int s1,
    int e1
);


//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector Rcpp_test(int k, int iGrid, arma::imat& rhb_t) {
  //
  Rcpp::IntegerVector output(32);
  output.fill(-1);
  std::uint32_t tmp(rhb_t(k, iGrid));
  for(int b = 0; b < 32; b++) {
    if (tmp & (1<<b)) {
      output(b) = 1;
    } else {
      output(b) = 0;
    }
  }
  return(output);
}

//' @export
// [[Rcpp::export]]
void Rcpp_make_gl_bound(arma::mat & gl, double minGLValue, Rcpp::IntegerVector& to_fix) {
    //
    double a, b;
    int n = to_fix.length();
    int i_col;
    for(int i = 0; i < n; i++) {
        i_col = to_fix(i);
        a = gl(0, i_col);
        b = gl(1, i_col);
        if (a > b) {
            b = b / a;
            a = 1;
            if (b < minGLValue) {
                b = minGLValue;
            }
        } else {
            a = a / b;
            b = 1;
            if (a < minGLValue) {
                a = minGLValue;
            }
        }
        gl(0, i_col) = a;
        gl(1, i_col) = b;
    }
    return;
}



//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_nth_partial_sort(
    Rcpp::NumericVector x,
    int nth
) {
    Rcpp::NumericVector y = clone(x);
    std::nth_element(y.begin(), y.begin()+nth, y.end());
    std::sort(y.begin(), y.begin()+nth);
    return(y);
}




// //' @export
// // [[Rcpp::export]]
// Rcpp::NumericVector Rcpp_get_top_K_or_more_matches(
//     Rcpp::NumericVector& x,
//     int K_top_matches
// ) {
//     val = -rcpp_nth_partial_sort(-x, as.integer(K_top_matches))[K_top_matches]
//     y <- which(x >= val) ## surprisingly slow
//     y <- y[order(-x[y])]
//     return(y)
// }



//' @export
// [[Rcpp::export]]
Rcpp::List Rcpp_get_top_K_or_more_matches_while_building_gamma(
    arma::mat& alphaHat_t,
    arma::colvec& betaHat_t_col,
    arma::colvec& gamma_t_col,
    int iGrid,
    int K,
    int K_top_matches,
    double special_multiplication_value = 1
) {
    // first pass, do gamma, and calculate minimal value for the get
    // second pass, get those values and store them
    //
    // note, possible that betaHat_t_col is off by a factor of not_jump_prob
    // this will uniformly scale the values. it won't affect order or rank
    //
    Rcpp::NumericVector top_K_values(K_top_matches); // reverse order i.e. first is smallest, last is highest
    top_K_values.fill(0);
    int beats_value = 0;
    int num_top_K = 1; //1-based for size
    for(int k = 0; k < K; k++) {
        gamma_t_col(k) = alphaHat_t(k, iGrid) * betaHat_t_col(k);
        if (gamma_t_col(k) == top_K_values(0)) {
            num_top_K++;
        } else if (gamma_t_col(k) > top_K_values(0)) {
            num_top_K++;
            // start storing procedure
            beats_value = 0;
            for(int j = 0; j < K_top_matches; j++) {
                if (gamma_t_col(k) > top_K_values(j)) {
                    beats_value = j;
                }
            }
            if (beats_value > 0) {
                // push back
                for(int i = 0; i < beats_value; i++) {
                    top_K_values(i) = top_K_values(i + 1);
                }
            }
            top_K_values(beats_value) = gamma_t_col(k);
        }
    }
    // second pass, get those that meet the threshold
    Rcpp::IntegerVector top_matches(num_top_K); // num_top_k is upper bound
    Rcpp::NumericVector top_matches_values(num_top_K);
    int count = -1;
    for(int k = 0; k < K; k++) {
        if (gamma_t_col(k) >= top_K_values(0)) {
            count++;
            top_matches(count) = k;
            top_matches_values(count) = gamma_t_col(k) * special_multiplication_value;
        }
    }
    if (count == -1) {
        std::cout << "problem with Rcpp_get_top_K_or_more_matches_while_building_gamma" << std::endl;
        std::cout << "iGrid = " << iGrid << std::endl;
        std::cout << "gamma_t_col(0) = " << gamma_t_col(0) << std::endl;
        std::cout << "alphaHat_t(0, iGrid) = " << alphaHat_t(0, iGrid) << std::endl;
        std::cout << "betaHat_t_col(0) = " << betaHat_t_col(0) << std::endl;
    }
    // this is good enough, save the rest for R next time, these should be much much smaller in all but rarest casests
    Rcpp::List to_return = Rcpp::List::create(
        Rcpp::Named("top_matches") = top_matches[Rcpp::Range(0, count)],
        Rcpp::Named("top_matches_values") = top_matches_values[Rcpp::Range(0, count)]
    );
    return to_return;
}




//' @export
// [[Rcpp::export]]
Rcpp::List Rcpp_get_top_K_or_more_matches_while_building_gamma_eigen(
    Eigen::Map<Eigen::MatrixXd> eigen_alphaHat_t,
    arma::colvec& betaHat_t_col,
    arma::colvec& gamma_t_col,
    int iGrid,
    int K,
    int K_top_matches,
    double special_multiplication_value = 1
) {
    // first pass, do gamma, and calculate minimal value for the get
    // second pass, get those values and store them
    //
    // note, possible that betaHat_t_col is off by a factor of not_jump_prob
    // this will uniformly scale the values. it won't affect order or rank
    //
    Rcpp::NumericVector top_K_values(K_top_matches); // reverse order i.e. first is smallest, last is highest
    top_K_values.fill(0);
    int beats_value = 0;
    int num_top_K = 1; //1-based for size
    for(int k = 0; k < K; k++) {
        gamma_t_col(k) = eigen_alphaHat_t(k, iGrid) * betaHat_t_col(k);
        if (gamma_t_col(k) == top_K_values(0)) {
            num_top_K++;
        } else if (gamma_t_col(k) > top_K_values(0)) {
            num_top_K++;
            // start storing procedure
            beats_value = 0;
            for(int j = 0; j < K_top_matches; j++) {
                if (gamma_t_col(k) > top_K_values(j)) {
                    beats_value = j;
                }
            }
            if (beats_value > 0) {
                // push back
                for(int i = 0; i < beats_value; i++) {
                    top_K_values(i) = top_K_values(i + 1);
                }
            }
            top_K_values(beats_value) = gamma_t_col(k);
        }
    }
    // second pass, get those that meet the threshold
    Rcpp::IntegerVector top_matches(num_top_K); // num_top_k is upper bound
    Rcpp::NumericVector top_matches_values(num_top_K);
    int count = -1;
    for(int k = 0; k < K; k++) {
        if (gamma_t_col(k) >= top_K_values(0)) {
            count++;
            top_matches(count) = k;
            top_matches_values(count) = gamma_t_col(k) * special_multiplication_value;
        }
    }
    if (count == -1) {
        std::cout << "problem with Rcpp_get_top_K_or_more_matches_while_building_gamma" << std::endl;
        std::cout << "iGrid = " << iGrid << std::endl;
        std::cout << "gamma_t_col(0) = " << gamma_t_col(0) << std::endl;
        std::cout << "eigen_alphaHat_t(0, iGrid) = " << eigen_alphaHat_t(0, iGrid) << std::endl;
        std::cout << "betaHat_t_col(0) = " << betaHat_t_col(0) << std::endl;
    }
    // this is good enough, save the rest for R next time, these should be much much smaller in all but rarest casests
    Rcpp::List to_return = Rcpp::List::create(
        Rcpp::Named("top_matches") = top_matches[Rcpp::Range(0, count)],
        Rcpp::Named("top_matches_values") = top_matches_values[Rcpp::Range(0, count)]
    );
    return to_return;
}



//' @export
// [[Rcpp::export]]
arma::mat Rcpp_build_eMatDH(
    arma::imat& distinctHapsB,
    const arma::mat& gl,
    const int nGrids,
    const int nSNPs,
    const double ref_error,
    const double ref_one_minus_error,
    const bool add_zero_row = false
) {
    const int nMaxDH = distinctHapsB.n_rows;
    arma::mat eMatDH;
    int kbump;
    if (add_zero_row) {
        eMatDH = arma::mat(nMaxDH + 1, nGrids);
        eMatDH.row(0).fill(1);
        kbump = 1;
    } else {
        eMatDH = arma::mat(nMaxDH, nGrids);
        kbump = 0;
    }
    arma::mat gl_local(2, 32);
    //
    for(int iGrid = 0; iGrid < nGrids; iGrid++) {
        //
        int s = 32 * iGrid; // 0-based here
        int e = 32 * (iGrid + 1) - 1;
        if (e > (nSNPs - 1)) {
            e = nSNPs - 1;
        }
        int nSNPsLocal = e - s + 1;
        // now subset
        for(int i = 0; i < nSNPsLocal; i++) {
            gl_local(0, i) = gl(0, i + s);
            gl_local(1, i) = gl(1, i + s);      
        }
        for(int k = 0; k < nMaxDH; k++) {
            //
            std::uint32_t tmp(distinctHapsB(k, iGrid));
            //
            double prob = 1;
            for(int b = 0; b < nSNPsLocal; b++) {
                double dR = gl_local(0, b);
                double dA = gl_local(1, b);
                if (tmp & (1<<b)) {
                    // alternate
                    prob *= (dR * ref_error + dA * ref_one_minus_error);
                } else {
                    prob *= (dR * ref_one_minus_error + dA * ref_error);
                }
            }
            eMatDH(kbump + k, iGrid) = prob;
        }
        if (add_zero_row) {
            eMatDH(0, iGrid) = arma::min(eMatDH.col(iGrid));
        }
    }
    return(eMatDH);
}


//' @export
// [[Rcpp::export]]
void rcpp_internal_make_eMatRead_t_using_binary(
    arma::mat & eMatRead_t,
    arma::imat& rhb,
    const int K,
    const int nSNPs,
    const Rcpp::NumericVector& u,
    const Rcpp::NumericMatrix& ps,
    const int nReads,
    const Rcpp::NumericVector& start,
    const Rcpp::NumericVector& end,
    const Rcpp::NumericVector& nr,
    const double ref_error,
    const int ceil_K_n,
    const int n
) {
    double ref_one_minus_error = 1 - ref_error;
    //
    arma::imat hap;
    int iRead, s, j;
    int e = 0;
    int KL, ik, sh, eh, i;
    double prob, pR, pA;
    Rcpp::IntegerVector to_pass_in(1);
    arma::imat haps;
    for(i = 0; i < ceil_K_n; i++) {
        sh = n * i; 
        eh = n * (i + 1) - 1;
        if (eh > (K - 1)) {
            eh = K - 1;
        }
        KL = eh - sh + 1; // 1-based, length
        Rcpp::IntegerVector haps_to_get(KL);
        for(ik = 0; ik < KL; ik++) {
            haps_to_get(ik) = ik + sh;
        }
        haps = inflate_fhb(rhb, haps_to_get, nSNPs);
        //
        for(ik = 0; ik < KL; ik++) {
            for(iRead = 0; iRead < nReads; iRead++) {
                s = start(iRead) - 1;
                e = end(iRead) - 1;
                prob = 1;
                for(j = 0; j < nr(iRead); j++) {
                    pR = ps(0, s + j);
                    pA = ps(1, s + j);
                    if (haps(u(s + j), ik) == 0) {
                        prob *= (pR * ref_one_minus_error + pA * ref_error);
                    } else {
                        prob *= (pA * ref_one_minus_error + pR * ref_error);
                    }
                }
                eMatRead_t((sh + ik), iRead) = prob;
            }
        }
    }
    return;
}



//' @export
// [[Rcpp::export]]
void Rcpp_haploid_reference_single_forward(
    Rcpp::IntegerVector& gammaSmall_cols_to_get,
    const arma::mat& gl,
    arma::mat& alphaHat_t,
    arma::rowvec& c,
    const arma::mat& transMatRate_t,
    const arma::imat& rhb_t,
    arma::imat& hapMatcher,
    arma::mat& eMatDH,
    const int& nGrids,
    const int& nSNPs,
    const int& K,
    const bool& use_eMatDH,
    double ref_error,
    const bool only_store_alpha_at_gamma_small,
    bool always_normalize,
    double min_emission_prob_normalization_threshold,
    const double maxEmissionMatrixDifference = 1e10,
    const bool normalize_emissions = false 
) {
    double jump_prob, not_jump_prob; //
    double one_minus_jump_prob = 1;
    int s, e, nSNPsLocal, i, iGrid, k;
    double ref_one_minus_error = 1 - ref_error;
    //const int n_which_hapMatcher_0 = which_hapMatcher_0.n_rows;
    //int hapMatcher_0_position = 0;
    //bool continue_hapMatcher_0 = true;
    double dR, dA, prob;
    double run_total = 0;
    arma::mat gl_local(2, 32);
    arma::icolvec dh_col(K);
    arma::colvec prob_col(K);
    arma::colvec eMatDH_col(K);
    arma::colvec alphaHat_t_col_prev(K);
    arma::colvec alphaHat_t_col(K);
    bool grid_has_variant;
    const double double_K = double(K);
    bool store_alpha_for_this_grid;
    double running_min_emission_prob = 1;
    double min_emission_prob = 1;
    double prev_alphaHat_t_col_sum = 1;
    double jump_prob_plus = 1;
    // double one_over_maxEmissionMatrixDifference = 1 / maxEmissionMatrixDifference;
    double emission_max = 1;
    //
    for(iGrid = 1; iGrid < nGrids; iGrid++) {
        jump_prob = transMatRate_t(1, iGrid - 1) / double_K;
        if (always_normalize) {
            jump_prob_plus = jump_prob;
        } else {
            jump_prob_plus = jump_prob * prev_alphaHat_t_col_sum;
            //std::cout << "prev_alphaHat_t_col_sum = " <<prev_alphaHat_t_col_sum << ", col sum version = " << arma::sum(alphaHat_t.col(iGrid - 1)) << std::endl;
        }
        not_jump_prob = transMatRate_t(0, iGrid - 1);
        one_minus_jump_prob = (1 - jump_prob);
        s = 32 * iGrid; // 0-based here
        e = 32 * (iGrid + 1) - 1;
        if (e > (nSNPs - 1)) {
            e = nSNPs - 1;
        }
        nSNPsLocal = e - s + 1;
        grid_has_variant = false;        
        for(i = 0; i < nSNPsLocal; i++) {
            gl_local(0, i) = gl(0, i + s);
            gl_local(1, i) = gl(1, i + s);
            if ((gl_local(0, i) != 1) | (gl_local(1, i) != 1)) {
                grid_has_variant = true;
            }
        }
        if (iGrid == 1) {
            grid_has_variant = true; // make sure re-scaling kicks on properly
	    alphaHat_t_col = alphaHat_t.col(0); // initialize with previous value
        }
	store_alpha_for_this_grid = true;
	if (only_store_alpha_at_gamma_small) {
            if (gammaSmall_cols_to_get(iGrid) < 0) {
	      store_alpha_for_this_grid = false;
	    }
	}
        if (use_eMatDH) {
            //
            // yes use eMatDH
            //
            if (grid_has_variant) {
                dh_col = hapMatcher.col(iGrid);
                eMatDH_col = eMatDH.col(iGrid);
                //
                // this chunk should be a function, but meh
                //
                if (normalize_emissions) {
                    emission_max = arma::max(eMatDH_col);
                    if (emission_max < 1) {
                        eMatDH_col *= (1 / emission_max);
                    }
                }
                //
                min_emission_prob = arma::min(eMatDH_col);
		run_total = 0;
                for(k = 0; k < K; k++) {
		    // if dh_col is 0 i.e. need to re-do this prob is 0 so we are OK
                    prob = eMatDH_col(dh_col(k));
                    if (dh_col(k) == 0) {
                        //
                        std::uint32_t tmp(rhb_t(k, iGrid));
                        //
                        prob = 1;
                        for(int b = 0; b < nSNPsLocal; b++) {
                            dR = gl_local(0, b);
                            dA = gl_local(1, b);
                            if (tmp & (1<<b)) {
                                // alternate
                                prob *= (dR * ref_error + dA * ref_one_minus_error);
                            } else {
                                prob *= (dR * ref_one_minus_error + dA * ref_error);
                            }
                        }
                        prob *= (1 / emission_max);
                        if (prob < min_emission_prob) {
                            min_emission_prob = prob;
                        }
                    } // NOTE - can I not multiply not_jump_prob, and do it at the end?
		    alphaHat_t_col(k) = (jump_prob_plus + not_jump_prob * alphaHat_t_col(k)) * prob;
                    run_total += alphaHat_t_col(k);
                }
                //
                if (always_normalize) {
                    c(iGrid) = 1 / run_total;
                    alphaHat_t_col *= c(iGrid);
                } else {
                    running_min_emission_prob *= min_emission_prob;
                    if (
                        (running_min_emission_prob < min_emission_prob_normalization_threshold) |
                        (iGrid == (nGrids - 1))
                        ) {
                        running_min_emission_prob = 1;
                        c(iGrid) = 1 / run_total;
                        alphaHat_t_col *= c(iGrid);
                        run_total = 1;
                    } else {
                        c(iGrid) = 1;
                    }
                }
            } else {
	        alphaHat_t_col = (jump_prob_plus + not_jump_prob * alphaHat_t_col);
	        c(iGrid) = 1;
                if (always_normalize) {
                    run_total = 1;
                } else {
                    run_total = K * jump_prob_plus + prev_alphaHat_t_col_sum * not_jump_prob;
                }
                if (always_normalize | (iGrid == (nGrids - 1))) {
                    c(iGrid) = 1 / run_total;
                    alphaHat_t_col *= c(iGrid);
                    run_total = 1;
                }
            }
            prev_alphaHat_t_col_sum = run_total;
	    if (store_alpha_for_this_grid) {
	        alphaHat_t.col(iGrid) = alphaHat_t_col;
	    }
        } else {
            //
            // no to eMatDH
            //
            for(k = 0; k < K; k++) {            
                std::uint32_t tmp(rhb_t(k, iGrid));                
                //
                prob = 1;
                for(int b = 0; b < nSNPsLocal; b++) {
                    dR = gl_local(0, b);
                    dA = gl_local(1, b);
                    if (tmp & (1<<b)) {
                        // alternate
                        prob *= (dR * ref_error + dA * ref_one_minus_error);
                    } else {
                        prob *= (dR * ref_one_minus_error + dA * ref_error);
                    }
                }
                alphaHat_t(k, iGrid) = (jump_prob + not_jump_prob * alphaHat_t(k, iGrid - 1)) * prob;
            }
            c(iGrid) = 1 / sum(alphaHat_t.col(iGrid));
            alphaHat_t.col(iGrid) = alphaHat_t.col(iGrid) * c(iGrid);
        }
        //
    }
    return;
}



//' @export
// [[Rcpp::export]]
void Rcpp_haploid_reference_single_forward_version2(
    Rcpp::IntegerVector& gammaSmall_cols_to_get,
    const arma::mat& gl,
    arma::mat& alphaHat_t,
    arma::rowvec& c,
    const arma::mat& transMatRate_t,
    const arma::imat& rhb_t,
    Rcpp::IntegerMatrix& eMatDH_special_matrix_helper,
    Rcpp::IntegerMatrix& eMatDH_special_matrix,
    const bool use_eMatDH_special_symbols,    
    arma::imat& hapMatcher,
    Rcpp::RawMatrix& hapMatcherR,
    bool use_hapMatcherR,
    arma::mat& eMatDH,
    const int& nGrids,
    const int& nSNPs,
    const int& K,
    const bool& use_eMatDH,
    double ref_error,
    const bool only_store_alpha_at_gamma_small,
    bool always_normalize,
    double min_emission_prob_normalization_threshold,
    const Rcpp::IntegerVector& eMatDH_special_grid_which,
    const Rcpp::List& eMatDH_special_values_list,
    const double maxEmissionMatrixDifference = 1e10,
    const bool normalize_emissions = false
) {
    double jump_prob, not_jump_prob; //
    double one_minus_jump_prob = 1;
    int s, e, nSNPsLocal, i, iGrid, k;
    int bvtd, s1, e1;
    double ref_one_minus_error = 1 - ref_error;
    //const int n_which_hapMatcher_0 = which_hapMatcher_0.n_rows;
    //int hapMatcher_0_position = 0;
    //bool continue_hapMatcher_0 = true;
    double dR, dA, prob;
    double run_total = 0;
    arma::mat gl_local(2, 32);
    arma::icolvec dh_col(K);
    Rcpp::RawVector dh_colR(K);    
    arma::colvec prob_col(K);
    arma::colvec eMatDH_col(K);
    arma::colvec temp_alphaHat_t_col(K);
    arma::colvec alphaHat_t_col(K);
    bool grid_has_variant;
    const double double_K = double(K);
    bool store_alpha_for_this_grid;
    double running_min_emission_prob = 1;
    double min_emission_prob = 1;
    double prev_alphaHat_t_col_sum = 1;
    double jump_prob_plus = 1;
    double jump_prob_plus_divided_by_not_jump = 1;
    double x = 0; // misc
    // double one_over_maxEmissionMatrixDifference = 1 / maxEmissionMatrixDifference;
    double emission_max = 1;
    Rcpp::IntegerVector vals_to_redo;    
    //
    for(iGrid = 1; iGrid < nGrids; iGrid++) {
        c(iGrid) = 1; // set, modify on the fly later
        jump_prob = transMatRate_t(1, iGrid - 1) / double_K;
        if (always_normalize) {
            jump_prob_plus = jump_prob;
        } else {
            jump_prob_plus = jump_prob * prev_alphaHat_t_col_sum;
            //std::cout << "prev_alphaHat_t_col_sum = " <<prev_alphaHat_t_col_sum << ", col sum version = " << arma::sum(alphaHat_t.col(iGrid - 1)) << std::endl;
        }
        not_jump_prob = transMatRate_t(0, iGrid - 1);
        jump_prob_plus_divided_by_not_jump = jump_prob_plus / not_jump_prob;
        one_minus_jump_prob = (1 - jump_prob);
        s = 32 * iGrid; // 0-based here
        e = 32 * (iGrid + 1) - 1;
        if (e > (nSNPs - 1)) {
            e = nSNPs - 1;
        }
        nSNPsLocal = e - s + 1;
        grid_has_variant = false;        
        for(i = 0; i < nSNPsLocal; i++) {
            gl_local(0, i) = gl(0, i + s);
            gl_local(1, i) = gl(1, i + s);
            if ((gl_local(0, i) != 1) | (gl_local(1, i) != 1)) {
                grid_has_variant = true;
            }
        }
        if (iGrid == 1) {
            grid_has_variant = true; // make sure re-scaling kicks on properly
	    alphaHat_t_col = alphaHat_t.col(0); // initialize with previous value
        }
	store_alpha_for_this_grid = true;
	if (only_store_alpha_at_gamma_small) {
            if (gammaSmall_cols_to_get(iGrid) < 0) {
	      store_alpha_for_this_grid = false;
	    }
	}
        if (use_eMatDH) {
            //
            // yes use eMatDH
            //
            if (grid_has_variant) {
                if (use_hapMatcherR) {
                    dh_colR = hapMatcherR(Rcpp::_, iGrid);
                } else {
                    dh_col = hapMatcher.col(iGrid);
                }
                eMatDH_col = eMatDH.col(iGrid);
                //
                if (normalize_emissions) {
                    emission_max = arma::max(eMatDH_col);
                    if (emission_max < 1) {
                        eMatDH_col *= (1 / emission_max);
                    }
                }
                min_emission_prob = arma::min(eMatDH_col);
                // no longer necessary, part of eMatDH
                //eMatDH_col(0) = min_emission_prob;                
                //
                // this setting of this to be the max here should protect against situations where many of the haplotypes are poor matches, e.g. if all of them are 1x-20, then multiplying here by 1, and then re-inserting later, would lead to underflow problems. unlikely to be a problem for dense reference panels
                //
                // double max_emission_prob = arma::max(eMatDH_col);
		run_total = 0;
                //
                // first, do special cases
                //
                if (eMatDH_special_grid_which(iGrid) > 0) {
                    if (use_eMatDH_special_symbols) {
                        s1 = eMatDH_special_matrix_helper(iGrid, 0);
                        e1 = eMatDH_special_matrix_helper(iGrid, 1);
                        // there is undoubtledly better code to do this
                        Rcpp::IntegerVector temp_vec(e1 - s1 + 1);
                        for(int ii = 0; ii < e1 - s1 + 1; ii++) {
                            temp_vec(ii) = eMatDH_special_matrix(s1 - 1 + ii, 0);
                        }
                        vals_to_redo = temp_vec;
                    } else {
                        vals_to_redo = Rcpp::as<Rcpp::IntegerVector>(eMatDH_special_values_list(eMatDH_special_grid_which(iGrid) - 1));
                    }
                    for(i = 0; i < vals_to_redo.size(); i++) {
                        k = vals_to_redo(i);
                        if (use_eMatDH_special_symbols) {
                            bvtd = rcpp_simple_binary_matrix_search(k, eMatDH_special_matrix, s1, e1);
                        } else {
                            bvtd = rhb_t(k, iGrid);
                        }
                        std::uint32_t tmp(bvtd);
                        prob = 1;
                        for(int b = 0; b < nSNPsLocal; b++) {
                            dR = gl_local(0, b);
                            dA = gl_local(1, b);
                            if (tmp & (1<<b)) {
                                // alternate
                                prob *= (dR * ref_error + dA * ref_one_minus_error);
                            } else {
                                prob *= (dR * ref_one_minus_error + dA * ref_error);
                            }
                        }
                        prob *= (1 / emission_max);
                        //store this here for now, re-put in after
                        temp_alphaHat_t_col(k) = (jump_prob_plus_divided_by_not_jump + alphaHat_t_col(k)) * prob;
                        run_total += temp_alphaHat_t_col(k);
                        if (prob < min_emission_prob) {
                            min_emission_prob = prob;
                        }
                    }
                }
                //
                eMatDH_col(0) = 0;
                //
                // next, do normal
                //
                //alphaHat_t_col(k) = (jump_prob_plus + not_jump_prob * alphaHat_t_col(k)) * eMatDH_col(dh_col(k));
                // if dh_col is 0 i.e. need to re-do this prob is 1 so we are OK                
                if (use_hapMatcherR) {
                    for(k = 0; k < K; k++) {
                        alphaHat_t_col(k) = (jump_prob_plus_divided_by_not_jump + alphaHat_t_col(k)) * eMatDH_col(dh_colR(k));
                        run_total += alphaHat_t_col(k);
                    }
                } else {
                    for(k = 0; k < K; k++) {
                        alphaHat_t_col(k) = (jump_prob_plus_divided_by_not_jump + alphaHat_t_col(k)) * eMatDH_col(dh_col(k));
                        run_total += alphaHat_t_col(k);
                    }
                }
                //
                // put back in
                //
                if (eMatDH_special_grid_which(iGrid) > 0) {
                    for(i = 0; i < vals_to_redo.size(); i++) {
                        k = vals_to_redo(i);
                       // remove influence from run_total
                        run_total -= alphaHat_t_col(k);
                        alphaHat_t_col(k) = temp_alphaHat_t_col(k);
                    }
                }
                //
                running_min_emission_prob *= min_emission_prob;
                //
                // now, the normal way, have new alphaHat_t_col, up to not_jump_prob
                //    and new run_total, accuracy on alphaHat_t_col
                //
            } else {
                //
                // here jump_prob_plus_divided_by_not_jump is A_{g-1} * \psi^{g-1} / \eta^{g-1}
                // 
	        alphaHat_t_col = (jump_prob_plus_divided_by_not_jump + alphaHat_t_col);
                //
                // then the sum is A_{g-1} / ega^{g-1}
                //
                run_total = prev_alphaHat_t_col_sum / not_jump_prob;                
                //run_total = K * jump_prob_plus_divided_by_not_jump + prev_alphaHat_t_col_sum;
            }
            //
            // capture the not_jump_prob
            //
            c(iGrid) /= not_jump_prob;
            //
            // then, possibly normalize
            //
            if (
                always_normalize |
                (running_min_emission_prob < min_emission_prob_normalization_threshold) |
                (iGrid == (nGrids - 1))
            ) {
                // normalize so that alphaHat_t_col has sum 1
                x = 1 / (run_total);
                alphaHat_t_col *= x;
                c(iGrid) /= run_total;
                run_total=1;
                running_min_emission_prob = 1;
            }
            prev_alphaHat_t_col_sum = run_total;
	    if (store_alpha_for_this_grid) {
	        alphaHat_t.col(iGrid) = alphaHat_t_col;
	    }
            //     x = arma::sum(alphaHat_t_col) - run_total;
            //     if (x > 1e-6 | x < (-1e-6)) {
            //         std::cout << "arma::sum(alphaHat_t_col) = " << arma::sum(alphaHat_t_col) << ", run_total = " << run_total << std::endl;
            //     }
	    //     //c(iGrid) = 1 / (not_jump_prob);
            //     if (iGrid == (nGrids - 1)) {
            //         c(iGrid) /= run_total;
            //         x = 1 / run_total;
            //         alphaHat_t_col *= x;
            //     }
            // }
            //
            //
        } else {
            //
            // no to eMatDH
            //
            for(k = 0; k < K; k++) {            
                std::uint32_t tmp(rhb_t(k, iGrid));                
                //
                prob = 1;
                for(int b = 0; b < nSNPsLocal; b++) {
                    dR = gl_local(0, b);
                    dA = gl_local(1, b);
                    if (tmp & (1<<b)) {
                        // alternate
                        prob *= (dR * ref_error + dA * ref_one_minus_error);
                    } else {
                        prob *= (dR * ref_one_minus_error + dA * ref_error);
                    }
                }
                alphaHat_t(k, iGrid) = (jump_prob + not_jump_prob * alphaHat_t(k, iGrid - 1)) * prob;
            }
            c(iGrid) = 1 / sum(alphaHat_t.col(iGrid));
            alphaHat_t.col(iGrid) = alphaHat_t.col(iGrid) * c(iGrid);
        }
        //
    }
    return;
}
















//' @export
// [[Rcpp::export]]
void Rcpp_haploid_reference_single_forward_version3(
    Rcpp::IntegerVector& gammaSmall_cols_to_get,
    const arma::mat& gl,
    Eigen::Map<Eigen::MatrixXd> alphaHat_t,
    arma::rowvec& c,
    const arma::mat& transMatRate_t,
    const arma::imat& rhb_t,
    Rcpp::IntegerMatrix& eMatDH_special_matrix_helper,
    Rcpp::IntegerMatrix& eMatDH_special_matrix,
    const bool use_eMatDH_special_symbols,    
    arma::imat& hapMatcher,
    Rcpp::RawMatrix& hapMatcherR,
    bool use_hapMatcherR,
    arma::mat& eMatDH,
    const int& nGrids,
    const int& nSNPs,
    const int& K,
    const bool& use_eMatDH,
    double ref_error,
    const bool only_store_alpha_at_gamma_small,
    bool always_normalize,
    double min_emission_prob_normalization_threshold,
    const Rcpp::IntegerVector& eMatDH_special_grid_which,
    const Rcpp::List& eMatDH_special_values_list,
    const double maxEmissionMatrixDifference = 1e10,
    const bool normalize_emissions = false
) {
    double jump_prob, not_jump_prob; //
    double one_minus_jump_prob = 1;
    int s, e, nSNPsLocal, i, iGrid, k;
    int bvtd, s1, e1;
    double ref_one_minus_error = 1 - ref_error;
    //const int n_which_hapMatcher_0 = which_hapMatcher_0.n_rows;
    //int hapMatcher_0_position = 0;
    //bool continue_hapMatcher_0 = true;
    double dR, dA, prob;
    double run_total = 0;
    arma::mat gl_local(2, 32);
    arma::icolvec dh_col(K);
    Rcpp::RawVector dh_colR(K);    
    arma::colvec prob_col(K);
    arma::colvec eMatDH_col(K);
    arma::colvec temp_alphaHat_t_col(K);
    arma::colvec alphaHat_t_col(K);
    bool grid_has_variant;
    const double double_K = double(K);
    bool store_alpha_for_this_grid;
    double running_min_emission_prob = 1;
    double min_emission_prob = 1;
    double prev_alphaHat_t_col_sum = 1;
    double jump_prob_plus = 1;
    double jump_prob_plus_divided_by_not_jump = 1;
    double x = 0; // misc
    // double one_over_maxEmissionMatrixDifference = 1 / maxEmissionMatrixDifference;
    double emission_max = 1;
    Rcpp::IntegerVector vals_to_redo;    
    //
    for(iGrid = 1; iGrid < nGrids; iGrid++) {
        c(iGrid) = 1; // set, modify on the fly later
        jump_prob = transMatRate_t(1, iGrid - 1) / double_K;
        if (always_normalize) {
            jump_prob_plus = jump_prob;
        } else {
            jump_prob_plus = jump_prob * prev_alphaHat_t_col_sum;
            //std::cout << "prev_alphaHat_t_col_sum = " <<prev_alphaHat_t_col_sum << ", col sum version = " << arma::sum(alphaHat_t.col(iGrid - 1)) << std::endl;
        }
        not_jump_prob = transMatRate_t(0, iGrid - 1);
        jump_prob_plus_divided_by_not_jump = jump_prob_plus / not_jump_prob;
        one_minus_jump_prob = (1 - jump_prob);
        s = 32 * iGrid; // 0-based here
        e = 32 * (iGrid + 1) - 1;
        if (e > (nSNPs - 1)) {
            e = nSNPs - 1;
        }
        nSNPsLocal = e - s + 1;
        grid_has_variant = false;        
        for(i = 0; i < nSNPsLocal; i++) {
            gl_local(0, i) = gl(0, i + s);
            gl_local(1, i) = gl(1, i + s);
            if ((gl_local(0, i) != 1) | (gl_local(1, i) != 1)) {
                grid_has_variant = true;
            }
        }
        if (iGrid == 1) {
            grid_has_variant = true; // make sure re-scaling kicks on properly
            for(k = 0; k < K; k++) {
                alphaHat_t_col(k) = alphaHat_t(k, 0);
            }
        }
	store_alpha_for_this_grid = true;
	if (only_store_alpha_at_gamma_small) {
            if (gammaSmall_cols_to_get(iGrid) < 0) {
	      store_alpha_for_this_grid = false;
	    }
	}
        if (use_eMatDH) {
            //
            // yes use eMatDH
            //
            if (grid_has_variant) {
                if (use_hapMatcherR) {
                    dh_colR = hapMatcherR(Rcpp::_, iGrid);
                } else {
                    dh_col = hapMatcher.col(iGrid);
                }
                eMatDH_col = eMatDH.col(iGrid);
                //
                if (normalize_emissions) {
                    emission_max = arma::max(eMatDH_col);
                    if (emission_max < 1) {
                        eMatDH_col *= (1 / emission_max);
                    }
                }
                min_emission_prob = arma::min(eMatDH_col);
                // no longer necessary, part of eMatDH
                //eMatDH_col(0) = min_emission_prob;                
                //
                // this setting of this to be the max here should protect against situations where many of the haplotypes are poor matches, e.g. if all of them are 1x-20, then multiplying here by 1, and then re-inserting later, would lead to underflow problems. unlikely to be a problem for dense reference panels
                //
                // double max_emission_prob = arma::max(eMatDH_col);
		run_total = 0;
                //
                // first, do special cases
                //
                if (eMatDH_special_grid_which(iGrid) > 0) {
                    if (use_eMatDH_special_symbols) {
                        s1 = eMatDH_special_matrix_helper(iGrid, 0);
                        e1 = eMatDH_special_matrix_helper(iGrid, 1);
                        // there is undoubtledly better code to do this
                        Rcpp::IntegerVector temp_vec(e1 - s1 + 1);
                        for(int ii = 0; ii < e1 - s1 + 1; ii++) {
                            temp_vec(ii) = eMatDH_special_matrix(s1 - 1 + ii, 0);
                        }
                        vals_to_redo = temp_vec;
                    } else {
                        vals_to_redo = Rcpp::as<Rcpp::IntegerVector>(eMatDH_special_values_list(eMatDH_special_grid_which(iGrid) - 1));
                    }
                    for(i = 0; i < vals_to_redo.size(); i++) {
                        k = vals_to_redo(i);
                        if (use_eMatDH_special_symbols) {
                            bvtd = rcpp_simple_binary_matrix_search(k, eMatDH_special_matrix, s1, e1);
                        } else {
                            bvtd = rhb_t(k, iGrid);
                        }
                        std::uint32_t tmp(bvtd);
                        prob = 1;
                        for(int b = 0; b < nSNPsLocal; b++) {
                            dR = gl_local(0, b);
                            dA = gl_local(1, b);
                            if (tmp & (1<<b)) {
                                // alternate
                                prob *= (dR * ref_error + dA * ref_one_minus_error);
                            } else {
                                prob *= (dR * ref_one_minus_error + dA * ref_error);
                            }
                        }
                        prob *= (1 / emission_max);
                        //store this here for now, re-put in after
                        temp_alphaHat_t_col(k) = (jump_prob_plus_divided_by_not_jump + alphaHat_t_col(k)) * prob;
                        run_total += temp_alphaHat_t_col(k);
                        if (prob < min_emission_prob) {
                            min_emission_prob = prob;
                        }
                    }
                }
                //
                eMatDH_col(0) = 0;
                //
                // next, do normal
                //
                //alphaHat_t_col(k) = (jump_prob_plus + not_jump_prob * alphaHat_t_col(k)) * eMatDH_col(dh_col(k));
                // if dh_col is 0 i.e. need to re-do this prob is 1 so we are OK                
                if (use_hapMatcherR) {
                    for(k = 0; k < K; k++) {
                        alphaHat_t_col(k) = (jump_prob_plus_divided_by_not_jump + alphaHat_t_col(k)) * eMatDH_col(dh_colR(k));
                        run_total += alphaHat_t_col(k);
                    }
                } else {
                    for(k = 0; k < K; k++) {
                        alphaHat_t_col(k) = (jump_prob_plus_divided_by_not_jump + alphaHat_t_col(k)) * eMatDH_col(dh_col(k));
                        run_total += alphaHat_t_col(k);
                    }
                }
                //
                // put back in
                //
                if (eMatDH_special_grid_which(iGrid) > 0) {
                    for(i = 0; i < vals_to_redo.size(); i++) {
                        k = vals_to_redo(i);
                       // remove influence from run_total
                        run_total -= alphaHat_t_col(k);
                        alphaHat_t_col(k) = temp_alphaHat_t_col(k);
                    }
                }
                //
                running_min_emission_prob *= min_emission_prob;
                //
                // now, the normal way, have new alphaHat_t_col, up to not_jump_prob
                //    and new run_total, accuracy on alphaHat_t_col
                //
            } else {
                //
                // here jump_prob_plus_divided_by_not_jump is A_{g-1} * \psi^{g-1} / \eta^{g-1}
                // 
	        alphaHat_t_col = (jump_prob_plus_divided_by_not_jump + alphaHat_t_col);
                //
                // then the sum is A_{g-1} / ega^{g-1}
                //
                run_total = prev_alphaHat_t_col_sum / not_jump_prob;                
                //run_total = K * jump_prob_plus_divided_by_not_jump + prev_alphaHat_t_col_sum;
            }
            //
            // capture the not_jump_prob
            //
            c(iGrid) /= not_jump_prob;
            //
            // then, possibly normalize
            //
            if (
                always_normalize |
                (running_min_emission_prob < min_emission_prob_normalization_threshold) |
                (iGrid == (nGrids - 1))
            ) {
                // normalize so that alphaHat_t_col has sum 1
                x = 1 / (run_total);
                alphaHat_t_col *= x;
                c(iGrid) /= run_total;
                run_total=1;
                running_min_emission_prob = 1;
            }
            prev_alphaHat_t_col_sum = run_total;
	    if (store_alpha_for_this_grid) {
                for(k = 0; k < K; k++) {
                    alphaHat_t(k, iGrid) = alphaHat_t_col(k);
                }
	    }
            //     x = arma::sum(alphaHat_t_col) - run_total;
            //     if (x > 1e-6 | x < (-1e-6)) {
            //         std::cout << "arma::sum(alphaHat_t_col) = " << arma::sum(alphaHat_t_col) << ", run_total = " << run_total << std::endl;
            //     }
	    //     //c(iGrid) = 1 / (not_jump_prob);
            //     if (iGrid == (nGrids - 1)) {
            //         c(iGrid) /= run_total;
            //         x = 1 / run_total;
            //         alphaHat_t_col *= x;
            //     }
            // }
            //
            //
        }
        //
    }
    return;
}


















//' @export
// [[Rcpp::export]]
void Rcpp_haploid_reference_single_backward(
    arma::mat& alphaHat_t,
    arma::mat& betaHat_t,
    arma::mat& gamma_t,
    arma::mat& gammaSmall_t,
    Rcpp::List& best_haps_stuff_list,
    Rcpp::IntegerVector& gammaSmall_cols_to_get,
    Rcpp::NumericVector& dosage,    
    const int& nGrids,
    const arma::mat& transMatRate_t,
    arma::mat& eMatDH,
    arma::imat& hapMatcher,
    const int& nSNPs,
    const int& K,
    const bool& use_eMatDH,
    const arma::imat& rhb_t,
    double ref_error,
    const arma::mat& gl,
    arma::rowvec& c,
    arma::mat& distinctHapsIE,
    bool return_betaHat_t,
    bool return_dosage,
    bool return_gamma_t,
    bool return_gammaSmall_t,
    bool get_best_haps_from_thinned_sites,
    const int nMaxDH,
    const int K_top_matches,
    const double maxEmissionMatrixDifference,    
    const bool normalize_emissions
) {
    // 
    //
    double ref_one_minus_error = 1 - ref_error;    
    double jump_prob, not_jump_prob;
    double one_minus_jump_prob = 1;
    int b, s, e, nSNPsLocal, i, iGrid, k, dh;
    const double double_K = double(K);
    arma::mat gl_local(2, 32);
    double dR, dA, prob, val, gk;
    arma::colvec ematcol(K);
    iGrid = nGrids - 1;
    arma::colvec alphaHat_t_col(K);    
    arma::colvec betaHat_t_col(K);
    arma::colvec gamma_t_col(K);
    bool calculate_small_gamma_t_col, grid_has_variant;
    arma::vec matched_gammas(nMaxDH + 1);
    arma::colvec e_times_b(K);
    arma::vec dosageL(32);
    dosageL.fill(0);
    arma::icolvec dh_col(K);
    arma::colvec eMatDH_col(nMaxDH + 1);
    arma::colvec prob_col(K);
    double sum_e_times_b = 0;
    double prev_val = -1; // argh
    double emission_max = 1;
    //
    //
    for(iGrid = nGrids - 1; iGrid >= 0; --iGrid) {
        if (iGrid == (nGrids - 1)) {
            betaHat_t_col.fill(1);
        } else {
            jump_prob = transMatRate_t(1, iGrid) / double_K;
            one_minus_jump_prob = (1 - jump_prob);
            not_jump_prob = transMatRate_t(0, iGrid);
            s = 32 * (iGrid + 1); // 0-based here
            e = 32 * (iGrid + 1 + 1) - 1;
            if (e > (nSNPs - 1)) {
                e = nSNPs - 1;
            }
	    grid_has_variant = false;
            nSNPsLocal = e - s + 1;
            for(i = 0; i < nSNPsLocal; i++) {
                gl_local(0, i) = gl(0, i + s);
                gl_local(1, i) = gl(1, i + s);
		if ((gl_local(0, i) != 1) | (gl_local(1, i) != 1)) {
		  grid_has_variant = true;
		}
            }
	    //
	    // use eMatDH
	    //
	    if (use_eMatDH) {
	        if (grid_has_variant) {
		    dh_col = hapMatcher.col(iGrid + 1);		  
		    eMatDH_col = eMatDH.col(iGrid + 1);
                    if (normalize_emissions) {
                        emission_max = arma::max(eMatDH_col);
                        if (emission_max < 1) {
                            eMatDH_col *= (1 / emission_max);
                        }
                    }
		    sum_e_times_b = 0;
		    for(k = 0; k < K; k++) {
                        prob = eMatDH_col(dh_col(k));
			if (dh_col(k) == 0) {
			    //
			    std::uint32_t tmp(rhb_t(k, iGrid + 1)); 
			    //
			    prob = 1;
			    for(int b = 0; b < nSNPsLocal; b++) {
			      dR = gl_local(0, b);
			      dA = gl_local(1, b);
			      if (tmp & (1<<b)) {
                                // alternate
                                prob *= (dR * ref_error + dA * ref_one_minus_error);
			      } else {
                                prob *= (dR * ref_one_minus_error + dA * ref_error);
			      }
			    }
                            prob *= (1 / emission_max);
			    //prob_col(k) = prob;
			}
			//ematcol(k) = prob;		    
			e_times_b(k) = betaHat_t_col(k) * prob;
			sum_e_times_b += e_times_b(k);
		    }
		    val = jump_prob * sum_e_times_b;
		    betaHat_t_col = ((not_jump_prob) * e_times_b + val);
		    prev_val = -1;
		} else {
		    // so here prob is uniformly 1, can then super simplify
		    // also, "val", can work out math, is the same if two consecutive no variants, so avoid another fullsum
		    if (prev_val == -1) {
		         val = jump_prob * sum(betaHat_t_col);
		    } else {
		     	double J1 = (transMatRate_t(1, iGrid + 1) / double_K);
		     	double N1 = transMatRate_t(0, iGrid + 1);
		     	val = prev_val * (N1 * jump_prob / J1 + double_K * jump_prob);
		    }
		    betaHat_t_col = ((not_jump_prob) * betaHat_t_col + val);
		    prev_val = val * c(iGrid); // c should not be relevant though, should be 1, has no variants
                }
	    } else {
	        for(k = 0; k < K; k++) {
		    //
		  std::uint32_t tmp(rhb_t(k, iGrid + 1));
		  //
		  prob = 1;
		  for(b = 0; b < nSNPsLocal; b++) {
		      dR = gl_local(0, b);
		      dA = gl_local(1, b);
		      if (tmp & (1<<b)) {
			// alternate
			prob *= (dR * ref_error + dA * ref_one_minus_error);
		      } else {
			prob *= (dR * ref_one_minus_error + dA * ref_error);
		      }
		  }
		  ematcol(k) = prob;
		}
		e_times_b = betaHat_t_col % ematcol;
		val = jump_prob * sum(e_times_b);
		betaHat_t_col = ((not_jump_prob) * e_times_b + val);
	    }
        }
	//
	// all this was to give us the unnormalized beta column that we can now work with
	// now with build gammas and dosages etc
        //
        calculate_small_gamma_t_col = false;
        if (return_gammaSmall_t) {
            if (gammaSmall_cols_to_get(iGrid) >= 0) {
                calculate_small_gamma_t_col = true;
            }
        }
        if (get_best_haps_from_thinned_sites && (gammaSmall_cols_to_get(iGrid) >= 0)) {
            best_haps_stuff_list(gammaSmall_cols_to_get(iGrid)) = Rcpp_get_top_K_or_more_matches_while_building_gamma(alphaHat_t, betaHat_t_col, gamma_t_col, iGrid, K, K_top_matches, 1);
            // special multiplication value is 1 because alpha and beta are normal here (up to c)
        } else {
            // otherwise, do per-entry, check for kth best value
            if (return_dosage | return_gamma_t | calculate_small_gamma_t_col) {
                gamma_t_col = alphaHat_t.col(iGrid) % betaHat_t_col; // betaHat_t_col includes effect of c, so this is accurate and does not need more c
            }
        }
        //
        if (return_dosage) {
            s = 32 * (iGrid); // 0-based here
            e = 32 * (iGrid + 1) - 1;
            if (e > (nSNPs - 1)) {
                e = nSNPs - 1;
            }
            nSNPsLocal = e - s + 1;
            matched_gammas.fill(0);
            dosageL.fill(0);
	    dh_col = hapMatcher.col(iGrid);
            if (use_eMatDH) {
                for(k = 0; k < K; k++) {
                    // some of the dh will be 0 and go to matched_gammas 0th entry, but that is OK we do not use that
		    // they are dealt with afterwards
		    matched_gammas(dh_col(k)) += gamma_t_col(k);
                    if (dh_col(k) == 0) {
		        gk = gamma_t_col(k);
                        std::uint32_t tmp(rhb_t(k, iGrid));
                        for(b = 0; b < nSNPsLocal; b++) {
                            if (tmp & (1<<b)) {
                                // alternate
                                dosageL(b) += gk * (ref_one_minus_error);
                            } else {
                                // reference
                                dosageL(b) += gk * (ref_error);	    
                            }
                        }
                    }
                }
                for(b = 0; b < nSNPsLocal; b++) {
                    for(dh = 0; dh < nMaxDH; dh++) {
                        dosageL(b) += distinctHapsIE(dh, s + b) * matched_gammas(dh + 1);
                    }
                    dosage(s + b) = dosageL(b);
                }
            } else {
                for(k = 0; k < K; k++) {
                    gk = gamma_t_col(k);
                    std::uint32_t tmp(rhb_t(k, iGrid));
                    for(b = 0; b < nSNPsLocal; b++) {
                        if (tmp & (1<<b)) {
                            // alternate
                            dosageL(b) += gk * (ref_one_minus_error);
                        } else {
                            // reference
                            dosageL(b) += gk * (ref_error);
                        }
                }
                }
                // put back
                for(b = 0; b < nSNPsLocal; b++) {
                    dosage(s + b) = dosageL(b);
                }
            }
        }
        // add in c here to betaHat_t_col (as long as not 1!)
	if (c(iGrid) != 1) {
	    betaHat_t_col *= c(iGrid);
	}
        if (return_betaHat_t) {
            betaHat_t.col(iGrid) = betaHat_t_col;
        }
        if (return_gamma_t) {
            gamma_t.col(iGrid) = gamma_t_col;
        }
        if (calculate_small_gamma_t_col) {
            gammaSmall_t.col(gammaSmall_cols_to_get(iGrid)) = gamma_t_col;
        }
    }
    return;
}







//' @export
// [[Rcpp::export]]
void Rcpp_haploid_reference_single_backward_version2(
    arma::mat& alphaHat_t,
    arma::mat& betaHat_t,
    arma::mat& gamma_t,
    arma::mat& gammaSmall_t,
    Rcpp::List& best_haps_stuff_list,
    Rcpp::IntegerVector& gammaSmall_cols_to_get,
    Rcpp::NumericVector& dosage,    
    const int& nGrids,
    const arma::mat& transMatRate_t,
    arma::mat& eMatDH,
    arma::imat& hapMatcher,
    Rcpp::RawMatrix& hapMatcherR,
    bool use_hapMatcherR,    
    const int& nSNPs,
    const int& K,
    const bool& use_eMatDH,
    const arma::imat& rhb_t,
    Rcpp::IntegerMatrix& eMatDH_special_matrix_helper,
    Rcpp::IntegerMatrix& eMatDH_special_matrix,
    const bool use_eMatDH_special_symbols,    
    double ref_error,
    const arma::mat& gl,
    arma::rowvec& c,
    arma::mat& distinctHapsIE,
    bool return_betaHat_t,
    bool return_dosage,
    bool return_gamma_t,
    bool return_gammaSmall_t,
    bool get_best_haps_from_thinned_sites,
    const int nMaxDH,
    const int K_top_matches,
    const Rcpp::IntegerVector& eMatDH_special_grid_which,
    const Rcpp::List& eMatDH_special_values_list,
    const double maxEmissionMatrixDifference,
    const bool normalize_emissions
) {
    // 
    //
    double ref_one_minus_error = 1 - ref_error;    
    double jump_prob; //
    double one_minus_jump_prob = 1;
    double not_jump_prob = 1;
    int b, s, e, nSNPsLocal, i, iGrid, k, dh;
    const double double_K = double(K);
    arma::mat gl_local(2, 32);
    double dR, dA, prob, val, gk;
    arma::colvec ematcol(K);
    iGrid = nGrids - 1;
    arma::colvec alphaHat_t_col(K);    
    arma::colvec betaHat_t_col(K);
    arma::colvec gamma_t_col(K);
    bool calculate_small_gamma_t_col, grid_has_variant;
    arma::vec matched_gammas(nMaxDH + 1);
    arma::colvec e_times_b(K);
    arma::colvec temp_e_times_b(K);    
    arma::vec dosageL(32);
    dosageL.fill(0);
    Rcpp::RawVector dh_colR(K);
    arma::icolvec dh_col(K);        
    arma::colvec eMatDH_col(nMaxDH + 1);
    arma::colvec prob_col(K);
    double sum_e_times_b = 0;
    //double prev_val = -1; // argh
    double min_emission_prob = 1;
    double max_emission_prob = 1;
    double B_prev = 1;
    double B_prev_star = 1;
    double emission_max = 1;
    Rcpp::IntegerVector vals_to_redo;
    int bvtd, s1, e1;    
    //
    //
    for(iGrid = nGrids - 1; iGrid >= 0; --iGrid) {
        if (iGrid == (nGrids - 1)) {
            betaHat_t_col.fill(1 / not_jump_prob);
            B_prev_star = K * c(iGrid) * not_jump_prob;
        } else {
            jump_prob = transMatRate_t(1, iGrid) / double_K;
            one_minus_jump_prob = (1 - jump_prob);
            not_jump_prob = transMatRate_t(0, iGrid);
            s = 32 * (iGrid + 1); // 0-based here
            e = 32 * (iGrid + 1 + 1) - 1;
            if (e > (nSNPs - 1)) {
                e = nSNPs - 1;
            }
	    grid_has_variant = false;
            nSNPsLocal = e - s + 1;
            for(i = 0; i < nSNPsLocal; i++) {
                gl_local(0, i) = gl(0, i + s);
                gl_local(1, i) = gl(1, i + s);
		if ((gl_local(0, i) != 1) | (gl_local(1, i) != 1)) {
		  grid_has_variant = true;
		}
            }
	    //
	    // use eMatDH
	    //
	    if (use_eMatDH) {
	        if (grid_has_variant) {
                    if (use_hapMatcherR) {
                        dh_colR = hapMatcherR(Rcpp::_, iGrid + 1);
                    } else {
                        dh_col = hapMatcher.col(iGrid + 1);
                    }
		    eMatDH_col = eMatDH.col(iGrid + 1);
                    //
                    if (normalize_emissions) {
                        emission_max = arma::max(eMatDH_col);
                        if (emission_max < 1) {
                            eMatDH_col *= (1 / emission_max);
                        }
                    }
                    //
                    min_emission_prob = arma::min(eMatDH_col);
                    max_emission_prob = arma::max(eMatDH_col);
                    //
		    sum_e_times_b = 0;
                    //
                    // first, do special cases
                    //
                    if (eMatDH_special_grid_which(iGrid + 1) > 0) {
                        if (use_eMatDH_special_symbols) {
                            s1 = eMatDH_special_matrix_helper(iGrid + 1, 0);
                            e1 = eMatDH_special_matrix_helper(iGrid + 1, 1);
                            // there is undoubtledly better code to do this
                            Rcpp::IntegerVector temp_vec(e1 - s1 + 1);
                            for(int ii = 0; ii < e1 - s1 + 1; ii++) {
                                temp_vec(ii) = eMatDH_special_matrix(s1 - 1 + ii, 0);
                            }
                            vals_to_redo = temp_vec;
                        } else {
                            vals_to_redo = Rcpp::as<Rcpp::IntegerVector>(eMatDH_special_values_list(eMatDH_special_grid_which(iGrid + 1) - 1));
                        }
                        for(i = 0; i < vals_to_redo.size(); i++) {
                            k = vals_to_redo(i);
                            if (use_eMatDH_special_symbols) {
                                bvtd = rcpp_simple_binary_matrix_search(k, eMatDH_special_matrix, s1, e1);
                            } else {
                                bvtd = rhb_t(k, iGrid + 1);
                            }
			    //
                            std::uint32_t tmp(bvtd);
			    //
			    prob = 1;
			    for(int b = 0; b < nSNPsLocal; b++) {
			      dR = gl_local(0, b);
			      dA = gl_local(1, b);
			      if (tmp & (1<<b)) {
                                // alternate
                                prob *= (dR * ref_error + dA * ref_one_minus_error);
			      } else {
                                prob *= (dR * ref_one_minus_error + dA * ref_error);
			      }
			    }
                            prob *= 1 / emission_max;
                            temp_e_times_b(k) = betaHat_t_col(k) * prob;
                            sum_e_times_b += temp_e_times_b(k);
                        }
                    }
                    eMatDH_col(0) = 0;
                    //
                    // now do normal
                    //
                    if (use_hapMatcherR) {
                        for(k = 0; k < K; k++) {
                            e_times_b(k) = betaHat_t_col(k) * eMatDH_col(dh_colR(k));
                            sum_e_times_b += e_times_b(k);
                        }
                    } else {
                        for(k = 0; k < K; k++) {
                            e_times_b(k) = betaHat_t_col(k) * eMatDH_col(dh_col(k));
                            sum_e_times_b += e_times_b(k);
                        }
                    }
                    //
                    // put in over-ridden values, if applicable
                    //
                    if (eMatDH_special_grid_which(iGrid + 1) > 0) {
                        for(i = 0; i < vals_to_redo.size(); i++) {
                            k = vals_to_redo(i);
                            e_times_b(k) = temp_e_times_b(k);
                        }
                    }
                    //
                    //
		    val = jump_prob / not_jump_prob * sum_e_times_b;
		    betaHat_t_col = (e_times_b + val);
                    //
                    B_prev = sum_e_times_b;
                    B_prev_star = c(iGrid) * B_prev;
                    //
		} else {
                    //
                    // previous version
                    //
                    val = jump_prob / not_jump_prob * B_prev_star;
		    betaHat_t_col = betaHat_t_col + val;
                    //
                    B_prev = B_prev_star;
                    B_prev_star = c(iGrid) * B_prev;
                }
                //
	    } else {
	        for(k = 0; k < K; k++) {
		  //
		  std::uint32_t tmp(rhb_t(k, iGrid + 1));
		  //
		  prob = 1;
		  for(b = 0; b < nSNPsLocal; b++) {
		      dR = gl_local(0, b);
		      dA = gl_local(1, b);
		      if (tmp & (1<<b)) {
			// alternate
			prob *= (dR * ref_error + dA * ref_one_minus_error);
		      } else {
			prob *= (dR * ref_one_minus_error + dA * ref_error);
		      }
		  }
		  ematcol(k) = prob;
		}
		e_times_b = betaHat_t_col % ematcol;
		//val = jump_prob * sum(e_times_b);
		//betaHat_t_col = ((not_jump_prob) * e_times_b + val);
		val = jump_prob  / not_jump_prob * sum(e_times_b);
		betaHat_t_col = (e_times_b + val);
	    }
            //
        }
        //
	// all this was to give us the unnormalized beta column that we can now work with
	// now with build gammas and dosages etc
        //
        calculate_small_gamma_t_col = false;
        if (return_gammaSmall_t) {
            if (gammaSmall_cols_to_get(iGrid) >= 0) {
                calculate_small_gamma_t_col = true;
            }
        }
        if (get_best_haps_from_thinned_sites && (gammaSmall_cols_to_get(iGrid) >= 0)) {
            best_haps_stuff_list(gammaSmall_cols_to_get(iGrid)) = Rcpp_get_top_K_or_more_matches_while_building_gamma(alphaHat_t, betaHat_t_col, gamma_t_col, iGrid, K, K_top_matches, not_jump_prob);
            // special multiplication value is 1 because alpha and beta are normal here (up to c)
            // note this is simply to simplify testing - it doesn't make any difference to algo
        } else {
            // otherwise, do per-entry, check for kth best value
            if (return_dosage | return_gamma_t | calculate_small_gamma_t_col) {
                // betaHat_t_col includes effect of c, so this is accurate and does not need more c
                // we do not multiply through by not_jump_prob here though
                gamma_t_col = alphaHat_t.col(iGrid) % betaHat_t_col;
            }
        }
        //
        if (return_dosage) {
            s = 32 * (iGrid); // 0-based here
            e = 32 * (iGrid + 1) - 1;
            if (e > (nSNPs - 1)) {
                e = nSNPs - 1;
            }
            nSNPsLocal = e - s + 1;
            matched_gammas.fill(0);
            dosageL.fill(0);
            if (use_hapMatcherR) {
                dh_colR = hapMatcherR(Rcpp::_, iGrid);
            } else {
                dh_col = hapMatcher.col(iGrid);
            }
            if (use_eMatDH) {
                // some of the dh will be 0 and go to matched_gammas 0th entry, but that is OK we do not use that
                // they are dealt with afterwards
                if (use_hapMatcherR) {
                    for(k = 0; k < K; k++) {
                        matched_gammas(dh_colR(k)) += gamma_t_col(k);
                    }
                } else {
                    for(k = 0; k < K; k++) {
                        matched_gammas(dh_col(k)) += gamma_t_col(k);
                    }
                }
                matched_gammas *= not_jump_prob;
                //
                // special cases here
                //
                if (eMatDH_special_grid_which(iGrid) > 0) {
                    if (use_eMatDH_special_symbols) {
                        s1 = eMatDH_special_matrix_helper(iGrid, 0);
                        e1 = eMatDH_special_matrix_helper(iGrid, 1);
                        // there is undoubtledly better code to do this
                        Rcpp::IntegerVector temp_vec(e1 - s1 + 1);
                        for(int ii = 0; ii < e1 - s1 + 1; ii++) {
                            temp_vec(ii) = eMatDH_special_matrix(s1 - 1 + ii, 0);
                        }
                        vals_to_redo = temp_vec;
                    } else {
                        vals_to_redo = Rcpp::as<Rcpp::IntegerVector>(eMatDH_special_values_list(eMatDH_special_grid_which(iGrid) - 1));
                    }
                    for(i = 0; i < vals_to_redo.size(); i++) {
                        k = vals_to_redo(i);
                        if (use_eMatDH_special_symbols) {
                            bvtd = rcpp_simple_binary_matrix_search(k, eMatDH_special_matrix, s1, e1);
                        } else {
                            bvtd = rhb_t(k, iGrid);
                        }
                        //
                        gk = gamma_t_col(k) * not_jump_prob;
                        std::uint32_t tmp(bvtd);
                        for(b = 0; b < nSNPsLocal; b++) {
                            if (tmp & (1<<b)) {
                                // alternate
                                dosageL(b) += gk * (ref_one_minus_error);
                            } else {
                                // reference
                                dosageL(b) += gk * (ref_error);	    
                            }
                        }
                    }
                }
                
                //
                // now do the sums
                //
                for(b = 0; b < nSNPsLocal; b++) {
                    for(dh = 0; dh < nMaxDH; dh++) {
                        dosageL(b) += distinctHapsIE(dh, s + b) * matched_gammas(dh + 1);
                    }
                    dosage(s + b) = dosageL(b);
                }
            } else {
                for(k = 0; k < K; k++) {
                    gk = gamma_t_col(k) * not_jump_prob;
                    std::uint32_t tmp(rhb_t(k, iGrid));
                    for(b = 0; b < nSNPsLocal; b++) {
                        if (tmp & (1<<b)) {
                            // alternate
                            dosageL(b) += gk * (ref_one_minus_error);
                        } else {
                            // reference
                            dosageL(b) += gk * (ref_error);
                        }
                    }
                }
                // put back
                for(b = 0; b < nSNPsLocal; b++) {
                    dosage(s + b) = dosageL(b);
                }
            }
        }
        //
        // add in c here to betaHat_t_col (as long as not 1!)
        //
        // I don't think I can ever skip these? they are off-by-1 on not-jump-prob?
        //std::cout << "c(iGrid) = " << c(iGrid) << ", 1 / not_jump_prob = " << 1 / not_jump_prob << std::endl;
        double x = c(iGrid) * not_jump_prob;
        betaHat_t_col *= x;
        if (return_betaHat_t) {
            betaHat_t.col(iGrid) = betaHat_t_col;
        }
        if (return_gamma_t) {
            gamma_t_col *= not_jump_prob;
            gamma_t.col(iGrid) = gamma_t_col;
        }
        if (calculate_small_gamma_t_col) {
            gammaSmall_t.col(gammaSmall_cols_to_get(iGrid)) = gamma_t_col;
        }
    }
    return;
}






//' @export
// [[Rcpp::export]]
void Rcpp_haploid_reference_single_backward_version3(
    Eigen::Map<Eigen::MatrixXd> alphaHat_t,
    arma::mat& betaHat_t,
    arma::mat& gamma_t,
    arma::mat& gammaSmall_t,
    Rcpp::List& best_haps_stuff_list,
    Rcpp::IntegerVector& gammaSmall_cols_to_get,
    Rcpp::NumericVector& dosage,    
    const int& nGrids,
    const arma::mat& transMatRate_t,
    arma::mat& eMatDH,
    arma::imat& hapMatcher,
    Rcpp::RawMatrix& hapMatcherR,
    bool use_hapMatcherR,    
    const int& nSNPs,
    const int& K,
    const bool& use_eMatDH,
    const arma::imat& rhb_t,
    Rcpp::IntegerMatrix& eMatDH_special_matrix_helper,
    Rcpp::IntegerMatrix& eMatDH_special_matrix,
    const bool use_eMatDH_special_symbols,    
    double ref_error,
    const arma::mat& gl,
    arma::rowvec& c,
    arma::mat& distinctHapsIE,
    bool return_betaHat_t,
    bool return_dosage,
    bool return_gamma_t,
    bool return_gammaSmall_t,
    bool get_best_haps_from_thinned_sites,
    const int nMaxDH,
    const int K_top_matches,
    const Rcpp::IntegerVector& eMatDH_special_grid_which,
    const Rcpp::List& eMatDH_special_values_list,
    const double maxEmissionMatrixDifference,
    const bool normalize_emissions
) {
    // 
    //
    double ref_one_minus_error = 1 - ref_error;    
    double jump_prob; //
    double one_minus_jump_prob = 1;
    double not_jump_prob = 1;
    int b, s, e, nSNPsLocal, i, iGrid, k, dh;
    const double double_K = double(K);
    arma::mat gl_local(2, 32);
    double dR, dA, prob, val, gk;
    arma::colvec ematcol(K);
    iGrid = nGrids - 1;
    arma::colvec alphaHat_t_col(K);    
    arma::colvec betaHat_t_col(K);
    arma::colvec gamma_t_col(K);
    bool calculate_small_gamma_t_col, grid_has_variant;
    arma::vec matched_gammas(nMaxDH + 1);
    arma::colvec e_times_b(K);
    arma::colvec temp_e_times_b(K);    
    arma::vec dosageL(32);
    dosageL.fill(0);
    Rcpp::RawVector dh_colR(K);
    arma::icolvec dh_col(K);        
    arma::colvec eMatDH_col(nMaxDH + 1);
    arma::colvec prob_col(K);
    double sum_e_times_b = 0;
    //double prev_val = -1; // argh
    double min_emission_prob = 1;
    double max_emission_prob = 1;
    double B_prev = 1;
    double B_prev_star = 1;
    double emission_max = 1;
    Rcpp::IntegerVector vals_to_redo;
    int bvtd, s1, e1;    
    //
    //
    for(iGrid = nGrids - 1; iGrid >= 0; --iGrid) {
        if (iGrid == (nGrids - 1)) {
            betaHat_t_col.fill(1 / not_jump_prob);
            B_prev_star = K * c(iGrid) * not_jump_prob;
        } else {
            jump_prob = transMatRate_t(1, iGrid) / double_K;
            one_minus_jump_prob = (1 - jump_prob);
            not_jump_prob = transMatRate_t(0, iGrid);
            s = 32 * (iGrid + 1); // 0-based here
            e = 32 * (iGrid + 1 + 1) - 1;
            if (e > (nSNPs - 1)) {
                e = nSNPs - 1;
            }
	    grid_has_variant = false;
            nSNPsLocal = e - s + 1;
            for(i = 0; i < nSNPsLocal; i++) {
                gl_local(0, i) = gl(0, i + s);
                gl_local(1, i) = gl(1, i + s);
		if ((gl_local(0, i) != 1) | (gl_local(1, i) != 1)) {
		  grid_has_variant = true;
		}
            }
	    //
	    // use eMatDH
	    //
	    if (use_eMatDH) {
	        if (grid_has_variant) {
                    if (use_hapMatcherR) {
                        dh_colR = hapMatcherR(Rcpp::_, iGrid + 1);
                    } else {
                        dh_col = hapMatcher.col(iGrid + 1);
                    }
		    eMatDH_col = eMatDH.col(iGrid + 1);
                    //
                    if (normalize_emissions) {
                        emission_max = arma::max(eMatDH_col);
                        if (emission_max < 1) {
                            eMatDH_col *= (1 / emission_max);
                        }
                    }
                    //
                    min_emission_prob = arma::min(eMatDH_col);
                    max_emission_prob = arma::max(eMatDH_col);
                    //
		    sum_e_times_b = 0;
                    //
                    // first, do special cases
                    //
                    if (eMatDH_special_grid_which(iGrid + 1) > 0) {
                        if (use_eMatDH_special_symbols) {
                            s1 = eMatDH_special_matrix_helper(iGrid + 1, 0);
                            e1 = eMatDH_special_matrix_helper(iGrid + 1, 1);
                            // there is undoubtledly better code to do this
                            Rcpp::IntegerVector temp_vec(e1 - s1 + 1);
                            for(int ii = 0; ii < e1 - s1 + 1; ii++) {
                                temp_vec(ii) = eMatDH_special_matrix(s1 - 1 + ii, 0);
                            }
                            vals_to_redo = temp_vec;
                        } else {
                            vals_to_redo = Rcpp::as<Rcpp::IntegerVector>(eMatDH_special_values_list(eMatDH_special_grid_which(iGrid + 1) - 1));
                        }
                        for(i = 0; i < vals_to_redo.size(); i++) {
                            k = vals_to_redo(i);
                            if (use_eMatDH_special_symbols) {
                                bvtd = rcpp_simple_binary_matrix_search(k, eMatDH_special_matrix, s1, e1);
                            } else {
                                bvtd = rhb_t(k, iGrid + 1);
                            }
			    //
                            std::uint32_t tmp(bvtd);
			    //
			    prob = 1;
			    for(int b = 0; b < nSNPsLocal; b++) {
			      dR = gl_local(0, b);
			      dA = gl_local(1, b);
			      if (tmp & (1<<b)) {
                                // alternate
                                prob *= (dR * ref_error + dA * ref_one_minus_error);
			      } else {
                                prob *= (dR * ref_one_minus_error + dA * ref_error);
			      }
			    }
                            prob *= 1 / emission_max;
                            temp_e_times_b(k) = betaHat_t_col(k) * prob;
                            sum_e_times_b += temp_e_times_b(k);
                        }
                    }
                    eMatDH_col(0) = 0;
                    //
                    // now do normal
                    //
                    if (use_hapMatcherR) {
                        for(k = 0; k < K; k++) {
                            e_times_b(k) = betaHat_t_col(k) * eMatDH_col(dh_colR(k));
                            sum_e_times_b += e_times_b(k);
                        }
                    } else {
                        for(k = 0; k < K; k++) {
                            e_times_b(k) = betaHat_t_col(k) * eMatDH_col(dh_col(k));
                            sum_e_times_b += e_times_b(k);
                        }
                    }
                    //
                    // put in over-ridden values, if applicable
                    //
                    if (eMatDH_special_grid_which(iGrid + 1) > 0) {
                        for(i = 0; i < vals_to_redo.size(); i++) {
                            k = vals_to_redo(i);
                            e_times_b(k) = temp_e_times_b(k);
                        }
                    }
                    //
                    //
		    val = jump_prob / not_jump_prob * sum_e_times_b;
		    betaHat_t_col = (e_times_b + val);
                    //
                    B_prev = sum_e_times_b;
                    B_prev_star = c(iGrid) * B_prev;
                    //
		} else {
                    //
                    // previous version
                    //
                    val = jump_prob / not_jump_prob * B_prev_star;
		    betaHat_t_col = betaHat_t_col + val;
                    //
                    B_prev = B_prev_star;
                    B_prev_star = c(iGrid) * B_prev;
                }
                //
	    } else {
	        for(k = 0; k < K; k++) {
		  //
		  std::uint32_t tmp(rhb_t(k, iGrid + 1));
		  //
		  prob = 1;
		  for(b = 0; b < nSNPsLocal; b++) {
		      dR = gl_local(0, b);
		      dA = gl_local(1, b);
		      if (tmp & (1<<b)) {
			// alternate
			prob *= (dR * ref_error + dA * ref_one_minus_error);
		      } else {
			prob *= (dR * ref_one_minus_error + dA * ref_error);
		      }
		  }
		  ematcol(k) = prob;
		}
		e_times_b = betaHat_t_col % ematcol;
		//val = jump_prob * sum(e_times_b);
		//betaHat_t_col = ((not_jump_prob) * e_times_b + val);
		val = jump_prob  / not_jump_prob * sum(e_times_b);
		betaHat_t_col = (e_times_b + val);
	    }
            //
        }
        //
	// all this was to give us the unnormalized beta column that we can now work with
	// now with build gammas and dosages etc
        //
        calculate_small_gamma_t_col = false;
        if (return_gammaSmall_t) {
            if (gammaSmall_cols_to_get(iGrid) >= 0) {
                calculate_small_gamma_t_col = true;
            }
        }
        if (get_best_haps_from_thinned_sites && (gammaSmall_cols_to_get(iGrid) >= 0)) {
            
            //std::cout << "in backward check, first one, iGrid = " << iGrid << std::endl;
        // for(k = 0; k < 1; k++) {
        //     std::cout << "before" << std::endl;
        //     std::cout << "alphaHat_t(k, iGrid)=" << alphaHat_t(k, iGrid) << std::endl;
        //     std::cout << "betaHat_t_col(k)=" << betaHat_t_col(k) << std::endl;                        
        //     std::cout << "gamma_t_col(k)=" << gamma_t_col(k) << std::endl;
        // }
        // std::cout << "sum(gamma_t_col)=" << sum(gamma_t_col) << std::endl;                
            
            best_haps_stuff_list(gammaSmall_cols_to_get(iGrid)) = Rcpp_get_top_K_or_more_matches_while_building_gamma_eigen(alphaHat_t, betaHat_t_col, gamma_t_col, iGrid, K, K_top_matches, not_jump_prob);
            // special multiplication value is 1 because alpha and beta are normal here (up to c)
            // note this is simply to simplify testing - it doesn't make any difference to algo

        // for(k = 0; k < 1; k++) {
        //     std::cout << "after" << std::endl;
        //     std::cout << "alphaHat_t(k, iGrid)=" << alphaHat_t(k, iGrid) << std::endl;
        //     std::cout << "betaHat_t_col(k)=" << betaHat_t_col(k) << std::endl;                        
        //     std::cout << "gamma_t_col(k)=" << gamma_t_col(k) << std::endl;
        // }
        // std::cout << "sum(gamma_t_col)=" << sum(gamma_t_col) << std::endl;        
            
        } else {
            // otherwise, do per-entry, check for kth best value
            if (return_dosage | return_gamma_t | calculate_small_gamma_t_col) {
                // betaHat_t_col includes effect of c, so this is accurate and does not need more c
                // we do not multiply through by not_jump_prob here though
                for(k = 0; k < K; k++) {
                    gamma_t_col(k) = alphaHat_t(k, iGrid) * betaHat_t_col(k);
                }
            }
        //     std::cout << "in backward check, second one, iGrid = " << iGrid << std::endl;
        // for(k = 0; k < 1; k++) {
        //     std::cout << "alphaHat_t(k, iGrid)=" << alphaHat_t(k, iGrid) << std::endl;
        //     std::cout << "betaHat_t_col(k)=" << betaHat_t_col(k) << std::endl;                        
        //     std::cout << "gamma_t_col(k)=" << gamma_t_col(k) << std::endl;
        // }
        // std::cout << "sum(gamma_t_col)=" << sum(gamma_t_col) << std::endl;        
            
        }


        
        
        //
        if (return_dosage) {
            s = 32 * (iGrid); // 0-based here
            e = 32 * (iGrid + 1) - 1;
            if (e > (nSNPs - 1)) {
                e = nSNPs - 1;
            }
            nSNPsLocal = e - s + 1;
            matched_gammas.fill(0);
            dosageL.fill(0);
            if (use_hapMatcherR) {
                dh_colR = hapMatcherR(Rcpp::_, iGrid);
            } else {
                dh_col = hapMatcher.col(iGrid);
            }
            if (use_eMatDH) {
                // some of the dh will be 0 and go to matched_gammas 0th entry, but that is OK we do not use that
                // they are dealt with afterwards
                if (use_hapMatcherR) {
                    for(k = 0; k < K; k++) {
                        matched_gammas(dh_colR(k)) += gamma_t_col(k);
                    }
                } else {
                    for(k = 0; k < K; k++) {
                        matched_gammas(dh_col(k)) += gamma_t_col(k);
                    }
                }
                matched_gammas *= not_jump_prob;
                //
                // special cases here
                //
                if (eMatDH_special_grid_which(iGrid) > 0) {
                    if (use_eMatDH_special_symbols) {
                        s1 = eMatDH_special_matrix_helper(iGrid, 0);
                        e1 = eMatDH_special_matrix_helper(iGrid, 1);
                        // there is undoubtledly better code to do this
                        Rcpp::IntegerVector temp_vec(e1 - s1 + 1);
                        for(int ii = 0; ii < e1 - s1 + 1; ii++) {
                            temp_vec(ii) = eMatDH_special_matrix(s1 - 1 + ii, 0);
                        }
                        vals_to_redo = temp_vec;
                    } else {
                        vals_to_redo = Rcpp::as<Rcpp::IntegerVector>(eMatDH_special_values_list(eMatDH_special_grid_which(iGrid) - 1));
                    }
                    for(i = 0; i < vals_to_redo.size(); i++) {
                        k = vals_to_redo(i);
                        if (use_eMatDH_special_symbols) {
                            bvtd = rcpp_simple_binary_matrix_search(k, eMatDH_special_matrix, s1, e1);
                        } else {
                            bvtd = rhb_t(k, iGrid);
                        }
                        //
                        gk = gamma_t_col(k) * not_jump_prob;
                        std::uint32_t tmp(bvtd);
                        for(b = 0; b < nSNPsLocal; b++) {
                            if (tmp & (1<<b)) {
                                // alternate
                                dosageL(b) += gk * (ref_one_minus_error);
                            } else {
                                // reference
                                dosageL(b) += gk * (ref_error);	    
                            }
                        }
                    }
                }
                
                //
                // now do the sums
                //
                for(b = 0; b < nSNPsLocal; b++) {
                    for(dh = 0; dh < nMaxDH; dh++) {
                        dosageL(b) += distinctHapsIE(dh, s + b) * matched_gammas(dh + 1);
                    }
                    dosage(s + b) = dosageL(b);
                }
            } else {
                for(k = 0; k < K; k++) {
                    gk = gamma_t_col(k) * not_jump_prob;
                    std::uint32_t tmp(rhb_t(k, iGrid));
                    for(b = 0; b < nSNPsLocal; b++) {
                        if (tmp & (1<<b)) {
                            // alternate
                            dosageL(b) += gk * (ref_one_minus_error);
                        } else {
                            // reference
                            dosageL(b) += gk * (ref_error);
                        }
                    }
                }
                // put back
                for(b = 0; b < nSNPsLocal; b++) {
                    dosage(s + b) = dosageL(b);
                }
            }
        }
        //
        // add in c here to betaHat_t_col (as long as not 1!)
        //
        // I don't think I can ever skip these? they are off-by-1 on not-jump-prob?
        //std::cout << "c(iGrid) = " << c(iGrid) << ", 1 / not_jump_prob = " << 1 / not_jump_prob << std::endl;
        double x = c(iGrid) * not_jump_prob;
        betaHat_t_col *= x;
        if (return_betaHat_t) {
            betaHat_t.col(iGrid) = betaHat_t_col;
        }
        if (return_gamma_t) {
            gamma_t_col *= not_jump_prob;
            gamma_t.col(iGrid) = gamma_t_col;
        }
        if (calculate_small_gamma_t_col) {
            gammaSmall_t.col(gammaSmall_cols_to_get(iGrid)) = gamma_t_col;
        }
    }
    return;
}







//' @export
// [[Rcpp::export]]
void Rcpp_haploid_dosage_versus_refs(
    const arma::mat& gl,
    arma::mat& arma_alphaHat_t,
    Eigen::Map<Eigen::MatrixXd> eigen_alphaHat_t,
    arma::mat& betaHat_t,
    arma::rowvec& c,
    arma::mat& gamma_t,
    arma::mat& gammaSmall_t,
    Rcpp::List& best_haps_stuff_list,
    Rcpp::NumericVector& dosage,
    const arma::mat& transMatRate_t,
    const arma::imat& rhb_t,
    double ref_error,
    const bool use_eMatDH,
    arma::imat& distinctHapsB,
    arma::mat& distinctHapsIE,
    Rcpp::IntegerMatrix& eMatDH_special_matrix_helper,
    Rcpp::IntegerMatrix& eMatDH_special_matrix,
    const bool use_eMatDH_special_symbols,    
    arma::imat& hapMatcher,
    Rcpp::RawMatrix& hapMatcherR,
    bool use_hapMatcherR,
    Rcpp::IntegerVector& gammaSmall_cols_to_get,
    const Rcpp::IntegerVector& eMatDH_special_grid_which,
    const Rcpp::List& eMatDH_special_values_list,
    const int K_top_matches,
    const int suppressOutput = 1,
    const double min_emission_prob_normalization_threshold = 1e-100,
    bool return_betaHat_t = true,
    bool return_dosage = true,
    bool return_gamma_t = true,
    bool return_gammaSmall_t = false,
    bool get_best_haps_from_thinned_sites = false,
    bool is_version_2 = true,
    bool is_version_3 = false,
    bool return_extra = false,
    bool always_normalize = true,
    bool use_eigen = false,
    const bool normalize_emissions = true
) {
    //
    // not used currently, the below
    // normalization, when invoked, just makes sure most likely has value 1
    const double maxEmissionMatrixDifference = 1e10;
    //
    double prev=clock();
    timeval atime;
    timeval btime;
    if (suppressOutput == 0) {
        gettimeofday(&atime, 0);
        //std::cout << "==== start at ==== " << std::endl;
        char buffer[30];
        struct timeval tv;
        time_t curtime;
        gettimeofday(&tv, NULL);
        curtime=tv.tv_sec;
        strftime(buffer,30,"%m-%d-%Y  %T.",localtime(&curtime));
        //printf("%s%ld\n",buffer,tv.tv_usec);
        //std::cout << "================== " << std::endl;        
    }
    std::string prev_section="Null";
    std::string next_section="Initialize variables";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    const int K = hapMatcher.n_rows;
    const int nGrids = transMatRate_t.n_cols + 1;
    const int nSNPs = gl.n_cols;
    double one_over_K = 1 / (double)K;
    double ref_one_minus_error = 1 - ref_error;
    //arma::rowvec c = arma::zeros(1, nGrids);
    arma::mat gl_local(2, 32);
    int i, k, dh;
    double prob;
    bool only_store_alpha_at_gamma_small;
    if ((get_best_haps_from_thinned_sites | return_gammaSmall_t) & (!return_gamma_t) & (!return_dosage) & (!return_betaHat_t)) {
        only_store_alpha_at_gamma_small = true;
    } else {
        only_store_alpha_at_gamma_small = false;
    }
    //
    //
    //
    arma::mat distinctHapsIE_local;
    arma::colvec e_times_b(K);
    arma::vec dosageL(32);
    dosageL.fill(0);
    const int nMaxDH = distinctHapsB.n_rows;
    arma::vec matched_gammas(nMaxDH);
    arma::mat eMatDH;
    if (use_eMatDH) {
        next_section="make eMatDH";
        prev=print_times(prev, suppressOutput, prev_section, next_section);
        prev_section=next_section;
        const bool add_zero_row_true = true;
        eMatDH = Rcpp_build_eMatDH(distinctHapsB, gl, nGrids, nSNPs, ref_error, ref_one_minus_error, add_zero_row_true);
    }
    //
    // initialize alphaHat_t
    //
    next_section="initialize alphaHat_t";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    int iGrid = 0;
    int s = 32 * iGrid; // 0-based here
    int e = 32 * (iGrid + 1) - 1;
    if (e > (nSNPs - 1)) {
      e = nSNPs - 1;
    }
    int nSNPsLocal = e - s + 1;
    // now subset
    for(i = 0; i < nSNPsLocal; i++) {
      gl_local(0, i) = gl(0, i + s);
      gl_local(1, i) = gl(1, i + s);      
    }
    int s1 = 0;
    int e1 = 0;
    int bvtd;
    if (use_eMatDH_special_symbols) {
        if (eMatDH_special_grid_which(iGrid) > 0) {
            s1 = eMatDH_special_matrix_helper(0, 0);
            e1 = eMatDH_special_matrix_helper(0, 1);
        }
    }
    arma::colvec alphaHat_t_col(K);
    for(k = 0; k < K; k++) {
        if (use_eMatDH & use_hapMatcherR) {
            dh = hapMatcherR(k, iGrid);
        } else if (use_eMatDH) {
            dh = hapMatcher(k, iGrid);
        } else {
            dh = 0;
        }
        if (dh > 0) {
            prob = eMatDH(dh, iGrid);
        } else {
            //
            if (use_eMatDH_special_symbols) {
                bvtd = rcpp_simple_binary_matrix_search(k, eMatDH_special_matrix, s1, e1);
            } else {
                bvtd = rhb_t(k, iGrid);
            }
            std::uint32_t tmp(bvtd);
            //
            prob = 1;
            for(int b = 0; b < nSNPsLocal; b++) {
                double dR = gl_local(0, b);
                double dA = gl_local(1, b);
                if (tmp & (1<<b)) {
                    // alternate
                    prob *= (dR * ref_error + dA * ref_one_minus_error);
                } else {
                    prob *= (dR * ref_one_minus_error + dA * ref_error);
                }
            }
        }
        alphaHat_t_col(k) = prob * one_over_K;
    }
    c(0) = 1 / sum(alphaHat_t_col);
    if (use_eigen) {
        for(k = 0; k < K; k++) {
            eigen_alphaHat_t(k, 0) = alphaHat_t_col(k) * c(0);
        }
    } else {
        arma_alphaHat_t.col(0) = alphaHat_t_col * c(0);
    }
    //
    
    // run forward
    //
    next_section="run forward";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    if (is_version_2) {
        Rcpp_haploid_reference_single_forward_version2(gammaSmall_cols_to_get, gl, arma_alphaHat_t, c, transMatRate_t, rhb_t, eMatDH_special_matrix_helper, eMatDH_special_matrix, use_eMatDH_special_symbols, hapMatcher, hapMatcherR, use_hapMatcherR, eMatDH, nGrids, nSNPs, K, use_eMatDH, ref_error, only_store_alpha_at_gamma_small, always_normalize, min_emission_prob_normalization_threshold, eMatDH_special_grid_which, eMatDH_special_values_list, maxEmissionMatrixDifference, normalize_emissions);
    } else if (use_eigen) {
        Rcpp_haploid_reference_single_forward_version3(gammaSmall_cols_to_get, gl, eigen_alphaHat_t, c, transMatRate_t, rhb_t, eMatDH_special_matrix_helper, eMatDH_special_matrix, use_eMatDH_special_symbols, hapMatcher, hapMatcherR, use_hapMatcherR, eMatDH, nGrids, nSNPs, K, use_eMatDH, ref_error, only_store_alpha_at_gamma_small, always_normalize, min_emission_prob_normalization_threshold, eMatDH_special_grid_which, eMatDH_special_values_list, maxEmissionMatrixDifference, normalize_emissions);
    } else {
        Rcpp_haploid_reference_single_forward(gammaSmall_cols_to_get, gl, arma_alphaHat_t, c, transMatRate_t, rhb_t, hapMatcher, eMatDH, nGrids, nSNPs, K, use_eMatDH, ref_error, only_store_alpha_at_gamma_small, always_normalize, min_emission_prob_normalization_threshold, maxEmissionMatrixDifference, normalize_emissions);
    }
    //
    // run backward algorithm
    //
    arma::colvec ematcol(K);
    iGrid = nGrids - 1;
    arma::colvec betaHat_t_col(K);
    arma::colvec gamma_t_col(K);
    //
    // run normal backward progression
    //
    
    next_section="run backward normal";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    if (is_version_2) {
        Rcpp_haploid_reference_single_backward_version2(arma_alphaHat_t, betaHat_t, gamma_t, gammaSmall_t, best_haps_stuff_list, gammaSmall_cols_to_get, dosage,     nGrids, transMatRate_t, eMatDH, hapMatcher, hapMatcherR, use_hapMatcherR, nSNPs, K, use_eMatDH, rhb_t, eMatDH_special_matrix_helper, eMatDH_special_matrix, use_eMatDH_special_symbols, ref_error, gl, c, distinctHapsIE, return_betaHat_t, return_dosage, return_gamma_t, return_gammaSmall_t, get_best_haps_from_thinned_sites, nMaxDH, K_top_matches, eMatDH_special_grid_which, eMatDH_special_values_list, maxEmissionMatrixDifference, normalize_emissions);
    } else if (use_eigen) {
        Rcpp_haploid_reference_single_backward_version3(eigen_alphaHat_t, betaHat_t, gamma_t, gammaSmall_t, best_haps_stuff_list, gammaSmall_cols_to_get, dosage,     nGrids, transMatRate_t, eMatDH, hapMatcher, hapMatcherR, use_hapMatcherR, nSNPs, K, use_eMatDH, rhb_t, eMatDH_special_matrix_helper, eMatDH_special_matrix, use_eMatDH_special_symbols, ref_error, gl, c, distinctHapsIE, return_betaHat_t, return_dosage, return_gamma_t, return_gammaSmall_t, get_best_haps_from_thinned_sites, nMaxDH, K_top_matches, eMatDH_special_grid_which, eMatDH_special_values_list, maxEmissionMatrixDifference, normalize_emissions);
    } else {
        Rcpp_haploid_reference_single_backward(arma_alphaHat_t, betaHat_t, gamma_t, gammaSmall_t, best_haps_stuff_list, gammaSmall_cols_to_get, dosage,     nGrids, transMatRate_t, eMatDH, hapMatcher, nSNPs, K, use_eMatDH, rhb_t, ref_error, gl, c, distinctHapsIE, return_betaHat_t, return_dosage, return_gamma_t, return_gammaSmall_t, get_best_haps_from_thinned_sites, nMaxDH, K_top_matches, maxEmissionMatrixDifference, normalize_emissions);
    }
    // 
    //
    next_section="done";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    if (suppressOutput == 0) {
        gettimeofday(&btime, 0);
        //std::cout << "==== end at ==== " << std::endl;
        //std::cout << btime << std::endl;
        long seconds  = btime.tv_sec  - atime.tv_sec;
        long useconds = btime.tv_usec - atime.tv_usec;
        long mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
        printf("Total elapsed time: %ld milliseconds\n", mtime);
        char buffer[30];
        struct timeval tv;
        time_t curtime;
        gettimeofday(&tv, NULL);
        curtime=tv.tv_sec;
        strftime(buffer,30,"%m-%d-%Y  %T.",localtime(&curtime));
        //printf("%s%ld\n",buffer,tv.tv_usec);
        //std::cout << "================ " << std::endl;        
    }
    return;
}





