// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

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
using namespace Rcpp;


double print_times(
    double prev,
    int suppressOutput,
    std::string past_text,
    std::string next_text
);

void rcpp_make_eMatRead_t(
    arma::mat& eMatRead_t,
    const Rcpp::List& sampleReads,
    const arma::cube& eHapsCurrent_tc,
    const int s,
    const double maxDifferenceBetweenReads,
    const int Jmax,
    arma::mat& eMatHapOri_t,
    const arma::vec& pRgivenH1,
    const arma::vec& pRgivenH2,
    double& prev,
    int suppressOutput,
    std::string& prev_section,
    std::string& next_section,
    const bool run_pseudo_haploid = false,
    const bool rescale_eMatRead_t = true
);



void rcpp_make_eMatGrid_t(
    arma::mat& eMatGrid_t,
    const arma::mat& eMatRead_t,
    const Rcpp::IntegerVector& H,
    const Rcpp::List sampleReads,
    const int hap,
    const int nGrids,
    double& prev,
    int suppressOutput,
    std::string& prev_section,
    std::string& next_section,
    const int run_fb_grid_offset = 0,
    const bool use_all_reads = false,
    const bool bound = false,
    const double maxEmissionMatrixDifference = 1000,
    const bool rescale = false
);


    

Rcpp::List rcpp_make_fb_snp_offsets(
    const arma::mat& alphaHat_t,
    const arma::mat& betaHat_t,
    const arma::mat& blocks_for_output
);


IntegerVector rcpp_order(Rcpp::NumericVector x) {
    //  if (is_true(any(duplicated(x)))) {
    //    Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
    //  }
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}

void Rcpp_run_forward_haploid(
    arma::mat& alphaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const arma::mat& priorCurrent_m,
    const int s,
    const Rcpp::NumericVector alphaStart = 0,
    bool run_fb_subset = false,
    const bool initialize_only = false
);


void Rcpp_run_backward_haploid(
    arma::mat& betaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const int s
);




Rcpp::List Rcpp_define_blocked_snps_using_gamma_on_the_fly(
    arma::mat& alphaHat_t1,
    arma::mat& alphaHat_t2,
    arma::mat& alphaHat_t3,
    arma::mat& betaHat_t1,
    arma::mat& betaHat_t2,
    arma::mat& betaHat_t3,    
    arma::rowvec& c1,
    arma::rowvec& c2,
    arma::rowvec& c3,
    arma::mat& eMatGrid_t1,
    arma::mat& eMatGrid_t2,
    arma::mat& eMatGrid_t3,
    Rcpp::NumericVector& smooth_cm,    
    const arma::cube& transMatRate_tc_H,
    const int shuffle_bin_radius,
    const Rcpp::IntegerVector& L_grid,
    const Rcpp::IntegerVector& grid,
    int s,
    const double block_gibbs_quantile_prob = 0.9,
    const bool verbose = false,
    const bool use_smooth_cm_in_block_gibbs = false
);
                            

Rcpp::List Rcpp_block_gibbs_resampler(
    arma::mat& alphaHat_t1,
    arma::mat& alphaHat_t2,
    arma::mat& alphaHat_t3,
    arma::mat& betaHat_t1,
    arma::mat& betaHat_t2,
    arma::mat& betaHat_t3,    
    arma::rowvec& c1,
    arma::rowvec& c2,
    arma::rowvec& c3,
    arma::mat& eMatGrid_t1,
    arma::mat& eMatGrid_t2,
    arma::mat& eMatGrid_t3,
    Rcpp::IntegerVector& H,
    Rcpp::IntegerVector& H_class,    
    const arma::mat& eMatRead_t,
    Rcpp::IntegerVector& blocked_snps,
    Rcpp::NumericVector& runif_block,    
    Rcpp::NumericVector& runif_total,
    Rcpp::NumericMatrix& runif_proposed,
    const Rcpp::IntegerVector& grid,
    Rcpp::IntegerVector& wif0,
    double ff,
    int s, // this is 0-based
    const arma::cube& alphaMatCurrent_tc,
    const arma::mat& priorCurrent_m,
    const arma::cube& transMatRate_tc_H,
    const int maxDifferenceBetweenReads,
    const int Jmax,
    bool do_checks = false,
    Rcpp::List initial_package = R_NilValue,
    bool verbose = false,
    Rcpp::List fpp_stuff = R_NilValue,
    bool use_cpp_bits_in_R = true,
    int block_approach = 4
);



//' @export
// [[Rcpp::export]]
void rcpp_initialize_gibbs_forward_backward(
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const arma::mat& priorCurrent_m,
    int s,
    arma::mat& alphaHat_t,
    arma::mat& betaHat_t,
    arma::rowvec& c,
    arma::mat& eMatGrid_t,
    const bool run_fb_subset = false,
    const Rcpp::NumericVector alphaStart = 0,
    const Rcpp::NumericVector betaEnd = 0
) {
    //
    //
    const int K = alphaMatCurrent_tc.n_rows;
    const int nGrids = alphaMatCurrent_tc.n_cols + 1;    
    //
    // alphaHat
    //
    const bool initialize_only = false;    
    Rcpp_run_forward_haploid(alphaHat_t, c, eMatGrid_t, alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s, alphaStart, run_fb_subset, initialize_only);
    //
    // betaHat
    //
    if (run_fb_subset == false) {    
        betaHat_t.col(nGrids-1).fill(c(nGrids-1));
    } else {
        for(int k = 0; k < K; k++) {
            betaHat_t(k, nGrids-1) = betaEnd(k);
        }
    }
    Rcpp_run_backward_haploid(betaHat_t, c, eMatGrid_t, alphaMatCurrent_tc, transMatRate_tc_H, s);
    return;
}

//' @export
// [[Rcpp::export]]
void rcpp_calculate_gn_genProbs_and_hapProbs(
    arma::mat& genProbsM_t,
    arma::mat& genProbsF_t,    
    arma::mat& hapProbs_t,
    int s,
    const arma::cube& eHapsCurrent_tc,
    const arma::mat& gammaMT_t,
    const arma::mat& gammaMU_t,
    const arma::mat& gammaP_t,
    const Rcpp::IntegerVector& grid,
    const int snp_start_1_based,
    const int snp_end_1_based,
    const int grid_offset = 0
) {
    // basically, copy and paste, either using grid or not
    const int nSNPs = snp_end_1_based - snp_start_1_based + 1;
    //    const int K = eHapsCurrent_t.n_rows;
    // new
    int iSNP, t, cur_grid;
    double g0, g1, g2;
    arma::colvec eHapsCurrent_t_col, gamma_t_col, one_minus_eHapsCurrent_t_col;
    int prev_grid = -1;
    // i_t is index from 0 to nSNPs + 1, controls where things go
    // t is the index in the whole set of SNPs
    // tt is the index in the grid
    for(iSNP = 0; iSNP < nSNPs; iSNP++) {
        g0 = 0;
        g1 = 0;
        g2 = 0;
        t = iSNP + snp_start_1_based - 1;
        cur_grid = grid(t) - grid_offset;
        if (cur_grid > prev_grid) {
            //gamma_t_col = gamma_t.col(cur_grid);
            prev_grid = cur_grid;
        }
        //eHapsCurrent_t_col = eHapsCurrent_t.col(t);
        //one_minus_eHapsCurrent_t_col = 1 - eHapsCurrent_t_col;
        g0 = sum(gammaMT_t.col(cur_grid) % (eHapsCurrent_tc.slice(s).col(t)));
        g1 = sum(gammaMU_t.col(cur_grid) % (eHapsCurrent_tc.slice(s).col(t)));
        g2 = sum(gammaP_t.col(cur_grid) % (eHapsCurrent_tc.slice(s).col(t)));
        //
        genProbsM_t(0, iSNP) = (1 - g0) * (1 - g1);
        genProbsM_t(1, iSNP) = (g0 * (1 - g1) + (1 - g0) * g1);
        genProbsM_t(2, iSNP) = g0 * g1;
        //
        genProbsF_t(0, iSNP) = (1 - g0) * (1 - g2);
        genProbsF_t(1, iSNP) = (g0 * (1 - g2) + (1 - g0) * g2);
        genProbsF_t(2, iSNP) = g0 * g2;
        //
        hapProbs_t(0, iSNP) = g0;
        hapProbs_t(1, iSNP) = g1;
        hapProbs_t(2, iSNP) = g2;
    }
    return;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_determine_label_probabilities(
    Rcpp::NumericVector rc,
    double ff,
    bool return_neutral = true
) {
    Rcpp::NumericMatrix prob_matrix(3, 3);
    if (return_neutral) {
        prob_matrix.fill(0);
        prob_matrix(0, 0) = 1;
        prob_matrix(1, 1) = 1;
        prob_matrix(2, 2) = 1;
        return(prob_matrix);
    }
    int i;
    Rcpp::NumericVector p(3);
    p(0) = 0.5;
    p(1) = (1 - ff) / 2;
    p(2) = ff / 2;
    Rcpp::NumericVector log_p(3);
    for(i = 0; i < 3; i++) {
        log_p(i) = std::log(p(i));
    }
    Rcpp::IntegerMatrix rc_mat(6, 3);
    rc_mat.row(0) = IntegerVector::create(1, 2, 3);
    rc_mat.row(1) = IntegerVector::create(1, 3, 2);
    rc_mat.row(2) = IntegerVector::create(2, 1, 3);
    rc_mat.row(3) = IntegerVector::create(2, 3, 1);
    rc_mat.row(4) = IntegerVector::create(3, 1, 2);
    rc_mat.row(5) = IntegerVector::create(3, 2, 1);
    //
    Rcpp::NumericVector options(6);
    double max = 0;
    for(i = 0; i < 6; i++) {
        options(i) = \
            rc(rc_mat(i, 0) - 1) * log_p(0) + \
            rc(rc_mat(i, 1) - 1) * log_p(1) + \
            rc(rc_mat(i, 2) - 1) * log_p(2);
        if (i == 0) {
            max = options(i);
        } else {
            if (options(i) > max) {
                max = options(i);
            }
        }
    }
    //
    Rcpp::NumericVector probs(6);
    double prob_sum = 0;
    for(i = 0; i < 6; i++) {    
        options(i) = options(i) - max;
        if (options(i) < (-100)) {
            options(i) = (-100);
        }
        probs(i) = std::exp(options(i));
        prob_sum += probs(i);
    }
    for(i = 0; i < 6; i++) {
        probs(i) /=prob_sum;
    }
    //
    //
    //
    // matrix is probability of that original labelling (row) being that option (column)
    for(int original_label = 0; original_label < 3; original_label++) {
        for(int true_label = 0; true_label < 3; true_label++) {
            for(i = 0; i < 6; i++) {                
                if (rc_mat(i, true_label) == (original_label + 1)) {
                    prob_matrix(original_label, true_label) += probs(i);
                }
            }
        }
    }
    return(prob_matrix);
}
    


    
//' @export
// [[Rcpp::export]]
void rcpp_alpha_forward_one(
    int s,
    const int iGrid,
    const int K,
    arma::mat& alphaHat_t,
    const arma::cube& transMatRate_tc_H,
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrent_tc,
    arma::rowvec& c,
    double& minus_log_c_sum,
    const bool normalize = false
) {
    // t is 1-based
    int iGrid_minus_1 = iGrid - 1;    
    double alphaConst = transMatRate_tc_H(1, iGrid_minus_1, s) * sum(alphaHat_t.col(iGrid_minus_1));
    double a;
    double x = transMatRate_tc_H(0, iGrid_minus_1, s);
    double c2 = c(iGrid);
    for(int k = 0; k < K; k++) {
        alphaHat_t(k, iGrid) = c2 * eMatGrid_t(k, iGrid) *       \
            ( x * alphaHat_t(k, iGrid_minus_1) +            \
              alphaConst * alphaMatCurrent_tc(k, iGrid_minus_1, s));
    }
    //alphaHat_t.col(iGrid) *= c(iGrid);
    if (normalize) {
        a = 1 / sum(alphaHat_t.col(iGrid));
        minus_log_c_sum -= std::log(a);
        c(iGrid) *= a;
        alphaHat_t.col(iGrid) *= a;
    }
    return;
}


//' @export
// [[Rcpp::export]]
void rcpp_reinitialize_in_iterations(
    int& s,
    arma::mat& alphaHat_t,
    arma::rowvec& c,
    const arma::mat& priorCurrent_m,
    const arma::mat& eMatGrid_t,
    const int& K
) {
    //
    for(int k = 0; k < K; k++) {
        alphaHat_t(k, 0) = priorCurrent_m(k, s) * eMatGrid_t(k, 0);
    }
    c(0) = 1 / sum(alphaHat_t.col(0));
    alphaHat_t.col(0) *= c(0);
    return;
}        



//' @export
// [[Rcpp::export]]
void sample_reads_in_grid(
    int& iRead,
    int& iGrid,
    bool& done_reads,
    int& read_wif_iRead,
    const bool& verbose,
    int& nReads,
    Rcpp::NumericVector& pC,
    Rcpp::NumericVector& pA1,
    Rcpp::NumericVector& pA2,
    arma::mat& alphaHat_m,
    arma::mat& betaHat_m,
    arma::mat& ab_m,
    const arma::mat& eMatRead_t,
    Rcpp::NumericVector& runif_reads,
    Rcpp::IntegerVector& H,
    arma::mat& eMatGrid_t1,
    arma::mat& eMatGrid_t2,
    arma::mat& eMatGrid_t3,
    arma::rowvec& c1,
    arma::rowvec& c2,
    arma::rowvec& c3,
    double& minus_log_c1_sum,
    double& minus_log_c2_sum,
    double& minus_log_c3_sum,    
    arma::mat& alphaHat_t1,
    arma::mat& alphaHat_t2,
    arma::mat& alphaHat_t3,
    arma::mat& betaHat_t1,
    arma::mat& betaHat_t2,
    arma::mat& betaHat_t3,    
    const Rcpp::List& sampleReads,
    const bool& return_p_store,
    int& iteration,
    Rcpp::NumericMatrix& p_store,
    const bool record_read_set,
    const Rcpp::NumericMatrix& rlc,
    Rcpp::IntegerVector& H_class,
    const double class_sum_cutoff,
    const int& i_gibbs_samplings,
    const int& n_gibbs_full_its,
    const Rcpp::NumericVector& prior_probs,
    const bool gibbs_initialize_iteratively = false,
    const int first_read_for_gibbs_initialization = 0
) {
    int h_rC = 0;
    int h_rA1 = 1;
    int h_rA2 = 2;
    int h_rN = 0;
    int w, j;
    double prod_pC, prod_pA1, prod_pA2;
    double norm_pC, norm_pA1, norm_pA2;
    Rcpp::NumericVector cumsum_flip_probs(3);
    double denom, chance, alphaConst;
    arma::colvec eMatRead_t_col;
    Rcpp::List readData;
    //
    //
    bool this_grid_has_at_least_one_read = false;
    bool at_least_one_read_has_changed = false;
    bool currently_doing_normal_progress = false;
    bool currently_doing_gibbs_initialization = false;
    bool currently_doing_pass_through = false;
    //
    //
    //
    while((done_reads == false) & ((read_wif_iRead) == iGrid)) {
        //
        // determine type of iteration
        //
        if (!gibbs_initialize_iteratively) {
            currently_doing_normal_progress = true;
        } else {
            if ((iRead < first_read_for_gibbs_initialization) & (iteration == 0)) { // iteration is 0-based here
                currently_doing_pass_through = true;
            } else if ((first_read_for_gibbs_initialization <= iRead) & (iteration == 0)) {
                currently_doing_pass_through = false;
                currently_doing_gibbs_initialization = true;
            } else if ((iRead < first_read_for_gibbs_initialization) & (iteration == 1)) {
                currently_doing_pass_through = false;                
                currently_doing_gibbs_initialization = true;
            } else {
                currently_doing_gibbs_initialization = false;
                currently_doing_normal_progress = true;
            }
        }
        //
        // re-set 
        //
        if (!this_grid_has_at_least_one_read) {
            //
            // build matrices to work with
            //
            pC.fill(0);
            pA1.fill(0);
            pA2.fill(0);
            alphaHat_m.col(0) = alphaHat_t1.col(iGrid);
            alphaHat_m.col(1) = alphaHat_t2.col(iGrid);
            alphaHat_m.col(2) = alphaHat_t3.col(iGrid);                    
            betaHat_m.col(0) = betaHat_t1.col(iGrid);
            betaHat_m.col(1) = betaHat_t2.col(iGrid);
            betaHat_m.col(2) = betaHat_t3.col(iGrid);
            ab_m = alphaHat_m % betaHat_m;
            pC(0) = sum(ab_m.col(0));
            pC(1) = sum(ab_m.col(1));                        
            pC(2) = sum(ab_m.col(2));
            // pC(0) = sum(alphaHat_m.col(0) % betaHat_m.col(0));
            // pC(1) = sum(alphaHat_m.col(1) % betaHat_m.col(1));
            // pC(2) = sum(alphaHat_m.col(2) % betaHat_m.col(2));
            this_grid_has_at_least_one_read = true;
            if (verbose) {
                std::cout << "set pC = " << pC << std::endl;
            }
        }
        //
        eMatRead_t_col = eMatRead_t.col(iRead);
        //
        if (verbose) {
            std::cout << "------------------- it = " << iteration << ", iRead = " << iRead << std::endl;
        }
        //
        //
        // these are made 0-based here
        if (currently_doing_normal_progress) {
            h_rC = H(iRead) - 1; // h_r current
            // ugh just do this manually
            if (h_rC == 0) {
                h_rA1 = 1; h_rA2 = 2;
            } else if (h_rC == 1) {
                h_rA1 = 0; h_rA2 = 2;
            }
            if (h_rC == 2) {
                h_rA1 = 0; h_rA2 = 1;                    
            }
            //
            if (verbose) {
                std::cout << "h_rC = " << h_rC << ", h_rA1 = " << h_rA1 << ", h_rA2 = " << h_rA2 << std::endl;
            }
            // so pC is current three probabilities
            // need same three probabilities for flip 1, flip 2
            pA1(h_rC) = 0;
            pA1(h_rA1) = 0;
            pA1(h_rA2) = pC(h_rA2); // stays the same
            //
            pA2(h_rC) = 0;
            pA2(h_rA1) = pC(h_rA1); // stays the same
            pA2(h_rA2) = 0;
            //
            // A1 - original hap loses
            pA1(h_rC) = sum(ab_m.col(h_rC) / eMatRead_t_col);
            // A2 - original hap loses
            pA2(h_rC) = pA1(h_rC);
            // A1 - new hap gains
            pA1(h_rA1) = sum(ab_m.col(h_rA1) % eMatRead_t_col);
            // A2 - new hap gains
            pA2(h_rA2) = sum(ab_m.col(h_rA2) % eMatRead_t_col);
        } else if (currently_doing_gibbs_initialization) {
            h_rC = 0; // wlog
            h_rA1 = 1;
            h_rA2 = 2;
            for(j=0; j < 3; j++) {
                pA1(j) = pC(j);
                pA2(j) = pC(j);                
            }
            // now, over-lay pC
            pC(h_rC) = sum(ab_m.col(h_rC) % eMatRead_t_col);
            pA1(h_rA1) = sum(ab_m.col(h_rA1) % eMatRead_t_col);
            pA2(h_rA2) = sum(ab_m.col(h_rA2) % eMatRead_t_col);
        } else {
            for(j=0; j < 3; j++) {
                pA1(j) = pC(j);
                pA2(j) = pC(j);                
            }
        }
        //
        // I think pC and pA1 and pA2 should always be in log space
        //
        // double log_prod_pC =   log(pC(0))  + log(pC(1))  + log(pC(2)) +  log(prior_probs(h_rC));
        // double log_prod_pA1 =  log(pA1(0)) + log(pA1(1)) + log(pA1(2)) + log(prior_probs(h_rA1));
        // double log_prod_pA2 =  log(pA2(0)) + log(pA2(1)) + log(pA2(2)) + log(prior_probs(h_rA2));
        // double max = std::max({log_prod_pC, log_prod_pA1, log_prod_pA2});
        // log_prod_pC -= max;
        // log_prod_pA1 -= max;
        // log_prod_pA2 -= max;
        // prod_pC=std::exp(log_prod_pC);
        // prod_pA1=std::exp(log_prod_pA1);
        // prod_pA2=std::exp(log_prod_pA2);
        //
        //
        prod_pC =  (pC(0) * pC(1) * pC(2)   ) * prior_probs(h_rC);
        prod_pA1 = (pA1(0) * pA1(1) * pA1(2)) * prior_probs(h_rA1);
        prod_pA2 = (pA2(0) * pA2(1) * pA2(2)) * prior_probs(h_rA2);
        //
        denom = prod_pC + prod_pA1 + prod_pA2;
        norm_pC = prod_pC / denom;
        norm_pA1 = prod_pA1 / denom;
        norm_pA2 = prod_pA2 / denom;
        if (verbose) {
            std::cout << "pC(0)=" << pC(0) << ", pC(1)=" << pC(1) << ", pC(2)=" << pC(2) << std::endl;
            std::cout << "pA1(0)=" << pA1(0) << ", pA1(1)=" << pA1(1) << ", pA1(2)=" << pA1(2) << std::endl;
            std::cout << "pA2(0)=" << pA2(0) << ", pA2(1)=" << pA2(1) << ", pA2(2)=" << pA2(2) << std::endl;            
        }
        //
        denom = prod_pC + prod_pA1 + prod_pA2;
        norm_pC = prod_pC / denom;
        norm_pA1 = prod_pA1 / denom;
        norm_pA2 = prod_pA2 / denom;
        //
        // mt = 0            -> 0.5
        // mu = 0.5          -> 0.5 + ff / 2
        // p  = 0.5 + ff / 2 -> 1
        //
        // which one is the winner - again, do this manually
        //
        // right - now need to choose one
        chance = runif_reads(nReads * (iteration) + iRead);                
        cumsum_flip_probs.fill(0);
        cumsum_flip_probs(h_rC) = norm_pC;
        cumsum_flip_probs(h_rA1) = norm_pA1;
        cumsum_flip_probs(h_rA2) = norm_pA2;
        //
        cumsum_flip_probs(1) += cumsum_flip_probs(0);
        cumsum_flip_probs(2) += cumsum_flip_probs(1);
        //
        h_rN = 0;
        for(int i = 2; i >= 0; i--) {
            if (chance < cumsum_flip_probs(i)) {
                h_rN = i;
            }
        }
        //
        if (verbose) {
            std::cout << "pC = " << pC << std::endl;
            std::cout << "pA1 = " << pA1 << std::endl;
            std::cout << "pA2 = " << pA2 << std::endl;                        
            std::cout << "denom = " << denom << std::endl;
        }
        //
        //
        //
        if (
            ((h_rN != h_rC) | currently_doing_gibbs_initialization) & (!currently_doing_pass_through)
        ) {
            //
            at_least_one_read_has_changed = true;
            H(iRead) = h_rN + 1; // store as 1-based                    
            // update alphas
            if (currently_doing_normal_progress) {
                alphaHat_m.col(h_rC) /= eMatRead_t_col;
                ab_m.col(h_rC) /= eMatRead_t_col;                
            }
            alphaHat_m.col(h_rN) %= eMatRead_t_col;
            // also, have to update ab
            ab_m.col(h_rN) %= eMatRead_t_col;
            //
            if (currently_doing_normal_progress) {            
                if (h_rC == 0) {eMatGrid_t1.col(iGrid) /= eMatRead_t_col;}
                if (h_rC == 1) {eMatGrid_t2.col(iGrid) /= eMatRead_t_col;}
                if (h_rC == 2) {eMatGrid_t3.col(iGrid) /= eMatRead_t_col;}
            }
            //
            if (h_rN == 0) {eMatGrid_t1.col(iGrid) %= eMatRead_t_col;}
            if (h_rN == 1) {eMatGrid_t2.col(iGrid) %= eMatRead_t_col;}
            if (h_rN == 2) {eMatGrid_t3.col(iGrid) %= eMatRead_t_col;}
            // reset pC - UGGGGGH - I wish I could name these (can I?)
            // see above where I define these
            // again this is terribble uggggggh
            if (currently_doing_normal_progress) {
                for(int i = 0; i < 3; i++) {
                    if (h_rC == 0) {
                        if (h_rN == 1) {pC(i) = pA1(i);}
                        if (h_rN == 2) {pC(i) = pA2(i);}
                    } else if (h_rC == 1) {
                        if (h_rN == 0) {pC(i) = pA1(i);}
                        if (h_rN == 2) {pC(i) = pA2(i);}
                    } else if (h_rC == 2) {
                        if (h_rN == 0) {pC(i) = pA1(i);}
                        if (h_rN == 1) {pC(i) = pA2(i);}
                    }
                }
            } else if (currently_doing_gibbs_initialization) {            
                //otherwise, if gibbs_initialize
                // reset no matter what, based on new sampled read
                for(int i = 0; i < 3; i++) {
                    // if h_rN is 0, we're good already
                    if (h_rN == 1) {pC(i) = pA1(i);}
                    if (h_rN == 2) {pC(i) = pA2(i);}
                }
            }
            //if (iRead <= 3) {
            //    std::cout << "the new pC is " << pC << std::endl;
            // }
        }
        //
        if (record_read_set) {
            Rcpp::NumericVector x(3);
            x(h_rC) = norm_pC;
            x(h_rA1) = norm_pA1;
            x(h_rA2) = norm_pA2;
            Rcpp::NumericVector y(7);
            double local_min = 2;
            int local_min_which = 8;
            for(int i = 0; i < 7; i++) {
                y(i) = std::abs(rlc(i, 0) - x(0)) + std::abs(rlc(i, 1) - x(1)) + std::abs(rlc(i, 2) - x(2));
                if (y(i) < local_min) {
                    local_min = y(i);
                    local_min_which = i;
                }
            }
            if (local_min < class_sum_cutoff) {
                H_class(iRead) = local_min_which + 1; // this is 0 for no match, 1 through 7 otherwise!
            } else {
                H_class(iRead) = 0;
            }
        }
        if (return_p_store) {
            // p_store_cols <- c("p_1", "p_2", "p_3", "chance", "h_rC", "h_rN", "c1", "c2", "c3", "p", "agreePer")
            w = \
                (i_gibbs_samplings) * n_gibbs_full_its * nReads + \
                nReads * (iteration) + iRead;
            // argh - can I access these by names? I think not?
            p_store(w, h_rC) = norm_pC;
            p_store(w, h_rA1) = norm_pA1;
            p_store(w, h_rA2) = norm_pA2;
            p_store(w, 3) = chance;
            p_store(w, 4) = h_rC + 1;
            p_store(w, 5) = h_rN + 1;
            p_store(w, 6) = 0; // agree percentage, requires true H, not done here
            if (currently_doing_normal_progress) {
                p_store(w, 7) = 2;
            } else if (currently_doing_gibbs_initialization) {
                p_store(w, 7) = 1;
            } else if (currently_doing_pass_through) {
                p_store(w, 7) = 0;
            }
            //
            if ((minus_log_c1_sum + std::log(sum(alphaHat_m.col(0)))) > 1) {
                 std::cout << "minus_log_c1_sum=" << minus_log_c1_sum << std::endl;
                 std::cout << "alphaHat_m.col(0) = " << alphaHat_m.col(0) << std::endl;                
                 std::cout << "std::log(sum(alphaHat_m.col(0))) = " << std::log(sum(alphaHat_m.col(0))) << std::endl;
            }
            p_store(w, 8) = (minus_log_c1_sum + std::log(sum(alphaHat_m.col(0))));
            p_store(w, 9) = (minus_log_c2_sum + std::log(sum(alphaHat_m.col(1))));
            p_store(w, 10) = (minus_log_c3_sum + std::log(sum(alphaHat_m.col(2))));
            p_store(w, 11) = (p_store(w, 8) + p_store(w, 9) + p_store(w, 10));
            // argh - fix in a minute
            double d = 0;
            for(int i = 0; i < nReads; i++) {
                d += std::log(prior_probs(H(i) - 1));
            }
            p_store(w, 12) = d; // argh - cannot figure out how to get/use character col names
            p_store(w, 13) = p_store(w, 11) + p_store(w, 12);
            //
            //p_store[w, "agreePer"] <- sum(abs(H ==  true_H)) / length(H) * 100
        }
        //
        iRead++;
        if ((nReads - 1) < iRead) {
            done_reads = true;
            read_wif_iRead = -1;
        } else {
            // update
            readData = as<Rcpp::List>(sampleReads[iRead]);
            read_wif_iRead = as<int>(readData[1]);
        }
    }
    //
    if (at_least_one_read_has_changed) {
        // arguably, only need to know which are lost, which are gained
        // can re-build here? but again, this only matters if lots of reads at same spot
        // only do if at least one of them changed?
        alphaHat_t1.col(iGrid) = alphaHat_m.col(0);
        alphaHat_t2.col(iGrid) = alphaHat_m.col(1);
        alphaHat_t3.col(iGrid) = alphaHat_m.col(2);
        //
        alphaConst = 1 / sum(alphaHat_m.col(0));
        c1(iGrid) *= alphaConst;
        minus_log_c1_sum -= std::log(alphaConst);
        alphaHat_t1.col(iGrid) *= alphaConst;
        //
        alphaConst = 1 / sum(alphaHat_m.col(1));
        c2(iGrid) *= alphaConst;
        minus_log_c2_sum -= std::log(alphaConst);
        alphaHat_t2.col(iGrid) *= alphaConst;
        //
        alphaConst = 1 / sum(alphaHat_m.col(2));
        c3(iGrid) *= alphaConst;
        minus_log_c3_sum -= std::log(alphaConst);
        alphaHat_t3.col(iGrid) *= alphaConst;
    }
    //
    return;
}


// note - all are passed by reference except test1v for either array or matrix in R
// array and matrix in R do not make a difference

// //' @export
// // [[Rcpp::export]]
// void test1r(
//     arma::mat& m1
// ) {
//     m1(0, 0) += 1.2;
//     return;
// }

// //' @export
// // [[Rcpp::export]]
// void test1v(
//     arma::mat m1
// ) {
//     m1(0, 0) += 1.2;
//     return;
// }

// //' @export
// // [[Rcpp::export]]
// void test2r(
//     Rcpp::NumericMatrix& m1
// ) {
//     m1(0, 0) += 1.2;
//     return;
// }

// //' @export
// // [[Rcpp::export]]
// void test2v(
//     Rcpp::NumericMatrix m1
// ) {
//     m1(0, 0) += 1.2;
//     return;
// }




//' @export
// [[Rcpp::export]]
void rcpp_apply_mat_relabel(
    arma::mat& m1,
    arma::mat& m2,
    arma::mat& m3,
    const int relabel
) {
    int iGrid;
    const int n = m1.n_cols;
    arma::colvec m_col;
    if (relabel == 2) {
        // 1, 3, 2  swap 2 and 3
        for(iGrid = 0; iGrid < n; iGrid++) {
            m_col = m2.col(iGrid) + 0.0;
            m2.col(iGrid) = m3.col(iGrid) + 0.0;
            m3.col(iGrid) = m_col + 0.0;
        }
    } else if (relabel == 3) {
        // 2 1 3   swap 1 and 2
        for(iGrid = 0; iGrid < n; iGrid++) {
            m_col = m2.col(iGrid);
            m2.col(iGrid) = m1.col(iGrid);
            m1.col(iGrid) = m_col;
        }
    } else if (relabel == 4) {
        // reorderX <- c(3, 1, 2)        
        for(iGrid = 0; iGrid < n; iGrid++) {
            m_col = m1.col(iGrid);
            m1.col(iGrid) = m2.col(iGrid);
            m2.col(iGrid) = m3.col(iGrid);            
            m3.col(iGrid) = m_col;
        }
    } else if (relabel == 5) {
        for(iGrid = 0; iGrid < n; iGrid++) {
            // reorderX <- c(2, 3, 1)
            m_col = m1.col(iGrid);
            m1.col(iGrid) = m3.col(iGrid);
            m3.col(iGrid) = m2.col(iGrid);            
            m2.col(iGrid) = m_col;
        }
    } else if (relabel == 6) {
        for(iGrid = 0; iGrid < n; iGrid++) {
            // reorder <- c(3, 2, 1)
            m_col = m3.col(iGrid);
            m3.col(iGrid) = m1.col(iGrid);
            m1.col(iGrid) = m_col;
        }
    }
    return;
}

//' @export
// [[Rcpp::export]]
void rcpp_apply_vec_relabel(
    arma::rowvec& m1,
    arma::rowvec& m2,
    arma::rowvec& m3,
    const int relabel
) {
    int iGrid;
    const int n = m1.n_cols;
    double m_col; // not really a col per-se
    if (relabel == 2) {
        // reorder <- c(1, 3, 2)
        // swap 2 and 3
        for(iGrid = 0; iGrid < n; iGrid++) {
            m_col = m2(iGrid);
            m2(iGrid) = m3(iGrid);
            m3(iGrid) = m_col;
        }
    } else if (relabel == 3) {
        for(iGrid = 0; iGrid < n; iGrid++) {
            m_col = m2(iGrid);
            m2(iGrid) = m1(iGrid);
            m1(iGrid) = m_col;
        }
    } else if (relabel == 4) {
        for(iGrid = 0; iGrid < n; iGrid++) {
            // reorderX <- c(3, 1, 2)
            m_col = m1(iGrid);
            m1(iGrid) = m2(iGrid);
            m2(iGrid) = m3(iGrid);            
            m3(iGrid) = m_col;
        }
    } else if (relabel == 5) {
        for(iGrid = 0; iGrid < n; iGrid++) {
            // reorderX <- c(2, 3, 1)
            m_col = m1(iGrid);
            m1(iGrid) = m3(iGrid);
            m3(iGrid) = m2(iGrid);            
            m2(iGrid) = m_col;
        }
    } else if (relabel == 6) {
        for(iGrid = 0; iGrid < n; iGrid++) {
            // reorder <- c(3, 2, 1)
            m_col = m3(iGrid);
            m3(iGrid) = m1(iGrid);
            m1(iGrid) = m_col;
        }
    }
    return;    
}

Rcpp::NumericVector calculate_rc(Rcpp::IntegerVector& H) {
    Rcpp::NumericVector rc(4); // cheeky - ignore 0th entry
    rc.fill(0);
    for(int iRead = 0; iRead < H.length(); iRead++) {
        rc(H(iRead)) += 1;
    }
    Rcpp::NumericVector rc_better(3);
    rc_better(0) = rc(1);
    rc_better(1) = rc(2);
    rc_better(2) = rc(3);    
    return(rc_better);
}

//' @export
// [[Rcpp::export]]
double rcpp_calc_prob_of_set_of_reads(
    const double ff,
    Rcpp::NumericVector rc,
    Rcpp::IntegerVector reorder = Rcpp::IntegerVector::create(0, 1, 2)
) {
    const Rcpp::NumericVector prior_probs = NumericVector::create(0.5, (1 - ff) / 2, (ff / 2));
    int n = rc(0) + rc(1) + rc(2);
    double r = std::lgamma(1.0 * (n + 1.0));
    for(int i=0; i < 3; i++) {
        if (prior_probs(i) > 0) {
            r += rc(reorder(i)) * std::log(prior_probs(i)) - std::lgamma(1.0 * (rc(reorder(i)) + 1.0));
        }
    }
    return(r);
}



//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calculate_likelihoods_values(
    arma::rowvec& c1,
    arma::rowvec& c2,
    arma::rowvec& c3,    
    Rcpp::IntegerVector& H,
    const int nGrids,
    const Rcpp::NumericVector prior_probs,
    const double ff
) {
    Rcpp::NumericVector to_out(7);
    double d1 = 0;
    double d2 = 0;
    double d3 = 0;
    double dH = 0;
    //
    for(int iGrid = 0; iGrid < nGrids; iGrid++) {
        d1 -= std::log(c1(iGrid));
        d2 -= std::log(c2(iGrid));
        d3 -= std::log(c3(iGrid));
    }
    //
    for(int iRead = 0; iRead < H.length(); iRead++) {
        dH += std::log(prior_probs(H(iRead) - 1));
    }
    //
    to_out(0) = d1;
    to_out(1) = d2;
    to_out(2) = d3;
    to_out(3) = d1 + d2 + d3;
    //
    to_out(4) = dH;
    to_out(5) = to_out(3) + to_out(4);
    Rcpp::NumericVector rc = calculate_rc(H);
    to_out(6) = rcpp_calc_prob_of_set_of_reads(ff, rc);
    return(to_out);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_get_weights_for_entire_relabelling(
    Rcpp::NumericVector rc,
    const double ff
) {
    //
    Rcpp::NumericVector reorder_log_probs = Rcpp::NumericVector::create(
        rcpp_calc_prob_of_set_of_reads(ff, rc, Rcpp::IntegerVector::create(0, 1, 2)),
        rcpp_calc_prob_of_set_of_reads(ff, rc, Rcpp::IntegerVector::create(0, 2, 1)),
        rcpp_calc_prob_of_set_of_reads(ff, rc, Rcpp::IntegerVector::create(1, 0, 2)),
        rcpp_calc_prob_of_set_of_reads(ff, rc, Rcpp::IntegerVector::create(1, 2, 0)),
        rcpp_calc_prob_of_set_of_reads(ff, rc, Rcpp::IntegerVector::create(2, 0, 1)),
        rcpp_calc_prob_of_set_of_reads(ff, rc, Rcpp::IntegerVector::create(2, 1, 0))
    );
    double maxval = Rcpp::max(reorder_log_probs);
    Rcpp::NumericVector weights(6);
    double sumval = 0;
    for(int i=0; i < 6; i++) {
        weights(i) = std::exp(reorder_log_probs(i) - maxval);
        sumval += weights(i);
    }
    for(int i=0; i < 6; i++) {
        weights(i) /= sumval;
    }
    return(weights);
}



//' @export
// [[Rcpp::export]]
int rcpp_consider_and_try_entire_relabelling(
    Rcpp::IntegerVector& H,
    const double ff,
    int relabel = -1
)  {
    Rcpp::NumericVector rc = calculate_rc(H);
    Rcpp::NumericVector weights = rcpp_get_weights_for_entire_relabelling(rc, ff);
    if (relabel == -1) {
        relabel = Rcpp::sample(6, 1, false, weights)(0);
    }
    Rcpp::IntegerVector reorderX;
    if (relabel > 1) {
        //
        if (relabel == 1) {reorderX = Rcpp::IntegerVector::create(0, 1, 2);}
        if (relabel == 2) {reorderX = Rcpp::IntegerVector::create(0, 2, 1);}
        if (relabel == 3) {reorderX = Rcpp::IntegerVector::create(1, 0, 2);}
        if (relabel == 4) {reorderX = Rcpp::IntegerVector::create(2, 0, 1);}
        if (relabel == 5) {reorderX = Rcpp::IntegerVector::create(1, 2, 0);}
        if (relabel == 6) {reorderX = Rcpp::IntegerVector::create(2, 1, 0);}
        for(int iRead = 0; iRead < H.length(); iRead ++ ) {
            H(iRead) = reorderX(H(iRead) - 1) + 1;
        }
    }
    return(relabel);
}





void add_to_per_it_likelihoods(
    int s,
    Rcpp::NumericMatrix& per_it_likelihoods,
    int i_gibbs_samplings,
    int iteration,
    int i_result_it,
    const int n_gibbs_full_its,
    Rcpp::IntegerVector& H,
    arma::rowvec& c1,
    arma::rowvec& c2,
    arma::rowvec& c3,
    const double ff,
    const Rcpp::NumericVector prior_probs,
    int relabel,
    int& i_per_it_likelihoods
) {
    // note - only ever done when not run_fb_subset
    int n = per_it_likelihoods.nrow();
    if (i_per_it_likelihoods > n) {
         std::cout << "PROBLEM - out of bounds on per_it_likelihoods" << std::endl;
         std::cout << "requested entry into 0-based row " << i_per_it_likelihoods << " but there are " << n << "rows" << std::endl;
         return;
    }
    const int nGrids = c1.n_cols;
    Rcpp::NumericVector temp = calculate_likelihoods_values(c1, c2, c3, H, nGrids, prior_probs, ff);
    per_it_likelihoods(i_per_it_likelihoods, 0) = s + 1; // 1-based
    per_it_likelihoods(i_per_it_likelihoods, 1) = i_gibbs_samplings + 1;
    per_it_likelihoods(i_per_it_likelihoods, 2) = iteration + 1;
    per_it_likelihoods(i_per_it_likelihoods, 3) = i_result_it + 1; // make these 1-based
    for(int j = 0; j < 7; j++) {
        per_it_likelihoods(i_per_it_likelihoods, 4 + j) = temp(j);
    }
    per_it_likelihoods(i_per_it_likelihoods, 11) = relabel;
    i_per_it_likelihoods++; // bump counter
    return;
};





//' @export
// [[Rcpp::export]]
void rcpp_gibbs_nipt_initialize(
    int s,
    std::string& prev_section,
    std::string& next_section,    
    const int suppressOutput,
    double& prev,
    const bool run_fb_subset,
    const Rcpp::List& alphaBetaBlocks_one,
    const int i_snp_block_for_alpha_beta,
    arma::mat& eMatRead_t,
    const Rcpp::List& sampleReads,
    Rcpp::IntegerVector& H,
    const int run_fb_grid_offset,
    const bool bound_eMatGrid_t,
    const bool rescale_eMatGrid_t,
    arma::mat& alphaHat_t1,
    arma::mat& betaHat_t1,
    arma::rowvec& c1,
    arma::mat& eMatGrid_t1,
    arma::mat& alphaHat_t2,
    arma::mat& betaHat_t2,
    arma::rowvec& c2,
    arma::mat& eMatGrid_t2,
    arma::mat& alphaHat_t3,
    arma::mat& betaHat_t3,
    arma::rowvec& c3,
    arma::mat& eMatGrid_t3,
    const arma::cube& transMatRate_tc_H,
    const bool gibbs_initialize_iteratively,
    const arma::mat& priorCurrent_m,
    const arma::cube& alphaMatCurrent_tc,
    const double maxEmissionMatrixDifference
) {
    int nGrids = alphaHat_t1.n_cols;
    Rcpp::NumericVector alphaStart1, betaEnd1, alphaStart2, betaEnd2, alphaStart3, betaEnd3;
    Rcpp::List alphaBetaBlocks1, alphaBetaBlocks2, alphaBetaBlocks3;
    //
    //
    //
    //
    //
    if (run_fb_subset) {
        next_section="Initialize from previous";
        prev=print_times(prev, suppressOutput, prev_section, next_section);
        prev_section=next_section;
        //
        alphaBetaBlocks1 = as<Rcpp::List>(alphaBetaBlocks_one["alphaBetaBlocks1"]);
        alphaStart1 = (as<Rcpp::NumericMatrix>(alphaBetaBlocks1["alphaHatBlocks_t"])).column(i_snp_block_for_alpha_beta - 1);
        betaEnd1 = (as<Rcpp::NumericMatrix>(alphaBetaBlocks1["betaHatBlocks_t"])).column(i_snp_block_for_alpha_beta - 1);
        //
        alphaBetaBlocks2 = as<Rcpp::List>(alphaBetaBlocks_one["alphaBetaBlocks2"]);
        alphaStart2 = (as<Rcpp::NumericMatrix>(alphaBetaBlocks2["alphaHatBlocks_t"])).column(i_snp_block_for_alpha_beta - 1);
        betaEnd2 = (as<Rcpp::NumericMatrix>(alphaBetaBlocks2["betaHatBlocks_t"])).column(i_snp_block_for_alpha_beta - 1);
        //
        alphaBetaBlocks3 = as<Rcpp::List>(alphaBetaBlocks_one["alphaBetaBlocks3"]);
        alphaStart3 = (as<Rcpp::NumericMatrix>(alphaBetaBlocks3["alphaHatBlocks_t"])).column(i_snp_block_for_alpha_beta - 1);
        betaEnd3 = (as<Rcpp::NumericMatrix>(alphaBetaBlocks3["betaHatBlocks_t"])).column(i_snp_block_for_alpha_beta - 1);
    }
    // 
    //
    //
    //
    next_section="Initialize eMatGrid_t";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    eMatGrid_t1.fill(1);
    eMatGrid_t2.fill(1);
    eMatGrid_t3.fill(1);
    if (!gibbs_initialize_iteratively) {
        const bool use_all_reads = false;
        rcpp_make_eMatGrid_t(eMatGrid_t1, eMatRead_t, H, sampleReads, 1, nGrids, prev, suppressOutput, prev_section, next_section, run_fb_grid_offset, use_all_reads, bound_eMatGrid_t, maxEmissionMatrixDifference, rescale_eMatGrid_t);
        rcpp_make_eMatGrid_t(eMatGrid_t2, eMatRead_t, H, sampleReads, 2, nGrids, prev, suppressOutput, prev_section, next_section, run_fb_grid_offset, use_all_reads, bound_eMatGrid_t, maxEmissionMatrixDifference, rescale_eMatGrid_t);
        rcpp_make_eMatGrid_t(eMatGrid_t3, eMatRead_t, H, sampleReads, 3, nGrids, prev, suppressOutput, prev_section, next_section, run_fb_grid_offset, use_all_reads, bound_eMatGrid_t, maxEmissionMatrixDifference, rescale_eMatGrid_t);
    }
    //
    // initialize forward backward
    //
    //
    next_section="Initialize forward backward";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    if (gibbs_initialize_iteratively) {
        // name these local variables so it is easier to read against function if function arguments change
        const bool initialize_only_local = true;
        const Rcpp::NumericVector alphaStart_local = 0;
        const bool run_fb_subset_local = false;
        //
        alphaHat_t1.fill(1); betaHat_t1.fill(1); c1.fill(1);
        Rcpp_run_forward_haploid(alphaHat_t1, c1, eMatGrid_t1, alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s, alphaStart_local, run_fb_subset_local, initialize_only_local);                
        //
        alphaHat_t2.fill(1); betaHat_t2.fill(1); c2.fill(1);
        Rcpp_run_forward_haploid(alphaHat_t2, c2, eMatGrid_t2, alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s, alphaStart_local, run_fb_subset_local, initialize_only_local);        
        //
        alphaHat_t3.fill(1); betaHat_t3.fill(1); c3.fill(1);        
        Rcpp_run_forward_haploid(alphaHat_t3, c3, eMatGrid_t3, alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s, alphaStart_local, run_fb_subset_local, initialize_only_local);
    } else {
        // recall s needs to be 0based in c++
        rcpp_initialize_gibbs_forward_backward(alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s, alphaHat_t1, betaHat_t1, c1, eMatGrid_t1, run_fb_subset, alphaStart1, betaEnd1);
        rcpp_initialize_gibbs_forward_backward(alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s, alphaHat_t2, betaHat_t2, c2, eMatGrid_t2, run_fb_subset, alphaStart2, betaEnd2);
        rcpp_initialize_gibbs_forward_backward(alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s, alphaHat_t3, betaHat_t3, c3, eMatGrid_t3, run_fb_subset, alphaStart3, betaEnd3);
    }
    return;
}



//' @export
// [[Rcpp::export]]
void rcpp_gibbs_nipt_iterate(
    int s,
    std::string& prev_section,
    std::string& next_section,    
    const int suppressOutput,
    double& prev,
    int iteration,
    const Rcpp::List& sampleReads,
    const arma::mat& priorCurrent_m,
    const arma::cube& transMatRate_tc_H,
    const arma::cube& alphaMatCurrent_tc,
    const int n_gibbs_full_its,
    const int nGrids,
    const int K,
    Rcpp::IntegerVector& H,
    const arma::mat& eMatRead_t,
    Rcpp::NumericVector& runif_reads,
    arma::mat& alphaHat_t1,
    arma::mat& betaHat_t1,
    arma::rowvec& c1,
    arma::mat& eMatGrid_t1,
    arma::mat& alphaHat_t2,
    arma::mat& betaHat_t2,
    arma::rowvec& c2,
    arma::mat& eMatGrid_t2,
    arma::mat& alphaHat_t3,
    arma::mat& betaHat_t3,
    arma::rowvec& c3,
    arma::mat& eMatGrid_t3,
    Rcpp::NumericMatrix& p_store,
    int& i_per_it_likelihoods,
    const bool verbose,
    const bool return_p_store,
    const bool run_fb_subset,
    const int i_gibbs_samplings,
    const int n_gibbs_starts,
    int i_result_it,
    const double ff,
    const bool record_read_set,
    const Rcpp::NumericMatrix& rlc,
    Rcpp::IntegerVector& H_class,
    const double class_sum_cutoff,
    const int run_fb_grid_offset,
    const Rcpp::NumericVector prior_probs,
    Rcpp::NumericMatrix & per_it_likelihoods,
    const bool gibbs_initialize_iteratively = false,
    const int first_read_for_gibbs_initialization = 0,
    const bool do_block_resampling = false,
    const int artificial_relabel = -1
) {
    //
    //
    //
    next_section="Begin an iteration";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    bool done_reads = false;
    int iRead = -1; // this is 0-based
    int iGrid = 0; // here this is 0-based
    double minus_log_c1_sum, minus_log_c2_sum, minus_log_c3_sum;
    int read_wif_iRead = -1;
    int relabel = 1;
    int iReadStart;
    int nReads = sampleReads.size();
    Rcpp::List readData;    
    arma::mat alphaHat_m(K, 3); // tested, this is the better orientation
    arma::mat betaHat_m(K, 3);
    arma::mat ab_m(K, 3);    
    Rcpp::NumericVector pC(3), pA1(3), pA2(3);
    //
    if (return_p_store | (relabel > 1)) {
        minus_log_c1_sum = 0; minus_log_c2_sum = 0; minus_log_c3_sum = 0;
        for(int iGridLocal = 0; iGridLocal < nGrids; iGridLocal++) {
            minus_log_c1_sum += -std::log(c1(iGridLocal));
            minus_log_c2_sum += -std::log(c2(iGridLocal));
            minus_log_c3_sum += -std::log(c3(iGridLocal));                
        }
    }
    for(iGrid = 0; iGrid < nGrids; iGrid++) {        
        if (verbose) {
            std::cout << "iGrid = " << iGrid << std::endl;
        }
        //
        if (iGrid > 0) {
            // move and normalize
            rcpp_alpha_forward_one(s, iGrid, K, alphaHat_t1, transMatRate_tc_H, eMatGrid_t1, alphaMatCurrent_tc, c1, minus_log_c1_sum, true);
            rcpp_alpha_forward_one(s, iGrid, K, alphaHat_t2, transMatRate_tc_H, eMatGrid_t2, alphaMatCurrent_tc, c2, minus_log_c2_sum, true);
            rcpp_alpha_forward_one(s, iGrid, K, alphaHat_t3, transMatRate_tc_H, eMatGrid_t3, alphaMatCurrent_tc, c3, minus_log_c3_sum, true);
        } else {
            rcpp_reinitialize_in_iterations(s, alphaHat_t1, c1, priorCurrent_m, eMatGrid_t1, K);
            rcpp_reinitialize_in_iterations(s, alphaHat_t2, c2, priorCurrent_m, eMatGrid_t2, K);
            rcpp_reinitialize_in_iterations(s, alphaHat_t3, c3, priorCurrent_m, eMatGrid_t3, K);            
        }
        //
        iRead++;
        iReadStart = iRead;
        if (!done_reads) {
            readData = as<Rcpp::List>(sampleReads[iRead]);
            read_wif_iRead = as<int>(readData[1]); // leading SNP from read
        } else {
            read_wif_iRead = -1;
        }
        if (verbose) {
            std::cout << "outside. read_wif_iRead = " << read_wif_iRead << std::endl;
        }
        //
        // for this grid
        // sample all the reads, one at a time
        // will be modifying alphaHat_m, betaHat_m, for re-injection below
        if (read_wif_iRead == iGrid) { // only start if makes sense!
            sample_reads_in_grid(
                iRead, iGrid, done_reads, read_wif_iRead,
                verbose, nReads, pC, pA1, pA2,
                alphaHat_m, betaHat_m, ab_m, eMatRead_t, runif_reads, H,
                eMatGrid_t1, eMatGrid_t2, eMatGrid_t3,
                c1, c2, c3,
                minus_log_c1_sum, minus_log_c2_sum, minus_log_c3_sum,
                alphaHat_t1, alphaHat_t2, alphaHat_t3,
                betaHat_t1, betaHat_t2, betaHat_t3,
                sampleReads, return_p_store, iteration, p_store,
                record_read_set, rlc, H_class, class_sum_cutoff,
                i_gibbs_samplings, n_gibbs_full_its, prior_probs,
                gibbs_initialize_iteratively, first_read_for_gibbs_initialization
            );
        }
        //
        iRead = iRead - 1;
        //
    }
    next_section="Finalize it use backward";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    betaHat_t1.col(nGrids - 1).fill(c1(nGrids-1));
    betaHat_t2.col(nGrids - 1).fill(c2(nGrids-1));
    betaHat_t3.col(nGrids - 1).fill(c3(nGrids-1));
    //
    Rcpp_run_backward_haploid(betaHat_t1, c1, eMatGrid_t1, alphaMatCurrent_tc, transMatRate_tc_H, s);
    Rcpp_run_backward_haploid(betaHat_t2, c2, eMatGrid_t2, alphaMatCurrent_tc, transMatRate_tc_H, s);
    Rcpp_run_backward_haploid(betaHat_t3, c3, eMatGrid_t3, alphaMatCurrent_tc, transMatRate_tc_H, s);
    //
    if (do_block_resampling) {
        next_section="Do block resampling";
        prev=print_times(prev, suppressOutput, prev_section, next_section);
        prev_section=next_section;
        Rcpp::IntegerVector cur_H(nReads);
        for(iRead = 0; iRead < nReads; iRead++) {
            cur_H(iRead) = H(iRead) + 0;
        }
        // artificial relabel is for testing purposes
        relabel = rcpp_consider_and_try_entire_relabelling(H, ff, artificial_relabel); // also will re-label H on the fly
        if (relabel > 1) { // re-label still 1-based, i.e. ranges from 1-6
            rcpp_apply_mat_relabel(alphaHat_t1, alphaHat_t2, alphaHat_t3, relabel);
            rcpp_apply_mat_relabel(betaHat_t1, betaHat_t2, betaHat_t3, relabel);
            rcpp_apply_vec_relabel(c1, c2, c3, relabel);
            rcpp_apply_mat_relabel(eMatGrid_t1, eMatGrid_t2, eMatGrid_t3, relabel);
        }
    }
    add_to_per_it_likelihoods(s, per_it_likelihoods, i_gibbs_samplings, iteration, i_result_it, n_gibbs_full_its, H, c1, c2, c3, ff, prior_probs, relabel, i_per_it_likelihoods);
    return;
}




Rcpp::IntegerVector random_gibbs_nipt_read_labels(
    const int nReads,
    const double ff
) {
    Rcpp::NumericVector x = Rcpp::runif(nReads);
    Rcpp::IntegerVector H(nReads);
    int iRead;
    for(iRead = 0; iRead < nReads; iRead++) {
        if (x(iRead) < 0.5) {
            H(iRead) = 1;
        } else if ((0.5 <= x(iRead)) & (x(iRead) < (0.5 + ff / 2))) {
            H(iRead) = 2;
        } else {
            H(iRead) = 3;
        }
    }
    return(H);
}






//
// old way of determining which haplotype is which
//
// 
// if ((rc(2) <= rc(0)) & (rc(1) <= rc(0))) { mt = 0;}
// if ((rc(2) <= rc(1)) & (rc(0) <= rc(1))) { mt = 1;}
// if ((rc(0) <= rc(2)) & (rc(1) <= rc(2))) { mt = 2;}
// if (mt == 0) {if (rc(2) < rc(1)) { mu = 1; p = 2;} else { mu = 2; p = 1;}}
// if (mt == 1) {if (rc(2) < rc(0)) { mu = 0; p = 2;} else { mu = 2; p = 0;}}
// if (mt == 2) {if (rc(1) < rc(0)) { mu = 0; p = 1;} else { mu = 1; p = 0;}}
//p = which_min(rc); // is this 0 or 1 based
//mt = which_max(rc);
//Rcpp::IntegerVector q(3), q2(2);
//q(0) = 0; q(1) = 1; q(2) = 2;
//q2(0) = mt; q2(1) = p;
//Rcpp::IntegerVector mTEMP = setdiff(q, q2);
//mu = mTEMP(0);
// if (verbose) {
//     std::cout << "haplotype designations are" << std::endl;
//     std::cout << "rc = " << rc << std::endl;            
//     std::cout << "mt = " << mt << std::endl;
//     std::cout << "mu = " << mu << std::endl;            
//     std::cout << "p = " << p << std::endl;
// }
//
// end of (old) code
//
// pre-normalize alpha - note - breaks normal alpha
// for(iGrid = 0; iGrid < nGrids; iGrid++) {
//     g_temp = 1 / c1(iGrid);
//     alphaHat_t1.col(iGrid) *= g_temp;
//     g_temp = 1 / c2(iGrid);
//     alphaHat_t2.col(iGrid) *= g_temp;
//     g_temp = 1 / c3(iGrid);
//     alphaHat_t3.col(iGrid) *= g_temp;
// }
// //
// uggggh - do manual re-naming now
// also uggggggggggggh - dividing by previous / updating
// is this efficient - does this prevent re-allocation
// could test eventually
// if (mt == 0) {gammaMT_t_local = alphaHat_t1 % betaHat_t1;}
// if (mt == 1) {gammaMT_t_local = alphaHat_t2 % betaHat_t2;}
// if (mt == 2) {gammaMT_t_local = alphaHat_t3 % betaHat_t3;}
// //
// if (mu == 0) {gammaMU_t_local = alphaHat_t1 % betaHat_t1;}
// if (mu == 1) {gammaMU_t_local = alphaHat_t2 % betaHat_t2;}
// if (mu == 2) {gammaMU_t_local = alphaHat_t3 % betaHat_t3;}
// //
// if (p == 0) {gammaP_t_local = alphaHat_t1 % betaHat_t1;}
// if (p == 1) {gammaP_t_local = alphaHat_t2 % betaHat_t2;}
// if (p == 2) {gammaP_t_local = alphaHat_t3 % betaHat_t3;}


// note - i is 1 based, which one of the weightings were performing, like if we're averaging 1-4, i would start at 1
//
//' @export
// [[Rcpp::export]]
void rcpp_fly_weighter(
    int i,
    const arma::mat & gCurrent,
    arma::mat & gAverage,
    double ll_current,
    Rcpp::NumericVector& log_mult,
    Rcpp::NumericVector& ll_rescaled,
    const arma::mat& gCurrent2,
    arma::mat& gAverage2,
    const arma::mat& gCurrent3,
    arma::mat& gAverage3,
    double& relative_difference,
    double& log_rescale,
    const int log_mult_max = 40,
    const bool equal_weighting = false
) {
    // initial values of relative_difference and log_rescale are irrelevant
    // they are written over anyway
    double relative_log_difference, c;
    int j;
    if (equal_weighting) {
        ll_current = 1;
    }
    if (i == 1) {
        log_mult(0) = ll_current;
        gAverage = gCurrent;
        gAverage2 = gCurrent2;
        gAverage3 = gCurrent3;
    } else {
        relative_log_difference = ll_current - log_mult(0);
        if (relative_log_difference > log_mult_max) {
            //
            log_rescale = log_mult(0) - ll_current;
            for(j = 0; j < (i - 1); j++) {
                ll_rescaled(j) += log_rescale;
            }
            //
            c = std::exp(log_rescale);
            gAverage *= c;
            gAverage2 *= c;
            gAverage3 *= c;            
            // this is the new value
            log_mult(0) = ll_current;
        }
        relative_difference = std::exp(ll_current - log_mult(0));
        gAverage += relative_difference * gCurrent;
        gAverage2 += relative_difference * gCurrent2;
        gAverage3 += relative_difference * gCurrent3;        
    }
    ll_rescaled(i - 1) = ll_current - log_mult(0);
    return;
}


void rcpp_single_fly_weighter(
    int i,
    arma::mat& gCurrent,
    arma::mat& gAverage,
    double& relative_difference,
    double& log_rescale
) {
    if (i == 1) {
        gAverage = gCurrent;
    } else {
        if (0 < log_rescale) {
            double c = std::exp(log_rescale);
            gAverage *= c;
        }
        gAverage += relative_difference * gCurrent;
    }
    return;
}





void set_seed(double seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(seed)));
}
//https://stackoverflow.com/questions/43221681/changing-rs-seed-from-rcpp-to-guarantee-reproducibility




Rcpp::NumericMatrix unpack_gammas(
    int s,
    std::string& prev_section,
    std::string& next_section,    
    const int suppressOutput,
    double& prev,
    const bool verbose,
    const int nReads,
    Rcpp::IntegerVector& H,
    const int nGrids,
    const arma::cube& eHapsCurrent_tc,
    const Rcpp::IntegerVector& grid,
    int snp_start_1_based,
    int snp_end_1_based,
    const int run_fb_grid_offset,
    arma::mat& gammaMT_t,
    arma::mat& gammaMU_t,
    arma::mat& gammaP_t,
    arma::mat& gammaMT_t_local,
    arma::mat& gammaMU_t_local,
    arma::mat& gammaP_t_local,
    arma::mat& genProbsM_t,
    arma::mat& genProbsF_t,
    arma::mat& genProbsM_t_local,
    arma::mat& genProbsF_t_local,
    arma::mat& hapProbs_t,
    arma::mat& hapProbs_t_local,
    arma::mat& alphaHat_t1,
    arma::mat& betaHat_t1,
    arma::rowvec& c1,
    arma::mat& eMatGrid_t1,
    arma::mat& alphaHat_t2,
    arma::mat& betaHat_t2,
    arma::rowvec& c2,
    arma::mat& eMatGrid_t2,
    arma::mat& alphaHat_t3,
    arma::mat& betaHat_t3,
    arma::rowvec& c3,
    arma::mat& eMatGrid_t3,
    const double ff,
    const bool run_fb_subset,
    const int i_snp_block_for_alpha_beta,    
    Rcpp::NumericMatrix previous_hap_label_prob_matrix,
    int i_gibbs_samplings,
    int i_result_it,
    const bool generate_fb_snp_offsets,
    const bool haploid_gibbs_equal_weighting,
    const bool return_gamma,
    const bool return_hapProbs,
    const Rcpp::NumericVector prior_probs,
    Rcpp::NumericVector& hg_log_mult,
    Rcpp::NumericVector& hg_ll_rescaled,
    const int log_mult_max = 40    
) {
    //
    next_section="unpack gammas";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    // get read counts
    //
    int iRead;    
    Rcpp::NumericVector rc(3);
    Rcpp::IntegerVector ro;
    int iGrid;
    double g_temp;
    arma::vec gamma_col;
    rc.fill(0);
    for(iRead = 0; iRead < nReads; iRead++) {
        rc(H(iRead) - 1) += 1;
    }
    Rcpp::NumericMatrix hap_label_prob_matrix;
    //
    // annoying need to keep separate version - hopefully not too slow
    //
    if (!run_fb_subset) {
        hap_label_prob_matrix = rcpp_determine_label_probabilities(rc, ff);
    } else {
        hap_label_prob_matrix = previous_hap_label_prob_matrix;
    }
    //
    next_section="add local gammas";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    // if diagonals are all close to 1, skip this step
    if (
        (hap_label_prob_matrix(0, 0) > 0.99) &&
        (hap_label_prob_matrix(1, 1) > 0.99) &&
        (hap_label_prob_matrix(2, 2) > 0.99)
    ) {
        //
        gammaMT_t_local = alphaHat_t1 % betaHat_t1;
        for(iGrid = 0; iGrid < nGrids; iGrid++) {
            g_temp = 1 / c1(iGrid);
            gammaMT_t_local.col(iGrid) *= g_temp;
        }
        gammaMU_t_local = alphaHat_t2 % betaHat_t2;
        for(iGrid = 0; iGrid < nGrids; iGrid++) {
            g_temp = 1 / c2(iGrid);
            gammaMU_t_local.col(iGrid) *= g_temp;
        }
        gammaP_t_local = alphaHat_t3 % betaHat_t3;
        for(iGrid = 0; iGrid < nGrids; iGrid++) {
            g_temp = 1 / c3(iGrid);
            gammaP_t_local.col(iGrid) *= g_temp;
        }
    } else {
        // ugh, see + vs += below...
        for(iGrid = 0; iGrid < nGrids; iGrid++) {
            g_temp = 1 / c1(iGrid);        
            gamma_col = (alphaHat_t1.col(iGrid) % betaHat_t1.col(iGrid)) * g_temp;
            gammaMT_t_local.col(iGrid) = hap_label_prob_matrix(0, 0) * gamma_col;
            gammaMU_t_local.col(iGrid) = hap_label_prob_matrix(1, 0) * gamma_col;
            gammaP_t_local.col(iGrid) = hap_label_prob_matrix(2, 0) * gamma_col;
        }
        for(iGrid = 0; iGrid < nGrids; iGrid++) {
            g_temp = 1 / c2(iGrid);        
            gamma_col = (alphaHat_t2.col(iGrid) % betaHat_t2.col(iGrid)) * g_temp;
            gammaMT_t_local.col(iGrid) += hap_label_prob_matrix(0, 1) * gamma_col;
            gammaMU_t_local.col(iGrid) += hap_label_prob_matrix(1, 1) * gamma_col;
            gammaP_t_local.col(iGrid) += hap_label_prob_matrix(2, 1) * gamma_col;
        }
        for(iGrid = 0; iGrid < nGrids; iGrid++) {
            g_temp = 1 / c3(iGrid);        
            gamma_col = (alphaHat_t3.col(iGrid) % betaHat_t3.col(iGrid)) * g_temp;
            gammaMT_t_local.col(iGrid) += hap_label_prob_matrix(0, 2) * gamma_col;
            gammaMU_t_local.col(iGrid) += hap_label_prob_matrix(1, 2) * gamma_col;
            gammaP_t_local.col(iGrid) += hap_label_prob_matrix(2, 2) * gamma_col;
        }
    }
    //
    next_section="calculate genProbs";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    // these are over-written entirely
    //genProbsM_t_local.fill(0); // ugh, this could be fixed I think
    //genProbsF_t_local.fill(0);
    //hapProbs_t_local.fill(0);
    // could do re-weighting on the build, but that might be exhausting
    rcpp_calculate_gn_genProbs_and_hapProbs(
        genProbsM_t_local,
        genProbsF_t_local,    
        hapProbs_t_local,
        s,
        eHapsCurrent_tc,
        gammaMT_t_local,
        gammaMU_t_local,
        gammaP_t_local,
        grid,
        snp_start_1_based,
        snp_end_1_based,
        run_fb_grid_offset
    );
    // duplicate sometimes, but OK, is easy
    next_section="fly weighter";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    // also when no full its
    Rcpp::NumericVector temp = calculate_likelihoods_values(c1, c2, c3, H, nGrids, prior_probs, ff);    
    double ll_current = temp(5); // "p_H_given_O_L_up_to_C"
    // if not run_fb_subset, this is the first time through
    // if it IS run_fb_subset, we already have the weights, and can use those
    double log_rescale, relative_difference;
    if (!run_fb_subset) {
        rcpp_fly_weighter(
            i_result_it + 1, // is 1-based
            genProbsM_t_local,
            genProbsM_t, // modified
            ll_current,
            hg_log_mult, // modified
            hg_ll_rescaled, // modified
            genProbsF_t_local,
            genProbsF_t, // modified
            hapProbs_t_local,
            hapProbs_t, // modified
            relative_difference,
            log_rescale,
            log_mult_max,
            haploid_gibbs_equal_weighting
        );
        // gammas (do I really need to save these? - make optional?)
        // yes, I think so - hapProbs
    } else {
        // re-scale to have 
        relative_difference = std::exp(hg_ll_rescaled(i_result_it));
        genProbsM_t += relative_difference * genProbsM_t_local;
        genProbsF_t += relative_difference * genProbsF_t_local;
        hapProbs_t += relative_difference * hapProbs_t_local;
    }
    //
    next_section="adding gamma";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    // store gamma if we need them later
    if (return_gamma) {
        rcpp_single_fly_weighter(i_result_it + 1, gammaMT_t_local, gammaMT_t, relative_difference, log_rescale);
        rcpp_single_fly_weighter(i_result_it + 1, gammaMU_t_local, gammaMU_t, relative_difference, log_rescale);
        rcpp_single_fly_weighter(i_result_it + 1, gammaP_t_local, gammaP_t, relative_difference, log_rescale);        
    }
    // convert names now
    next_section="done unpacking gamms";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    return(hap_label_prob_matrix);
}




//' @export
// [[Rcpp::export]]
Rcpp::List rcpp_forwardBackwardGibbsNIPT(
    const Rcpp::List& sampleReads,
    const arma::mat& priorCurrent_m,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& eHapsCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const double ff,
    const arma::mat& blocks_for_output,
    arma::mat& alphaHat_t1,
    arma::mat& betaHat_t1,
    arma::mat& alphaHat_t2,
    arma::mat& betaHat_t2, 
    arma::mat& alphaHat_t3, 
    arma::mat& betaHat_t3,
    arma::cube& hapSum_tc,
    Rcpp::IntegerVector& wif0,
    Rcpp::IntegerVector& L_grid,
    Rcpp::NumericVector& smooth_cm,
    const int Jmax_local = 100,
    const double maxDifferenceBetweenReads = 1000,
    const double maxEmissionMatrixDifference = 10000000000,
    const bool run_fb_subset = false,
    const int run_fb_grid_offset = 0,
    const bool return_genProbs = true,
    const bool return_hapProbs = true,
    bool return_gamma = true,
    const bool return_alpha = false,
    const bool return_p_store = false,
    const bool return_extra = false,
    const Rcpp::IntegerVector& grid = 0,
    int snp_start_1_based = -1,
    int snp_end_1_based = -1,
    const bool generate_fb_snp_offsets = false,
    const int suppressOutput = 1,
    int n_gibbs_starts = 1,    
    const int n_gibbs_sample_its = 1,
    const int n_gibbs_burn_in_its = 1,    
    const bool use_starting_read_labels = false,
    const bool verbose = false,
    const Rcpp::List& double_list_of_starting_read_labels = R_NilValue,
    Rcpp::IntegerVector seed_vector = -1,
    const Rcpp::List& prev_list_of_alphaBetaBlocks = R_NilValue,
    const int i_snp_block_for_alpha_beta = 1,
    const bool haploid_gibbs_equal_weighting = true,
    const bool gibbs_initialize_iteratively = false,
    const bool gibbs_initialize_at_first_read = true,
    const bool do_block_resampling = false,
    const int artificial_relabel = -1,
    const bool pass_in_alphaBeta = false,
    const bool update_in_place = false,
    const bool update_hapSum = false,
    const bool record_read_set = false,
    const double class_sum_cutoff = 0.06,
    const bool perform_block_gibbs = false,
    const int shuffle_bin_radius = 5000,
    const Rcpp::IntegerVector block_gibbs_iterations = Rcpp::IntegerVector::create(0),
    const bool return_gibbs_block_output = false,
    const bool rescale_eMatRead_t = true,
    const bool use_smooth_cm_in_block_gibbs = false,
    const double block_gibbs_quantile_prob = 0.9
) {
    // I think these break the gibbs-ness - disable for now!
    // rescale_eMatRead_t should be fine to reset - will be constant across reads - only the read not per-base input considered
    const bool bound_eMatGrid_t = false;
    const bool rescale_eMatGrid_t = false;
    double prev=clock();
    std::string prev_section="Null";
    std::string next_section="Initialize variables";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    //
    // initialize variables 
    //
    //
    int n_gibbs_full_its = n_gibbs_burn_in_its + n_gibbs_sample_its;    
    int nReads = sampleReads.size();
    int nGrids = transMatRate_tc_H.n_cols + 1;
    int K = eHapsCurrent_tc.n_rows;
    const int S = eHapsCurrent_tc.n_slices;
    if ((snp_start_1_based == -1) & (snp_end_1_based == -1)) {
        snp_start_1_based = 1;
        snp_end_1_based = grid.size();
    }
    int nSNPsLocal;
    if (snp_start_1_based == -1) {
        nSNPsLocal = grid.length();
    } else {
        nSNPsLocal = snp_end_1_based - snp_start_1_based + 1; // also 1-based
    }
    //
    const Rcpp::NumericVector prior_probs = NumericVector::create(0.5, (1 - ff) / 2, (ff / 2));    
    //
    Rcpp::List to_return;
    Rcpp::List list_of_alphaBetaBlocks, double_list_of_alphaBetaBlocks; // store all of them!
    Rcpp::List alphaBetaBlocks_one, alphaBetaBlocks1, alphaBetaBlocks2, alphaBetaBlocks3;
    //
    Rcpp::NumericMatrix previous_hap_label_prob_matrix, hap_label_prob_matrix;
    //
    //
    arma::mat gammaMT_t_master, gammaMU_t_master, gammaP_t_master;
    arma::mat gammaMT_t, gammaMU_t, gammaP_t;    
    if (return_gamma) {
        gammaMT_t_master = arma::zeros(K, nGrids);
        gammaMU_t_master = arma::zeros(K, nGrids);
        gammaP_t_master = arma::zeros(K, nGrids);
        gammaMT_t = arma::zeros(K, nGrids);
        gammaMU_t = arma::zeros(K, nGrids);
        gammaP_t = arma::zeros(K, nGrids);
    }
    arma::mat gammaMT_t_local = arma::zeros(K, nGrids);
    arma::mat gammaMU_t_local = arma::zeros(K, nGrids);
    arma::mat gammaP_t_local = arma::zeros(K, nGrids);
    arma::mat genProbsM_t_master = arma::zeros(3, nSNPsLocal);
    arma::mat genProbsF_t_master = arma::zeros(3, nSNPsLocal);
    arma::mat genProbsM_t = arma::zeros(3, nSNPsLocal);
    arma::mat genProbsF_t = arma::zeros(3, nSNPsLocal);
    arma::mat genProbsM_t_local = arma::zeros(3, nSNPsLocal);
    arma::mat genProbsF_t_local = arma::zeros(3, nSNPsLocal);
    arma::mat hapProbs_t_master = arma::zeros(3, nSNPsLocal);    
    arma::mat hapProbs_t = arma::zeros(3, nSNPsLocal);
    arma::mat hapProbs_t_local = arma::zeros(3, nSNPsLocal);
    //
    //
    // things for weightings
    //
    Rcpp::NumericVector hg_log_mult = NumericVector::create(-1, -1);
    Rcpp::NumericVector hg_ll_rescaled;
    if (n_gibbs_full_its > 0) {
        hg_ll_rescaled = Rcpp::NumericVector(n_gibbs_starts * n_gibbs_sample_its);
    } else {
        hg_ll_rescaled = Rcpp::NumericVector(n_gibbs_starts);
    }
    hg_ll_rescaled.fill(-1);
    //
    //
    //
    Rcpp::NumericMatrix p_store;
    Rcpp::CharacterVector p_store_cols;
    if (return_p_store) {
        p_store_cols = CharacterVector::create(
            "p_1", "p_2", "p_3",
            "chance", "h_rC", "h_rN",
            "agreePer",
            "itType",
            "p_O1_given_H1_L",
            "p_O2_given_H2_L",
            "p_O3_given_H3_L",
            "p_O_given_H_L",
            "p_H_given_L", 
            "p_H_given_O_L_up_to_C"
        );
        p_store = Rcpp::NumericMatrix(n_gibbs_starts * n_gibbs_full_its * nReads, 3 + 3 + 2 + 3 + 3);
        p_store.fill(0);
        colnames(p_store) = p_store_cols;
    }
    //
    //
    //
    //
    //
    //
    arma::mat eMatHapOri_t;
    arma::vec pRgivenH1;
    arma::vec pRgivenH2;
    arma::mat eMatRead_t = arma::ones(K, nReads);
    const bool run_pseudo_haploid = false;
    //
    //
    //
    next_section="Initialize containers";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    if (!pass_in_alphaBeta) {
        alphaHat_t1 = arma::zeros(K, nGrids);
        betaHat_t1 = arma::zeros(K, nGrids);
        alphaHat_t2 = arma::zeros(K, nGrids);
        betaHat_t2 = arma::zeros(K, nGrids);
        alphaHat_t3 = arma::zeros(K, nGrids);
        betaHat_t3 = arma::zeros(K, nGrids);
    }
    arma::rowvec c1 = arma::zeros(1, nGrids);
    arma::mat eMatGrid_t1 = arma::ones(K, nGrids);    
    arma::rowvec c2 = arma::zeros(1, nGrids);
    arma::mat eMatGrid_t2 = arma::ones(K, nGrids);    
    arma::rowvec c3 = arma::zeros(1, nGrids);
    arma::mat eMatGrid_t3 = arma::ones(K, nGrids);
    if (update_hapSum & !update_in_place) {
        hapSum_tc = arma::zeros(K, nGrids, S);
    }
    //
    //
    //
    //
    // loop comes next
    //
    //
    double gamma_temp;    
    int first_read_for_gibbs_initialization;
    Rcpp::NumericVector runif_reads;
    Rcpp::IntegerVector H, H_class;
    Rcpp::List list_of_starting_read_labels;
    Rcpp::List double_list_of_ending_read_labels, list_of_ending_read_labels;
    Rcpp::NumericMatrix rlc(7, 3);
    if (record_read_set) {
        H_class = Rcpp::IntegerVector(nReads);
        //
        rlc(0, 0) = 1; rlc(0, 1) = 0; rlc(0, 2) = 0;
        rlc(1, 0) = 0; rlc(1, 1) = 1; rlc(1, 2) = 0;
        rlc(2, 0) = 0; rlc(2, 1) = 0; rlc(2, 2) = 1;
        // longer ones
        rlc(3, 0) = prior_probs(0) / (prior_probs(0) + prior_probs(1));
        rlc(3, 1) = prior_probs(1) / (prior_probs(0) + prior_probs(1));
        rlc(3, 2) = 0;
        //
        rlc(4, 0) = prior_probs(0) / (prior_probs(0) + prior_probs(2));
        rlc(4, 1) = 0;
        rlc(4, 2) = prior_probs(2) / (prior_probs(0) + prior_probs(2));
        //
        rlc(5, 0) = 0;
        rlc(5, 1) = prior_probs(1) / (prior_probs(1) + prior_probs(2));
        rlc(5, 2) = prior_probs(2) / (prior_probs(1) + prior_probs(2));
        //
        rlc(6, 0) = prior_probs(0);
        rlc(6, 1) = prior_probs(1);
        rlc(6, 2) = prior_probs(2);
    } else {
        H_class = Rcpp::IntegerVector(1);
    }
    //
    // this outer loop is:
    // normal: gibbs starts
    // run_fb_subset: each one represents a sampled iteration
    //
    int iteration = 0;
    int i_gibbs_samplings, i_result_it, i_ever_it, n_outer, n_results, n_per_it_likelihoods;
    //
    // outer loop in normal (i.e. !run_fb_subset) is the gibbs starts
    // outer loop in run_fb_subset is over each result we care about (a sample of H)
    //
    // n_results = how many total samplings are we considering here e.g. adding to genProbsM_t
    //
    // i_outer = outer loop, useful for run_fb_subset
    // i_gibbs_sampling = which gibbs start this is
    // i_result_it = among all results (n_gibbs_starts * n_sample_its), which is this
    // i_ever_it = among all iterations (n_gibbs_starts * n_full_its), which is this
    // 
    if (!run_fb_subset) {
        // i.e. normal
        n_outer = n_gibbs_starts;
        n_results = n_gibbs_starts * n_gibbs_sample_its;
        n_per_it_likelihoods = S * n_gibbs_starts * (n_gibbs_sample_its + n_gibbs_burn_in_its);
        if (n_gibbs_sample_its == 0) {
            n_per_it_likelihoods = S * n_gibbs_starts;
        }
    } else {
        n_outer = n_gibbs_starts * n_gibbs_sample_its;
        n_results = n_outer;
        n_per_it_likelihoods = S * n_outer;
    }
    //defined original in R
    // column names are P(O^1 | \lambda) for (1, 2, 3), then P(H | \lambda), then P(O, H | \lambda) aka the likelihood (for H) ll = log10(P(H | O, \lambda))
    int i_per_it_likelihoods = 0; // easier to just keep counter!
    Rcpp::NumericMatrix per_it_likelihoods = Rcpp::NumericMatrix(n_per_it_likelihoods, 12);
    colnames(per_it_likelihoods) = CharacterVector::create("s", "i_samp", "i_it", "i_result_it", "p_O1_given_H1_L", "p_O2_given_H2_L", "p_O3_given_H3_L", "p_O_given_H_L", "p_H_given_L", "p_H_given_O_L_up_to_C", "p_set_H_given_L", "relabel");
    //
    //
    //
    for(int s = 0; s < S; s++) {
        //
        list_of_ending_read_labels = R_NilValue; // null-ify
        if (s > 0) {
            eMatRead_t.fill(1);
            genProbsM_t.fill(0);
            genProbsF_t.fill(0);
            hapProbs_t.fill(0);            
            if (return_gamma) {
                gammaMT_t.fill(0);
                gammaMU_t.fill(0);
                gammaP_t.fill(0);
            }
            list_of_alphaBetaBlocks = R_NilValue;
        }
        if (run_fb_subset) {
            // so previously worked on double list
            // now gets back to the single (previous) list
            list_of_alphaBetaBlocks = as<Rcpp::List>(prev_list_of_alphaBetaBlocks[s]);
        }
        //
        //
        //
        next_section="Initialize eMatRead_t";
        prev=print_times(prev, suppressOutput, prev_section, next_section);
        prev_section=next_section;
        rcpp_make_eMatRead_t(eMatRead_t, sampleReads, eHapsCurrent_tc, s, maxDifferenceBetweenReads, Jmax_local, eMatHapOri_t, pRgivenH1, pRgivenH2, prev, suppressOutput, prev_section, next_section, run_pseudo_haploid, rescale_eMatRead_t);
        
        //
        //
        //
        for(int i_outer = 0; i_outer < n_outer; i_outer++) {
            //
            //
            next_section="Start outer loop";
            prev=print_times(prev, suppressOutput, prev_section, next_section);
            prev_section=next_section;
            //
            //
            if (!run_fb_subset) {
                i_gibbs_samplings = i_outer;
            } else {
                i_gibbs_samplings = floor(i_outer / n_gibbs_sample_its); // 0-based
            }
            //
            if (seed_vector(0) > 0) {
                set_seed(seed_vector(i_gibbs_samplings + S * s));
            }
            //
            runif_reads = Rcpp::runif(nReads * n_gibbs_full_its);
            if (!gibbs_initialize_at_first_read & !run_fb_subset) {
                if (nReads > 0) {
                    first_read_for_gibbs_initialization = Rcpp::sample(nReads, 1)(0) - 1;
                } else {
                    first_read_for_gibbs_initialization = 0; // or will this just break anyway
                }
            } else {
                first_read_for_gibbs_initialization = 0; // note - this is 0-based in cpp
            }
            //
            if (use_starting_read_labels) {
                // holy fuck Rcpp would copy this by reference
                // despite the const above. argh
                list_of_starting_read_labels = as<Rcpp::List>(double_list_of_starting_read_labels(s));
                H = Rcpp::clone(as<Rcpp::IntegerVector>(list_of_starting_read_labels(i_outer)));
            } else {
                H = random_gibbs_nipt_read_labels(nReads, ff);
            }
            //
            if (run_fb_subset) {
                hg_ll_rescaled = as<Rcpp::NumericVector>(list_of_alphaBetaBlocks["hg_ll_rescaled"]);
                alphaBetaBlocks_one = as<Rcpp::List>(list_of_alphaBetaBlocks[i_outer]);
                previous_hap_label_prob_matrix = as<Rcpp::NumericMatrix>(alphaBetaBlocks_one["hap_label_prob_matrix"]);
            }
            //
            //
            //
            next_section="Initialize";
            prev=print_times(prev, suppressOutput, prev_section, next_section);
            prev_section=next_section;
            //
            //
            rcpp_gibbs_nipt_initialize(
                s, prev_section, next_section, suppressOutput, prev,
                run_fb_subset, alphaBetaBlocks_one, i_snp_block_for_alpha_beta,
                eMatRead_t, sampleReads, H, run_fb_grid_offset, bound_eMatGrid_t, rescale_eMatGrid_t,
                alphaHat_t1, betaHat_t1, c1, eMatGrid_t1,
                alphaHat_t2, betaHat_t2, c2, eMatGrid_t2,
                alphaHat_t3, betaHat_t3, c3, eMatGrid_t3,
                transMatRate_tc_H, gibbs_initialize_iteratively,
                priorCurrent_m, alphaMatCurrent_tc, maxEmissionMatrixDifference
            );
            //
            //
            // if running fb subset, do this here. also do if no its
            if ((run_fb_subset) | (n_gibbs_full_its == 0)) {
                // if no iterations, still need to update gamma
                i_result_it = i_outer;
                i_ever_it = i_outer;
                hap_label_prob_matrix = unpack_gammas(
                    s, prev_section, next_section, suppressOutput, prev,
                    verbose, nReads, H, nGrids,
                    eHapsCurrent_tc, grid, snp_start_1_based, snp_end_1_based, run_fb_grid_offset,
                    gammaMT_t, gammaMU_t, gammaP_t,
                    gammaMT_t_local, gammaMU_t_local, gammaP_t_local,
                    genProbsM_t, genProbsF_t, genProbsM_t_local, genProbsF_t_local,
                    hapProbs_t, hapProbs_t_local,                
                    alphaHat_t1, betaHat_t1, c1, eMatGrid_t1,
                    alphaHat_t2, betaHat_t2, c2, eMatGrid_t2,
                    alphaHat_t3, betaHat_t3, c3, eMatGrid_t3,
                    ff, run_fb_subset, i_snp_block_for_alpha_beta,
                    previous_hap_label_prob_matrix, i_gibbs_samplings, i_result_it,
                    generate_fb_snp_offsets, haploid_gibbs_equal_weighting,
                    return_gamma, return_hapProbs,
                    prior_probs, hg_log_mult, hg_ll_rescaled
                );
                //
                if (!run_fb_subset) {
                    // want this to go into i_result_it which is i_outer
                    add_to_per_it_likelihoods(s, per_it_likelihoods, i_gibbs_samplings, iteration, i_result_it, n_gibbs_full_its, H, c1, c2, c3, ff, prior_probs, false, i_per_it_likelihoods);
                }
            } else {
                for(iteration = 0; iteration < n_gibbs_full_its; iteration++) {
                    if ((iteration + 1) > n_gibbs_burn_in_its) {
                        i_result_it = n_gibbs_sample_its * (i_gibbs_samplings) + iteration - n_gibbs_burn_in_its;
                    } else {
                        i_result_it = -2;
                    }
                    i_ever_it = n_gibbs_full_its * (i_gibbs_samplings) + iteration;
                    //
                    rcpp_gibbs_nipt_iterate(
                        s, prev_section, next_section, suppressOutput, prev,
                        iteration, sampleReads, priorCurrent_m, transMatRate_tc_H, alphaMatCurrent_tc, n_gibbs_full_its,
                        nGrids, K, H, eMatRead_t, runif_reads,
                        alphaHat_t1, betaHat_t1, c1, eMatGrid_t1,
                        alphaHat_t2, betaHat_t2, c2, eMatGrid_t2,
                        alphaHat_t3, betaHat_t3, c3, eMatGrid_t3,
                        p_store, i_per_it_likelihoods, verbose, return_p_store, run_fb_subset,
                        i_gibbs_samplings, n_gibbs_starts, i_result_it, ff,
                        record_read_set, rlc, H_class, class_sum_cutoff,
                        run_fb_grid_offset, prior_probs, per_it_likelihoods,
                        gibbs_initialize_iteratively, first_read_for_gibbs_initialization,
                        do_block_resampling, artificial_relabel
                    );
                    // check
                    bool to_block_gibbs = false;
                    if (perform_block_gibbs) {
                        for(int i = 0; i < block_gibbs_iterations.length(); i++) {
                            if (iteration == block_gibbs_iterations(i)) {
                                to_block_gibbs = true;
                            }
                        }
                    }
                    // on burn in its, try every 10th iteration, within the first 100
                    if (perform_block_gibbs & to_block_gibbs) {
                        // define blocks
                        Rcpp::NumericVector temp; //
                        next_section="Block gibbs - define sites";
                        prev=print_times(prev, suppressOutput, prev_section, next_section);
                        prev_section=next_section;
                        Rcpp::List out = Rcpp_define_blocked_snps_using_gamma_on_the_fly(alphaHat_t1, alphaHat_t2, alphaHat_t3, betaHat_t1, betaHat_t2, betaHat_t3, c1, c2, c3, eMatGrid_t1, eMatGrid_t2, eMatGrid_t3, smooth_cm, transMatRate_tc_H, shuffle_bin_radius, L_grid, grid, s, block_gibbs_quantile_prob, verbose, use_smooth_cm_in_block_gibbs); // fix me verbose
                        Rcpp::IntegerVector blocked_snps = as<Rcpp::IntegerVector>(out["blocked_snps"]);
                        int n_blocks = blocked_snps(nSNPsLocal - 1) + 1;
                        //
                        Rcpp::NumericMatrix runif_proposed(6, nReads);
                        for(int j = 0; j < 6; j++) {
                            runif_proposed.row(j) = Rcpp::runif(nReads);
                        }
                        Rcpp::NumericVector runif_block = Rcpp::runif(nReads);
                        Rcpp::NumericVector runif_total = Rcpp::runif(nReads);
                        //
                        next_section="Block gibbs - sample";
                        prev=print_times(prev, suppressOutput, prev_section, next_section);
                        prev_section=next_section;
                        Rcpp::List out2 = Rcpp_block_gibbs_resampler(alphaHat_t1, alphaHat_t2, alphaHat_t3, betaHat_t1, betaHat_t2, betaHat_t3, c1,c2,c3, eMatGrid_t1, eMatGrid_t2, eMatGrid_t3, H, H_class, eMatRead_t, blocked_snps, runif_block, runif_total, runif_proposed, grid, wif0, ff, s, alphaMatCurrent_tc, priorCurrent_m, transMatRate_tc_H, maxDifferenceBetweenReads, Jmax_local);
                        // , false, R_NilValue, true
                        //std::cout << "remove last three values above me" << std::endl;
                        //
                        if (return_gibbs_block_output) {
                            to_return.push_back(out2, "gibbs_block_output");
                            to_return.push_back(out, "block_stuff");
                        }
                    }
                    //
                    // save things when either: sampling iteration (yay!) or run_fb_subset (always!)
                    // 
                    if ((iteration + 1) > n_gibbs_burn_in_its) {
                        hap_label_prob_matrix = unpack_gammas(
                            s, prev_section, next_section, suppressOutput, prev,
                            verbose, nReads, H, nGrids,
                            eHapsCurrent_tc, grid, snp_start_1_based, snp_end_1_based, run_fb_grid_offset,
                            gammaMT_t, gammaMU_t, gammaP_t,
                            gammaMT_t_local, gammaMU_t_local, gammaP_t_local,
                            genProbsM_t, genProbsF_t, genProbsM_t_local, genProbsF_t_local,
                            hapProbs_t, hapProbs_t_local,
                            alphaHat_t1, betaHat_t1, c1, eMatGrid_t1,
                            alphaHat_t2, betaHat_t2, c2, eMatGrid_t2,
                            alphaHat_t3, betaHat_t3, c3, eMatGrid_t3,
                            ff, run_fb_subset, i_snp_block_for_alpha_beta,
                            previous_hap_label_prob_matrix, i_gibbs_samplings, i_result_it,
                            generate_fb_snp_offsets, haploid_gibbs_equal_weighting,
                            return_gamma, return_hapProbs,
                            prior_probs, hg_log_mult, hg_ll_rescaled
                        );
                        list_of_ending_read_labels.push_back(Rcpp::clone(H), "H") ;
                        //
                        if (generate_fb_snp_offsets) {
                            alphaBetaBlocks_one = R_NilValue; // null-ify
                            alphaBetaBlocks1 = rcpp_make_fb_snp_offsets(alphaHat_t1, betaHat_t1, blocks_for_output);
                            alphaBetaBlocks2 = rcpp_make_fb_snp_offsets(alphaHat_t2, betaHat_t2, blocks_for_output);
                            alphaBetaBlocks3 = rcpp_make_fb_snp_offsets(alphaHat_t3, betaHat_t3, blocks_for_output);
                            alphaBetaBlocks_one.push_back(alphaBetaBlocks1, "alphaBetaBlocks1");
                            alphaBetaBlocks_one.push_back(alphaBetaBlocks2, "alphaBetaBlocks2");
                            alphaBetaBlocks_one.push_back(alphaBetaBlocks3, "alphaBetaBlocks3");
                            alphaBetaBlocks_one.push_back(hap_label_prob_matrix, "hap_label_prob_matrix");
                            list_of_alphaBetaBlocks.push_back(alphaBetaBlocks_one);
                        }
                    }
                }
            }
            //
            next_section="Done all iterations";
            prev=print_times(prev, suppressOutput, prev_section, next_section);
            prev_section=next_section;
            //
            //
            //
        }
        next_section="Done all samplings";
        prev=print_times(prev, suppressOutput, prev_section, next_section);
        prev_section=next_section;
        //
        // do end of S bits here
        //
        //
        next_section="Add gammas and genProbs to out";
        prev=print_times(prev, suppressOutput, prev_section, next_section);
        prev_section=next_section;
        //
        // undo weighting here. first, determine weighting, through gamma_temp, then apply
        //
        if (!haploid_gibbs_equal_weighting) {
            gamma_temp = 0;
            for(int i =0 ; i < n_gibbs_starts * n_gibbs_sample_its; i++) {
                gamma_temp += std::exp(hg_ll_rescaled(i));
            }
            gamma_temp = 1 / gamma_temp;
        } else {
            // if sampling and burn in are both 0, so a pass through
            if (n_results > 0) {
                gamma_temp = 1 / double(n_results);
            } else {
                // i.e. no burn in or sample - just pass through
                if (n_gibbs_starts > 1) {
                    gamma_temp = 1 / double(n_gibbs_starts);
                } else {
                    gamma_temp = 1;
                }
            }
        }
        if (return_gamma) {
            //
            if (gamma_temp != 1) {
                gammaMT_t *= gamma_temp;
                gammaMU_t *= gamma_temp;
                gammaP_t *= gamma_temp;
            }
            gammaMT_t_master += gammaMT_t;
            gammaMU_t_master += gammaMU_t;
            gammaP_t_master += gammaP_t;            
        }
        if (return_genProbs) {
            if (gamma_temp != 1) {
                genProbsM_t *= gamma_temp;
                genProbsF_t *= gamma_temp;
            }
            genProbsM_t_master += genProbsM_t;
            genProbsF_t_master += genProbsF_t;
        }
        if (return_hapProbs) {
            if (gamma_temp != 1) {
                hapProbs_t *= gamma_temp;
            }
            hapProbs_t_master += hapProbs_t;
        }
        if (!((run_fb_subset) | (n_gibbs_full_its == 0))) {
            double_list_of_ending_read_labels.push_back(list_of_ending_read_labels, "list_of_ending_read_labels");
        }
        if (generate_fb_snp_offsets) {
            list_of_alphaBetaBlocks.push_back(hg_ll_rescaled, "hg_ll_rescaled");
            double_list_of_alphaBetaBlocks.push_back(list_of_alphaBetaBlocks);
        }
        if (update_hapSum & !run_fb_subset) {
            hapSum_tc.slice(s) += 0.5 * (gammaMT_t + gammaMU_t); // close enough
        }
    }
    //
    //
    next_section="Done all S";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    //
    //
    gamma_temp = 1 / double(S);
    if (return_gamma) {
        if (1 < S) {
            gammaMT_t_master *= gamma_temp;
            gammaMU_t_master *= gamma_temp;
            gammaP_t_master *= gamma_temp;
        }
        to_return.push_back(gammaMT_t_master, "gammaMT_t");
        to_return.push_back(gammaMU_t_master, "gammaMU_t");
        to_return.push_back(gammaP_t_master, "gammaP_t");            
    }
    if (return_genProbs) {
        if (1 < S) {
            genProbsM_t_master *= gamma_temp;
            genProbsF_t_master *= gamma_temp;
        }
        to_return.push_back(genProbsM_t_master, "genProbsM_t");
        to_return.push_back(genProbsF_t_master, "genProbsF_t");
    }
    if (return_hapProbs) {
        if (1 < S) {
            hapProbs_t_master *= gamma_temp;
        }
        to_return.push_back(hapProbs_t_master, "hapProbs_t");
    }
    //
    // only bother with these on final if there are repeats
    //
    next_section="Finish up, p1, p2, H";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    // now, at the end, push back the normalized version
    //
    to_return.push_back(H, "H");
    to_return.push_back(double_list_of_ending_read_labels, "double_list_of_ending_read_labels");  // now they've been reset
    if (return_alpha) {
        to_return.push_back(alphaHat_t1, "alphaHat_t1");
        to_return.push_back(alphaHat_t2, "alphaHat_t2");
        to_return.push_back(alphaHat_t3, "alphaHat_t3");
        to_return.push_back(betaHat_t1, "betaHat_t1");
        to_return.push_back(betaHat_t2, "betaHat_t2");
        to_return.push_back(betaHat_t3, "betaHat_t3");
        to_return.push_back(c1, "c1");
        to_return.push_back(c2, "c2");
        to_return.push_back(c3, "c3");
        to_return.push_back(eMatGrid_t1, "eMatGrid_t1");
        to_return.push_back(eMatGrid_t2, "eMatGrid_t2");
        to_return.push_back(eMatGrid_t3, "eMatGrid_t3");
    }
    to_return.push_back(per_it_likelihoods, "per_it_likelihoods");        
    //
    if (return_p_store) {
        to_return.push_back(p_store, "p_store");
    }
    if (record_read_set) {
        to_return.push_back(H_class, "H_class");
    }
    if (return_extra) {
        to_return.push_back(eMatRead_t, "eMatRead_t");
        to_return.push_back(eMatGrid_t1, "eMatGrid_t1");
        to_return.push_back(eMatGrid_t2, "eMatGrid_t2");
        to_return.push_back(eMatGrid_t3, "eMatGrid_t3");
    }
    if (generate_fb_snp_offsets) {
        next_section="offsets";
        prev=print_times(prev, suppressOutput, prev_section, next_section);
        prev_section=next_section;
        //
        // re-cast with original (in R) name
        to_return.push_back(double_list_of_alphaBetaBlocks, "list_of_alphaBetaBlocks"); // different format!
    }
    if (update_hapSum & !update_in_place) {
        to_return.push_back(hapSum_tc, "hapSum_tc");
    }
    //
    next_section="Done";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    return(to_return);
}