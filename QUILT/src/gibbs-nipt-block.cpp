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


void Rcpp_run_backward_haploid(
    arma::mat& betaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const int s
);

void Rcpp_run_backward_haploid_QUILT_faster(
    arma::mat& betaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& transMatRate_tc_H,
    const Rcpp::LogicalVector& grid_has_read,
    const int s
);

Rcpp::NumericVector rcpp_make_smoothed_rate(
    const Rcpp::NumericVector & sigma_rate,
    const Rcpp::IntegerVector & L_grid,
    const int shuffle_bin_radius,
    const bool verbose = false
);

int rcpp_determine_where_to_stop(
    const Rcpp::NumericVector& smoothed_rate,
    const Rcpp::LogicalVector& available,
    int& snp_best, // 0-based here
    double& thresh,
    int& nGrids,
    bool is_left
);



void rcpp_apply_vec_relabel(
    arma::rowvec& m1,
    arma::rowvec& m2,
    arma::rowvec& m3,
    const int relabel
);

void rcpp_apply_mat_relabel(
    arma::mat& m1,
    arma::mat& m2,
    arma::mat& m3,
    const int relabel
);


//' @export
// [[Rcpp::export]]
double rcpp_simple_quantile(arma::vec x, double q) {
    int v = int(x.size() * q);
    arma::uvec a = arma::sort_index(x, "ascend");
    return x(a(v));
}


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
);


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


//' @export
// [[Rcpp::export]]
double rcpp_get_log_p_H_class(
    Rcpp::IntegerVector& H_class,
    double ff
) {
    double out = 0;
    Rcpp::NumericVector vals(8);
    vals(0) = 0;
    vals(1) = log(0.5);
    vals(2) = log(0.5 - ff * 0.5);
    vals(3) = log(ff * 0.5);
    vals(4) = log(1.0 - ff * 0.5);
    vals(5) = log(0.5 + ff * 0.5);
    vals(6) = log(0.5);
    vals(7) = 0;
    for(int i = 0; i < H_class.length(); i++) {
        out += vals[H_class[i]];
    }
    return(out);
}


//' @export
// [[Rcpp::export]]
double rcpp_get_log_p_H_class2(
    int n1,
    int n2,
    int n3,
    int n4,
    int n5,
    int n6,
    double ff
) {
    // so if ff == 0, or ff == 1, some options should not be possible
    // so we should penalize them heavily!
    double out = 0;
    if (ff == 0) {
        // problem is that ff / 2 is now 0 so log of 0
        out = 0 +
            n1 * log(0.5) +
            n2 * log(0.5 - ff * 0.5) +
            n3 * log(0.001) +
            n4 * log(1 - ff * 0.5) + 
            n5 * log(1 * 0.5 + ff * 0.5) + 
            n6 * log(1 * 0.5);
    } else if (ff == 1) {
        out = 0 +
            n1 * log(0.5) +
            n2 * log(0.001) +
            n3 * log(ff * 0.5) +
            n4 * log(1 - ff * 0.5) + 
            n5 * log(1 * 0.5 + ff * 0.5) + 
            n6 * log(1 * 0.5);
    } else {
        out = 0 +
            n1 * log(0.5) +
            n2 * log(0.5 - ff * 0.5) +
            n3 * log(ff * 0.5) +
            n4 * log(1 - ff * 0.5) + 
            n5 * log(1 * 0.5 + ff * 0.5) + 
            n6 * log(1 * 0.5);
    }
    return(out);
}


//' @export
// [[Rcpp::export]]
void rcpp_sample_H_using_H_class(
    Rcpp::IntegerVector & H_class,
    Rcpp::IntegerVector & H,
    double ff
) {
    Rcpp::IntegerVector one_through_three = Rcpp::IntegerVector::create(1,2,3);
    //
    Rcpp::NumericVector probs07 = Rcpp::NumericVector::create(0.5, 0.5 - ff * 0.5, ff * 0.5);
    Rcpp::NumericVector probs4 = Rcpp::NumericVector::create(0.5, 0.5 - 0.5 * ff, 0);
    Rcpp::NumericVector probs5 = Rcpp::NumericVector::create(0.5, 0, 0.5 * ff);
    Rcpp::NumericVector probs6 = Rcpp::NumericVector::create(0, 0.5 - ff * 0.5, ff * 0.5);
    //
    Rcpp::IntegerVector i_temp(1);
    for(int iRead = 0; iRead < H.length(); iRead++) {
        int hc = H_class[iRead];
        if (hc == 0 | hc == 7) {
            i_temp = Rcpp::sample(one_through_three, 1, false, probs07);
        } else if (hc == 1) {
            i_temp(0) = 1;
        } else if (hc == 2) {
            i_temp(0) = 2;
        } else if (hc == 3) {
            i_temp(0) = 3;
        } else if (hc == 4) {
            i_temp = Rcpp::sample(one_through_three, 1, false, probs4);
        } else if (hc == 5) {
            i_temp = Rcpp::sample(one_through_three, 1, false, probs5);            
        } else if (hc == 6) {
            i_temp = Rcpp::sample(one_through_three, 1, false, probs6);                        
        }
        H[iRead] = i_temp(0);
    }
    return;
}


//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_calculate_block_read_label_probabilities_using_H_class(
    int read_start_0_based,
    int read_end_0_based,
    Rcpp::IntegerVector & H_class,
    Rcpp::NumericVector log_prior_probs,
    const arma::imat & rr,
    double ff
) {
    Rcpp::NumericVector choice_log_probs_H = Rcpp::NumericVector(6);
    Rcpp::IntegerVector ns = Rcpp::IntegerVector(8);
    for(int iRead = read_start_0_based; iRead <= read_end_0_based; iRead++) {
        ns[H_class[iRead]]++;
    }
    for(int ir = 0; ir < 6; ir++) {
        //std::cout << "ir = " << ir << std::endl;
        //std::cout << "n1 = " << ns[rr(ir, 0)] << std::endl;
        //std::cout << "n4 = " << ns[7 - rr(ir, 2)] << std::endl;
        choice_log_probs_H[ir] = rcpp_get_log_p_H_class2(
            ns[rr(ir, 0)],
            ns[rr(ir, 1)],
            ns[rr(ir, 2)],
            ns[7 - rr(ir, 2)],
            ns[7 - rr(ir, 1)],
            ns[7 - rr(ir, 0)],
            ff
        );
    }
    return(choice_log_probs_H);
}


//' @export
// [[Rcpp::export]]
void rcpp_test_one_based_swap(
    const arma::imat & rx,
    int ir_chosen,
    Rcpp::IntegerVector & H,
    Rcpp::IntegerVector & H_class
) {
    // recall rx is 1-based
    Rcpp::IntegerVector zero_based_swap(8);
    zero_based_swap(0) = 0;
    // ir_chosen is 1-based here i.e. values of 1 through 6
    zero_based_swap(1) = rx(ir_chosen - 1, 0);
    zero_based_swap(2) = rx(ir_chosen - 1, 1);
    zero_based_swap(3) = rx(ir_chosen - 1, 2);
    zero_based_swap(4) = 7 - rx(ir_chosen - 1, 2);
    zero_based_swap(5) = 7 - rx(ir_chosen - 1, 1);
    zero_based_swap(6) = 7 - rx(ir_chosen - 1, 0);
    zero_based_swap(7) = 7;
    //
    for(int iRead = 0; iRead < H.length(); iRead++) {
        H_class[iRead] = zero_based_swap[H_class[iRead]];
        H[iRead] = zero_based_swap[H[iRead]];
    }
    return;
}

//' @export
// [[Rcpp::export]]
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
    const bool use_smooth_cm_in_block_gibbs = false,
    double ff = 0
) {
    //
    //
    //
    const int nSNPs = grid.length();
    int nGrids = transMatRate_tc_H.n_cols + 1;
    Rcpp::List to_return;
    //
    //
    if (verbose) {
        std::cout << "calculate rates" << std::endl;
    }
    //
    Rcpp::NumericVector rate2(nGrids - 1);
    double d;
    int iGrid;
    // probably best to run these as three - keep more close to RAM at same time
    for(iGrid = 0; iGrid < (nGrids - 2); iGrid++) {
        d = transMatRate_tc_H(0, iGrid, s);
        rate2(iGrid) += 1 - d * sum(alphaHat_t1.col(iGrid) % betaHat_t1.col(iGrid + 1) % eMatGrid_t1.col(iGrid + 1));
    }
    for(iGrid = 0; iGrid < (nGrids - 2); iGrid++) {
        d = transMatRate_tc_H(0, iGrid, s);
        rate2(iGrid) += 1 - d * sum(alphaHat_t2.col(iGrid) % betaHat_t2.col(iGrid + 1) % eMatGrid_t2.col(iGrid + 1));
    }
    if (ff > 0) {
        for(iGrid = 0; iGrid < (nGrids - 2); iGrid++) {
            d = transMatRate_tc_H(0, iGrid, s);
            rate2(iGrid) += 1 - d * sum(alphaHat_t3.col(iGrid) % betaHat_t3.col(iGrid + 1) % eMatGrid_t3.col(iGrid + 1));
        }
    }
    //
    if (verbose) {
        std::cout << "calculate smoothed rate" << std::endl;
    }
    //
    Rcpp::NumericVector smoothed_rate = rcpp_make_smoothed_rate(rate2, L_grid, shuffle_bin_radius, verbose);
    //
    //
    //
    if (use_smooth_cm_in_block_gibbs) {
        if (verbose) {
            std::cout << "multiply by recombination rate" << std::endl;
        }
        for(iGrid = 0; iGrid < (nGrids - 2); iGrid++) {
            rate2(iGrid) *= smooth_cm(iGrid);
        }
    }
    //
    //
    if (verbose) {
        std::cout << "make simple quintile" << std::endl;
    }
    double break_thresh = 1;

    
    d = rcpp_simple_quantile(smoothed_rate, block_gibbs_quantile_prob);
    if (d < break_thresh) {
        break_thresh = d;
    }
    //
    if (verbose) {
        std::cout << "make available" << std::endl;
    }
    //
    //int n = smoothed_rate.size(); // this has size nGrids - 1 by definition
    Rcpp::LogicalVector available(nGrids - 1); // this is of 
    available.fill(false);
    for(int i = 0; i < (nGrids - 1); i++) {
        if (smoothed_rate(i) == NA_REAL) {
            available(i) = false;
        } else if (smoothed_rate(i) < 0.01) {
            available(i) = false;
        }
        if (break_thresh < smoothed_rate(i)) {
            available(i) = true;
        }
    }
    //
    // cheap out
    //
    Rcpp::IntegerVector blocked_snps(nSNPs);
    if (sum(available) == 0) {
        to_return.push_back(blocked_snps, "blocked_snps");
        return(to_return);
    }
    //
    if (verbose) {
        std::cout << "do best sorting" << std::endl;
    }
    // I hope
    arma::uvec best2 = arma::sort_index( as<arma::vec>(smoothed_rate) , "descend");
    double nAvailable = sum(available);
    Rcpp::IntegerVector best(nAvailable);
    for(int i = 0; i < nAvailable; i++) {
        best(i) = int(best2(i));
    }
    //
    Rcpp::IntegerVector to_keep(0);
    //
    if (verbose) {
        std::cout << "start peak finding" << std::endl;
    }
    int a, b, j, snp_left, snp_right;
    for(int iBest = 0; iBest < nAvailable; iBest++) {
        if (verbose) {
            std::cout << "iBest = " << iBest << std::endl;
        }
        if (available(best(iBest))) {
            if (verbose) {
                std::cout << "insides" << std::endl;
            }
            int snp_best = best(iBest); // 0-based
            a = std::max(snp_best - 1, 0);
            b = std::min(snp_best + 1, nGrids - 1 - 1); // nGrids - 1 is size, another - 1 as below sum is inclusive
            if (verbose) {
                std::cout << "snp_best = " << snp_best << ", a = " << a << ", b = " << b << std::endl;
            }
            d = 0;
            for(j = a; j <= b; j++) {
                if (available(j)) {
                    d += 1;
                }
            }
            if (d == 3) {
                if (verbose) {
                    std::cout << "left" << std::endl;
                }
                snp_left = rcpp_determine_where_to_stop(smoothed_rate, available, snp_best, break_thresh, nGrids, true);
                if (verbose) {
                    std::cout << "right" << std::endl;
                }
                snp_right = rcpp_determine_where_to_stop(smoothed_rate, available, snp_best, break_thresh, nGrids, false);
                for(j = snp_left; j <= snp_right; j++) {
                    available(j) = false;
                }
                if (verbose) {
                    std::cout << "ends" << std::endl;
                }
                //std::cout << "snp_left = " << snp_left << ", snp_right = " << snp_right << std::endl;
            } else {
                for(j = a; j <= b; j++) {                
                    available(j) = false;
                }
            }
            to_keep.push_back(snp_best + 1); // store this with a +1, matches R code
        }
    }
    if (verbose) {
        std::cout << "finish off" << std::endl;
    }
    // add 0 and nGrids - 1 if not already in
    if (min(to_keep) != 0) {
        to_keep.push_back(0);
    }
    if (max(to_keep) != (nGrids - 1)) {
        to_keep.push_back(nGrids - 1);
    }
    to_keep.sort(); // hmm?
    int n = to_keep.size();    
    //Rcpp::IntegerVector blocks_to_consider = std::sort(to_keep);
    //
    Rcpp::IntegerVector blocked_grid(nGrids);
    //
    for(int i = 0; i < (n - 1); i++) {
        // a = blocks_to_consider(i);
        // b = blocks_to_consider(i + 1);
        a = to_keep(i);
        b = to_keep(i + 1);
        for(j = a; j <= b; j++) {
            blocked_grid(j) = i;
        }
    }
    // 
    for(int iSNP = 0; iSNP < nSNPs; iSNP++) {
        blocked_snps(iSNP) = blocked_grid(grid(iSNP));
    }
    //
    // save and quit
    //
    to_return.push_back(blocked_snps, "blocked_snps");
    to_return.push_back(break_thresh, "break_thresh");
    to_return.push_back(smoothed_rate, "smoothed_rate");    
    //
    if (verbose) {
        std::cout << "done making blocked snps" << std::endl;
    }
    //
    return(to_return);
}




Rcpp::NumericVector rcpp_calculate_block_read_label_probabilities(
    const double& read_start_0_based,
    const double& read_end_0_based,
    const Rcpp::IntegerVector& H,
    const arma::vec& log_prior_probs,
    const arma::imat& rr0
) {
  //
  Rcpp::NumericVector choice_log_probs_H(6);
  for(int iRead = read_start_0_based; iRead <= read_end_0_based; iRead++) {
    int h = H(iRead) - 1;
    for(int ir = 0; ir < 6; ir++) {
      choice_log_probs_H(ir) += log_prior_probs(rr0(ir, h));
    }
  }
  return(choice_log_probs_H);
}

Rcpp::NumericVector rcpp_calculate_block_read_label_probabilities_using_read_informativeness(
    const double& read_start_0_based,
    const double& read_end_0_based,
    const Rcpp::IntegerVector& H,
    const arma::vec& log_prior_probs,
    const arma::imat& rr0,
    const Rcpp::LogicalVector& read_is_uninformative    
) {
  //
  Rcpp::NumericVector choice_log_probs_H(6);
  for(int iRead = read_start_0_based; iRead <= read_end_0_based; iRead++) {
      if (!read_is_uninformative(iRead)) {
          int h = H(iRead) - 1;
          for(int ir = 0; ir < 6; ir++) {
              choice_log_probs_H(ir) += log_prior_probs(rr0(ir, h));
          }
      }
  }
  return(choice_log_probs_H);
}


Rcpp::NumericVector rcpp_calculate_block_read_label_probabilities_using_proposed_H(
    const double& read_start_0_based,
    const double& read_end_0_based,
    arma::imat& proposed_H,
    const arma::vec& log_prior_probs,
    const arma::imat& rr0
) {
  //
  Rcpp::NumericVector choice_log_probs_H(6);
  for(int iRead = read_start_0_based; iRead <= read_end_0_based; iRead++) {
      for(int ir = 0; ir < 6; ir++) {
          int h = rr0(ir, proposed_H(ir, iRead) - 1);
          choice_log_probs_H(ir) += log_prior_probs(h);
      }
  }
  return(choice_log_probs_H);
}



//' @export
// [[Rcpp::export]]
void Rcpp_consider_block_relabelling(
    const int iBlock,
    const Rcpp::NumericVector& runif_block,
    Rcpp::NumericVector& sum_H,
    const int s,
    const arma::imat& rr,    
    const arma::imat& rr0,
    double ff,
    const Rcpp::NumericVector& log_prior_probs,
    Rcpp::NumericVector& logC_before,
    Rcpp::NumericVector& logC_after,
    const bool verbose,
    const Rcpp::List& swap_list,
    arma::mat& eMatGridLocal,
    arma::mat& betaHatLocal,
    const int& iGrid,
    const int& grid_start_0_based,
    const int& grid_end_0_based,
    const int& read_start_0_based,
    const int& read_end_0_based,
    Rcpp::IntegerVector& wif0,
    arma::cube& log_cStore,
    const arma::cube& alphaStore,
    const Rcpp::LogicalVector& read_is_uninformative,
    const int block_approach,
    const bool& do_checks,
    Rcpp::List& all_packages,
    Rcpp::NumericMatrix& block_results,
    Rcpp::NumericVector& ever_changed,
    const arma::cube& transMatRate_tc_H,
    const arma::cube& alphaMatCurrent_tc,
    const arma::mat& priorCurrent_m,
    Rcpp::List& fpp_stuff,
    Rcpp::IntegerVector& H,
    Rcpp::IntegerVector& H_class,    
    arma::imat& proposed_H,
    int nReads,
    const arma::mat& eMatRead_t,
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
    const bool sample_is_diploid
) {
  //
  if (verbose) {
    std::cout << "in Rcpp_consider_block_relabelling" << std::endl;
  }
  int i, ir, iGrid2;
  //
  //
  //
  betaHatLocal.col(0) = betaHat_t1.col(iGrid);
  betaHatLocal.col(1) = betaHat_t2.col(iGrid);
  if(!sample_is_diploid) betaHatLocal.col(2) = betaHat_t3.col(iGrid);
  //
  //
  if (verbose) {
    std::cout << "calculate choice log probs" << std::endl;
  }
  //
  Rcpp::NumericMatrix choice_log_probs_Pm(6, 3);
  Rcpp::NumericVector choice_log_probs_P(6);
  for(ir = 0; ir < 6; ir++) {
    for(i = 0; i < 3; i++) {
        double logC_inside = 0;
        for(iGrid2 = grid_start_0_based; iGrid2 <= grid_end_0_based; iGrid2++) {
            logC_inside += log_cStore(iGrid2, i, ir);
        }
        // go
        choice_log_probs_Pm(ir, i) =                                    \
            log(sum(alphaStore.slice(ir).col(i) % betaHatLocal.col(i))) +
            - logC_before(i) +                  \
            - logC_inside +                     \
            - logC_after(i);
        choice_log_probs_P(ir) += choice_log_probs_Pm(ir, i);
        //
    }
  }
  if (verbose) {
    std::cout << "done calculate choice log probs" << std::endl;
    std::cout << "block_approach = " << block_approach << std::endl;    
  }
  //
  Rcpp::NumericVector choice_log_probs_H;
  if (block_approach == 1) {
      choice_log_probs_H = rcpp_calculate_block_read_label_probabilities(read_start_0_based, read_end_0_based, H, log_prior_probs, rr0);
  } else if (block_approach == 2) {
      choice_log_probs_H = rcpp_calculate_block_read_label_probabilities_using_read_informativeness(read_start_0_based, read_end_0_based, H, log_prior_probs, rr0, read_is_uninformative);
  } else if (block_approach == 4) {
      if (verbose) {
          std::cout << "realculate block read label probs" << std::endl;
      }
      //
      // problem looks to be with this, but why
      //
      choice_log_probs_H = rcpp_calculate_block_read_label_probabilities_using_proposed_H(read_start_0_based, read_end_0_based, proposed_H, log_prior_probs, rr0);
      if (verbose) {
          std::cout << "done realculate block read label probs" << std::endl;
      }
  } else if ( block_approach == 6) {
      choice_log_probs_H = rcpp_calculate_block_read_label_probabilities_using_H_class(read_start_0_based, read_end_0_based, H_class, log_prior_probs, rr, ff);
  }
  //
  if (verbose) {
    std::cout << "add them together" << std::endl;
  }
  // add them together
  Rcpp::NumericVector choice_log_probs(6);
  Rcpp::NumericVector choice_probs(6);  
  for(ir = 0; ir < 6; ir++) {
      choice_log_probs(ir) = choice_log_probs_H(ir) + choice_log_probs_P(ir);
  }
  // remove maximum
  double a = -(Rcpp::max(choice_log_probs));
  for(ir = 0; ir < 6; ir++) {
    choice_log_probs(ir) += a;
  }
  // bound and exponentiate
  for(ir = 0; ir < 6; ir++) {
    if (choice_log_probs(ir) < (-100)) {
      choice_log_probs(ir) = -100;
    }
    choice_probs(ir) = std::exp(choice_log_probs(ir));
  }
  if (ff == 0) {
      // only allow the two options
      choice_probs(1) = 0;
      choice_probs(3) = 0;
      choice_probs(4) = 0;
      choice_probs(5) = 0;      
  }
  // normalize to have sum 1
  double d = (1 /sum(choice_probs));
  for(ir = 0; ir < 6; ir++) {  
    choice_probs(ir) *= d;
  }
  //
  if (verbose) {
    std::cout << "now sample / choose" << std::endl;
  }
  //
  // ir_chosen <- sample(1:6, size = 1, prob = choice_probs)  
  double chance = runif_block(iBlock);
  int ir_chosen = 0;
  Rcpp::NumericVector cumsum_probs(6);
  cumsum_probs(0) = choice_probs(0);
  for(ir = 1; ir < 6; ir++) {
    cumsum_probs(ir) += choice_probs(ir) + cumsum_probs(ir - 1);
  }
  for(int ir = 5; ir >= 0; ir--) {
    if (chance < cumsum_probs(ir)) {
      ir_chosen = ir;
    }
  }
  //int ir_chosen1 = ir_chosen + 1;
  if (verbose) {
    std::cout << "In block" << iBlock << ", see the following probabilities" << std::endl;
    std::cout << choice_probs << std::endl;
    std::cout << "Have selected block relabelling (1-based):" << ir_chosen + 1 << std::endl;
  }
  // in general I need this later
  //
  if (verbose) {
      std::cout << "set up things needed for later" << std::endl;
  }
  // one_based_swap <- as.integer(c(1, 1 + rx[ir_chosen, 1:3], 8 - rx[ir_chosen, 3:1], 8))
  // recall rx is 1-based
    arma::imat rx(6, 3);
    // ugh
    rx(0, 0) = 1; rx(0, 1) = 2; rx(0, 2) = 3;
    rx(1, 0) = 1; rx(1, 1) = 3; rx(1, 2) = 2;
    rx(2, 0) = 2; rx(2, 1) = 1; rx(2, 2) = 3;
    rx(3, 0) = 3; rx(3, 1) = 1; rx(3, 2) = 2;
    rx(4, 0) = 2; rx(4, 1) = 3; rx(4, 2) = 1;
    rx(5, 0) = 3; rx(5, 1) = 2; rx(5, 2) = 1;
  Rcpp::IntegerVector zero_based_swap(8);
  zero_based_swap(0) = 0;
  // be clear about which ir_chosen this is, 1 or 0 based
  // ir_chosen is 1-based here i.e. values of 1 through 6
  zero_based_swap(1) = rx(ir_chosen, 0);
  zero_based_swap(2) = rx(ir_chosen, 1);
  zero_based_swap(3) = rx(ir_chosen, 2);
  zero_based_swap(4) = 7 - rx(ir_chosen, 2);
  zero_based_swap(5) = 7 - rx(ir_chosen, 1);
  zero_based_swap(6) = 7 - rx(ir_chosen, 0);
  zero_based_swap(7) = 7;
  // see test with rcpp_test_one_based_swap in it in gibbs-block-nipt to confirm this code
  //for(int iRead = 0; iRead < H.length(); iRead++) {
  //H_class[iRead] = zero_based_swap[H_class[iRead]];
  //H[iRead] = zero_based_swap[H[iRead]];
  //  }
  //
  if (verbose) {
      std::cout << "save some things for block_results" << std::endl;
  }
  //
  // now have ir_chosen, save results
  //
  int ibr = 2 * iBlock;
  block_results(ibr, 0) = iBlock;
  block_results(ibr, 1) = 0;
  for(ir = 0; ir < 6; ir++) {
    block_results(ibr, 2 + ir) = choice_probs(ir);
  }
  block_results(ibr, 8) = ir_chosen + 1; // internally 0-based, store 1-based, same as R
  //
  block_results(ibr, 9) = choice_log_probs_Pm(ir_chosen, 0);
  block_results(ibr, 10) = choice_log_probs_Pm(ir_chosen, 1);
  block_results(ibr, 11) = choice_log_probs_Pm(ir_chosen, 2);
  for(i = 0; i < 3; i++) {
    block_results(ibr, 12) += choice_log_probs_Pm(ir_chosen, i);
  }
  // so have rest of the H
  // then also change for those H that are changed!
  double x = 0;
  for(i = 0; i < 3; i++) {
      if (sum_H(i) > 0) {
          x += log_prior_probs(i) * sum_H(i);
      }
  }
  int h1, h2;
  for(int iRead = read_start_0_based; iRead <= read_end_0_based; iRead++) {
      h1 = int(H(iRead) - 1);
      if ((block_approach == 1) | (block_approach == 2)) {
          h2 = h1;
      } else if (block_approach == 4) {
          h2 = proposed_H(ir_chosen, iRead) - 1;
      } else if (block_approach == 6) {
          h2 = zero_based_swap[H[iRead]] - 1;
      }
      x -= log_prior_probs(h1);
      x += log_prior_probs(rr0(ir_chosen, h2));
  }
  block_results(ibr, 13) = x;
  block_results(ibr, 14) = block_results(ibr, 12) + block_results(ibr, 13);
  //
  //
  // 
  //
  if (!((ever_changed(0) == 1) | (ir_chosen != 0))) {
    if (verbose) {
      std::cout << "No change warranted" << std::endl;
    }
  } else {
      if (ir_chosen != 0) {
          if (verbose) {
	    std::cout << "Apply block relabelling with ir_chosen=" << ir_chosen << std::endl;
	    std::cout << "i.e. the following relabelling occurs:" << std::endl;
	    for(i = 0; i < 3; i++) {
	      std::cout << "previous label " << i << " -> " << rr0(ir_chosen, i) << std::endl;
	    }
	  }
      }
      ever_changed(0) = 1;
      //
      int iRead = read_start_0_based;
      int wif_read = wif0(iRead);
      int Hl = -1;
      int h = -1;
      for(iGrid2 = grid_start_0_based; iGrid2 <= grid_end_0_based; iGrid2++) {
          //if (verbose) {
          //    std::cout << "changing iGrid2 (0-based) = " <<  iGrid2 << std::endl;
	  //}
          if ((block_approach == 1) | (block_approach == 2)) {          
              eMatGridLocal.col(rr0(ir_chosen, 0)) = eMatGrid_t1.col(iGrid2);
              eMatGridLocal.col(rr0(ir_chosen, 1)) = eMatGrid_t2.col(iGrid2);
              if(!sample_is_diploid) eMatGridLocal.col(rr0(ir_chosen, 2)) = eMatGrid_t3.col(iGrid2);
          } else if ((block_approach == 4) | (block_approach == 6)) {
              eMatGridLocal.fill(1);
              // if necessary, move
              while ((iRead <= (nReads - 1)) & (wif_read < iGrid2)) {
                  iRead += 1;
                  if (iRead < (nReads - 1)) {
                      wif_read = wif0(iRead);
                  }
              }
              while ((iRead <= (nReads - 1)) & (wif_read == iGrid2)) {
                  if (block_approach == 4) {
                      Hl = proposed_H(ir_chosen, iRead) - 1;
                      h = rr0(ir_chosen, Hl);
                  } else {
                      h = zero_based_swap[H[iRead]] - 1; // I think
                  }
                  eMatGridLocal.col(h) %= eMatRead_t.col(iRead);
                  //
                  iRead += 1;
                  if (iRead <= (nReads - 1)) {
                      wif_read = wif0(iRead);
                  }
              }
          }
	  eMatGrid_t1.col(iGrid2) = eMatGridLocal.col(0);
	  eMatGrid_t2.col(iGrid2) = eMatGridLocal.col(1);
	  if(!sample_is_diploid) eMatGrid_t3.col(iGrid2) = eMatGridLocal.col(2);
	  //
	  if (iGrid2 == 0) {
	      // first grid!
	      alphaHat_t1.col(iGrid2) = priorCurrent_m.col(s) % eMatGrid_t1.col(iGrid2);
	      alphaHat_t2.col(iGrid2) = priorCurrent_m.col(s) % eMatGrid_t2.col(iGrid2);
	      if(!sample_is_diploid) alphaHat_t3.col(iGrid2) = priorCurrent_m.col(s) % eMatGrid_t3.col(iGrid2);
	  } else {
	      // normal!
  	      alphaHat_t1.col(iGrid2) = eMatGrid_t1.col(iGrid2) % (
	          transMatRate_tc_H(0, iGrid2 - 1, s) * alphaHat_t1.col(iGrid2 - 1) + \
		  transMatRate_tc_H(1, iGrid2 - 1, s) * alphaMatCurrent_tc.slice(s).col(iGrid2 - 1)
	      );
	      alphaHat_t2.col(iGrid2) = eMatGrid_t2.col(iGrid2) % (
                  transMatRate_tc_H(0, iGrid2 - 1, s) * alphaHat_t2.col(iGrid2 - 1) + \
	          transMatRate_tc_H(1, iGrid2 - 1, s) * alphaMatCurrent_tc.slice(s).col(iGrid2 - 1)
	      );
	      if(!sample_is_diploid) alphaHat_t3.col(iGrid2) = eMatGrid_t3.col(iGrid2) % (
	        transMatRate_tc_H(0, iGrid2 - 1, s) * alphaHat_t3.col(iGrid2 - 1) + \
		transMatRate_tc_H(1, iGrid2 - 1, s) * alphaMatCurrent_tc.slice(s).col(iGrid2 - 1)
	      );
	  }
	  // re-normalize
	  c1(iGrid2) = 1 / arma::sum(alphaHat_t1.col(iGrid2));
	  alphaHat_t1.col(iGrid2) *= c1(iGrid2);
	  c2(iGrid2) = 1 / arma::sum(alphaHat_t2.col(iGrid2));
	  alphaHat_t2.col(iGrid2) *= c2(iGrid2);
      if(!sample_is_diploid){
          c3(iGrid2) = 1 / arma::sum(alphaHat_t3.col(iGrid2));
          alphaHat_t3.col(iGrid2) *= c3(iGrid2);
      }
      }
      // finally, re-do read labels
      int gained = -1;
      for(int iRead = read_start_0_based; iRead <= read_end_0_based; iRead++) {
	int lost = H(iRead) - 1;
        if ((block_approach == 1) | (block_approach == 2)) {
            gained = rr0(ir_chosen, lost);
        } else if (block_approach == 4) {
            gained = rr0(ir_chosen, proposed_H(ir_chosen, iRead) - 1);
        } else if (block_approach == 6) {
            gained = zero_based_swap[H[iRead]] - 1; // very weird?
            H_class[iRead] = zero_based_swap[H_class[iRead]]; // should be OK
        }
	H(iRead) = gained + 1;
	sum_H(gained) += 1.0;
	sum_H(lost) -= 1.0;
      }
  }
  // add H_class stuff
  // eventually, make this faster! will be pretty slow
  // werwer-H_class
  //
  block_results(ibr, 15) = rcpp_get_log_p_H_class(H_class, ff);
  block_results(ibr, 16) = block_results(ibr, 12) + block_results(ibr, 15);
  
  return;
};



//' @export
// [[Rcpp::export]]
void Rcpp_consider_total_relabelling(
    const int iBlock,
    const arma::imat& rr,
    const arma::imat& rr0,
    double ff,
    const Rcpp::NumericVector& log_prior_probs,
    Rcpp::NumericVector& logC_before,
    Rcpp::NumericVector& logC_after,
    const bool& verbose,
    const Rcpp::List& swap_list,
    Rcpp::NumericMatrix& block_results,
    Rcpp::NumericVector& runif_total,
    Rcpp::NumericVector& sum_H,
    Rcpp::IntegerVector& H,
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
    const bool do_checks = false
) {
  //
    if (verbose) {
        std::cout << "Am considering entire relabelling here" << std::endl;
    }
  //
  const int nReads = H.length();
  int iRead, ir, i, h;
  //
  Rcpp::IntegerVector counts = Rcpp::IntegerVector(3);
  Rcpp::NumericMatrix choice_log_probs_Pm(6, 3);
  Rcpp::NumericVector choice_log_probs_P(6);
  Rcpp::NumericVector choice_probs(6);  
  for(ir = 0; ir < 6; ir++) {
    for(i = 0; i < 3; i++) {
        int x = sum_H(rr0(ir, i));
        if (x > 0) {
            choice_log_probs_Pm(ir, i) = log_prior_probs(i) * x;
        }
        choice_log_probs_P(ir) += choice_log_probs_Pm(ir, i);
    }
  }
  double a = -(Rcpp::max(choice_log_probs_P));
  for(ir = 0; ir < 6; ir++) {
    choice_log_probs_P(ir) += a;
  }
  for(ir = 0; ir < 6; ir++) {
    if (choice_log_probs_P(ir) < (-100)) {
      choice_log_probs_P(ir) = -100;
    }
    choice_probs(ir) = std::exp(choice_log_probs_P(ir));
  }
  if (ff == 0) {
      // only allow the two options
      choice_probs(1) = 0;
      choice_probs(3) = 0;
      choice_probs(4) = 0;
      choice_probs(5) = 0;      
  }
  a = (1 /sum(choice_probs));
  for(ir = 0; ir < 6; ir++) {  
    choice_probs(ir) *= a;
  }
  //
  // make selection of ir_chosen
  //
  double chance = runif_total(iBlock);
  int ir_chosen = 0;
  Rcpp::NumericVector cumsum_probs(6);
  cumsum_probs(0) = choice_probs(0);
  for(ir = 1; ir < 6; ir++) {
    cumsum_probs(ir) += choice_probs(ir) + cumsum_probs(ir - 1);
  }
  for(int ir = 5; ir >= 0; ir--) {
    if (chance < cumsum_probs(ir)) {
      ir_chosen = ir;
    }
  }
  //
  int ibr = 2 * iBlock + 1;
  block_results(ibr, 0) = iBlock;
  block_results(ibr, 1) = 1;
  for(ir = 0; ir < 6; ir++) {
    block_results(ibr, 2 + ir) = choice_probs(ir);
  }
  block_results(ibr, 8) = ir_chosen + 1; // internally 0-based, store 1-based
  // define based on previous ones
  block_results(ibr, 12) = 0;
  for(i = 0; i < 3; i++) {
    h = rr0(ir_chosen, i);
    block_results(ibr, 9 + h) = block_results(ibr - 1, 9 + i);
    block_results(ibr, 12) += block_results(ibr, 9 + h);
  }
  // block_results(ibr, 12) done above as sum
  // record this efficiently
  // Htemp = rr[ir_chosen, H]; sum(log_prior_probs[Htemp])
  for(i = 0; i < 3; i++) {
    h = rr0(ir_chosen, i);
    block_results(ibr, 13) += (log_prior_probs(h) * double(sum_H(i)));
  }
  block_results(ibr, 14) = block_results(ibr, 12) + block_results(ibr, 13);
  //
  if (verbose) {
    std::cout << "In block" << iBlock << ", see the following probabilities" << std::endl;
    std::cout << choice_probs << std::endl;
    std::cout << "Have selected block relabelling:" << ir_chosen << std::endl;
  }
  if (ir_chosen == 0) {
    if (verbose) {
      std::cout << "No total relabelling to apply" << std::endl;
      return;
    }
  } else {
      if (verbose) {
	std::cout << "Apply total relabelling" << std::endl;
	std::cout << "i.e. the following relabelling applies " << std::endl;
	for(i = 0; i < 3; i++) {
	  std::cout << "previous label " << i << " -> " << rr0(ir_chosen, i) << std::endl;
	}
      }
      // swap H here
      for(iRead = 0; iRead < nReads; iRead++) {
          H(iRead) = rr0(ir_chosen, H(iRead) - 1) + 1;
      }
      //
      // hope I got this right!
      // cannot seem to get pointer swap to work
      // oh well this is veeeeeery rare on real scale data
      //
      if (ir_chosen > 0) { // re-label still 1-based, i.e. ranges from 1-6
          int relabel = ir_chosen + 1;
          rcpp_apply_mat_relabel(alphaHat_t1, alphaHat_t2, alphaHat_t3, relabel);
          rcpp_apply_mat_relabel(betaHat_t1, betaHat_t2, betaHat_t3, relabel);
          rcpp_apply_vec_relabel(c1, c2, c3, relabel);
          rcpp_apply_mat_relabel(eMatGrid_t1, eMatGrid_t2, eMatGrid_t3, relabel);
      }
    // 
    Rcpp::NumericVector a(3);
    Rcpp::NumericVector b(3);
    Rcpp::NumericVector c(3);
    for(i = 0; i < 3; i++) {
      a(i) = logC_before(i);
      b(i) = logC_after(i);
      c(i) = sum_H(i);
    }
    // easy
    for(i = 0; i < 3; i++) {
      int h = rr0(ir_chosen, i);
      logC_before(h) = a(i);
      logC_after(h) = b(i);
      sum_H(h) = c(i);
    }
  }
  return;
}



//' @export
// [[Rcpp::export]]
void Rcpp_gibbs_block_forward_one(
    Rcpp::IntegerVector& approach2_iRead,
    const int iGrid,
    const int s,
    double ff,    
    arma::cube& alphaStore,
    arma::cube& log_cStore,
    const arma::imat& rr,    
    const arma::imat& rr0,
    arma::mat& eMatGridLocal,
    arma::cube& eMatGridLocalc,
    const arma::cube& transMatRate_tc_H,
    const arma::cube& alphaMatCurrent_tc,
    const arma::mat& priorCurrent_m,
    const Rcpp::LogicalVector& read_is_uninformative,
    const int block_approach,
    Rcpp::IntegerVector& wif0,
    const arma::mat& eMatRead_t,
    const int nReads,
    Rcpp::IntegerVector& H,
    arma::imat& proposed_H,
    Rcpp::IntegerVector& H_class,
    arma::mat& rlc,
    arma::cube& rlcM,
    const Rcpp::NumericMatrix& runif_proposed
) {
    //
    //eMatGridLocal.col(0) = eMatGrid_t1.col(iGrid);
    //eMatGridLocal.col(1) = eMatGrid_t2.col(iGrid);
    //eMatGridLocal.col(2) = eMatGrid_t3.col(iGrid);
    int wif_read = -1;
    int Hl, Hori, Hc;
    const double one_over_K = 1 / double(eMatRead_t.n_rows);
    arma::colvec el, eMatRead_t_col;
    if ((block_approach == 2) | (block_approach == 4)) {
        // so here, just work directly on eMatGridLocal
        if (approach2_iRead(0) <= (nReads - 1)) {        
            wif_read = wif0(approach2_iRead(0));
        }
        while (approach2_iRead(0) <= (nReads - 1) & wif_read < iGrid) {
            approach2_iRead(0) += 1;
            if (approach2_iRead(0) < (nReads - 1)) {
                wif_read = wif0(approach2_iRead(0));
            }
        }
    }
    if (block_approach == 2) {
        while (approach2_iRead(0) <= (nReads - 1) & wif_read == iGrid) {
            el = eMatRead_t.col(approach2_iRead(0));
            Hl = H(approach2_iRead(0)) - 1;
            if (read_is_uninformative(approach2_iRead(0))) {
                eMatGridLocal.col(Hl) /= el;
            }
            approach2_iRead(0) += 1;
            if (approach2_iRead(0) < (nReads - 1)) {
                wif_read = wif0(approach2_iRead(0));
            }
        }
    } else if (block_approach == 4) {
        eMatGridLocalc.fill(1);
        if (approach2_iRead(0) <= (nReads - 1)) {
            wif_read = wif0(approach2_iRead(0));
        }
        while ((approach2_iRead(0) <= (nReads - 1)) & (wif_read == iGrid)) {
            Hori = int(H(approach2_iRead(0)) - 1);
            Hc = H_class(approach2_iRead(0)); // this is its own base, 0 through 7, invariant
            eMatRead_t_col = eMatRead_t.col(approach2_iRead(0));
            for(int ir = 0; ir < 6; ir++) {
                if (Hc == 0) {
                    Hl = Hori;
                } else {
                    //
                    double chance_prob = runif_proposed(ir, approach2_iRead(0));
                    if (chance_prob < rlcM(0, ir, Hc - 1)) {
                        Hl = 0;
                    } else if ( \
                        (rlcM(0, ir, Hc - 1) < chance_prob) &
                        (chance_prob < (rlcM(0, ir, Hc - 1) + rlcM(1, ir, Hc - 1))) \
                    ) {
                        Hl = 1;
                    } else {
                        Hl = 2; // right?
                    }
                }
                proposed_H(ir, approach2_iRead(0)) = Hl + 1; // store as 1-based ALWAYS
                eMatGridLocalc.slice(ir).col(Hl) %= eMatRead_t_col;
            }
            approach2_iRead(0) += 1;
            if (approach2_iRead(0) <= (nReads - 1)) {
                wif_read = wif0(approach2_iRead(0));
            }
        }
    }
    //
    //
    //
    if (iGrid == 0) {
        for(int ir = 0; ir < 6; ir++) {
            if (block_approach == 4) {
                eMatGridLocal = eMatGridLocalc.slice(ir);
            }
            for(int i = 0; i < 3; i++) {
                int h = rr0(ir, i);
                alphaStore.slice(ir).col(h) = priorCurrent_m.col(s) % eMatGridLocal.col(i);
                double d = 1 / sum(alphaStore.slice(ir).col(h));
                log_cStore(iGrid, h, ir) = log(d);
                alphaStore.slice(ir).col(h) = d * alphaStore.slice(ir).col(h);
            }
        }
    } else {
        for(int ir = 0; ir < 6; ir++) {
            // this is OK for ff > 0 OR (ff == 0 & ir == 0 | (ir == 2))
            if ((ff > 0) | ((ff == 0) & ((ir == 0) | (ir == 2)))) {
                if (block_approach == 4) {
                    eMatGridLocal = eMatGridLocalc.slice(ir);
                }
                for(int i = 0; i < 3; i++) {
                    int h = rr0(ir, i);
                    alphaStore.slice(ir).col(h) = eMatGridLocal.col(i) % ( \
                        transMatRate_tc_H(0, iGrid - 1, s) * alphaStore.slice(ir).col(h) + \
                        transMatRate_tc_H(1, iGrid - 1, s) * one_over_K \
                    );
                    // transMatRate_tc_H(1, iGrid - 1, s) * alphaMatCurrent_tc.slice(s).col(iGrid - 1)
                    double d = 1 / sum(alphaStore.slice(ir).col(h));
                    log_cStore(iGrid, h, ir) = log(d);
                    alphaStore.slice(ir).col(h) = d * alphaStore.slice(ir).col(h);
                }
            }
        }
    }
    return;
}

//' @export
// [[Rcpp::export]]
void Rcpp_reset_local_variables(
    const int iGrid,
    const bool verbose,
    arma::mat& alphaHatLocal,    
    arma::cube& alphaStore,
    arma::mat& alphaHat_t1,
    arma::mat& alphaHat_t2,
    arma::mat& alphaHat_t3,
    arma::rowvec& c1,
    arma::rowvec& c2,
    arma::rowvec& c3,
    arma::cube& log_cStore,
    const bool sample_is_diploid
) {
    if (verbose) {
        std::cout << "Reset alphaHatLocal, cStore" << std::endl;
    }
    // reset alphaStore, cStore
    // they are all using the current data
    alphaHatLocal.col(0) = alphaHat_t1.col(iGrid);
    alphaHatLocal.col(1) = alphaHat_t2.col(iGrid);
    if(!sample_is_diploid) alphaHatLocal.col(2) = alphaHat_t3.col(iGrid);
    else alphaHatLocal.col(2).fill(0);
    Rcpp::NumericVector cLocal(3);
    cLocal(0) = c1(iGrid);
    cLocal(1) = c2(iGrid);
    cLocal(2) = c3(iGrid);    
    for(int ir = 0; ir < 6; ir++) {
        for(int i = 0; i < 3; i++) {
            int h = i;
            alphaStore.slice(ir).col(h) = alphaHatLocal.col(i);
            log_cStore(iGrid, h, ir) = log(cLocal(i));
        }
    }
    return;
}

               
double ceiling_point5(double x) {
    if (double(int(x)) < x) {
        return(x + 0.5);
    } else {
        return(x);
    }
}



//' @export
// [[Rcpp::export]]
Rcpp::List Rcpp_make_gibbs_considers(
    Rcpp::IntegerVector& blocked_snps,
    const Rcpp::IntegerVector& grid,
    Rcpp::IntegerVector& wif0,
    const int nGrids,
    const bool do_removal = true,
    const bool verbose = false
) {
    const int nSNPs = blocked_snps.length();
    int n_blocks = blocked_snps(nSNPs - 1) + 1;
    const int nReads = wif0.length();
    //
    if (verbose) {
        std::cout << "check assumption" << std::endl;
    }
    //
    Rcpp::List to_return;
    int iSNP, iBlock, a, b;
    for(iSNP = 0; iSNP < (nSNPs - 1); iSNP++) {
        if ((blocked_snps(iSNP + 1) - blocked_snps(iSNP)) > 1) {
            std::cout << "failure of assumption on blocked snps" << std::endl;
            return to_return;
        }
    }
    //
    if (verbose) {
        std::cout << "first bit" << std::endl;
    }
    //
    Rcpp::IntegerVector consider_snp_start_0_based(n_blocks);
    Rcpp::IntegerVector consider_snp_end_0_based(n_blocks);
    int cur_block = -1;
    int start = 0;
    bool record = false;
    iBlock = 0;
    for(iSNP = 0; iSNP < nSNPs; iSNP++) {
        cur_block = blocked_snps(iSNP); // 0-based
        if (iSNP == (nSNPs - 1)) {
            record = true;
        } else if (blocked_snps(iSNP) < blocked_snps(iSNP + 1)) {
            record = true;
        } else {
            record = false;
        }
        if (record) {
            consider_snp_start_0_based(iBlock) = start;
            consider_snp_end_0_based(iBlock) = iSNP;
            start = iSNP + 1;
            iBlock = iBlock + 1;
        }
    }
    //
    if (verbose) {
        std::cout << "grid" << std::endl;
    }
    //
    Rcpp::IntegerVector consider_grid_start_0_based (n_blocks);
    Rcpp::IntegerVector consider_grid_end_0_based(n_blocks);
    for(iBlock = 0; iBlock < n_blocks; iBlock++) {
        consider_grid_start_0_based(iBlock) = grid(consider_snp_start_0_based(iBlock));
        consider_grid_end_0_based(iBlock) = grid(consider_snp_end_0_based(iBlock));
    }
    //
    if (verbose) {
        std::cout << "reads" << std::endl;
    }
    //
    Rcpp::IntegerVector blocked_grid(nGrids);
    for(iBlock = 0; iBlock < n_blocks; iBlock++) {
        a = consider_grid_start_0_based(iBlock);
        b = consider_grid_end_0_based(iBlock);
        for(int i = a; i <= b; i++) {
            blocked_grid(i) = iBlock;
        }
    }
    //
    if (verbose) {
        std::cout << "setup" << std::endl;
    }
    //
    Rcpp::IntegerVector consider_reads_start_0_based(n_blocks);
    consider_reads_start_0_based.fill(-1);
    Rcpp::IntegerVector consider_reads_end_0_based(n_blocks);
    consider_reads_end_0_based.fill(-1);
    int previous_block_first_iRead = 0;
    int previous_grid = wif0(previous_block_first_iRead);
    int previous_block = blocked_grid(previous_grid);
    int this_grid, this_block;
    for(int this_iRead = 1; this_iRead < nReads; this_iRead++) {
        this_grid = wif0(this_iRead);
        this_block = blocked_grid(this_grid);
        if (this_iRead == (nReads - 1)) {
            consider_reads_start_0_based(this_block) = previous_block_first_iRead;
            consider_reads_end_0_based(this_block) = this_iRead;
        } else if (previous_block < this_block) {
            consider_reads_start_0_based(previous_block) = previous_block_first_iRead;
            consider_reads_end_0_based(previous_block) = this_iRead - 1; // this minus 1 is because it is before
            // reset
            previous_block_first_iRead = this_iRead;
            previous_block = blocked_grid(wif0(this_iRead));
            previous_grid = this_grid;
        }
    }
    //
    //
    if (do_removal) {
        if (verbose) {
            std::cout << "removal" << std::endl;
        }
        Rcpp::LogicalVector remove(n_blocks);
        for(iBlock = 0; iBlock < n_blocks; iBlock++) {
            if (consider_reads_start_0_based(iBlock) == -1) {
                remove(iBlock) = true;
            } else {
                remove(iBlock) = false;
            }
        }
        int n_to_remove = Rcpp::sum(remove);
        if (n_to_remove > 0) {
            if (verbose) {
                std::cout << "something to remove" << std::endl;
            }
            int n_new_blocks = n_blocks - n_to_remove;
            Rcpp::IntegerVector w(n_to_remove);
            a = 0;
            for(iBlock = 0; iBlock < n_blocks; iBlock++) {
                if (remove(iBlock)) {
                    w(a) = iBlock;
                    a += 1;
                }
            }
            //
            if (verbose) {
                std::cout << "loop bit" << std::endl;
                std::cout << "also, w is as follows" << std::endl;
                std::cout << w << std::endl;
            }
            //
            int jBefore = 0;
            bool todo = false;
            for(int jNow = 0; jNow < n_to_remove; jNow++) {
                if (verbose) {
                    std::cout << "jNow = " << jNow << std::endl;
                }
                if (jNow == (n_to_remove - 1)) {
                    todo = true;
                } else {
                    if ((w(jNow + 1) - w(jNow)) == 1) {
                        todo = false;
                        jBefore -= 1;
                    } else {
                        todo = true;
                    }
                }
                if (todo) {
                    if (verbose) {
                        std::cout << "todo, jBefore = " << jBefore << ", jNow = " << jNow << std::endl;
                    }
                    int s1 = w(jBefore);
                    int e1 = w(jNow);
                    if (verbose) {
                        std::cout << "s1 = " << s1 << ", e1 = " << e1 << std::endl;
                    }
                    double x = ceiling_point5(0.5 * double(consider_grid_start_0_based(s1) + consider_grid_end_0_based(e1)));
                    double y = ceiling_point5(0.5 * double(consider_snp_start_0_based(s1) + consider_snp_end_0_based(e1)));
                    //
                    if (s1 == 0) {
                        // if need to remove the first one, then set this to 0
                        s1 = 1;
                        x = 0;
                        y = 0;
                    }
                    if (e1 == (n_blocks - 1)) {
                        e1 = e1 - 1;
                        x = consider_grid_end_0_based(n_blocks - 1);
                        y = consider_snp_end_0_based(n_blocks - 1);
                    }
                    if (verbose) {
                        std::cout << "x = " << x << ", y = " << y << std::endl;
                    }
                    consider_grid_start_0_based(e1 + 1) = x;
                    consider_grid_end_0_based(s1 - 1) = x - 1;
                    //
                    consider_snp_start_0_based(e1 + 1) = y;
                    consider_snp_end_0_based(s1 - 1) = y - 1;
                    jBefore = jNow;
                }
                jBefore += 1;
            }
            //
            if (verbose) {
                std::cout << "rebuild" << std::endl;
            }
            Rcpp::IntegerVector new_consider_reads_start_0_based(n_new_blocks);
            Rcpp::IntegerVector new_consider_reads_end_0_based(n_new_blocks);
            Rcpp::IntegerVector new_consider_grid_start_0_based(n_new_blocks);
            Rcpp::IntegerVector new_consider_grid_end_0_based(n_new_blocks);
            Rcpp::IntegerVector new_consider_snp_start_0_based(n_new_blocks);
            Rcpp::IntegerVector new_consider_snp_end_0_based(n_new_blocks);
            //
            int i_prev_block = -1;
            for(iBlock = 0; iBlock < n_blocks; iBlock++) {
                if (!remove(iBlock)) {
                    i_prev_block += 1;
                    //
                    new_consider_reads_start_0_based(i_prev_block) = consider_reads_start_0_based(iBlock);
                    new_consider_reads_end_0_based(i_prev_block) = consider_reads_end_0_based(iBlock);
                    //
                    new_consider_grid_start_0_based(i_prev_block) = consider_grid_start_0_based(iBlock);
                    new_consider_grid_end_0_based(i_prev_block) = consider_grid_end_0_based(iBlock);
                    //
                    new_consider_snp_start_0_based(i_prev_block) = consider_snp_start_0_based(iBlock);
                    new_consider_snp_end_0_based(i_prev_block) = consider_snp_end_0_based(iBlock);
                }
            }
            //
            consider_reads_start_0_based = new_consider_reads_start_0_based;
            consider_reads_end_0_based = new_consider_reads_end_0_based;
            consider_snp_start_0_based = new_consider_snp_start_0_based;
            consider_snp_end_0_based = new_consider_snp_end_0_based;            
            consider_grid_start_0_based = new_consider_grid_start_0_based;
            consider_grid_end_0_based = new_consider_grid_end_0_based;            
        }
    }
    //
    // last bit
    //
    n_blocks = consider_snp_end_0_based.length();
    Rcpp::IntegerVector consider_grid_where_0_based(nGrids);
    consider_grid_where_0_based.fill(-1);
    for(iBlock = 0; iBlock < n_blocks; iBlock++) {
        int b = consider_grid_end_0_based(iBlock);
        consider_grid_where_0_based(b) = iBlock;
    }
    //
    // return
    //
    to_return.push_back(consider_snp_start_0_based, "consider_snp_start_0_based");
    to_return.push_back(consider_snp_end_0_based, "consider_snp_end_0_based");
    to_return.push_back(consider_grid_start_0_based, "consider_grid_start_0_based");
    to_return.push_back(consider_grid_end_0_based, "consider_grid_end_0_based");
    to_return.push_back(consider_reads_start_0_based, "consider_reads_start_0_based");
    to_return.push_back(consider_reads_end_0_based, "consider_reads_end_0_based");
    to_return.push_back(consider_grid_where_0_based, "consider_grid_where_0_based");
    to_return.push_back(n_blocks, "n_blocks");
    return to_return;
}



//' @export
// [[Rcpp::export]]
void Rcpp_fill_rlcM(
    arma::cube& rlcM,
    const arma::mat& rlc,
    const arma::imat& rr0
) {
    // manually define rlcI
    arma::imat rlcI(7, 3);    
    rlcI(0, 0) = 1;
    rlcI(0, 1) = 0;
    rlcI(0, 2) = 0;
    //
    rlcI(1, 0) = 0;
    rlcI(1, 1) = 1;
    rlcI(1, 2) = 0;
    //
    rlcI(2, 0) = 0;
    rlcI(2, 1) = 0;
    rlcI(2, 2) = 1;
    //
    rlcI(3, 0) = 1;
    rlcI(3, 1) = 1;
    rlcI(3, 2) = 0;
    //
    rlcI(4, 0) = 1;
    rlcI(4, 1) = 0;
    rlcI(4, 2) = 1;
    //
    rlcI(5, 0) = 0;
    rlcI(5, 1) = 1;
    rlcI(5, 2) = 1;
    //
    rlcI(6, 0) = 1;
    rlcI(6, 1) = 1;
    rlcI(6, 2) = 1;
    // make rlcI,
    Rcpp::NumericVector x(3);
    Rcpp::IntegerVector future_labels(3);
    int hcF, i, j, d;
    for(int hcC = 0; hcC < 7; hcC++) {
        for(int ir = 0; ir < 6; ir++) {
            x.fill(0);
            for(i = 0; i < 3; i++) {
                if (rlc(hcC, i) > 0) {
                    x(rr0(ir, i)) = 1;
                }
            }
            // future one
            hcF = -1;
            for(i = 0; i < 7; i++) {
                d = 0;
                for(j = 0; j < 3; j++) {
                    if ((rlcI(i, j) == 1) == (x(j) == 1)) {
                        d += 1;
                    }
                }
                if (d == 3) {
                    if ((-1) < hcF) {
                        std::cout << "Rcpp: Faulty assumption in determining hcF while filling rlcM" << std::endl;
                    }
                    hcF = i;
                }
            }
            if (hcF == -1) {
                std::cout << "Rcpp: Faulty assumption in setting hfC while filling rlcM" << std::endl;
            }
            for(j = 0; j < 3; j++) {
                rlcM(j, ir, hcC) = rlc(hcF, rr0(ir, j));
            }
        }
    }
    return;
}



//' @export
// [[Rcpp::export]]
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
    Rcpp::LogicalVector& grid_has_read,
    double ff,
    int s, // this is 0-based
    const double maxEmissionMatrixDifference,
    const Rcpp::List& sampleReads,    
    const arma::cube& alphaMatCurrent_tc,
    const arma::mat& priorCurrent_m,
    const arma::cube& transMatRate_tc_H,
    const int maxDifferenceBetweenReads,
    const int Jmax,
    std::string& prev_section,
    std::string& next_section, 
    const int suppressOutput,
    double& prev,
    const bool sample_is_diploid,
    bool do_checks = false,
    Rcpp::List initial_package = R_NilValue,
    bool verbose = false,
    Rcpp::List fpp_stuff = R_NilValue,
    bool use_cpp_bits_in_R = true,
    int block_approach = 6,
    bool consider_total_relabelling = false,
    const bool resample_H_using_H_class = true
) {
    //
    // note for resample_H_using_Hclass
    // turn this ON for default use
    // but for testing, keep it off
    //
    next_section="block gibbs - begin";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    Rcpp::List to_return;
    Rcpp::List all_packages;
    const Rcpp::LogicalVector read_is_uninformative(1);
    //
    //
    //
    const int K = alphaMatCurrent_tc.n_rows;
    // const int nSNPs = grid.length();
    const int nGrids = alphaMatCurrent_tc.n_cols + 1;
    // const int S = alphaMatCurrent_tc.n_slices + 1;
    //
    const Rcpp::NumericVector prior_probs = NumericVector::create(0.5, (1 - ff) / 2, (ff / 2)); 
    Rcpp::NumericVector log_prior_probs = NumericVector::create(log(0.5), log((1 - ff) / 2), log((ff / 2)));
    //
    next_section="block gibbs - make considers";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    Rcpp::List out = Rcpp_make_gibbs_considers(
        blocked_snps,
        grid,
        wif0,
        nGrids
    );
    // unset
    Rcpp::IntegerVector consider_reads_start_0_based = as<Rcpp::IntegerVector>(out["consider_reads_start_0_based"]);
    Rcpp::IntegerVector consider_reads_end_0_based = as<Rcpp::IntegerVector>(out["consider_reads_end_0_based"]);
    Rcpp::IntegerVector consider_grid_start_0_based = as<Rcpp::IntegerVector>(out["consider_grid_start_0_based"]);
    Rcpp::IntegerVector consider_grid_end_0_based = as<Rcpp::IntegerVector>(out["consider_grid_end_0_based"]);
    Rcpp::IntegerVector consider_snp_start_0_based = as<Rcpp::IntegerVector>(out["consider_snp_start_0_based"]);
    Rcpp::IntegerVector consider_snp_end_0_based = as<Rcpp::IntegerVector>(out["consider_snp_end_0_based"]);
    Rcpp::IntegerVector consider_grid_where_0_based = as<Rcpp::IntegerVector>(out["consider_grid_where_0_based"]);
    const int n_blocks = as<int>(out["n_blocks"]);
    //
    //
    next_section="block gibbs - check considers";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    if (verbose) {
        if (consider_grid_start_0_based(0) != 0) {
            std::cout << "problem with consider grid start entry" << std::endl;
            return(to_return);
        }
        if (consider_grid_end_0_based(n_blocks - 1) != (nGrids - 1)) {
            std::cout << "problem with consider grid end entry" << std::endl;
            return(to_return);            
        }
        for(int i = 0; i < (n_blocks - 1); i++) {
            if (1 != (consider_grid_start_0_based(i + 1) - consider_grid_end_0_based(i))) {
                std::cout << "problem with difference with consider grids" << std::endl;
                return(to_return); 
            }
        }
    }
    //
    next_section="block gibbs - initialize containers";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    Rcpp::NumericMatrix block_results;
    Rcpp::CharacterVector block_results_columns;
    block_results_columns = CharacterVector::create("iBlock", "total", "p1", "p2", "p3", "p4", "p5", "p6", "ir_chosen", "p_O1_given_H1_L", "p_O2_given_H2_L", "p_O3_given_H3_L", "p_O_given_H_L", "p_H_given_L", "p_H_given_O_L_up_to_C", "p_H_class_given_L", "p_H_class_given_O_L_up_to_C");
    block_results = Rcpp::NumericMatrix(2 * n_blocks, block_results_columns.length());
    block_results.fill(0);
    colnames(block_results) = block_results_columns; // neat
    //
    if (verbose) {
        std::cout << "make rr" << std::endl;
    }
    //
    arma::imat rr(6, 3);
    // ugh
    rr(0, 0) = 1; rr(0, 1) = 2; rr(0, 2) = 3;
    rr(1, 0) = 1; rr(1, 1) = 3; rr(1, 2) = 2;
    rr(2, 0) = 2; rr(2, 1) = 1; rr(2, 2) = 3;
    rr(3, 0) = 2; rr(3, 1) = 3; rr(3, 2) = 1;
    rr(4, 0) = 3; rr(4, 1) = 1; rr(4, 2) = 2;
    rr(5, 0) = 3; rr(5, 1) = 2; rr(5, 2) = 1;
    arma::imat rr0 = rr - 1; 
    Rcpp::List swap_list; // not used in cpp
    //
    if (verbose) {
        std::cout << "make rlc" << std::endl;
    }
    // for convenience also make "p"
    const Rcpp::NumericVector p = NumericVector::create(0.5, (1 - ff) / 2, (ff / 2));
    arma::mat rlc(7, 3);
    //
    rlc(0, 0) = 1;
    rlc(0, 1) = 0;
    rlc(0, 2) = 0;
    //
    rlc(1, 0) = 0;
    rlc(1, 1) = 1;
    rlc(1, 2) = 0;
    //
    rlc(2, 0) = 0;
    rlc(2, 1) = 0;
    rlc(2, 2) = 1;
    //
    rlc(3, 0) = p(0) / (p(0) + p(1));
    rlc(3, 1) = p(1) / (p(0) + p(1));
    rlc(3, 2) = 0;
    //
    rlc(4, 0) = p(0) / (p(0) + p(2));
    rlc(4, 1) = 0;
    rlc(4, 2) = p(2) / (p(0) + p(2));
    //
    rlc(5, 0) = 0;
    rlc(5, 1) = p(1) / (p(1) + p(2));
    rlc(5, 2) = p(2) / (p(1) + p(2));
    //
    rlc(6, 0) = p(0);
    rlc(6, 1) = p(1);
    rlc(6, 2) = p(2);
    //
    // am here, make rlcM somehow - does this need to be tested, a function, etc
    if (verbose) {
        std::cout << "make rlcm" << std::endl;
    }
    arma::cube rlcM(3, 6, 7);
    Rcpp_fill_rlcM(rlcM, rlc, rr0);
    //
    if (verbose) {
        std::cout << "more initializing" << std::endl;
    }
    //
    Rcpp::NumericVector logC_before(3);
    Rcpp::NumericVector logC_after(3);    
    logC_after(0) = sum(log(c1));
    logC_after(1) = sum(log(c2));
    logC_after(2) = sum(log(c3));
    //
    int iBlock = 0;
    Rcpp::NumericVector ever_changed(1);
    ever_changed(0) = 0; // argh
    //
    Rcpp::NumericVector sum_H(3);
    const int nReads = H.length();
    for(int iRead = 0; iRead < nReads; iRead++) {
        sum_H(H(iRead) - 1) += 1; // H here is 0-based
    }
    arma::imat proposed_H(6, nReads);
    proposed_H.fill(-1);
    Rcpp::IntegerVector approach2_iRead(1);
    approach2_iRead = 0; // here 0-based
    next_section="block gibbs - go forwards";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    // more initializations
    //
    arma::cube alphaStore(K, 3, 6);
    arma::mat alphaHatLocal(K, 3);
    arma::mat betaHatLocal(K, 3);
    arma::mat eMatGridLocal(K, 3);
    arma::cube eMatGridLocalc(K, 3, 6);    
    arma::cube log_cStore(nGrids, 3, 6);
    for(int iGrid = 0; iGrid < nGrids; iGrid++) {    
        //if (verbose) {
        //    std::cout << "iBlock = " << iBlock << ", iGrid = " << iGrid << ", nGrids = " << nGrids << std::endl;
        //}
        eMatGridLocal.col(0) = eMatGrid_t1.col(iGrid);
        eMatGridLocal.col(1) = eMatGrid_t2.col(iGrid);
        if(!sample_is_diploid) eMatGridLocal.col(2) = eMatGrid_t3.col(iGrid);
        //
        // go forward one
        //
        Rcpp_gibbs_block_forward_one(approach2_iRead, iGrid, s, ff, alphaStore, log_cStore, rr, rr0, eMatGridLocal, eMatGridLocalc, transMatRate_tc_H, alphaMatCurrent_tc, priorCurrent_m, read_is_uninformative, block_approach, wif0, eMatRead_t, nReads, H, proposed_H, H_class, rlc, rlcM, runif_proposed);
        //
        if ((-1) < consider_grid_where_0_based(iGrid)) {
            iBlock = consider_grid_where_0_based(iGrid);
            int grid_start_0_based = consider_grid_start_0_based(iBlock);
            int grid_end_0_based = consider_grid_end_0_based(iBlock);
            int read_start_0_based = consider_reads_start_0_based(iBlock);
            int read_end_0_based = consider_reads_end_0_based(iBlock);
            //
            // check block relabellings
            //
            Rcpp_consider_block_relabelling(iBlock, runif_block, sum_H, s, rr, rr0, ff, log_prior_probs, logC_before, logC_after, verbose, swap_list, eMatGridLocal, betaHatLocal, iGrid, grid_start_0_based, grid_end_0_based, read_start_0_based, read_end_0_based, wif0, log_cStore, alphaStore, read_is_uninformative, block_approach, do_checks, all_packages, block_results, ever_changed, transMatRate_tc_H, alphaMatCurrent_tc, priorCurrent_m, fpp_stuff, H, H_class, proposed_H, nReads, eMatRead_t, alphaHat_t1, betaHat_t1, c1, eMatGrid_t1, alphaHat_t2, betaHat_t2, c2, eMatGrid_t2, alphaHat_t3, betaHat_t3, c3, eMatGrid_t3, sample_is_diploid);
            //
            // total relabelling 
            //
            if ((ff > 0) & consider_total_relabelling) {
                // do not bother for ff = 0, this means diploid, so switches are 50-50 and pointless
                Rcpp_consider_total_relabelling(iBlock, rr, rr0, ff, log_prior_probs, logC_before, logC_after, verbose, swap_list, block_results, runif_total, sum_H, H, alphaHat_t1, betaHat_t1, c1, eMatGrid_t1, alphaHat_t2, betaHat_t2, c2, eMatGrid_t2, alphaHat_t3, betaHat_t3, c3, eMatGrid_t3);
            }
            //
            // if this is NOT the last block, reset, and checks
            //
            if ((iBlock + 1) < n_blocks) {
                // note above - iBlock + 1, first one, is 0 vs 1-based. second one is need to check future
                Rcpp_reset_local_variables(iGrid, verbose, alphaHatLocal, alphaStore, alphaHat_t1, alphaHat_t2, alphaHat_t3, c1, c2, c3, log_cStore, sample_is_diploid);
            }
            //
            //
            //
            //if (verbose) {
            //    std::cout << "add to logC from previous run" << std::endl;
            //}
            for(int iGrid2 = grid_start_0_based; iGrid2 <= grid_end_0_based; iGrid2++) {
                logC_before(0) += log(c1(iGrid2));
                logC_before(1) += log(c2(iGrid2));
                logC_before(2) += log(c3(iGrid2));
            }
        }
        logC_after(0) -= log(c1(iGrid));
        logC_after(1) -= log(c2(iGrid));
        logC_after(2) -= log(c3(iGrid));        
    }
    if ((ff > 0) && block_approach == 6 && resample_H_using_H_class) {
        // now have new H_class
        // easy to sample H using that
        // then easy to make eMatGrid_t1 etc
        // then can do forwards
        // then can sample new read labels based on that 
        // sample new read labels based on that
        // then re-do forwards
        // I can probably do this in-line but meh just do this for now
        next_section="block gibbs - approach 6 - re-sample H, then re-do eMatGrid and forward";
        prev=print_times(prev, suppressOutput, prev_section, next_section);
        prev_section=next_section;
        //
        //
        rcpp_sample_H_using_H_class(H_class, H, ff);
        //
        eMatGrid_t1.fill(1);
        eMatGrid_t2.fill(1);
        eMatGrid_t3.fill(1);
        const bool use_all_reads = false;
        const Rcpp::NumericVector alphaStart_local = 0;        
        const bool initialize_only_local = false; // initialize the whole thing!!!
        bool run_fb_subset_local = false;
        //
        const int run_fb_grid_offset = 0;  // I htink
        const bool rescale_eMatGrid_t = true; // done safely I think, just multiplied (see main gibbs-nipt)
        const bool bound_eMatGrid_t = false; // similarly, see main bit of script
        //
        rcpp_make_eMatGrid_t(eMatGrid_t1, eMatRead_t, H, sampleReads, 1, nGrids, prev, suppressOutput, prev_section, next_section, run_fb_grid_offset, use_all_reads, bound_eMatGrid_t, maxEmissionMatrixDifference, rescale_eMatGrid_t);
        rcpp_make_eMatGrid_t(eMatGrid_t2, eMatRead_t, H, sampleReads, 2, nGrids, prev, suppressOutput, prev_section, next_section, run_fb_grid_offset, use_all_reads, bound_eMatGrid_t, maxEmissionMatrixDifference, rescale_eMatGrid_t);
        rcpp_make_eMatGrid_t(eMatGrid_t3, eMatRead_t, H, sampleReads, 3, nGrids, prev, suppressOutput, prev_section, next_section, run_fb_grid_offset, use_all_reads, bound_eMatGrid_t, maxEmissionMatrixDifference, rescale_eMatGrid_t);
        //
        alphaHat_t1.fill(1); betaHat_t1.fill(1); c1.fill(1);
        Rcpp_run_forward_haploid(alphaHat_t1, c1, eMatGrid_t1, alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s, alphaStart_local, run_fb_subset_local, initialize_only_local);                
        //
        alphaHat_t2.fill(1); betaHat_t2.fill(1); c2.fill(1);
        Rcpp_run_forward_haploid(alphaHat_t2, c2, eMatGrid_t2, alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s, alphaStart_local, run_fb_subset_local, initialize_only_local);        
        //
        alphaHat_t3.fill(1); betaHat_t3.fill(1); c3.fill(1);        
        Rcpp_run_forward_haploid(alphaHat_t3, c3, eMatGrid_t3, alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s, alphaStart_local, run_fb_subset_local, initialize_only_local);
    }
    //
    // re-run backward
    //
    next_section="block gibbs - finalize by re-running backwards";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    betaHat_t1.col(nGrids - 1).fill(c1(nGrids-1));
    betaHat_t2.col(nGrids - 1).fill(c2(nGrids-1));
    //
    Rcpp_run_backward_haploid(betaHat_t1, c1, eMatGrid_t1, alphaMatCurrent_tc, transMatRate_tc_H, s);
    Rcpp_run_backward_haploid(betaHat_t2, c2, eMatGrid_t2, alphaMatCurrent_tc, transMatRate_tc_H, s);
    // yes for the long term, not sure what I will ultimately do with this though
    Rcpp_run_backward_haploid_QUILT_faster(betaHat_t1, c1, eMatGrid_t1, transMatRate_tc_H, grid_has_read, s);    
    Rcpp_run_backward_haploid_QUILT_faster(betaHat_t2, c2, eMatGrid_t2, transMatRate_tc_H, grid_has_read, s);
    if(!sample_is_diploid) {
        betaHat_t3.col(nGrids - 1).fill(c3(nGrids-1));
        Rcpp_run_backward_haploid(betaHat_t3, c3, eMatGrid_t3, alphaMatCurrent_tc, transMatRate_tc_H, s);
        Rcpp_run_backward_haploid_QUILT_faster(betaHat_t3, c3, eMatGrid_t3, transMatRate_tc_H, grid_has_read, s);
    }
    if (verbose) {
        std::cout << "done rcpp block gibbs sampling" << std::endl;
    }
    to_return.push_back(block_results, "block_results");
    return(to_return);


}





//' @export
// [[Rcpp::export]]
Rcpp::List Rcpp_shard_block_gibbs_resampler(
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
    double ff,
    const arma::mat& eMatRead_t,
    Rcpp::IntegerVector& blocked_snps,
    const Rcpp::IntegerVector& grid,
    Rcpp::IntegerVector& wif0,
    int s, // this is 0-based
    const arma::cube& alphaMatCurrent_tc,
    const arma::mat& priorCurrent_m,
    const arma::cube& transMatRate_tc_H,
    bool do_checks = false,
    Rcpp::List initial_package = R_NilValue,
    bool verbose = false,
    Rcpp::List fpp_stuff = R_NilValue,
    bool shard_check_every_pair = false
) {

    const int K = alphaMatCurrent_tc.n_rows;
    const int nGrids = alphaMatCurrent_tc.n_cols + 1;
    const int nReads = H.length();
    Rcpp::List to_return;
    Rcpp::IntegerVector consider_reads_start_0_based;
    Rcpp::IntegerVector consider_reads_end_0_based;
    Rcpp::IntegerVector consider_grid_start_0_based;
    Rcpp::IntegerVector consider_grid_end_0_based;
    Rcpp::IntegerVector consider_snp_start_0_based;
    Rcpp::IntegerVector consider_snp_end_0_based;
    Rcpp::IntegerVector consider_grid_where_0_based;
    int n_blocks;
    if (!shard_check_every_pair) {
        if (verbose) {
            std::cout << "make considers" << std::endl;
        }
        //
        Rcpp::List out = Rcpp_make_gibbs_considers(
            blocked_snps,
            grid,
            wif0,
            nGrids
        );
        // unset
        consider_reads_start_0_based = as<Rcpp::IntegerVector>(out["consider_reads_start_0_based"]);
        consider_reads_end_0_based = as<Rcpp::IntegerVector>(out["consider_reads_end_0_based"]);
        consider_grid_start_0_based = as<Rcpp::IntegerVector>(out["consider_grid_start_0_based"]);
        consider_grid_end_0_based = as<Rcpp::IntegerVector>(out["consider_grid_end_0_based"]);
        consider_snp_start_0_based = as<Rcpp::IntegerVector>(out["consider_snp_start_0_based"]);
        consider_snp_end_0_based = as<Rcpp::IntegerVector>(out["consider_snp_end_0_based"]);
        consider_grid_where_0_based = as<Rcpp::IntegerVector>(out["consider_grid_where_0_based"]);
        n_blocks = as<int>(out["n_blocks"]);
    } else {
        n_blocks = nGrids;
    }
    //
    // make output container
    //
    Rcpp::NumericMatrix shard_block_results;
    Rcpp::CharacterVector shard_block_columns;
    if (ff == 0) {
        shard_block_columns = CharacterVector::create("iBlock", "p_stay", "p_flip", "pA1", "pA2", "pB1", "pB2", "flip_mode", "p_O_stay", "p_O_flip");
    } else {
        shard_block_columns = CharacterVector::create("iBlock", "ir_chosen", "p1", "p2", "p3", "p4", "p5", "p6");
    }
    shard_block_results = Rcpp::NumericMatrix(n_blocks - 1, shard_block_columns.length());
    shard_block_results.fill(0);
    colnames(shard_block_results) = shard_block_columns; // neat
    Rcpp::NumericVector runif_block = Rcpp::runif(n_blocks - 1);

    double minus_log_c_sum = 0; // forward_one does something different here
    double minus_log_c1_sum = 0;
    double minus_log_c2_sum = 0;
    double minus_log_c3_sum = 0;    
    double minus_log_original_c1_sum = 0;
    double minus_log_original_c2_sum = 0;
    double minus_log_original_c3_sum = 0;    
    double original_c1_this_grid, original_c2_this_grid, original_c3_this_grid;
    double pA1, pA2, pB1, pB2;
    for(int iGrid2 = 0; iGrid2 < nGrids; iGrid2++) {
        minus_log_original_c1_sum -= log(c1(iGrid2));
        minus_log_original_c2_sum -= log(c2(iGrid2));
        if(ff>0) minus_log_original_c3_sum -= log(c3(iGrid2));        
    }
    
    bool in_flip_mode = false;
    int iRead = 0;
    int iGrid, iGridConsider, k;
    int split_grid = -1;
    arma::colvec emat_temp_col;
    arma::colvec e1, e2, e3;
    bool done_reads = false;
    double calculated_difference, probs1, probs2, probs_sum, x1, x2;

    // do the ff stuff here
    arma::imat rr(6, 3);
    // ugh
    rr(0, 0) = 1; rr(0, 1) = 2; rr(0, 2) = 3;
    rr(1, 0) = 1; rr(1, 1) = 3; rr(1, 2) = 2;
    rr(2, 0) = 2; rr(2, 1) = 1; rr(2, 2) = 3;
    rr(3, 0) = 2; rr(3, 1) = 3; rr(3, 2) = 1;
    rr(4, 0) = 3; rr(4, 1) = 1; rr(4, 2) = 2;
    rr(5, 0) = 3; rr(5, 1) = 2; rr(5, 2) = 1;
    arma::imat rr0 = rr - 1;
    // similar
    arma::imat rx(6, 3);
    // ugh
    rx(0, 0) = 1; rx(0, 1) = 2; rx(0, 2) = 3;
    rx(1, 0) = 1; rx(1, 1) = 3; rx(1, 2) = 2;
    rx(2, 0) = 2; rx(2, 1) = 1; rx(2, 2) = 3;
    rx(3, 0) = 3; rx(3, 1) = 1; rx(3, 2) = 2;
    rx(4, 0) = 2; rx(4, 1) = 3; rx(4, 2) = 1;
    rx(5, 0) = 3; rx(5, 1) = 2; rx(5, 2) = 1;
    int ir_applied = 1;
    Rcpp::IntegerVector one_through_6 = Rcpp::IntegerVector::create(1,2,3,4,5,6);
    Rcpp::NumericVector rr_probs(6);
    double A, B, C;
    int ir, ir_chosen;
            
    if (verbose) {
        std::cout << "start" << std::endl;
    }
    
    for(int iGrid = 0; iGrid < nGrids; iGrid++) {
        if (verbose) {
            std::cout << "iGrid = " << iGrid << ", nGrids = " << nGrids << std::endl;
            std::cout << "normal forward one" << std::endl;
        }
        //
        original_c1_this_grid = c1(iGrid);
        original_c2_this_grid = c2(iGrid);
        if(ff>0) original_c3_this_grid = c3(iGrid);
        //
        // normal forward one (includes initialization)
        //
        if (iGrid == 0) {
            alphaHat_t1.col(iGrid) = priorCurrent_m.col(s) % eMatGrid_t1.col(iGrid);
            c1(iGrid) = 1 / sum(alphaHat_t1.col(iGrid));
            alphaHat_t1.col(iGrid) *= c1(iGrid);
            //
            alphaHat_t2.col(iGrid) = priorCurrent_m.col(s) % eMatGrid_t2.col(iGrid);
            c2(iGrid) = 1 / sum(alphaHat_t2.col(iGrid));
            alphaHat_t2.col(iGrid) *= c2(iGrid);
            if (ff > 0) {
                alphaHat_t3.col(iGrid) = priorCurrent_m.col(s) % eMatGrid_t3.col(iGrid);
                c3(iGrid) = 1 / sum(alphaHat_t3.col(iGrid));
                alphaHat_t3.col(iGrid) *= c3(iGrid);
            }
        } else {
            if (ff == 0 && in_flip_mode) {
                emat_temp_col = eMatGrid_t1.col(iGrid);
                eMatGrid_t1.col(iGrid) = eMatGrid_t2.col(iGrid);
                eMatGrid_t2.col(iGrid) = emat_temp_col;
            }
            if (ff > 0 && ir_applied != 1) {
                e1 = eMatGrid_t1.col(iGrid);
                e2 = eMatGrid_t2.col(iGrid);
                e3 = eMatGrid_t3.col(iGrid);
                if (ir_applied == 2) {
                    eMatGrid_t2.col(iGrid) = e3;
                    eMatGrid_t3.col(iGrid) = e2;                  
                } else if (ir_applied == 3) {
                    eMatGrid_t1.col(iGrid) = e2;
                    eMatGrid_t2.col(iGrid) = e1;
                } else if (ir_applied == 4) {
                    eMatGrid_t1.col(iGrid) = e2;
                    eMatGrid_t2.col(iGrid) = e3;
                    eMatGrid_t3.col(iGrid) = e1;                  
                } else if (ir_applied == 5) {
                    eMatGrid_t1.col(iGrid) = e3;
                    eMatGrid_t2.col(iGrid) = e1;
                    eMatGrid_t3.col(iGrid) = e2;
                } else if (ir_applied == 6) {
                    eMatGrid_t1.col(iGrid) = e3;
                    eMatGrid_t2.col(iGrid) = e2;
                    eMatGrid_t3.col(iGrid) = e1;                  
                }
            }
            rcpp_alpha_forward_one(s, iGrid, K, alphaHat_t1, transMatRate_tc_H, eMatGrid_t1, alphaMatCurrent_tc, c1, minus_log_c_sum, true);
            rcpp_alpha_forward_one(s, iGrid, K, alphaHat_t2, transMatRate_tc_H, eMatGrid_t2, alphaMatCurrent_tc, c2, minus_log_c_sum, true);
            if (ff > 0) {
                rcpp_alpha_forward_one(s, iGrid, K, alphaHat_t3, transMatRate_tc_H, eMatGrid_t3, alphaMatCurrent_tc, c3, minus_log_c_sum, true);
            }
        }
        // note, forward_one does something different, this is what we want here, we are over-writing c
        minus_log_c1_sum -= log(c1(iGrid));
        minus_log_c2_sum -= log(c2(iGrid));
        if(ff>0) minus_log_c3_sum -= log(c3(iGrid));
        //
        //  go over read labels too, maybe flip them
        //
        if (verbose) {
            std::cout << "check reads, in_flip_mode = " << in_flip_mode << ", ir_applied = " << ir_applied << std::endl;
        }
        done_reads = false;
        while(!done_reads) {
            if (iRead > (nReads - 1)) {
                done_reads = true;
            } else {
                if (wif0(iRead) == iGrid) {
                    if (ff == 0 && in_flip_mode) {
                        H(iRead) = 3 - H(iRead);
                    } else if (ff > 0 && ir_applied != 1) {
                        H(iRead) = rx(ir_applied - 1, H(iRead) - 1);
                    }
                    iRead++;
                }
                if (iRead > (nReads - 1)) {
                    done_reads = true;
                } else {
                    if (wif0(iRead) > iGrid) {
                        done_reads = true;
                    }
                }
            }
        }
        //
        // now do this bit
        //
        if (verbose) {
            std::cout << "decide if to go in" << std::endl;
        }
        bool check = false;
        if (!shard_check_every_pair) {
            // the standard one
            iGridConsider = consider_grid_where_0_based(iGrid);
            if ((-1 < iGridConsider) & (iGridConsider < (n_blocks - 1))) {
                check = true;
            }
        } else {
            //do all but the last one
            if (iGrid < (nGrids - 1)) {
                check = true;
            }
            iGridConsider = iGrid;
        }
        // do not do last one
        if (check) {
            if (!shard_check_every_pair) {            
                split_grid = consider_grid_end_0_based(iGridConsider);
            } 
            if (verbose) {
                std::cout << "Considering split_grid = " << split_grid << std::endl;
                std::cout << "iGridConsider = " << iGridConsider << std::endl;                
            }
            //
            // on the fly version
            //
            if (ff == 0) {
                pA1 = minus_log_c1_sum + minus_log_original_c1_sum + log(sum(alphaHat_t1.col(iGrid) % betaHat_t1.col(iGrid)));
                pA2 = minus_log_c2_sum + minus_log_original_c2_sum + log(sum(alphaHat_t2.col(iGrid) % betaHat_t2.col(iGrid)));
                pB1 = minus_log_c2_sum + minus_log_original_c1_sum + log(sum(alphaHat_t2.col(iGrid) % betaHat_t1.col(iGrid)));
                pB2 = minus_log_c1_sum + minus_log_original_c2_sum + log(sum(alphaHat_t1.col(iGrid) % betaHat_t2.col(iGrid)));
                //
                calculated_difference = pB1 + pB2 - pA1 - pA2;
                probs1 = 1;
                probs2 = exp(calculated_difference);
                probs_sum = probs1 + probs2;
                probs1 /= probs_sum;
                probs2 /= probs_sum;
                in_flip_mode = runif_block(iGridConsider) > probs1;
                // record stuff now yyyyyyyyyeeeeeeeeeeeeeeeeessssssssssssss
                shard_block_results(iGridConsider, 0) = iGridConsider;
                shard_block_results(iGridConsider, 1) = probs1;
                shard_block_results(iGridConsider, 2) = probs2;
                shard_block_results(iGridConsider, 3) = pA1;
                shard_block_results(iGridConsider, 4) = pA2;
                shard_block_results(iGridConsider, 5) = pB1;
                shard_block_results(iGridConsider, 6) = pB2;
                if (in_flip_mode) {
                    shard_block_results(iGridConsider, 7) = 1;
                } else {
                    shard_block_results(iGridConsider, 7) = 0;
                }
                shard_block_results(iGridConsider, 8) = pA1 + pA2;
                shard_block_results(iGridConsider, 9) = pB1 + pB2;
            } else {
                for(ir = 0; ir < 6; ir++) {
                    if (ir == 0) {
                        A = log(sum(alphaHat_t1.col(iGrid) % betaHat_t1.col(iGrid))) + minus_log_c1_sum;
                        B = log(sum(alphaHat_t2.col(iGrid) % betaHat_t2.col(iGrid))) + minus_log_c2_sum;
                        C = log(sum(alphaHat_t3.col(iGrid) % betaHat_t3.col(iGrid))) + minus_log_c3_sum;
                    } else if (ir == 1) {
                        A = log(sum(alphaHat_t1.col(iGrid) % betaHat_t1.col(iGrid))) + minus_log_c1_sum;
                        B = log(sum(alphaHat_t2.col(iGrid) % betaHat_t3.col(iGrid))) + minus_log_c3_sum;
                        C = log(sum(alphaHat_t3.col(iGrid) % betaHat_t2.col(iGrid))) + minus_log_c2_sum;
                    } else if (ir == 2) {
                        A = log(sum(alphaHat_t1.col(iGrid) % betaHat_t2.col(iGrid))) + minus_log_c2_sum;
                        B = log(sum(alphaHat_t2.col(iGrid) % betaHat_t1.col(iGrid))) + minus_log_c1_sum;
                        C = log(sum(alphaHat_t3.col(iGrid) % betaHat_t3.col(iGrid))) + minus_log_c3_sum;
                    } else if (ir == 3) {
                        A = log(sum(alphaHat_t1.col(iGrid) % betaHat_t2.col(iGrid))) + minus_log_c2_sum;
                        B = log(sum(alphaHat_t2.col(iGrid) % betaHat_t3.col(iGrid))) + minus_log_c3_sum;
                        C = log(sum(alphaHat_t3.col(iGrid) % betaHat_t1.col(iGrid))) + minus_log_c1_sum;
                    } else if (ir == 4) {
                        A = log(sum(alphaHat_t1.col(iGrid) % betaHat_t3.col(iGrid))) + minus_log_c3_sum;
                        B = log(sum(alphaHat_t2.col(iGrid) % betaHat_t1.col(iGrid))) + minus_log_c1_sum;
                        C = log(sum(alphaHat_t3.col(iGrid) % betaHat_t2.col(iGrid))) + minus_log_c2_sum;
                    } else if (ir == 5) {
                        A = log(sum(alphaHat_t1.col(iGrid) % betaHat_t3.col(iGrid))) + minus_log_c3_sum;
                        B = log(sum(alphaHat_t2.col(iGrid) % betaHat_t2.col(iGrid))) + minus_log_c2_sum;
                        C = log(sum(alphaHat_t3.col(iGrid) % betaHat_t1.col(iGrid))) + minus_log_c1_sum;
                    }
                    double pA1 = A + minus_log_original_c1_sum;
                    double pA2 = B + minus_log_original_c2_sum;
                    double pA3 = C + minus_log_original_c3_sum;
                    double p = pA1 + pA2 + pA3;
                    rr_probs(ir) = p;
                }
                A = rr_probs(0);
                for(ir = 0; ir < 6; ir++) {
                    rr_probs(ir) -= A;
                    rr_probs(ir) = exp(rr_probs(ir));
                }
                A = sum(rr_probs);
                for(ir = 0; ir < 6; ir++) {
                    rr_probs(ir) /= A;
                }
                // should be the same now
                Rcpp::IntegerVector ir_temp = Rcpp::sample(one_through_6, 1, false, rr_probs);
                ir_chosen = ir_temp(0);
                ir_applied = ir_chosen; // just take the same one
                shard_block_results(iGridConsider, 0) = iGridConsider;
                shard_block_results(iGridConsider, 1) = ir_chosen;
                shard_block_results(iGridConsider, 2) = rr_probs(0);
                shard_block_results(iGridConsider, 3) = rr_probs(1);
                shard_block_results(iGridConsider, 4) = rr_probs(2);
                shard_block_results(iGridConsider, 5) = rr_probs(3);
                shard_block_results(iGridConsider, 6) = rr_probs(4);
                shard_block_results(iGridConsider, 7) = rr_probs(5);            
                //
            }
            //
            // at the end 
            //
            if (verbose) {
                if (ff == 0 && in_flip_mode) {
                    std::cout << "FLIP ME UP BRO" << std::endl;
                } else if (ff > 0 && (ir_applied > 1)) {
                    std::cout << "FLIP ME UP BRO, ff more than 0" << std::endl;
                } else {
                    std::cout << "no no no flip thanks" << std::endl;
                }
            }
        }
        minus_log_original_c1_sum += log(original_c1_this_grid);
        minus_log_original_c2_sum += log(original_c2_this_grid);
        if (ff > 0) {
            minus_log_original_c3_sum += log(original_c3_this_grid);
        }
    }
    //
    // re-run backward
    //
    if (verbose) {
        std::cout << "finalize by re-running backward" << std::endl;
    }
    betaHat_t1.col(nGrids - 1).fill(c1(nGrids-1));
    betaHat_t2.col(nGrids - 1).fill(c2(nGrids-1));
    if(ff>0) betaHat_t3.col(nGrids - 1).fill(c3(nGrids-1));
    //
    Rcpp_run_backward_haploid(betaHat_t1, c1, eMatGrid_t1, alphaMatCurrent_tc, transMatRate_tc_H, s);
    Rcpp_run_backward_haploid(betaHat_t2, c2, eMatGrid_t2, alphaMatCurrent_tc, transMatRate_tc_H, s);
    if(ff>0) Rcpp_run_backward_haploid(betaHat_t3, c3, eMatGrid_t3, alphaMatCurrent_tc, transMatRate_tc_H, s);
    if (verbose) {
        std::cout << "done rcpp block gibbs sampling" << std::endl;
    }
    to_return.push_back(shard_block_results, "shard_block_results");
    return(to_return);
}
