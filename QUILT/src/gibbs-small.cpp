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


Rcpp::IntegerVector rcpp_int_expand(arma::ivec& hapc, const int nSNPs);



//' @export
// [[Rcpp::export]]
int rcpp_simple_binary_search(
    int val,
    Rcpp::IntegerVector vec
) {
    //
    // note that this returns the 0-based value
    //
    int nori = vec.length();
    if (nori == 1) {
        return(0);
    }
    int n = nori;
    int i = n / 2;
    n = n / 4;
    while(true) {
        if (vec(i) == val) {
            return(i);
        } else if (vec(i) < val) {
            i += n;
        } else {
            i -= n;
        }
        n = n / 2;
        if (n < 1) {
            n = 1;
        }
        if (i < 0) {
            i = 0;
        }
        if (i > (nori - 1)) {
            i = nori - 1;
        }
    }
}







//' @export
// [[Rcpp::export]]
int rcpp_simple_binary_matrix_search(
    int val,
    Rcpp::IntegerMatrix mat,
    int s1,
    int e1
) {
    int nori = e1 - s1 + 1; // width (1-based)
    if (nori == 1) {
        return(0);
    }
    int n = nori;
    int i = n / 2; // 0-based here
    n = n / 4;
    int c = 0;
    while(c < 100) {
        c++;
        if (mat(s1 - 1 + i, 0) == val) {
            return(mat(s1 - 1 + i, 1));
        } else if (mat(s1 - 1 + i, 0) < val) {
            i += n;
        } else {
            i -= n;
        }
        n = n / 2;
        if (n < 1) {
            n = 1;
        }
        if (i < 0) {
            i = 0;
        }
        if (i > (nori - 1)) {
            i = nori - 1;
        }
    }
    std::cout << "Something has gone wrong with binary matrix search, val = " << val << ", s1 = " << s1 << ", e1 = " << e1 << std::endl;
    return(mat(s1, 1));
}








//' @export
// [[Rcpp::export]]
void Rcpp_make_eMatRead_t_for_gibbs_using_objects(
    arma::mat& eMatRead_t,
    const Rcpp::List& sampleReads,
    const arma::imat& hapMatcher,
    const Rcpp::RawMatrix hapMatcherR,
    const bool use_hapMatcherR,
    const Rcpp::IntegerVector& grid,
    const arma::imat& rhb_t,
    const arma::mat& distinctHapsIE,
    const Rcpp::IntegerMatrix& eMatDH_special_matrix_helper,
    const Rcpp::IntegerMatrix& eMatDH_special_matrix,
    const double ref_error,
    const Rcpp::IntegerVector& which_haps_to_use,
    const bool rescale_eMatRead_t,
    const int Jmax,
    const double maxDifferenceBetweenReads,
    const bool use_eMatDH_special_symbols
) {
    const int nReads = sampleReads.size();
    const int K = which_haps_to_use.length();
    int iRead, iGrid0, iGrid0_prev, k, j;
    double eps, d1;
    double pR = 1;
    double pA = 1;
    double e;
    Rcpp::RawVector haps_at_gridR(K);
    Rcpp::IntegerVector haps_at_grid(K);
    Rcpp::IntegerVector hap(32);
    double d2 = 1 / maxDifferenceBetweenReads;
    int bvtd; // binary value to decompose
    int s1 = 0;
    int e1 = 0;
    for(iRead = 0; iRead < nReads; iRead++) {
        Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
        int J = as<int>(readData[0]); // number of Unique SNPs on read
        arma::ivec bq = as<arma::ivec>(readData[2]); // bq for each SNP
        arma::ivec u = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
	iGrid0 = grid(u(0));
        if (use_hapMatcherR) {
            for(k = 0; k < K; k++) {
                haps_at_gridR(k) = hapMatcherR(which_haps_to_use(k) - 1, iGrid0);
            }
        } else {
            for(k = 0; k < K; k++) {
                haps_at_grid(k) = hapMatcher(which_haps_to_use(k) - 1, iGrid0);
            }
        }
        if (use_eMatDH_special_symbols) {
            s1 = eMatDH_special_matrix_helper(iGrid0, 0);
            e1 = eMatDH_special_matrix_helper(iGrid0, 1);
        }
	iGrid0_prev = iGrid0;	
        if(J >= Jmax) {
            J = Jmax;
        }
        for(j = 0; j <= J; j++) {
            if(bq(j) < 0) {
                eps = pow(10,(double(bq(j)) / 10));
                pR = 1 - eps;
                pA = eps / 3;
            }
            if(bq(j) > 0) {
                eps = pow(10, (-double(bq(j)) / 10));
                pR = eps / 3;
                pA = 1 - eps;
            }
	    iGrid0 = grid(u(j));
            if (iGrid0 != iGrid0_prev) {
                if (use_hapMatcherR) {
                    for(k = 0; k < K; k++) {
                        haps_at_gridR(k) = hapMatcherR(which_haps_to_use(k) - 1, iGrid0);
                    }
                } else {
                    for(k = 0; k < K; k++) {
                        haps_at_grid(k) = hapMatcher(which_haps_to_use(k) - 1, iGrid0);
                    }
                }
                if (use_eMatDH_special_symbols) {
                    s1 = eMatDH_special_matrix_helper(iGrid0, 0);
                    e1 = eMatDH_special_matrix_helper(iGrid0, 1);
                }
	    }
	    iGrid0_prev = iGrid0;
	    for(k = 0; k < K; k++) {

                // if common, do something
                // if rare, do something else, IF necessary
                
                if (use_hapMatcherR & (haps_at_gridR(k) > 0)) {
                    e = distinctHapsIE(haps_at_gridR(k) - 1, u(j));
                } else if (haps_at_grid(k) > 0) {
                    e = distinctHapsIE(haps_at_grid(k) - 1, u(j));
		} else {
                    if (use_eMatDH_special_symbols) {
                        bvtd = rcpp_simple_binary_matrix_search(
                            which_haps_to_use(k) - 1,
                            eMatDH_special_matrix,
                            s1,
                            e1
                        );
                    } else {
                        bvtd = rhb_t(which_haps_to_use(k) - 1, iGrid0);
                    }
                    std::uint32_t tmp(bvtd);
                    int j2 = 0;
                    for (int i = 0; i < 32; i++, tmp >>= 1) {
                        hap(j2++) = tmp & 0x1;
                    }
                    if (hap(u(j) % 32) == 1) {
                        e = 1 - ref_error;
                    } else {
                        e = ref_error;
                    }
		}
                
                eMatRead_t(k, iRead) *= (e * pA + (1 - e) * pR);
	    }
	    //
        }
        if (rescale_eMatRead_t) {
	  //
	  // this is copied from rcpp_make_eMatRead_t, informed by real data!
	  // 
	  double x=0;
	  for(k=0; k < K; k++) {
	    if(eMatRead_t(k,iRead)>x) {
	      x = eMatRead_t(k,iRead);
	    }
	  }
	  // of course I hit an error where x was defined but d1 was not for a read with x as 4.51168e-311
	  d1 = 1 / x;
	  if ((x == R_NaN) | (x == R_PosInf) | (x == R_NegInf) | (x == 0) | (d1 == R_PosInf) | (d1 == R_NegInf) | (d1 == R_NaN)) {
	    // e.g. with extremely long molecule, just ignore, hopefully not frequent
	    // note - could probably get around by doing this on the fly periodically
	    for(k=0; k < K; k++) {                
	        eMatRead_t(k, iRead) = 1;
	    }
	  } else {
	      // x is the maximum now
	      for(k=0; k < K; k++) {
	          eMatRead_t(k,iRead) *= d1;
	          if(eMatRead_t(k,iRead) < d2) {
		      eMatRead_t(k,iRead) = d2;
		  }
	      }
	  }
	}
    }
    return;
}







//' @export
// [[Rcpp::export]]
void Rcpp_make_eMatRead_t_for_final_rare_common_gibbs_using_objects(
    arma::mat& eMatRead_t,
    const Rcpp::List& rare_per_hap_info,
    const Rcpp::IntegerVector& common_snp_index,
    const Rcpp::LogicalVector& snp_is_common,
    const Rcpp::List& sampleReads,
    const Rcpp::RawMatrix hapMatcherR,
    const Rcpp::IntegerVector& grid,
    const arma::mat& distinctHapsIE,
    const Rcpp::IntegerMatrix& eMatDH_special_matrix_helper,
    const Rcpp::IntegerMatrix& eMatDH_special_matrix,
    const double ref_error,
    const Rcpp::IntegerVector& which_haps_to_use,
    const bool rescale_eMatRead_t,
    const int Jmax,
    const double maxDifferenceBetweenReads,
    const Rcpp::List& rare_per_snp_info
) {
    const int nReads = sampleReads.size();
    const int K = which_haps_to_use.length();
    int iRead, iGrid0, iGrid0_prev, k, j;
    double eps, d1;
    double pR = 1;
    double pA = 1;
    double e;
    Rcpp::RawVector haps_at_gridR(K);
    Rcpp::IntegerVector haps_at_grid(K);
    Rcpp::IntegerVector hap(32);
    double d2 = 1 / maxDifferenceBetweenReads;
    int bvtd; // binary value to decompose
    int s1 = 0;
    int e1 = 0;
    int u_j_common;
    int ik;
    arma::colvec ori_eMatRead_t_col(K);
    bool check_val;
    Rcpp::IntegerVector vals;
    const double one_minus_ref_error = 1 - ref_error;
    Rcpp::IntegerVector k_with_alt;
    double xe1 = 1;
    double xe2 = 1;
    //
    for(iRead = 0; iRead < nReads; iRead++) {
        // for(k = 0; k < K; k++) {        
        //     ori_eMatRead_t_col(k) = eMatRead_t(k, iRead);
        // }
        // // reset this to 1
        // eMatRead_t.col(iRead).fill(1);
        Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
        int J = as<int>(readData[0]); // number of Unique SNPs on read
        arma::ivec bq = as<arma::ivec>(readData[2]); // bq for each SNP
        arma::ivec u = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
        iGrid0_prev = -1;
        if(J >= Jmax) {
            J = Jmax;
        }
        for(j = 0; j <= J; j++) {
            if(bq(j) < 0) {
                eps = pow(10,(double(bq(j)) / 10));
                pR = 1 - eps;
                pA = eps / 3;
            }
            if(bq(j) > 0) {
                eps = pow(10, (-double(bq(j)) / 10));
                pR = eps / 3;
                pA = 1 - eps;
            }
            //
            // remember rare thing now MUCH more common
            //
            if (snp_is_common(u(j))) {
                //
                // potentially update this
                //
                //iGrid0 = grid(u(j));
                u_j_common = common_snp_index(u(j)) - 1; // common_snp_index is 1-based
                iGrid0 = floor(u_j_common / 32); // grid among common SNPs
                if (iGrid0 != iGrid0_prev) {
                    for(k = 0; k < K; k++) {
                        haps_at_gridR(k) = hapMatcherR(which_haps_to_use(k) - 1, iGrid0);
                    }
                    s1 = eMatDH_special_matrix_helper(iGrid0, 0);
                    e1 = eMatDH_special_matrix_helper(iGrid0, 1);
                }
                iGrid0_prev = iGrid0;
                // this is the position of u(j) (which is common) among all SNPs
                for(k = 0; k < K; k++) {
                    if (haps_at_gridR(k) > 0) {
                        e = distinctHapsIE(haps_at_gridR(k) - 1, u_j_common);
                    } else {
                        // this implementation might be very expensive, but oh well
                        bvtd = rcpp_simple_binary_matrix_search(
                            which_haps_to_use(k) - 1,
                            eMatDH_special_matrix,
                            s1,
                            e1
                        );
                        std::uint32_t tmp(bvtd);
                        int j2 = 0;
                        for (int i = 0; i < 32; i++, tmp >>= 1) {
                            hap(j2++) = tmp & 0x1;
                        }
                        // this could be a bug, heeds to be grid of specific location
                        if (hap(u_j_common % 32) == 1) {
                            e = 1 - ref_error;
                        } else {
                            e = ref_error;
                        }
                    }
                    eMatRead_t(k, iRead) *= (e * pA + (1 - e) * pR);
                }
            } else {
                // snp is rare so go for it
                k_with_alt = rare_per_snp_info[u(j)];
                // if all a 1, not much to do
                if (k_with_alt.length() == 1) {
                    if (!rescale_eMatRead_t) {
                        xe1 = ref_error * pA + one_minus_ref_error * pR;
                        eMatRead_t.col(iRead) *= xe1;
                    }
                } else {
                    // add in assuming all ref
                    xe1 =           ref_error * pA + one_minus_ref_error * pR;
                    xe2 = one_minus_ref_error * pA +           ref_error * pR;
                    eMatRead_t.col(iRead) *= xe1;                
                    // un-do and re-do because these are alt so should only have alts in them
                    for(ik = 1; ik < k_with_alt.length(); ik++) {
                        k = k_with_alt(ik) - 1; // these are 1-based                    
                        eMatRead_t(k, iRead) *= xe2 / xe1;
                    }
                }
                // int u_j_plus_1 = u(j) + 1;
                // for(k = 0; k < K; k++) {
                //     if (ori_eMatRead_t_col(k) == 0) {
                //         if (iRead == 12 && k < 6) {
                //             std::cout << "k = " << k << " and j = " << j << " is worth study" << std::endl;
                //         }
                //         // need to check this one
                //         vals = rare_per_hap_info[which_haps_to_use(k) - 1];
                //         check_val = false;
                //         for(int v = 0; v < vals.length(); v++) {
                //             if (u_j_plus_1 == (vals(v))) {
                //                 check_val = true;
                //             }
                //         }
                //         // 
                //         if (check_val) {
                //             e = 1 - ref_error;
                //         } else {
                //             e = ref_error;
                //         }
                //     } else {
                //         // assume ref
                //         e = ref_error;
                //     }
                //     eMatRead_t(k, iRead) *= (e * pA + (1 - e) * pR);
                // }
            }
	    //
        }
        if (rescale_eMatRead_t) {
	  //
	  // this is copied from rcpp_make_eMatRead_t, informed by real data!
	  // 
	  double x=0;
	  for(k=0; k < K; k++) {
	    if(eMatRead_t(k,iRead)>x) {
	      x = eMatRead_t(k,iRead);
	    }
	  }
	  // of course I hit an error where x was defined but d1 was not for a read with x as 4.51168e-311
	  d1 = 1 / x;
	  if ((x == R_NaN) | (x == R_PosInf) | (x == R_NegInf) | (x == 0) | (d1 == R_PosInf) | (d1 == R_NegInf) | (d1 == R_NaN)) {
	    // e.g. with extremely long molecule, just ignore, hopefully not frequent
	    // note - could probably get around by doing this on the fly periodically
	    for(k=0; k < K; k++) {                
	        eMatRead_t(k, iRead) = 1;
	    }
	  } else {
	      // x is the maximum now
	      for(k=0; k < K; k++) {
	          eMatRead_t(k,iRead) *= d1;
	          if(eMatRead_t(k,iRead) < d2) {
		      eMatRead_t(k,iRead) = d2;
		  }
	      }
	  }
	}
    }
    return;
}




//' @export
// [[Rcpp::export]]
void rcpp_calculate_gibbs_small_genProbs_and_hapProbs_using_binary_objects(
    arma::mat& alphaHat_t1,
    arma::mat& alphaHat_t2,
    arma::mat& alphaHat_t3,     
    arma::mat& betaHat_t1,
    arma::mat& betaHat_t2, 
    arma::mat& betaHat_t3,
    arma::rowvec& c1,
    arma::rowvec& c2,
    arma::rowvec& c3,
    arma::mat& genProbsM_t,
    arma::mat& genProbsF_t,    
    arma::mat& hapProbs_t,
    const arma::mat& gammaMT_t,
    const arma::mat& gammaMU_t,
    const arma::mat& gammaP_t,
    const arma::imat& hapMatcher,
    const Rcpp::RawMatrix& hapMatcherR,
    bool use_hapMatcherR,
    const arma::imat& distinctHapsB,
    const arma::mat& distinctHapsIE,
    const Rcpp::IntegerMatrix& eMatDH_special_matrix_helper,
    const Rcpp::IntegerMatrix& eMatDH_special_matrix,
    const Rcpp::IntegerVector& which_haps_to_use,
    const double ref_error,
    const arma::imat& rhb_t,
    const bool use_eMatDH_special_symbols,
    const bool calculate_gamma_on_the_fly,
    const bool sample_is_diploid
) {
    // loop over grids
    int k, iGrid, s, e, b, dh, nSNPsLocal, bvtd, kk;
    int K;
    int nGrids;
    if (calculate_gamma_on_the_fly) {
        nGrids = alphaHat_t1.n_cols;
        K = alphaHat_t1.n_rows;
    } else {
        nGrids = gammaMT_t.n_cols;
        K = gammaMT_t.n_rows;        
    }
    const int nSNPs = genProbsM_t.n_cols;
    const int nMaxDH = distinctHapsIE.n_rows;    
    const double ref_one_minus_error = 1 - ref_error;
    arma::icolvec dh_col(K);
    arma::vec g0(32);
    g0.fill(0);
    arma::vec g1(32);
    g1.fill(0);
    arma::vec g2(32);
    g2.fill(0);
    arma::vec h0(32);
    h0.fill(0);
    arma::vec h1(32);
    h1.fill(0);
    arma::vec h2(32);
    h2.fill(0);
    double gk0, gk1, gk2;
    double x1, x2, x3;
    arma::vec mg0(nMaxDH + 1);
    arma::vec mg1(nMaxDH + 1);
    arma::vec mg2(nMaxDH + 1);
    //
    arma::colvec gammaMT_t_local(K);
    arma::colvec gammaMU_t_local(K);
    arma::colvec gammaP_t_local(K);
    if(sample_is_diploid) gammaP_t_local.fill(0);
    //
    for(iGrid = 0; iGrid < nGrids; iGrid++) {
        s = 32 * (iGrid); // 0-based here
        e = 32 * (iGrid + 1) - 1;
        if (e > (nSNPs - 1)) {
            e = nSNPs - 1;
        }
        nSNPsLocal = e - s + 1;
        g0.fill(0);
        g1.fill(0);
        g2.fill(0);
        h0.fill(0);
        h1.fill(0);
        h2.fill(0);
        if (calculate_gamma_on_the_fly) {
            // hope this is right / don't get this wrong with labelling
            x1 = 1 / c1(iGrid);        
            gammaMT_t_local = (alphaHat_t1.col(iGrid) % betaHat_t1.col(iGrid)) * x1;
            x2 = 1 / c2(iGrid);        
            gammaMU_t_local = (alphaHat_t2.col(iGrid) % betaHat_t2.col(iGrid)) * x2;
            if(!sample_is_diploid) {
                x3 = 1 / c3(iGrid);        
                gammaP_t_local = (alphaHat_t3.col(iGrid) % betaHat_t3.col(iGrid)) * x3;   
            }
            //
        } else {
            gammaMT_t_local = gammaMT_t.col(iGrid);
            gammaMU_t_local = gammaMU_t.col(iGrid);
            if(!sample_is_diploid) gammaP_t_local = gammaP_t.col(iGrid);
        }
        for(k = 0; k < K; k++) {
            gk0 = gammaMT_t_local(k);
            gk1 = gammaMU_t_local(k);
            gk2 = gammaP_t_local(k);
            //gk0 = gammaMT_t(k, iGrid);
            //gk1 = gammaMU_t(k, iGrid);
            //gk2 = gammaP_t(k, iGrid);
            //
            // get the binary value to decompose
            //
            if (use_hapMatcherR) {
                kk = hapMatcherR(which_haps_to_use(k) - 1, iGrid);
            } else {
                kk = hapMatcher(which_haps_to_use(k) - 1, iGrid);
            }
            if (kk > 0) {
                bvtd = distinctHapsB(kk - 1, iGrid);
            } else {
                if (use_eMatDH_special_symbols) {
                    bvtd = rcpp_simple_binary_matrix_search(
                        which_haps_to_use(k) - 1,
                        eMatDH_special_matrix,
                        eMatDH_special_matrix_helper(iGrid, 0),
                        eMatDH_special_matrix_helper(iGrid, 1)
                    );
                } else {
                    bvtd = rhb_t(which_haps_to_use(k) - 1, iGrid);
                }
            }
            std::uint32_t tmp(bvtd);
            for(b = 0; b < nSNPsLocal; b++, tmp >>= 1) {
                if ((tmp & 0x1) == 0) {
                //if (tmp & (1<<b)) {                    
                    // alternate
                    h0(b) += gk0;
                    h1(b) += gk1;
                    h2(b) += gk2;
                } else {
                    g0(b) += gk0;
                    g1(b) += gk1;
                    g2(b) += gk2;
                }
            }
        }
        for(b = 0; b < nSNPsLocal; b++) {
            g0(b) = g0(b) * ref_one_minus_error + h0(b) * ref_error;
            g1(b) = g1(b) * ref_one_minus_error + h1(b) * ref_error;
            g2(b) = g2(b) * ref_one_minus_error + h2(b) * ref_error;            
        }
        //
        // add back in properly
        //
        for(b = 0; b < nSNPsLocal; b++) {
            genProbsM_t(0, s + b) = (1 - g0(b)) * (1 - g1(b));
            genProbsM_t(1, s + b) = (g0(b) * (1 - g1(b)) + (1 - g0(b)) * g1(b));
            genProbsM_t(2, s + b) = g0(b) * g1(b);
            //
            genProbsF_t(0, s + b) = (1 - g0(b)) * (1 - g2(b));
            genProbsF_t(1, s + b) = (g0(b) * (1 - g2(b)) + (1 - g0(b)) * g2(b));
            genProbsF_t(2, s + b) = g0(b) * g2(b);
            //
            hapProbs_t(0, s + b) = g0(b);
            hapProbs_t(1, s + b) = g1(b);
            hapProbs_t(2, s + b) = g2(b);
        }
    }
    return;
    //     s = 32 * (iGrid); // 0-based here
    //     e = 32 * (iGrid + 1) - 1;
    //     if (e > (nSNPs - 1)) {
    //         e = nSNPs - 1;
    //     }
    //     nSNPsLocal = e - s + 1;
    //     //
    //     g0.fill(0);
    //     g1.fill(0);
    //     g2.fill(0);
    //     mg0.fill(0);
    //     mg1.fill(0);
    //     mg2.fill(0);                
    //     //
    //     // loop
    //     //
    //     for(k = 0; k < K; k++) {
    //         dh = hapMatcher(which_haps_to_use(k) - 1, iGrid);
    //         mg0(dh) += gammaMT_t(k, iGrid);
    //         mg1(dh) += gammaMU_t(k, iGrid);
    //         mg2(dh) += gammaP_t(k, iGrid);            
    //         if (dh == 0) {
    //             // need to build this one
    //             gk0 = gammaMT_t(k, iGrid);
    //             gk1 = gammaMU_t(k, iGrid);
    //             gk2 = gammaP_t(k, iGrid);                
    //             std::uint32_t tmp(rhb_t(which_haps_to_use(k) - 1, iGrid));
    //             for(b = 0; b < nSNPsLocal; b++) {
    //                 if (tmp & (1<<b)) {
    //                     // alternate
    //                     g0(b) += gk0 * (ref_one_minus_error);
    //                     g1(b) += gk1 * (ref_one_minus_error);
    //                     g2(b) += gk2 * (ref_one_minus_error);                        
    //                 } else {
    //                     // reference
    //                     g0(b) += gk0 * (ref_error);
    //                     g1(b) += gk1 * (ref_error);
    //                     g2(b) += gk2 * (ref_error);
    //                 }
    //             }
    //         }
    //     }
    //     //
    //     // now do the rest!
    //     //
    //     for(b = 0; b < nSNPsLocal; b++) {
    //         for(dh = 0; dh < nMaxDH; dh++) {
    //             g0(b) += distinctHapsIE(dh, s + b) * mg0(dh + 1);
    //             g1(b) += distinctHapsIE(dh, s + b) * mg1(dh + 1);
    //             g2(b) += distinctHapsIE(dh, s + b) * mg2(dh + 1);
    //         }
    //     }
    //     // add back in properly
    //     for(b = 0; b < nSNPsLocal; b++) {
    //         genProbsM_t(0, s + b) = (1 - g0(b)) * (1 - g1(b));
    //         genProbsM_t(1, s + b) = (g0(b) * (1 - g1(b)) + (1 - g0(b)) * g1(b));
    //         genProbsM_t(2, s + b) = g0(b) * g1(b);
    //         //
    //         genProbsF_t(0, s + b) = (1 - g0(b)) * (1 - g2(b));
    //         genProbsF_t(1, s + b) = (g0(b) * (1 - g2(b)) + (1 - g0(b)) * g2(b));
    //         genProbsF_t(2, s + b) = g0(b) * g2(b);
    //         //
    //         hapProbs_t(0, s + b) = g0(b);
    //         hapProbs_t(1, s + b) = g1(b);
    //         hapProbs_t(2, s + b) = g2(b);
    //     }
    // }
}





//' @export
// [[Rcpp::export]]
void rcpp_calculate_genProbs_and_hapProbs_final_rare_common(
    arma::mat& alphaHat_t1,
    arma::mat& alphaHat_t2,
    arma::mat& alphaHat_t3,     
    arma::mat& betaHat_t1,
    arma::mat& betaHat_t2, 
    arma::mat& betaHat_t3,
    arma::rowvec& c1,
    arma::rowvec& c2,
    arma::rowvec& c3,
    arma::mat& genProbsM_t,
    arma::mat& genProbsF_t,    
    arma::mat& hapProbs_t,
    const arma::mat& gammaMT_t,
    const arma::mat& gammaMU_t,
    const arma::mat& gammaP_t,
    const Rcpp::RawMatrix& hapMatcherR,
    const arma::imat& distinctHapsB,
    const arma::mat& distinctHapsIE,
    const Rcpp::IntegerMatrix& eMatDH_special_matrix_helper,
    const Rcpp::IntegerMatrix& eMatDH_special_matrix,
    const Rcpp::IntegerVector& which_haps_to_use,
    const double ref_error,
    const arma::imat& rhb_t,
    const bool use_eMatDH_special_symbols,
    const Rcpp::List& rare_per_hap_info,
    const Rcpp::IntegerVector& common_snp_index,
    const Rcpp::LogicalVector& snp_is_common,
    const Rcpp::List& rare_per_snp_info,
    const bool sample_is_diploid
) {
    // loop over grids
    int k, bvtd, kk;
    int K;
    int nGrids;
    nGrids = alphaHat_t1.n_cols;
    K = alphaHat_t1.n_rows;
    const int nSNPs = genProbsM_t.n_cols;
    const int nMaxDH = distinctHapsIE.n_rows;    
    arma::icolvec dh_col(K);
    double x1, x2, x3;
    arma::vec mg0(nMaxDH + 1);
    arma::vec mg1(nMaxDH + 1);
    arma::vec mg2(nMaxDH + 1);
    //
    arma::colvec gammaMT_t_local;
    arma::colvec gammaMU_t_local;
    arma::colvec gammaP_t_local;
    //
    // first do common SNPs
    //
    int iFullGrid_prev = -1;
    int iFullGrid;
    Rcpp::IntegerVector k_with_alt;
    double h1, h2;
    int common_snp, common_grid;
    double d;
    int ik;
    Rcpp::IntegerVector hap(32);
    double one_minus_2_times_ref_error = 1 - 2 * ref_error;
    //
    //
    //
    for(int iFullSNP = 0; iFullSNP < nSNPs; iFullSNP++) {
        // iFullSNP is the index of all the SNPs
        iFullGrid = floor(iFullSNP / 32); // for gamma
        //
        // get gammas
        //
        if (iFullGrid != iFullGrid_prev) {
            // not sure about this re labelling
            // note here assuming will always do on the fly
            x1 = 1 / c1(iFullGrid);        
            gammaMT_t_local = (alphaHat_t1.col(iFullGrid) % betaHat_t1.col(iFullGrid)) * x1;
            x2 = 1 / c2(iFullGrid);        
            gammaMU_t_local = (alphaHat_t2.col(iFullGrid) % betaHat_t2.col(iFullGrid)) * x2;
            if(!sample_is_diploid){
                x3 = 1 / c3(iFullGrid);        
                gammaP_t_local = (alphaHat_t3.col(iFullGrid) % betaHat_t3.col(iFullGrid)) * x3;  
            } 
            iFullGrid_prev = iFullGrid;
        }
        //
        // check, if common, can do distinctHapsIE
        //
        if (snp_is_common(iFullSNP)) {
            common_snp = common_snp_index(iFullSNP) - 1; // the 0-based index of the common SNP
            common_grid = floor(common_snp / 32);
            for(k = 0; k < K; k++) {
                kk = hapMatcherR(which_haps_to_use(k) - 1, common_grid);
                if (kk > 0) {
                    d = distinctHapsIE(kk - 1, common_snp);
                } else {
                    // need to re-build, argh
                    bvtd = rcpp_simple_binary_matrix_search(
                        which_haps_to_use(k) - 1,
                        eMatDH_special_matrix,
                        eMatDH_special_matrix_helper(common_grid, 0),
                        eMatDH_special_matrix_helper(common_grid, 1)
                    );
                    std::uint32_t tmp(bvtd);
                    int j2 = 0;
                    for (int i = 0; i < 32; i++, tmp >>= 1) {
                        hap(j2++) = tmp & 0x1;
                    }
                    // this could be a bug, heeds to be grid of specific location
                    if (hap(common_snp % 32) == 1) {
                        d = 1 - ref_error;
                    } else {
                        d = ref_error;
                    }
                }
                hapProbs_t(0, iFullSNP) += gammaMT_t_local(k) * d;
                hapProbs_t(1, iFullSNP) += gammaMU_t_local(k) * d;
                if(!sample_is_diploid) hapProbs_t(2, iFullSNP) += gammaP_t_local(k) * d;
            }
        } else {
            // use rare info, slightly re-jigged (small so OK?)
            k_with_alt = rare_per_snp_info[iFullSNP];
            if (k_with_alt.length() == 1) {
                // nothing here!
                hapProbs_t(0, iFullSNP) = ref_error;
                hapProbs_t(1, iFullSNP) = ref_error;
                if(!sample_is_diploid) hapProbs_t(2, iFullSNP) = ref_error;
            } else {
                // first add all in assuming ref
                for(k = 0; k < K; k++) {
                    hapProbs_t(0, iFullSNP) += gammaMT_t_local(k) * ref_error;
                    hapProbs_t(1, iFullSNP) += gammaMU_t_local(k) * ref_error;
                    if(!sample_is_diploid) hapProbs_t(2, iFullSNP) += gammaP_t_local(k) * ref_error;
                }
                // now undo and redo
                for(ik = 1; ik < k_with_alt.length(); ik++) {
                    k = k_with_alt(ik) - 1; // these are 1-based
                    hapProbs_t(0, iFullSNP) += gammaMT_t_local(k) * one_minus_2_times_ref_error;
                    hapProbs_t(1, iFullSNP) += gammaMU_t_local(k) * one_minus_2_times_ref_error;
                    if(!sample_is_diploid) hapProbs_t(2, iFullSNP) += gammaP_t_local(k) * one_minus_2_times_ref_error;
                }
            }
        }
        //"maternal"
        h1 = hapProbs_t(0, iFullSNP);
        h2 = hapProbs_t(1, iFullSNP);
        genProbsM_t(0, iFullSNP) = (1 - h1) * (1 - h2);
        genProbsM_t(1, iFullSNP) = h1 * (1 - h2) + h2 * (1 - h1);
        genProbsM_t(2, iFullSNP) = h1 * h2;
        // "fetal"
        if(!sample_is_diploid) {
            h1 = hapProbs_t(0, iFullSNP);
            h2 = hapProbs_t(2, iFullSNP);
            genProbsF_t(0, iFullSNP) = (1 - h1) * (1 - h2);
            genProbsF_t(1, iFullSNP) = h1 * (1 - h2) + h2 * (1 - h1);
            genProbsF_t(2, iFullSNP) = h1 * h2;
        }
    }
    return;
}
