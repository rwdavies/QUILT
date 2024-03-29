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
void Rcpp_make_eMatRead_t_for_gibbs_using_objects(
    arma::mat& eMatRead_t,
    const Rcpp::List& sampleReads,
    const arma::imat& hapMatcher,
    const Rcpp::IntegerVector& grid,
    const arma::imat& rhb_t,
    const arma::mat& distinctHapsIE,
    const double ref_error,
    const Rcpp::IntegerVector& which_haps_to_use,
    const bool rescale_eMatRead_t,
    const int Jmax,
    const double maxDifferenceBetweenReads
) {
    const int nReads = sampleReads.size();
    const int K = which_haps_to_use.length();
    int iRead, iGrid0, iGrid0_prev, k, j;
    double eps, d1;
    double pR = 1;
    double pA = 1;
    double e;
    Rcpp::IntegerVector haps_at_grid(K);
    Rcpp::IntegerVector hap(32);
    double d2 = 1 / maxDifferenceBetweenReads;    
    for(iRead = 0; iRead < nReads; iRead++) {
        Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
        int J = as<int>(readData[0]); // number of Unique SNPs on read
        arma::ivec bq = as<arma::ivec>(readData[2]); // bq for each SNP
        arma::ivec u = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
	iGrid0 = grid(u(0));	
	for(k = 0; k < K; k++) {
	  haps_at_grid(k) = hapMatcher(which_haps_to_use(k) - 1, iGrid0);
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
	        for(k = 0; k < K; k++) {
		    haps_at_grid(k) = hapMatcher(which_haps_to_use(k) - 1, iGrid0);
	        }
	    }
	    iGrid0_prev = iGrid0;
	    //
	    for(k = 0; k < K; k++) {
	        if (haps_at_grid(k) > 0) {
		  e = distinctHapsIE(haps_at_grid(k) - 1, u(j));
		} else {
                    // pretty efficient, not perfect, but meh
                    std::uint32_t tmp(rhb_t(which_haps_to_use(k) - 1, iGrid0));
		  int j2 = 0;
		  for (int i = 0; i < 32; i++, tmp >>= 1) {
		      hap(j2++) = tmp & 0x1;
		  }
		  if (hap(u(j) - iGrid0 * 32) == 1) {
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
void rcpp_calculate_gibbs_small_genProbs_and_hapProbs_using_binary_objects(
    arma::mat& genProbsM_t,
    arma::mat& genProbsF_t,    
    arma::mat& hapProbs_t,
    const arma::mat& gammaMT_t,
    const arma::mat& gammaMU_t,
    const arma::mat& gammaP_t,
    const arma::imat& hapMatcher,
    const arma::mat& distinctHapsIE,
    const Rcpp::IntegerVector& which_haps_to_use,
    const double ref_error,
    const arma::imat& rhb_t
) {
    // loop over grids
    int k, iGrid, s, e, b, dh, nSNPsLocal;
    const int K = gammaMT_t.n_rows;
    const int nGrids = gammaMT_t.n_cols;
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
    arma::vec mg0(nMaxDH + 1);
    arma::vec mg1(nMaxDH + 1);
    arma::vec mg2(nMaxDH + 1);        
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
        for(k = 0; k < K; k++) {
            gk0 = gammaMT_t(k, iGrid);
            gk1 = gammaMU_t(k, iGrid);
            gk2 = gammaP_t(k, iGrid);                
            std::uint32_t tmp(rhb_t(which_haps_to_use(k) - 1, iGrid));
            //std::uint32_t tmp(rhb_t_subset(k, iGrid));
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
