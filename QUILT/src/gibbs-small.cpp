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
