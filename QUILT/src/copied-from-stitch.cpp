// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

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



//
//
// not yet sure of the proper R-ish way to do this
// for now copy these in from R. is OK because of both GPL license
//
//

double print_times(
    double prev,
    int suppressOutput,
    std::string past_text,
    std::string next_text
) {
    if( suppressOutput == 0 ) {    
        double cur=clock();
        std::cout << std::setw (40) << past_text;
        printf ("- %.6f cpu sec -", ((double)cur - (double)prev)* 1.0e-6);
        std::cout << next_text << std::endl;
        prev=cur;
    }
    return prev;
}




Rcpp::IntegerVector rcpp_int_expand(arma::ivec& hapc, const int nSNPs) {
  const int nbSNPs = hapc.size();
  //const int nSNPs = nbSNPs * 32;
  Rcpp::IntegerVector hap(nSNPs);
  int j = 0;
  int imax;
  for(int bs = 0; bs < nbSNPs; bs++) {
    if (bs < (nbSNPs - 1)) {
      imax = 32;
    } else {
      // final one!
      imax = nSNPs - 32 * bs;
    }
    std::uint32_t tmp(hapc(bs));
    for (int i = 0; i < imax; i++, tmp >>= 1) {
      hap(j++) = tmp & 0x1;
    }
  }
  return(hap);
}



arma::imat inflate_fhb(
    arma::imat& rhb,
    Rcpp::IntegerVector& haps_to_get,
    const int nSNPs
) {
    const int K = rhb.n_cols;
    const int nbSNPs = rhb.n_rows;
    // i think this function might work without kmax
    // but probably safer / simpler to keep it in
    int imax;
    int n_haps_to_get = haps_to_get.size();    
    arma::imat rhi_subset(nSNPs, n_haps_to_get);
    // outer loop is on k, as the columns are "K", i.e. haps
    for(int ik = 0; ik < n_haps_to_get; ik++) {
        int k = haps_to_get(ik);
        for(int bs = 0; bs < nbSNPs; bs++) {
            int d32_times_bs = 32 * bs;
            if (bs < (nbSNPs - 1)) {
                imax = 32;
            } else {
                // final one!
                imax = nSNPs - d32_times_bs;
            }
	    std::uint32_t tmp(rhb(bs, k));
	    // this weird looking code taken largely from R c code
	    // e.g. search for this in R github
	    // SEXP attribute_hidden do_intToBits(SEXP call, SEXP op, SEXP args, SEXP env)
	    for (int i = 0; i < imax; i++, tmp >>= 1) {
	      // might be inefficient
	      // might need to work with temporary matrix?
	      // revisit later if actually slow!
                rhi_subset(d32_times_bs + i, ik) = tmp & 0x1;
	    }
	}
    }
    return(rhi_subset);
}

