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
using namespace Rcpp;


//
//
// not yet sure of the proper R-ish way to do this
// for now copy these in from STITCH. is OK because of both GPL license. but man oh man this will be hard to maintain
// maybe I can copy in STITCH like I do seqlib, probably the right move here
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
        printf ("- %.6f cpu sec *1000 -", ((double)cur - (double)prev)* 1.0e-3);
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




//' @export
// [[Rcpp::export]]
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
) {
    next_section="make eMatRead_t";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    // constants
    //
    const int K = eHapsCurrent_tc.n_rows; // traditional K for haplotypes
    const int nReads = sampleReads.size();
    //
    // new variables
    //
    double eps, x, d1;
    double pR = 1;
    double pA = 1; 
    double d2 = 1 / maxDifferenceBetweenReads;
    int j, k, J, jj, iRead;
    //arma::mat eMatRead_t = arma::ones(K,nReads);
    //
    // now build
    //
    for(iRead = 0; iRead < nReads; iRead++) {
        Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
        // recal that below is what is used to set each element of sampleRead
        // note - this is no longer quite accurate
        // sampleReads.push_back(Rcpp::List::create(nU,d,phiU,pRU));
        J = as<int>(readData[0]); // number of Unique SNPs on read
        //int readSNP = as<int>(readData[1]); // leading SNP from read
        arma::ivec bqU = as<arma::ivec>(readData[2]); // bq for each SNP
        arma::ivec pRU = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
        // once each SNP is done, have P(read | k), can multiply to get P(read|(k1,k2))
        if(J >= Jmax) {
            J = Jmax;
        }
        for(j = 0; j <= J; j++) {
            if(bqU(j) < 0) {
                eps = pow(10,(double(bqU(j)) / 10));
                pR = 1 - eps;
                pA = eps / 3;
            }
            if(bqU(j) > 0) {
                eps = pow(10, (-double(bqU(j)) / 10));
                pR = eps / 3;
                pA = 1 - eps;
            }
            jj=pRU(j);
            eMatRead_t.col(iRead) %= ( eHapsCurrent_tc.slice(s).col(jj) * pA + (1 - eHapsCurrent_tc.slice(s).col(jj)) * pR);
            //
            if (run_pseudo_haploid == true) {
                x = pRgivenH1(iRead) / (pRgivenH1(iRead) + pRgivenH2(iRead));
                //
                for(k = 0; k < K; k++) {
                    // inefficient?
                    eMatHapOri_t(k, iRead) = eMatRead_t(k, iRead);
                    eMatRead_t(k,iRead) = x * eMatRead_t(k,iRead) + (1-x) * pRgivenH2(iRead);
                }
            }
        }
        //
        // cap P(read|k) to be within maxDifferenceBetweenReads orders of magnitude
        // also, if doing this, reset to be 1 at maximum
        //
        if (rescale_eMatRead_t) {
            x=0;
            for(k=0; k < K; k++)
                if(eMatRead_t(k,iRead)>x)
                    x=eMatRead_t(k,iRead);
            if (run_pseudo_haploid) {
                // original behaviour, ignore x == 0, probably unlikely
                x = x / maxDifferenceBetweenReads;
                // x is the maximum now
                for(k=0; k < K; k++)
                    if(eMatRead_t(k,iRead)<x)
                        eMatRead_t(k,iRead) = x;
                //
            } else {
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
    }
    return;
}


//' @export
// [[Rcpp::export]]
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
) {
    //
    next_section="Make eMat";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //    
    int nReads = sampleReads.size(); //
    const int K = eMatRead_t.n_rows; // traditional K for haplotypes        
    // arma::mat eMatGrid_t = arma::ones(K, nGrids); // why is this called SNP? 
    int iRead, k, w, readSNP;
    bool proceed;
    //
    for(iRead = 0; iRead < nReads; iRead++) {
        proceed = false;
        if (use_all_reads) {
            proceed = true;
        } else {
            if (H(iRead) == hap) {
                proceed = true;
            }
        }
        if (proceed) {
            //w = wif(iRead) - run_fb_grid_offset;            
            // unclear how slow this is, but carry around less in RAM
            Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
            readSNP = as<int>(readData[1]); // leading SNP from read
            w = readSNP - run_fb_grid_offset;
            eMatGrid_t.col(w) %= eMatRead_t.col(iRead);
            //std::cout << "eMatRead_t.col(iRead) = " << eMatRead_t(0, iRead) << ", " << eMatRead_t(1, iRead) << ", " << eMatRead_t(2, iRead) << ", " << eMatRead_t(3, iRead);
            //                std::cout << std::endl;                
        }
    }
    // now - afterward - cap eMatHapSNP
    double x, rescale_val, d2;
    int t;
    if (bound) {
        for(t = 0; t < nGrids; t++) {
            // if this is less than exactly 1 (i.e. there are results here), proceed
            if (eMatGrid_t(0, t) < 1) {
                x = 0;
                for (k = 0; k < K; k++) {
                    if (eMatGrid_t(k, t) > x) {
                        x = eMatGrid_t(k, t);
                    }
                }
                // x is the maximum now
                rescale_val = 1 / x;        
                for (k = 0; k < K; k++) {
                    if (rescale) {                    
                        eMatGrid_t(k, t) *= rescale_val;
                    }
                    d2 = 1 / maxEmissionMatrixDifference;
                    if(eMatGrid_t(k, t) < (d2)) {
                        eMatGrid_t(k, t) = d2;
                    }
                }
            }
        }
    }
    return;
}


//' @export
// [[Rcpp::export]]
Rcpp::List rcpp_make_fb_snp_offsets(
    const arma::mat& alphaHat_t,
    const arma::mat& betaHat_t,
    const arma::mat& blocks_for_output
) {
    int s, e;
    arma::mat alphaHatBlocks_t = arma::zeros(alphaHat_t.n_rows, blocks_for_output.n_rows);
    arma::mat betaHatBlocks_t = arma::zeros(betaHat_t.n_rows, blocks_for_output.n_rows);
    const int blocks_for_output_n_rows = blocks_for_output.n_rows;
    for(int i_output=0; i_output < blocks_for_output_n_rows; i_output++) {
        s = blocks_for_output(i_output, 2); // these are 0-based. these are the grid entries
        e = blocks_for_output(i_output, 3);
        alphaHatBlocks_t.col(i_output) = alphaHat_t.col(s);
        betaHatBlocks_t.col(i_output) = betaHat_t.col(e);
    }
    return(wrap(Rcpp::List::create(
                                   Rcpp::Named("alphaHatBlocks_t") = alphaHatBlocks_t,
                                   Rcpp::Named("betaHatBlocks_t") = betaHatBlocks_t
                                   )));
}



//' @export
// [[Rcpp::export]]
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
) {
    const int K = alphaMatCurrent_tc.n_rows;
    const int nGrids = alphaMatCurrent_tc.n_cols + 1;    
    //
    // initialize
    //
    int k;
    if (run_fb_subset == false) {
        for(k = 0; k < K; k++) {
            alphaHat_t(k, 0) = priorCurrent_m(k, s) * eMatGrid_t(k, 0);
        }
    } else {
        for(k=0; k < K; k++) {
            alphaHat_t(k, 0) = alphaStart(k);
        }
    }
    c(0) = 1 / sum(alphaHat_t.col(0));
    alphaHat_t.col(0) = alphaHat_t.col(0) * c(0);
    if (initialize_only) {
        return;
    }
    //
    //
    for(int iGrid = 1; iGrid < nGrids; iGrid++) {
        // NOTE - previously used this code, which is mathematically right
        // BUT since scaling is being used here, arma::sum(alphaHat_t.col(t - 1) is equal to 1 by definition
        // so can use the below (uncommented) code to simplify
        // alphaConst = transMatRate_t_H(1, t-1) * arma::sum(alphaHat_t.col(t - 1));
        //
        alphaHat_t.col(iGrid) = eMatGrid_t.col(iGrid) % (		   \
            transMatRate_tc_H(0, iGrid - 1, s) * alphaHat_t.col(iGrid - 1) + \
            transMatRate_tc_H(1, iGrid - 1, s) * alphaMatCurrent_tc.slice(s).col(iGrid - 1) );
        c(iGrid) = 1 / arma::sum(alphaHat_t.col(iGrid));
        alphaHat_t.col(iGrid) *= c(iGrid);
    }
    return ;
}


//' @export
// [[Rcpp::export]]
void Rcpp_run_backward_haploid(
    arma::mat& betaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const int s
) {
    const int nGrids = eMatGrid_t.n_cols;
    double x;
    arma::colvec e_times_b;
    for(int iGrid = nGrids - 2; iGrid >= 0; --iGrid) {
        e_times_b = eMatGrid_t.col(iGrid + 1) % betaHat_t.col(iGrid + 1);
        x = transMatRate_tc_H(1, iGrid, s) * sum(alphaMatCurrent_tc.slice(s).col(iGrid) % e_times_b);
        betaHat_t.col(iGrid) = c(iGrid) * (x + transMatRate_tc_H(0, iGrid, s) * e_times_b);
    }
    return;
}


// modified from above for when alphaMatCurrent_tc is a constant and is 1/K
// also does something simpler when nothing in that grid

//' @export
// [[Rcpp::export]]
void Rcpp_run_backward_haploid_QUILT_faster(
    arma::mat& betaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& transMatRate_tc_H,
    const Rcpp::LogicalVector& grid_has_read,
    const int s
) {
    const int nGrids = eMatGrid_t.n_cols;
    const double one_over_K = 1 / double(eMatGrid_t.n_rows);
    double x;
    arma::colvec e_times_b;
    for(int iGrid = nGrids - 2; iGrid >= 0; --iGrid) {
        if (grid_has_read(iGrid + 1)) {
            e_times_b = eMatGrid_t.col(iGrid + 1) % betaHat_t.col(iGrid + 1);
            x = transMatRate_tc_H(1, iGrid, s) * sum(e_times_b) * one_over_K;
            betaHat_t.col(iGrid) = c(iGrid) * (x + transMatRate_tc_H(0, iGrid, s) * e_times_b);
        } else {
            x = transMatRate_tc_H(1, iGrid, s) * sum(betaHat_t.col(iGrid + 1)) * one_over_K;
            betaHat_t.col(iGrid) = c(iGrid) * (x + transMatRate_tc_H(0, iGrid, s) * betaHat_t.col(iGrid + 1));
        }
    }
    return;
}



//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_make_smoothed_rate(
    const Rcpp::NumericVector & sigma_rate,
    const Rcpp::IntegerVector & L_grid,
    const int shuffle_bin_radius,
    const bool verbose = false
) {
    const int nGrids = L_grid.length();
    int iGrid, iGrid_left, iGrid_right, bp_remaining, bp_to_add;
    int focal_point, bp_prev;
    double total_bp_added;
    //int min_L_grid = Rcpp::min(L_grid);
    //int max_L_grid = Rcpp::max(L_grid);
    Rcpp::NumericVector smoothed_rate(nGrids - 1);
    for(iGrid = 0; iGrid < (nGrids - 1); iGrid++) {
        if (verbose) {
            std::cout << "iGrid=" << iGrid << std::endl;
        }
        focal_point = (L_grid(iGrid) + L_grid(iGrid + 1)) / 2; // automatically rounded        
        if (verbose) {
            std::cout << "focal_point=" << focal_point << std::endl;
        }
        //
        // left
        //
        iGrid_left=iGrid;
        bp_remaining = shuffle_bin_radius;
        bp_prev = focal_point;
        total_bp_added = 0;
        while ((0 < bp_remaining) & (0 <= iGrid_left)) {
            bp_to_add = (bp_prev - L_grid(iGrid_left));
            if ((bp_remaining - bp_to_add) < 0) {
                bp_to_add = bp_remaining;
                bp_remaining = 0;
            } else {
                bp_remaining = bp_remaining - bp_to_add;
            }
            // add bit
            smoothed_rate(iGrid) = smoothed_rate(iGrid) + \
                bp_to_add * sigma_rate(iGrid_left);
            // move
            total_bp_added += bp_to_add;            
            bp_prev = L_grid(iGrid_left);
            iGrid_left = iGrid_left - 1;                
        }
        //
        // right
        //
        iGrid_right = iGrid + 1;
        bp_remaining = shuffle_bin_radius;
        bp_prev = focal_point;            
        while ((0 < bp_remaining) & (iGrid_right < nGrids)) {
            bp_to_add = (L_grid(iGrid_right) - bp_prev);
            if ((bp_remaining - bp_to_add) < 0) {
                bp_to_add = bp_remaining;
                bp_remaining = 0;
            } else {
                bp_remaining = bp_remaining - bp_to_add;
            }
            // add bit
            smoothed_rate(iGrid) = smoothed_rate(iGrid) +       \
                bp_to_add * sigma_rate(iGrid_right - 1);
            // move
            total_bp_added += bp_to_add;
            bp_prev = L_grid(iGrid_right);
            iGrid_right = iGrid_right + 1;                
        }
        //
        // normalize
        //
        smoothed_rate(iGrid) /= total_bp_added;
    }
    return(smoothed_rate);
}

//' @export
// [[Rcpp::export]]
int rcpp_determine_where_to_stop(
    const Rcpp::NumericVector& smoothed_rate,
    const Rcpp::LogicalVector& available,
    int& snp_best, // 0-based here
    double& thresh,
    int& nGrids,
    bool is_left
) {
    double mult;
    if (is_left) {
        mult = 1;
    } else {
        mult = -1;
    }
    int snp_consider = snp_best;
    double val_cur = smoothed_rate(snp_consider);
    double val_prev = smoothed_rate(snp_best);
    //
    int snp_min = snp_consider;
    double val_min = smoothed_rate(snp_min);
    int c = 1;
    bool are_done = false;
    while(!are_done) {
        snp_consider = snp_consider + (-1) * mult;
        val_cur = smoothed_rate(snp_consider);
        if (5 <= c) {
            val_prev = smoothed_rate(snp_consider + 5 * mult);
        }
        c += 1;
        if (val_cur < val_min) {
            snp_min = snp_consider;
            val_min = val_cur;
        }
        // do not continue if would go out of bound
        if ((snp_consider <= 2) | ((nGrids - 3) <= snp_consider)) {
            are_done = true;
        } else if (available(snp_consider + (-1) * mult) == false) {
            are_done = true;
        } else if ((3 * val_min) < val_cur) {
            are_done = true;
        } else if ((val_cur < thresh) & (val_prev < val_cur)) {
            are_done = true;
        }
    }
    return(snp_min);
}



//' @export
// [[Rcpp::export]]
Rcpp::NumericVector increment2N(int yT, int xT, Rcpp::NumericVector y, Rcpp::NumericVector z) {
  Rcpp::NumericVector x(xT+1);
  int t;
  for(t=0; t<=yT-1; t++)
    x[z[t]]=x[z[t]]+y[t];
  return(x);
}
