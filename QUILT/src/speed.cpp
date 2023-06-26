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


//' @export
// [[Rcpp::export]]
void rcpp_test1(
    arma::mat& alphaHat_t1,
    arma::mat& alphaHat_t2,
    arma::mat& betaHat_t1,
    arma::mat& betaHat_t2,
    const Rcpp::List& sampleReads,
    arma::mat& eMatRead_t    
) {
    int K = alphaHat_t1.n_rows;
    int nGrids = alphaHat_t1.n_cols;
    int nReads = eMatRead_t.n_cols;
    arma::mat ab_m(K, 2);
    Rcpp::NumericVector pC(3), pA1(3), pA2(3);
    arma::mat alphaHat_m(K, 3); // tested, this is the better orientation
    arma::mat betaHat_m(K, 3);
    arma::colvec eMatRead_t_col;    
    //
    // loop over grids
    //
    int iRead = 0;
    bool done_reads = false;
    int h_rC = 0;
    int h_rA1 = 1;
    for(int iGrid = 0; iGrid < nGrids; iGrid++) {
      //
      // reset stuff for this grid
      //
      alphaHat_m.col(0) = alphaHat_t1.col(iGrid);
      alphaHat_m.col(1) = alphaHat_t2.col(iGrid);
      betaHat_m.col(0) = betaHat_t1.col(iGrid);
      betaHat_m.col(1) = betaHat_t2.col(iGrid);
      ab_m = alphaHat_m % betaHat_m;
      pC(0) = sum(ab_m.col(0));
      pC(1) = sum(ab_m.col(1));
      //
      // check reads
      //
      Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
      int read_wif_iRead = as<int>(readData[1]);
      while(!done_reads && (read_wif_iRead == iGrid)) {
	//
	// key bit here
	//
        eMatRead_t_col = eMatRead_t.col(iRead);          
	pA1(h_rC) = arma::sum(ab_m.col(h_rC) / eMatRead_t_col);
	pA1(h_rA1) = arma::sum(ab_m.col(h_rA1) % eMatRead_t_col);
	//
	// boring stuff afterwards with updating probabilities etc
	//
        iRead++;
	if ((nReads - 1) < iRead) {
	  done_reads = true;
	  read_wif_iRead = -1;
	} else {
            readData = as<Rcpp::List>(sampleReads[iRead]);
            read_wif_iRead = as<int>(readData[1]);
	}
      }
      //
      // here and in the boring stuff, things would be updated, etc
      //
    }
}



//' @export
// [[Rcpp::export]]
void rcpp_test2(
    arma::mat& alphaHat_t1,
    arma::mat& alphaHat_t2,
    arma::mat& betaHat_t1,
    arma::mat& betaHat_t2,
    const Rcpp::List& sampleReads,
    arma::mat& eMatRead_t    
) {
    int K = alphaHat_t1.n_rows;
    int nGrids = alphaHat_t1.n_cols;
    int nReads = eMatRead_t.n_cols;
    arma::mat ab_m(K, 2);
    Rcpp::NumericVector pC(3), pA1(3), pA2(3);
    arma::mat alphaHat_m(K, 3); // tested, this is the better orientation
    arma::mat betaHat_m(K, 3);
    arma::colvec eMatRead_t_col;    
    //
    // loop over grids
    //
    int iRead = 0;
    bool done_reads = false;
    int h_rC = 0;
    int h_rA1 = 1;
    for(int iGrid = 0; iGrid < nGrids; iGrid++) {
      //
      // reset stuff for this grid
      //
      alphaHat_m.col(0) = alphaHat_t1.col(iGrid);
      alphaHat_m.col(1) = alphaHat_t2.col(iGrid);
      betaHat_m.col(0) = betaHat_t1.col(iGrid);
      betaHat_m.col(1) = betaHat_t2.col(iGrid);
      ab_m = alphaHat_m % betaHat_m;
      pC(0) = sum(ab_m.col(0));
      pC(1) = sum(ab_m.col(1));
      //
      // check reads
      //
      Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
      int read_wif_iRead = as<int>(readData[1]);
      while(!done_reads && (read_wif_iRead == iGrid)) {
	//
	// key bit here
	//
        eMatRead_t_col = eMatRead_t.col(iRead);
        pA1.fill(0);
        for(int k = 0; k < K; k++) {
            pA1(h_rC) += ab_m(k, h_rC) / eMatRead_t_col(k);
            pA1(h_rA1) += ab_m(k, h_rC) * eMatRead_t_col(k);
        }
	//
	// boring stuff afterwards with updating probabilities etc
	//
        iRead++;
	if ((nReads - 1) < iRead) {
	  done_reads = true;
	  read_wif_iRead = -1;
	} else {
            readData = as<Rcpp::List>(sampleReads[iRead]);
            read_wif_iRead = as<int>(readData[1]);
	}
      }
      //
      // here and in the boring stuff, things would be updated, etc
      //
    }
}


//' @export
// [[Rcpp::export]]
void rcpp_test3(
    arma::mat& alphaHat_t1,
    arma::mat& alphaHat_t2,
    arma::mat& betaHat_t1,
    arma::mat& betaHat_t2,
    const Rcpp::List& sampleReads,
    arma::mat& eMatRead_t    
) {
    int K = alphaHat_t1.n_rows;
    int nGrids = alphaHat_t1.n_cols;
    int nReads = eMatRead_t.n_cols;
    arma::mat ab_m(K, 2);
    Rcpp::NumericVector pC(3), pA1(3), pA2(3);
    arma::mat alphaHat_m(K, 3); // tested, this is the better orientation
    arma::mat betaHat_m(K, 3);
    arma::colvec eMatRead_t_col;    
    //
    // loop over grids
    //
    int iRead = 0;
    bool done_reads = false;
    int h_rC = 0;
    int h_rA1 = 1;
    for(int iGrid = 0; iGrid < nGrids; iGrid++) {
      //
      // reset stuff for this grid
      //
      alphaHat_m.col(0) = alphaHat_t1.col(iGrid);
      alphaHat_m.col(1) = alphaHat_t2.col(iGrid);
      betaHat_m.col(0) = betaHat_t1.col(iGrid);
      betaHat_m.col(1) = betaHat_t2.col(iGrid);
      ab_m = alphaHat_m % betaHat_m;
      pC(0) = sum(ab_m.col(0));
      pC(1) = sum(ab_m.col(1));
      //
      // check reads
      //
      Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
      int read_wif_iRead = as<int>(readData[1]);
      while(!done_reads && (read_wif_iRead == iGrid)) {
	//
	// key bit here
	//
        if (h_rC == 0) {
            pA1(h_rC) = arma::sum(alphaHat_t1(iGrid) / eMatRead_t.col(iRead));
            pA1(h_rA1) = arma::sum(alphaHat_t2(iGrid) / eMatRead_t.col(iRead));
        } else {
            pA1(h_rA1) = arma::sum(alphaHat_t1(iGrid) / eMatRead_t.col(iRead));
            pA1(h_rC) = arma::sum(alphaHat_t2(iGrid) / eMatRead_t.col(iRead));
        }                                  
	//
	// boring stuff afterwards with updating probabilities etc
	//
        iRead++;
	if ((nReads - 1) < iRead) {
	  done_reads = true;
	  read_wif_iRead = -1;
	} else {
            readData = as<Rcpp::List>(sampleReads[iRead]);
            read_wif_iRead = as<int>(readData[1]);
	}
      }
      //
      // here and in the boring stuff, things would be updated, etc
      //
    }
}
