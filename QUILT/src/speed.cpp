// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

#include <algorithm>
#include <cmath>
#include <ctype.h>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
void rcpp_test1(arma::mat& alphaHat_t1, arma::mat& alphaHat_t2, arma::mat& betaHat_t1, arma::mat& betaHat_t2,
                const Rcpp::List& sampleReads, arma::mat& eMatRead_t)
{
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
    for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
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
        while (!done_reads && (read_wif_iRead == iGrid))
        {
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
            if ((nReads - 1) < iRead)
            {
                done_reads = true;
                read_wif_iRead = -1;
            }
            else
            {
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
void rcpp_test2(arma::mat& alphaHat_t1, arma::mat& alphaHat_t2, arma::mat& betaHat_t1, arma::mat& betaHat_t2,
                const Rcpp::List& sampleReads, arma::mat& eMatRead_t)
{
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
    for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
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
        while (!done_reads && (read_wif_iRead == iGrid))
        {
            //
            // key bit here
            //
            eMatRead_t_col = eMatRead_t.col(iRead);
            pA1.fill(0);
            for (int k = 0; k < K; k++)
            {
                pA1(h_rC) += ab_m(k, h_rC) / eMatRead_t_col(k);
                pA1(h_rA1) += ab_m(k, h_rC) * eMatRead_t_col(k);
            }
            //
            // boring stuff afterwards with updating probabilities etc
            //
            iRead++;
            if ((nReads - 1) < iRead)
            {
                done_reads = true;
                read_wif_iRead = -1;
            }
            else
            {
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
void rcpp_test3(arma::mat& alphaHat_t1, arma::mat& alphaHat_t2, arma::mat& betaHat_t1, arma::mat& betaHat_t2,
                const Rcpp::List& sampleReads, arma::mat& eMatRead_t)
{
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
    for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
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
        while (!done_reads && (read_wif_iRead == iGrid))
        {
            //
            // key bit here
            //
            eMatRead_t_col = eMatRead_t.col(iRead);
            // A1 - original hap loses
            pA1(h_rC) = sum(ab_m.col(h_rC) / eMatRead_t_col);
            // A1 - new hap gains
            pA1(h_rA1) = sum(ab_m.col(h_rA1) % eMatRead_t_col);
            //
            // boring stuff afterwards with updating probabilities etc
            //
            iRead++;
            if ((nReads - 1) < iRead)
            {
                done_reads = true;
                read_wif_iRead = -1;
            }
            else
            {
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
void rcpp_test4(Eigen::Map<Eigen::ArrayXXd> alphaHat_t1, Eigen::Map<Eigen::ArrayXXd> alphaHat_t2,
                Eigen::Map<Eigen::ArrayXXd> betaHat_t1, Eigen::Map<Eigen::ArrayXXd> betaHat_t2,
                const Rcpp::List& sampleReads, Eigen::Map<Eigen::ArrayXXd> eMatRead_t)
{
    const int K = alphaHat_t1.rows();
    const int nGrids = alphaHat_t1.cols();
    const int nReads = eMatRead_t.cols();
    Rcpp::NumericVector pC(3), pA1(3), pA2(3);
    Eigen::ArrayXXd alphaHat_m(K, 2); // tested, this is the better orientation
    Eigen::ArrayXXd betaHat_m(K, 2);
    Eigen::ArrayXXd ab_m(K, 2);
    //
    // loop over grids
    //
    int iRead = 0;
    bool done_reads = false;
    int h_rC = 0;
    int h_rA1 = 1;
    for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
        //
        // reset stuff for this grid
        //
        alphaHat_m.col(0) = alphaHat_t1.col(iGrid);
        alphaHat_m.col(1) = alphaHat_t2.col(iGrid);
        betaHat_m.col(0) = betaHat_t1.col(iGrid);
        betaHat_m.col(1) = betaHat_t2.col(iGrid);
        ab_m = alphaHat_m * betaHat_m;
        pC(0) = ab_m.col(0).sum();
        pC(1) = ab_m.col(1).sum();
        //
        // check reads
        //
        Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
        int read_wif_iRead = as<int>(readData[1]);
        while (!done_reads && (read_wif_iRead == iGrid))
        {
            //
            // key bit here
            //
            // A1 - original hap loses
            pA1(h_rC) = (ab_m.col(h_rC) / eMatRead_t.col(iRead)).sum();
            // A1 - new hap gains
            pA1(h_rA1) = (ab_m.col(h_rA1) * eMatRead_t.col(iRead)).sum();
            //
            // boring stuff afterwards with updating probabilities etc
            //
            iRead++;
            if ((nReads - 1) < iRead)
            {
                done_reads = true;
                read_wif_iRead = -1;
            }
            else
            {
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
void rcpp_test5(arma::mat& alphaHat_t1, arma::mat& alphaHat_t2, arma::mat& betaHat_t1, arma::mat& betaHat_t2,
                const Rcpp::List& sampleReads, arma::mat& eMatRead_t)
{
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
    for (int iGrid = 0; iGrid < nGrids; iGrid++)
    {
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
        while (!done_reads && (read_wif_iRead == iGrid))
        {
            //
            // key bit here
            //
	  pA1(h_rC) = arma::sum(ab_m.col(h_rC) / eMatRead_t.col(iRead));
	  pA1(h_rA1) = arma::sum(ab_m.col(h_rA1) % eMatRead_t.col(iRead));
            //
            // boring stuff afterwards with updating probabilities etc
            //
            iRead++;
            if ((nReads - 1) < iRead)
            {
                done_reads = true;
                read_wif_iRead = -1;
            }
            else
            {
                readData = as<Rcpp::List>(sampleReads[iRead]);
                read_wif_iRead = as<int>(readData[1]);
            }
        }
        //
        // here and in the boring stuff, things would be updated, etc
        //
    }
}
