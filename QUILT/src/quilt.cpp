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
#include <math.h>
using namespace Rcpp;



//' @export
// [[Rcpp::export]]
double Rcpp_quilt_test_doubler(
    double a
) {
  return(a * a);
}


//' @export
// [[Rcpp::export]]
int Rcpp_raw_test(
    Rcpp::RawVector a,
    Rcpp::IntegerVector b
) {
    int c = 0;
    for(int i = 0; i < a.length(); i++) {
        c += b(a(i));
    }
    return(c);
}

//' @export
// [[Rcpp::export]]
int Rcpp_raw_test_int(
    Rcpp::IntegerVector a,
    Rcpp::IntegerVector b
) {
    int c = 0;
    for(int i = 0; i < a.length(); i++) {
        c += b(a(i));
    }
    return(c);
}
