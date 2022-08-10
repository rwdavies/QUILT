// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; -*-

#include <Rcpp.h>
#include "vcfpp.h"

using namespace Rcpp;
using namespace vcfpp;


//' @export
// [[Rcpp::export]]
List pbwt_build(String vcf, String samples, String region, int N, int M) {
    // N is number of SNPs, M is number of Haps
    IntegerMatrix a(N, M), d(N, M), u(N, M + 1), v(N, M + 1);
    IntegerVector a0(M), a1(M), d0(M), d1(M);
    // get genotype from vcf/bcf
    BcfReader br(vcf, samples, region);
    BcfRecord var(br.header);
    std::vector<char> gt; // can be bool, char, int
    int i =0, k = 0, u_= 0, v_= 0, p= 0, q= 0;
    while (br.getNextVariant(var)) {
      if (!var.isSNP()) continue;
      var.getGenotypes(gt);
      assert(gt.size() == M);
      u_ = 0, v_ = 0, p = k + 1, q = k + 1;
      for (i = 0; i < M; i++)
      {
          int d_ = k > 0 ? d(k - 1, i) : 0;
          int a_ = k > 0 ? a(k - 1, i) : i;
          p = std::max(p, d_);
          q = std::max(q, d_);
          u(k, i) = u_;
          v(k, i) = v_;
          if (gt[a_])
          {
              a1[v_] = a_;
              d1[v_] = q;
              v_++;
              q = 0;
          }
          else
          {
              a0[u_] = a_;
              d0[u_] = p;
              u_++;
              p = 0;
          }
      }
      u(k, M) = u_;
      v(k, M) = M;
      for (i = 0; i < M; i++)
      {
          v(k, i) += u_;
          if (i < u_)
          {
              a(k, i) = a0[i];
              d(k, i) = d0[i];
          }
          else
          {
              a(k, i) = a1[i - u_];
              d(k, i) = d1[i - u_];
          }
      }
      k++;
    }
    assert(k == N);

    Rcpp::List out = List::create(Rcpp::Named("a") = a,
                        Rcpp::Named("u") = u,
                        Rcpp::Named("v") = v,
                        Rcpp::Named("n") = N,
                        Rcpp::Named("m") = M
                       );
    out.attr("class") = "pbwt";

    return out;
}

//' @export
// [[Rcpp::export]]
std::vector<int> find_neighour_haps(List p, IntegerVector z, int L = 1, int Step = 2) {
    if (!p.inherits("pbwt")) stop("Input must contain pbwt struct!");
    int N = p["n"];
    int M = p["m"];
    if (z.size() != N) stop("the query z must has the same number of sites N as the panel!");
    IntegerVector t(N);
    t[0] = z[0] ? M : 0;
    IntegerMatrix u = as<IntegerMatrix>(p["u"]);
    IntegerMatrix v = as<IntegerMatrix>(p["v"]);
    IntegerMatrix a = as<IntegerMatrix>(p["a"]);
    int i = 0, j =0;
    for (int i = 1; i < N; i++)
    {
        if (z[i])
            t[i] = v(i, t[i - 1]);
        else
            t[i] = u(i, t[i - 1]);
    }
    NumericVector rng = runif(1); // flip the coin to make a decision;
    int s = floor(rng[0] * Step); // the first start of selection;
    std::vector<int> out;
    for (i = s; i < N; i += Step) {
        if (t[i] == M) {
            // where we hit the bottom
            for (j = 1; j <= L; j++)
                out.push_back(a(i, M-j));
        } else if (t[i] == 0) {
            // where we hit the top;
            for (j = 0; j < L; j++)
                out.push_back(a(i, 0+j));
        } else {
            // L haps before z;
            for (j = 1; j <= L; j++)
                out.push_back(a(i, std::max(t[i]-j, 0)));
            // L haps after z;
            for (j = 0; j < L; j++)
                out.push_back(a(i, std::min(t[i]+j, M-1)));
        }
    }
    return out;
}
