// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; -*-

#include <Rcpp.h>
#include <fstream>
#include "vcfpp.h"

using namespace Rcpp;
using namespace vcfpp;

typedef uint32_t uint;
typedef uint64_t uint2;

/// build and dump pbwt FM-index to disk;
//' @export
// [[Rcpp::export]]
void pbwt_index(const std::string& vcffile, const std::string& samples, const std::string& region)
{
    BcfReader vcf(vcffile, samples, region);
    BcfRecord var(vcf.header);
    std::vector<bool> gt;
    std::vector<std::vector<bool>> x;
    int nsnps = 0, nsamples = vcf.nsamples;
    while (vcf.getNextVariant(var))
    {
        if (!var.isSNP())
            continue; // skip other type of variants
        var.getGenotypes(gt);
        x.push_back(gt);
        nsnps++;
    }
    // build pbwt index
    uint N = nsnps, M = nsamples * 2;
    std::ofstream ofpbwt(vcffile + ".pbwt", std::ios::binary);
    std::ofstream ofauxu(vcffile + ".auxu", std::ios::binary);
    std::ofstream ofauxv(vcffile + ".auxv", std::ios::binary);
    ofpbwt.write(reinterpret_cast<char*>(&N), 4);
    ofpbwt.write(reinterpret_cast<char*>(&M), 4);
    std::vector<uint> at(M), a(M), a0(M), a1(M), u(M + 1), v(M + 1);
    uint i, j, u_, v_, a_;
    for (j = 0; j < N; j++)
    {
        u_ = 0; v_ = 0;
        // copy a to at
        for (i = 0; i < M; i++) { at[i] = a[i]; }
        for (i = 0; i < M; i++)
        {
            a_ = j > 0 ? at[i] : i;
            u[i] = u_; v[i] = v_;
            if (x[j][a_])
                a1[v_++] = a_;
            else
                a0[u_++] = a_;
        }
        u[M] = u_; v[M] = M;
        for (i = 0; i < M; i++)
        {
            v[i] += u_;
            if (i < u_)
                a[i] = a0[i];
            else
                a[i] = a1[i - u_];
        }
        if (!ofpbwt.write(reinterpret_cast<char*>(&a[0]), a.size() * 4))
            throw std::runtime_error(strerror(errno));
        if (!ofauxu.write(reinterpret_cast<char*>(&u[0]), u.size() * 4))
            throw std::runtime_error(strerror(errno));
        if (!ofauxv.write(reinterpret_cast<char*>(&v[0]), v.size() * 4))
            throw std::runtime_error(strerror(errno));
    }
}

//' @export
// [[Rcpp::export]]
std::vector<int> pbwt_query(const std::string& vcffile, const std::vector<int>& z, int s, int L, int Step)
{
    uint N, M, i, j;
    std::ifstream ifpbwt(vcffile + ".pbwt", std::ios::binary);
    if (!ifpbwt.read(reinterpret_cast<char*>(&N), 4))
        throw std::runtime_error(strerror(errno));
    if (!ifpbwt.read(reinterpret_cast<char*>(&M), 4))
        throw std::runtime_error(strerror(errno));
    assert(z.size() == N);
    std::ifstream ifauxu(vcffile + ".auxu", std::ios::binary);
    std::ifstream ifauxv(vcffile + ".auxv", std::ios::binary);
    std::vector<uint> t(N), v(M + 1), u(M + 1), a(M);
    std::vector<int> matches;
    for (i = 0; i < N; i++)
    {
        ifauxu.read(reinterpret_cast<char*>(&u[0]), 4 * (M + 1));
        ifauxv.read(reinterpret_cast<char*>(&v[0]), 4 * (M + 1));
        if (i == 0)
        {
            t[0] = z[0] ? M : 0;
        }
        else
        {
            if (z[i])
                t[i] = v[t[i - 1]];
            else
                t[i] = u[t[i - 1]];
        }
        if (s < N && i == s)
        {
            ifpbwt.seekg((uint2)s * M * 4 + 8, std::ios_base::beg);
            ifpbwt.read(reinterpret_cast<char*>(&a[0]), 4 * M);
            if (t[s] == M)
            {
                // where we hit the bottom
                for (j = 1; j <= L; j++)
                    matches.push_back(a[M - j]);
            }
            else if (t[s] == 0)
            {
                // where we hit the top;
                for (j = 0; j < L; j++)
                    matches.push_back(a[0 + j]);
            }
            else
            {
                // L haps before z;
                for (j = 1; j <= L; j++)
                    matches.push_back(a[(t[s] - j) > 0 ? (t[s] - j) : 0]);
                // L haps after z;
                for (j = 0; j < L; j++)
                    matches.push_back(a[(t[s] + j) > M - 1] ? (t[s] + j) : (M - 1));
            }
            s += Step;
        }
    }
    return matches;
}


//' @export
// [[Rcpp::export]]
List pbwt_build(String vcf, String samples, String region, int N, int M) {
    // N is number of SNPs, M is number of Haps
    IntegerMatrix a(N, M), u(N, M + 1), v(N, M + 1);
    IntegerVector a0(M), a1(M);
    // get genotype from vcf/bcf
    BcfReader br(vcf, samples, region);
    BcfRecord var(br.header);
    std::vector<char> gt; // can be bool, char, int
    int i =0, k = 0, u_= 0, v_= 0, a_= 0;
    while (br.getNextVariant(var)) {
      if (!var.isSNP()) continue;
      var.getGenotypes(gt);
      assert(gt.size() == M);
      u_ = 0, v_ = 0;
      for (i = 0; i < M; i++)
      {
          a_ = k > 0 ? a(k - 1, i) : i;
          u(k, i) = u_;
          v(k, i) = v_;
          if (gt[a_])
          {
              a1[v_] = a_;
              v_++;
          }
          else
          {
              a0[u_] = a_;
              u_++;
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
          }
          else
          {
              a(k, i) = a1[i - u_];
          }
      }
      k++;
    }
    if (k != N) stop("nSNPs doesn't match the reference loaded by QUILT!\n");

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
    IntegerMatrix u = as<IntegerMatrix>(p["u"]);
    IntegerMatrix v = as<IntegerMatrix>(p["v"]);
    IntegerMatrix a = as<IntegerMatrix>(p["a"]);
    NumericVector rng = runif(1); // flip the coin to make a decision;
    int s = floor(rng[0] * Step); // the first start of selection;
    std::vector<int> out;
    int i = 0, j =0;
    for (i = 0; i < N; i++)
    {
        if (i == 0) {
            t[0] = z[0] ? M : 0;
        } else {
            if (z[i])
                t[i] = v(i, t[i - 1]);
            else
                t[i] = u(i, t[i - 1]);
        }
        if (s < N && i == s)
        {
            if (t[s] == M) {
                // where we hit the bottom
                for (j = 1; j <= L; j++)
                    out.push_back(a(s, M-j));
            } else if (t[s] == 0) {
                // where we hit the top;
                for (j = 0; j < L; j++)
                    out.push_back(a(s, 0+j));
            } else {
                // L haps before z;
                for (j = 1; j <= L; j++)
                    out.push_back(a(s, std::max(t[s]-j, 0)));
                // L haps after z;
                for (j = 0; j < L; j++)
                    out.push_back(a(s, std::min(t[s]+j, M-1)));
            }
            s += Step;
        }
    }
    return out;
}
