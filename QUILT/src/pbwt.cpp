// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; -*-

#include <Rcpp.h>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "vcfpp.h"

using namespace Rcpp;
using namespace vcfpp;


//' @export
// [[Rcpp::export]]
Rcpp::List pbwt_index(String vcf, String samples, String region) {
    // get genotype from vcf/bcf
    BcfReader br(vcf, samples, region);
    BcfRecord var(br.header);
    std::vector<bool> gt; // can be bool, char, int
    std::vector<std::vector<bool>> x;
    int nsnps = 0, nsamples = br.nsamples;
    while (br.getNextVariant(var))
    {
        if (!var.isSNP())
            continue; // skip other type of variants
        var.getGenotypes(gt);
        x.push_back(gt);
        nsnps++;
    }
    // build pbwt index
    int N = nsnps, M = nsamples * 2;
    // N is number of SNPs, M is number of Haps
    IntegerMatrix a(N, M), u(N, M + 1), v(N, M + 1);
    IntegerVector a0(M), a1(M);
    int i, j, u_, v_, a_;
    for (j = 0; j < N; j++)
    {
        u_ = 0; v_ = 0;
        for (i = 0; i < M; i++)
        {
            a_ = j > 0 ? a(j - 1, i) : i;
            u(j, i) = u_;
            if (x[j][a_])
                a1[v_++] = a_;
            else
                a0[u_++] = a_;
        }
        u(j, M) = u_;
        for (i = 0; i < M; i++)
        {
            if (i < u_)
                a(j, i) = a0[i];
            else
                a(j, i) = a1[i - u_];
        }
    }

    Rcpp::List out = List::create(Rcpp::Named("a") = a,
                        Rcpp::Named("u") = u,
                        Rcpp::Named("n") = N,
                        Rcpp::Named("m") = M
                       );
    out.attr("class") = "pbwt";

    return out;
}

//' @export
// [[Rcpp::export]]
std::vector<int> pbwt_query(List p, IntegerVector z, int L = 2, int Step = 8) {
    if (!p.inherits("pbwt")) stop("Input must contain pbwt struct!");
    int N = p["n"];
    int M = p["m"];
    if (z.size() != N) stop("the query z must has the same number of sites N as the panel!");
    IntegerVector t(N);
    IntegerMatrix u = as<IntegerMatrix>(p["u"]);
    IntegerMatrix a = as<IntegerMatrix>(p["a"]);
    std::vector<int> matches;
    std::vector<int> selects; // random pick a selection point at each Step;
    int i, j, s, k = 0;
    NumericVector rng = runif( N / Step); // flip the coin to make a decision;
    for (i = 1; i < rng.size(); i++)
        selects.push_back(i * Step + floor(rng[i] * Step));
    for (i = 0; i < N; i++)
    {
        if (i == 0) {
            t[0] = z[0] ? M : 0;
        } else {
            if (z[i])
                t[i] = u(i, M) + t[i - 1] - u(i, t[i - 1]);
            else
                t[i] = u(i, t[i - 1]);
        }
        if (selects[k] < N && selects[k] == i)
        {
            s = selects[k];
            if (t[s] == M) {
                // where we hit the bottom
                for (j = 1; j <= L; j++)
                    matches.push_back(a(s, M-j));
            } else if (t[s] == 0) {
                // where we hit the top;
                for (j = 0; j < L; j++)
                    matches.push_back(a(s, 0+j));
            } else {
                // L haps before z;
                for (j = 1; j <= L; j++)
                    matches.push_back(a(s, std::max(t[s]-j, 0)));
                // L haps after z;
                for (j = 0; j < L; j++)
                    matches.push_back(a(s, std::min(t[s]+j, M-1)));
            }
            k = k + 1;
        }
    }
    return matches;
}

// obslote
std::vector<int> old_pbwt_query(const std::string& pbwtfile, const std::vector<int>& z, int L, int Step)
{
    int N, M, i, j, s, k = 0;
    std::ifstream ifpbwt(pbwtfile + ".pbwt", std::ios::binary);
    if (!ifpbwt.read(reinterpret_cast<char*>(&N), 4))
        throw std::runtime_error(strerror(errno));
    if (!ifpbwt.read(reinterpret_cast<char*>(&M), 4))
        throw std::runtime_error(strerror(errno));
    assert(z.size() == N);
    std::ifstream ifauxu(pbwtfile + ".auxu", std::ios::binary);
    std::vector<int> t(N), u(M + 1), a(M);
    std::vector<int> matches;
    std::vector<int> selects; // random pick a selection point at each Step;
    NumericVector rng = runif( N / Step); // flip the coin to make a decision;
    // for (i = 1; i < std::ceil(N / Step); i++)
    for (i = 1; i < rng.size(); i++)
        selects.push_back(i * Step + floor(rng[i] * Step));
    for (i = 0; i < N; i++)
    {
        ifauxu.read(reinterpret_cast<char*>(&u[0]), 4 * (M + 1));
        if (i == 0)
        {
            t[0] = z[0] ? M : 0;
        }
        else
        {
            if (z[i])
                t[i] = u[M] + t[i - 1] - u[t[i - 1]];
            else
                t[i] = u[t[i - 1]];
        }
        if (selects[k] < N && selects[k] == i)
        {
            s = selects[k];
            ifpbwt.seekg((uint64_t)s * M * 4 + 8, std::ios_base::beg);
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
                    matches.push_back(a[(t[s] + j) < M - 1 ? (t[s] + j) : (M - 1)]);
            }
            k = k + 1;
        }
    }
    return matches;
}
