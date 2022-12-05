// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; -*-

#include <Rcpp.h>
#include "vcfpp.h"

#include <algorithm>
#include <bitset>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <unordered_map>

using namespace Rcpp;
using namespace vcfpp;
using namespace std;

using IntGridVec = vector<uint32_t>;
using IntSet = set<uint32_t, less<uint32_t>>;
using IntMap = unordered_map<double, int32_t>; // use double to replace uint32_t as symbol
using IntVecMap = unordered_map<uint32_t, std::vector<uint32_t>>;
using SymbolIdxMap = map<double, int32_t>;

template <typename T>
T reverseBits(T n, size_t B = sizeof(T) * 8)
{
    assert(B <= std::numeric_limits<T>::digits);
    T rv = 0;
    for (size_t i = 0; i < B; ++i, n >>= 1)
        rv = (rv << 1) | (n & 0x01);
    return rv;
}

IntMap build_C(const IntGridVec& x, const IntSet& s)
{
    uint32_t c{0}, n{0}, cap{0};
    IntMap C;
    for (const auto& si : s)
    {
        c = 0;
        for (const auto& xi : x)
        {
            if (xi < si)
                c++;
        }
        C[si] = c;
        if (++n == cap && cap > 0)
            break;
    }

    return C;
}

Rcpp::List save_W(const IntGridVec& x, const IntSet& s)
{
    size_t i{0}, n{0}, k{0};
    Rcpp::List W(s.size());
    Rcpp::CharacterVector Wn(s.size());
    for (const auto& si : s)
    {
        k = 0;
        SymbolIdxMap Ws;
        for (i = 0; i < x.size(); i++)
        {
            if (x[i] == si)
                Ws[i] = k++; // k: how many occurrences before i
        }
        W[n] = wrap(Ws);
        Wn[n] = std::to_string(si);
        n++;
    }
    W.attr("names") = Wn;

    return W;
}

IntVecMap build_W(const IntGridVec& x, const IntSet& s, const IntMap& C)
{
    uint32_t c{0}, i{0}, n{0}, cap{0};
    IntVecMap W;
    for (const auto& si : s)
    {
        W[si] = vector<uint32_t>(x.size());
        c = 0;
        for (i = 0; i < x.size(); i++)
        {
            if (x[i] == si)
                c++;
            W[si][i] = c + C.at(si);
        }
        if (++n == cap && cap > 0)
            break;
    }

    return W;
}

//' @export
// [[Rcpp::export]]
Rcpp::List mspbwt_index(const std::string& vcfpanel, const std::string& samples, const std::string& region)
{
    const int B = 32;
    BcfReader vcf(vcfpanel, samples, region);
    BcfRecord var(vcf.header);
    uint64_t N{0}, M{0}, G{0}, k{0}, m{0}, i{0}; // N haplotypes, M SNPs, G Grids, k Current grid, m Current SNP
    N = vcf.nsamples * 2;
    M = vcf.get_region_records(region);
    G = (M + B - 1) / B;
    cerr << "Haps(N):" << N << "\tSNPs(M):" << M << "\tGrids(G):" << G << "\tInt(B):" << B << endl;
    vector<IntGridVec> X; // Grids x Haps
    X.resize(G, IntGridVec(N));
    vcf.setRegion(region); // seek back to region
    vector<bool> gt;
    while (vcf.getNextVariant(var))
    {
        var.getGenotypes(gt);
        // update current grids
        for (i = 0; i < N; i++)
            X[k][i] = (X[k][i] << 1) | (gt[i] != 0);
        m++;
        if (m % B == 0)
        {
            X[k][i] = reverseBits(X[k][i]); // reverset bits
            k++;                            // update next grid
        }
    }
    if (G == k + 1)
    {
        int pad = G * B - M;
        for (i = 0; i < N; i++)
            X[k][i] <<= pad;
        X[k][i] = reverseBits(X[k][i]); // reverset bits
    }
    else if (G == k)
    {
        cerr << "no need padding\n";
    }
    else
    {
        throw std::runtime_error("something wrong\n");
    }

    IntGridVec y0(N);
    vector<vector<int>> A; // (Grids+1) x Haps, we use int here so that R can take it
    A.resize(G + 1, vector<int>(N));
    vector<int> a0(N);
    Rcpp::List C(G);
    Rcpp::List W(G);
    Rcpp::List Symbols(G);
    for (k = 0; k < G; k++)
    {
        if (k == 0)
        {
            std::iota(a0.begin(), a0.end(), 0);
            A[k] = a0;
        }
        for (i = 0; i < N; i++)
            y0[i] = X[k][a0[i]];
        IntSet s(y0.begin(), y0.end()); // convert to set which unique sorted
        IntMap Cg = build_C(y0, s);
        C[k] = wrap(Cg);
        auto Wg = build_W(y0, s, Cg); // here Wg is S x N
        for (i = 0; i < N; i++)
            A[k + 1][Wg[y0[i]][i] - 1] = a0[i];
        // next run
        a0 = A[k + 1];
        // here save current W, which differs from the previous complete table
        W[k] = save_W(y0, s);
        Symbols[k] = wrap(s);
    }

    Rcpp::List out = List::create(Named("A",A),
                                  Named("C",C),
                                  Named("W", W),
                                  Named("Symbols", Symbols),
                                  Named("G", G),
                                  Named("N", N),
                                  Named("M", M));

    out.attr("class") = "mspbwt";

    return out;
}

vector<uint32_t> encodeZgrid(IntegerVector z, int G)
{
    vector<uint32_t> zg(G);
    const int B = 32;
    R_xlen_t m{0}, k{0}, M{z.size()};
    for (m = 0; m < M; m++)
    {
        zg[k] = (zg[k] << 1) | (z[m] != 0);
        if ((m + 1) % B == 0)
        {
            zg[k] = reverseBits(zg[k]);
            k++;
        }
    }
    if (G == k + 1)
    {
        zg[k] <<= G * B - M; // padding 0s
        zg[k] = reverseBits(zg[k]);
    }
    else if (G == k)
    {
        cerr << "no need padding\n";
    }
    else
    {
        throw std::runtime_error("something wrong\n");
    }
    return zg;
}

//' @export
// [[Rcpp::export]]
NumericVector mspbwt_query(List p, IntegerVector z, int L = 2, int Step = 1)
{
    if (!p.inherits("mspbwt")) stop("Input must contain mspbwt struct!");
    int G = p["G"];
    Rcpp::List A = as<Rcpp::List>(p["A"]);
    Rcpp::List W = as<Rcpp::List>(p["W"]);
    Rcpp::List C = as<Rcpp::List>(p["C"]);
    Rcpp::List Symbols = as<Rcpp::List>(p["Symbols"]);
    vector<uint32_t> zg = encodeZgrid(z, G);
    NumericVector az(G); // use int for index to be compatibable to R
    int k = 0;
    NumericVector ks = Symbols[k];
    auto kzus = *std::prev(std::upper_bound(ks.begin(), ks.end(), zg[k]));
    Rcpp::List Wk = W[k];
    Rcpp::List wkz = Wk[std::to_string(kzus)];
    az[k] = wkz[wkz.length()-1];
    for (k = 1; k < G; k++)
    {
        ks = Symbols[k];
        kzus = *std::prev(std::upper_bound(ks.begin(), ks.end(), zg[k]));
        wkz = Wk[to_string(kzus)];
        az[k] = C[k][to_string(kzus)];
        kzus = *std::prev(std::upper_bound(wkz.begin(), wkz.end(), az[k-1]));
        az[k] = az[k] + kzus;
    }

    return az;
}
