// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; -*-

#include <Rcpp.h>
#include "vcfpp.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <unordered_map>

using namespace Rcpp;
using namespace vcfpp;
using namespace std;

using grid_t = uint16_t; // options: uint8_t, uint16_t

using IntGridVec = vector<grid_t>;
using IntSet = set<grid_t, less<grid_t>>;
using IntMap = unordered_map<grid_t, int32_t>; // use double to replace uint32_t as symbol
using IntVecMap = unordered_map<grid_t, vector<uint32_t>>;
using SymbolIdxMap = map<int32_t, int32_t>;

class Timer
{
protected:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_timing_clock, prev_timing_clock;

public:
    Timer();
    ~Timer();
    void clock();
    unsigned int reltime();
    unsigned int abstime();
};

Timer::Timer()
{
    start_timing_clock = std::chrono::high_resolution_clock::now();
}

Timer::~Timer()
{
}

void Timer::clock()
{
    prev_timing_clock = std::chrono::high_resolution_clock::now();
}

unsigned int Timer::reltime()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - prev_timing_clock).count();
}

unsigned int Timer::abstime()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start_timing_clock).count();
}
template <typename T>
T reverseBits(T n, size_t B = sizeof(T) * 8)
{
    assert(B <= std::numeric_limits<T>::digits);
    T rv = 0;
    for (size_t i = 0; i < B; ++i, n >>= 1)
        rv = (rv << 1) | (n & 0x01);
    return rv;
}

IntGridVec encodeZgrid(const IntegerVector& z, int G)
{
    IntGridVec zg(G);
    const int B = sizeof(grid_t) * 8;
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
        throw std::runtime_error("something wrong! G:" + to_string(G) + ", k: " + to_string(k) + "\n");
    }
    return zg;
}

IntVecMap build_W(const IntGridVec& x, const IntSet& s, const IntMap& C)
{
    uint32_t c{0}, i{0};
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
    }

    return W;
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

IntegerMatrix save_W(NumericVector x, const IntSet& s)
{
    size_t i{0}, n{0}, k{0};
    IntegerMatrix Wg(2, s.size()); // first row is symbol index, second is symbol occ
    for (const auto& si : s)
    {
        k = 0;
        for (i = 0; i < x.size(); i++)
        {
            if (x[i] == si)
            {
                Wg(0, n) = i;
                Wg(1, n) = k;
                k++;
            }
        }
        n++;
    }
    return Wg;
}

IntMap build_Symbols(const IntSet& s)
{
    uint32_t n{0};
    IntMap symbol;
    for (const auto& si : s)
    {
        symbol[si] = n++;
    }
    return symbol;
}

//' @export
// [[Rcpp::export]]
Rcpp::List mspbwt_index(const std::string& vcfpanel, const std::string& samples, const std::string& region, bool fast = 1)
{
    Timer tm;
    tm.clock();
    const int B = sizeof(grid_t) * 8;
    BcfReader vcf(vcfpanel, samples, region);
    BcfRecord var(vcf.header);
    // N haplotypes, M SNPs, G Grids, k Current grid, m Current SNP
    uint64_t N{0}, M{0}, G{0}, k{0}, m{0}, i{0}, j{0};
    N = vcf.nsamples * 2;
    while (vcf.getNextVariant(var))
        M++;
    G = (M + B - 1) / B;
    Rcout << "Haps(N):" << N << "\tSNPs(M):" << M << "\tGrids(G):" << G << "\tInt(B):" << B << endl;
    vector<IntGridVec> X; // Grids x Haps
    IntegerMatrix XG(N, G); // (Grids+1) x Haps
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
            for (i = 0; i < N; i++)
            {
                X[k][i] = reverseBits(X[k][i]); // reverset bits
                XG(i, k) = X[k][i];
            }
            k++;                                // update next grid
        }
    }
    if (G == k + 1)
    {
        for (i = 0; i < N; i++)
        {
            X[k][i] <<= G * B - M;
            X[k][i] = reverseBits(X[k][i]); // reverset bits
            XG(i, k) = X[k][i];
        }
    }
    else if (G == k)
    {
        cerr << "no need padding\n";
    }
    else
    {
        throw std::runtime_error("something wrong! G:" + to_string(G) + ", k: " + to_string(k) + "\n");
    }
    Rcout << "elapsed time of compressing panel: " << tm.reltime() << " milliseconds" << endl;


    tm.clock();
    IntGridVec y0(N);
    IntegerMatrix A(G+1, N); // (Grids+1) x Haps
    IntegerVector a0(N);
    Rcpp::List C(G);
    Rcpp::List W(G);
    Rcpp::List Symbols(G);
    IntSet s; // convert to set which unique sorted
    IntMap sm, Ct;
    IntVecMap Wt;
    vector<vector<int32_t>> Wg;
    vector<int32_t> idx, Cg;
    uint32_t c{0}, w{0}, occ{0};
    for (k = 0; k < G; k++)
    {
        if (k == 0)
        {
            for (j = 0; j < N; j++) {
                a0(j) = j;
                A(k, j) = a0(j);
            }
        }
        for (i = 0; i < N; i++)
        {
            y0[i] = X[k][a0[i]];
            s.insert(y0[i]);
        }
        // build symbols
        sm = build_Symbols(s);
        // build C
        Cg.resize(s.size());
        j=0;
        for (const auto& si : s)
        {
            c = 0;
            for (i = 0; i < N; i++)
            {
                if (y0[i] < si)
                    c++;
                if (y0[i] == si)
                    idx.push_back(i);
            }
            Cg[j] = c;
            Wg.push_back(idx);
            idx.clear();
            j++;
        }
        C[k] = Cg;
        W[k] = Wg;
        Symbols[k] = wrap(s);
        // build A from W
        if ((s.size() < 256) || fast)
        {
            // small number of symbols then build_W or forcing fast mode
            Ct = build_C(y0, s);
            Wt = build_W(y0, s, Ct);
            for (i = 0; i < N; i++)
                A(k + 1, Wt[y0[i]][i] - 1) = a0[i];
        }
        else
        {
            // too symbols then do not build_W directly
            for (i = 0; i < N; i++)
            {
                // map y0 to index first
                j = sm[y0[i]];
                w = Cg[j];
                if (i >= Wg[j][0]) {
                    for (occ = 0; occ < Wg[j].size(); occ++) {
                        if(i == Wg[j][occ])
                        {
                            w += occ + 1;
                            break;
                        }
                    }
                }
                A(k + 1, w - 1) = a0[i];
            }
        }
        // next run
        a0 = A(k+1, _);
        s.clear();
        Wg.clear();
    }
    Rcout << "elapsed time of building indices: " << tm.reltime() << " milliseconds" << endl;

    Rcpp::List out = List::create(Named("A",A),
                                  Named("C",C),
                                  Named("W", W),
                                  Named("X", XG),
                                  Named("Symbols", Symbols),
                                  Named("G", G),
                                  Named("N", N),
                                  Named("M", M));

    out.attr("class") = "mspbwt";

    return out;
}

//' @export
// [[Rcpp::export]]
IntegerVector mspbwt_query(const IntegerMatrix& XG,const IntegerMatrix& A, const List& C, const List& W, const List& S, int G, int M, int N, const IntegerVector& z, int L = 1)
{
    assert(M == z.size());
    Timer tm;
    tm.clock();
    Rcout << "RefHaps(N):" << N << "\tSNPs(M):" << M << "\tGrids(G):" << G << endl;
    int k = 0, n = 0, s;
    IntGridVec zg = encodeZgrid(z, G);
    IntegerVector az(G);
    vector<vector<double>> Symbols = as< vector<vector<double>> >(S);
    auto kzus = upper_bound(Symbols[k].begin(), Symbols[k].end(), zg[k]) == Symbols[k].begin() ? Symbols[k].begin() : prev(upper_bound(Symbols[k].begin(), Symbols[k].end(), zg[k])) ;
    int j = std::distance(Symbols[k].begin(), kzus); // symbol rank
    List Wk = as<List>(W[k]);
    IntegerVector wkz = Wk[j];
    az[k] = wkz[wkz.size()-1];
    IntegerVector selects;
    // L haps before z;
    for (s = 0; s < L; s++)
    {
        j = A(k+1, std::max(az(k)-s-1, 0));
        if (XG(j, k) == zg[k])
            selects.push_back(j);
        j = A(k+1, std::min(az(k)+s, M-1));
        if (XG(j, k) == zg[k])
            selects.push_back(j);
    }
    Rcout << "elapsed time of processing first grid of z: " << tm.reltime() << " milliseconds" << endl;

    tm.clock();
    for (k = 1; k < G; k++)
    {
        auto kzus = upper_bound(Symbols[k].begin(), Symbols[k].end(), zg[k]) == Symbols[k].begin() ? Symbols[k].begin() : prev(upper_bound(Symbols[k].begin(), Symbols[k].end(), zg[k])) ;
        j = std::distance(Symbols[k].begin(), kzus); // symbol rank. expensive. maybe use List lookup
        az[k] = as<IntegerVector>(C[k])[j];
        if (az[k-1] >= az[k])
        {
            wkz = as<IntegerVector>(as<List>(W[k])[j]);
            auto kzui = upper_bound(wkz.begin(), wkz.end(), az[k - 1]) == wkz.begin() ? wkz.begin() : prev(upper_bound(wkz.begin(), wkz.end(), az[k - 1])) ;
            az[k] = az[k] + std::distance(wkz.begin(), kzui);
        } else {
            n++;
        }
        for (s = 0; s < L; s++)
        {
            j = A(k+1, std::max(az(k)-s-1, 0));
            if (XG(j, k) == zg[k])
                selects.push_back(j);
            j = A(k+1, std::min(az(k)+s, M-1));
            if (XG(j, k) == zg[k])
                selects.push_back(j);
        }
    }
    Rcout << "elapsed time of processing all grids of z: " << tm.reltime() << " milliseconds" << endl;
    Rcout << "elapsed time of mspbwt query: " << tm.abstime() << " milliseconds" << endl;
    Rcout << "skip binary search for " << n << "/" << G << " Grids\n";

    return selects;
}
