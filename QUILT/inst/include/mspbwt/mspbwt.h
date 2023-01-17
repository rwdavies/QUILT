#ifndef MSPBWT_H_
#define MSPBWT_H_

#include "vcfpp/vcfpp.h"

#include <algorithm>
#include <bitset>
#include <fstream>
#include <iterator>
#include <map>
#include <numeric>
#include <random>
#include <unordered_map>
#include <unordered_set>


using namespace std;
using namespace vcfpp;

using IntMapU = unordered_map<int, int>;
using IntVec = vector<int>;
using IntVec2D = vector<IntVec>;
using IntVec3D = vector<IntVec2D>;

template <typename T>
T reverseBits(T n, size_t B = sizeof(T) * 8)
{
    assert(B <= std::numeric_limits<T>::digits);
    T rv = 0;
    for (size_t i = 0; i < B; ++i, n >>= 1)
        rv = (rv << 1) | (n & 0x01);
    return rv;
}

// 0-based
vector<int> inline seq_by(int start, int end, int by)
{
    int n = (end - start + 1) % by == 0 ? (end - start + 1) / by : ((end - start + 1) / by + 1);
    vector<int> seq(n);
    int i{0}, x;
    for (x = start; x <= end; x = x + by)
        seq[i++] = x;
    return seq;
}

// C++11 compatible
template <typename T, typename = typename std::enable_if<std::is_unsigned<T>::value>::type>
class msPBWT
{
private:
    using grid_t = T;
    using GridVec = vector<grid_t>;
    using GridVec2D = vector<GridVec>;
    using GridVec3D = vector<GridVec2D>;
    using GridSetU = unordered_set<grid_t>;
    using GridMapU = unordered_map<grid_t, int>;                 // {symbol : index}
    using GridVecMapU = unordered_map<grid_t, std::vector<int>>; // {symbol : vec(index)}
    using SymbolIdxMap = map<int, int, less<int>>;               // {index: rank}
    using WgSymbolMap = map<grid_t, SymbolIdxMap, less<grid_t>>; // {symbol:{index:rank}}

    int B{sizeof(T) * 8}, N{0}, M{0}, G{0}, G1{0}, G2{0}, nindices{4};
    bool is_save_X{1}, is_save_D{0};
    GridVec2D X; //  Grids x Haps
    GridVec2D S; // Grids x Sorted and Unique symbols
    IntVec3D W;
    IntVec2D C;
    IntVec2D A;          // nindices x Grids x Haps
    IntVec2D D;          // nindices x Grids x Haps
    vector<int> reorder; // (M) reorder SNPs or just subset SNPs

public:
    msPBWT(int nindices_ = 4) : nindices(nindices_)
    {
    }

    virtual ~msPBWT(){};

    bool verbose{1};

    vector<int> randhapz()
    {
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> dist(0, 101);
        vector<int> z(M);
        for (int i = 0; i < M; i++)
        {
            z[i] = dist(gen) > 95;
        }
        return z;
    }

    GridVec encodezg(const vector<int>& z_)
    {
        assert(z_.size() == M);
        int k{0};
        GridVec zg(G);
        vector<bool> z(M);
        for (k = 0; k < M; k++)
            z[k] = z_[reorder[k]] != 0;
        k = 0;
        for (size_t m = 0; m < z.size(); m++)
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
            zg[k] <<= G * B - M; // padding zeros
            zg[k] = reverseBits(zg[k]);
        }
        else if (G == k)
        {
            if (verbose)
                cerr << "no need padding\n";
        }
        else
        {
            throw std::runtime_error("something wrong\n");
        }
        return zg;
    }

    GridVec randzg()
    {
        auto z = randhapz();
        return encodezg(z);
    }

    GridMapU build_C(const GridVec& y, const GridVec& s)
    {
        int c{0};
        GridMapU C;
        for (const auto& si : s)
        {
            c = 0;
            for (const auto& xi : y)
            {
                if (xi < si)
                    c++;
            }
            C[si] = c;
        }

        return C;
    }

    WgSymbolMap save_W(const GridVec& y, const GridVec& s)
    {
        int k{0};
        WgSymbolMap W;
        for (const auto& si : s)
        {
            k = 0;
            for (size_t i = 0; i < y.size(); i++)
            {
                if (y[i] == si)
                    W[si][i] = k++; // k: how many occurrences before i
            }
        }

        return W;
    }

    vector<int> save_C(const GridVec& y, const GridVec& s)
    {
        int c{0}, i{0};
        vector<int> C(s.size());
        for (const auto& si : s)
        {
            c = 0;
            for (const auto& yi : y)
            {
                if (yi < si)
                    c++;
            }
            C[i++] = c;
        }

        return C;
    }

    vector<vector<int>> save_Occ(const GridVec& y, const GridVec& s)
    {
        vector<vector<int>> Occ(s.size());
        vector<int> idx;
        size_t i{0};
        int k{0};
        for (const auto& si : s)
        {
            for (i = 0; i < y.size(); i++)
            {
                if (y[i] == si)
                    idx.push_back(i);
            }
            Occ[k++] = idx;
            idx.clear();
        }

        return Occ;
    }

    GridVecMapU build_W(const GridVec& y, const GridVec& s)
    {
        int c{0};
        size_t i{0};
        auto C = build_C(y, s);
        GridVecMapU W;
        for (const auto& si : s)
        {
            W[si] = vector<int>(y.size());
            c = 0;
            for (i = 0; i < y.size(); i++)
            {
                if (y[i] == si)
                    c++;
                W[si][i] = c + C.at(si);
            }
        }

        return W;
    }

    void build_indices_k(int ni, int ki, const GridVec& xk, GridVec& yk, GridVec& sk,
                         vector<int>& ck, vector<vector<int>>& wk, vector<int>& Occ,
                         vector<vector<int>>& kas, vector<vector<int>>& kds, vector<int>& sqp,
                         vector<int>& a0, vector<int>& d0, GridSetU& symbols)
    {
        int s, n, c, e, i;
        for (n = 0; n < N; n++)
        {
            yk[n] = xk[a0[n]];
            symbols.insert(yk[n]);
        }
        sk = GridVec(symbols.begin(), symbols.end());
        std::sort(sk.begin(), sk.end());
        symbols.clear();
        // C[k] = save_C(y1, S[k]);
        // W[k] = save_Occ(y1, S[k]);
        // auto Wg = build_W(y1, S[k]); // here Wg is S x N
        // for (n = 0; n < N; n++)
        //     A[k][Wg[y1[n]][n] - 1] = a0[n];
        // next run
        // a0 = A[k];

        int sn = sk.size();
        ck.resize(sn);
        wk.resize(sn);
        for (s = 0; s < sn; s++)
        {
            c = 0;
            for (n = 0; n < N; n++)
            {
                if (yk[n] < sk[s])
                    c++;
                if (yk[n] == sk[s])
                    Occ.push_back(n);
            }
            ck[s] = c;
            wk[s] = Occ;
            Occ.clear();
        }

        // save A and D
        kas.resize(sn);
        if (is_save_D)
        {
            kds.resize(sn);
            sqp.resize(sn, ki + 1);
        }
        // iterate over all haplotyps in reverese prefix sorted order
        for (n = 0; n < N; n++)
        {
            i = a0[n]; // index ihap at Y[k]
            if (is_save_D)
                e = d0[n];
            for (s = 0; s < sn; s++)
            {
                if (is_save_D)
                {
                    if (e > sqp[s])
                        sqp[s] = e;
                }
                if (xk[i] == sk[s])
                {
                    kas[s].push_back(i);
                    if (is_save_D)
                    {
                        kds[s].push_back(sqp[s]);
                        sqp[s] = 0;
                    }
                }
            }
        }
        // make sure A[k], D[k] are empty
        for (s = 0; s < sn; s++)
        {
            A[ni].insert(A[ni].end(), kas[s].begin(), kas[s].end());
            if (is_save_D)
                D[ni].insert(D[ni].end(), kds[s].begin(), kds[s].end());
        }
        kas.clear();
        a0 = A[ni];
        if (is_save_D)
        {
            kds.clear();
            sqp.clear();
            d0 = D[ni];
        }
    }

    void build(const std::string& vcfpanel, const std::string& samples, const std::string& region,
               double maf = 0.0002)
    {
        size_t i{0};
        int k{0}, m{0};
        BcfReader vcf(vcfpanel, samples, region);
        BcfRecord var(vcf.header);
        N = vcf.nsamples * 2;
        M = 0;
        vector<bool> gt;
        vector<vector<bool>> allgts, gt_rares, gt_commons;
        vector<int> snp_rares, snp_commons;
        double af;
        while (vcf.getNextVariant(var))
        {
            var.getGenotypes(gt);
            if (!var.isNoneMissing() || !var.allPhased())
                continue;
            // keep track of snp index with AF < minaf
            af = 0;
            for (auto g : gt)
                af += g;
            af /= N;
            if (af < maf)
            {
                snp_rares.push_back(M);
                gt_rares.push_back(gt);
            }
            else
            {
                snp_commons.push_back(M);
                gt_commons.push_back(gt);
            }
            if ((M + 1) % B == 0)
            {
                reorder.insert(reorder.end(), snp_rares.begin(), snp_rares.end());
                snp_rares.clear();
                reorder.insert(reorder.end(), snp_commons.begin(), snp_commons.end());
                snp_commons.clear();
                allgts.insert(allgts.cend(), gt_rares.begin(), gt_rares.end());
                gt_rares.clear();
                allgts.insert(allgts.cend(), gt_commons.begin(), gt_commons.end());
                gt_commons.clear();
            }
            M++;
        }
        reorder.insert(reorder.end(), snp_rares.begin(), snp_rares.end());
        snp_rares.clear();
        reorder.insert(reorder.end(), snp_commons.begin(), snp_commons.end());
        snp_commons.clear();
        allgts.insert(allgts.cend(), gt_rares.begin(), gt_rares.end());
        gt_rares.clear();
        allgts.insert(allgts.cend(), gt_commons.begin(), gt_commons.end());
        gt_commons.clear();
        G = (M + B - 1) / B;
        if (verbose)
            cerr << "N: " << N << ",M: " << M << ",G: " << G << ",B: " << B << endl;
        X.resize(G, GridVec(N));
        k = 0;
        for (m = 0; m < M; m++)
        {
            for (i = 0; i < N; i++)
                X[k][i] = (X[k][i] << 1) | (allgts[m][i] != 0);
            if ((m + 1) % B == 0)
            {
                for (i = 0; i < N; i++)
                    X[k][i] = reverseBits(X[k][i]); // reverset bits
                k++;                                // update next grid
            }
        }
        if (G == k + 1)
        {
            for (i = 0; i < N; i++)
            {
                X[k][i] <<= G * B - M;
                X[k][i] = reverseBits(X[k][i]); // reverset bits
            }
        }
        else if (G == k)
        {
            if (verbose)
                cerr << "no need padding\n";
        }
        else
        {
            throw std::runtime_error("something wrong\n");
        }

        // building A and save indices now
        A.resize(G);
        S.resize(G);
        C.resize(G);
        W.resize(G);
        if (is_save_D)
            D.resize(G);
        // create some intermediate varibales
        vector<int> Occ;
        vector<vector<int>> kas; // for building A
        vector<vector<int>> kds; // for building D
        vector<int> sqp;         // for building D
        GridSetU symbols;        // keep track of unique symbols at k
        GridVec y1(N);
        vector<int> a0(N); // initilize a[k]
        vector<int> d0;    // initilize d[k]
        int Gi, ki, ni{-1};
        for (i = 0; i < nindices; i++)
        {
            auto Gv = seq_by(i, G - 1, nindices);
            Gi = Gv.size();
            // initial a and d at ki=-1
            std::iota(a0.begin(), a0.end(), 0);
            if (is_save_D)
                d0.resize(N, 0);
            for (ki = 0; ki < Gi; ki++)
            {
                ni++;
                k = Gv[ki];
                build_indices_k(ni, ki, X[k], y1, S[ni], C[ni], W[ni], Occ, kas, kds, sqp, a0, d0,
                                symbols);
            }
        }
    }

    int save(const std::string& filename)
    {
        ofstream out(filename, ios::out | ios::binary);
        if (out.fail())
            return 1;
        out.write((char*)&B, sizeof(B));
        out.write((char*)&M, sizeof(M));
        out.write((char*)&N, sizeof(N));
        out.write((char*)&G, sizeof(G));
        out.write((char*)&nindices, sizeof(nindices));
        // write reorder (M)
        for (int m = 0; m < M; m++)
            out.write((char*)&reorder[m], sizeof(int));
        if (is_save_X)
        {
            // write X
            for (int k = 0; k < G; k++)
            {
                for (int n = 0; n < N; n++)
                    out.write((char*)&X[k][n], sizeof(T));
            }
        }
        // write S
        for (int k = 0; k < G; k++)
        {
            size_t sz = S[k].size();
            out.write((char*)&sz, sizeof(sz));
            for (size_t i = 0; i < sz; i++)
                out.write((char*)&S[k][i], sizeof(T));
        }
        // write A
        for (int k = 0; k < G; k++)
        {
            for (int n = 0; n < N; n++)
                out.write((char*)&A[k][n], sizeof(int));
        }
        if (is_save_D)
        {
            // write D
            for (int k = 0; k < G; k++)
            {
                for (int n = 0; n < N; n++)
                    out.write((char*)&D[k][n], sizeof(int));
            }
        }
        // write C
        for (int k = 0; k < G; k++)
        {
            size_t sz = C[k].size();
            out.write((char*)&sz, sizeof(sz));
            for (size_t i = 0; i < sz; i++)
                out.write((char*)&C[k][i], sizeof(int));
        }
        // write W
        for (int k = 0; k < G; k++)
        {
            size_t sz = W[k].size();
            out.write((char*)&sz, sizeof(sz));
            for (size_t i = 0; i < sz; i++)
            {
                size_t sz1 = W[k][i].size();
                out.write((char*)&sz1, sizeof(sz1));
                for (size_t j = 0; j < sz1; j++)
                    out.write((char*)&W[k][i][j], sizeof(int));
            }
        }
        return 0;
    }

    int load(const std::string& filename)
    {
        ifstream in(filename, ios::in | ios::binary);
        if (in.fail())
            return 1;
        if (!in.read((char*)&B, sizeof(B)))
            return 2;
        if (B != sizeof(T) * 8)
            throw invalid_argument("the binary file may be created with different B!\n");
        if (!in.read((char*)&M, sizeof(M)))
            return 2;
        if (!in.read((char*)&N, sizeof(N)))
            return 2;
        if (!in.read((char*)&G, sizeof(G)))
            return 2;
        if (!in.read((char*)&nindices, sizeof(nindices)))
            return 2;
        cerr << "N: " << N << ",M: " << M << ",G: " << G << ",B: " << B << ",nindices " << nindices
             << endl;
        // read reorder (M)
        reorder.resize(M);
        for (int m = 0; m < M; m++)
            in.read((char*)&reorder[m], sizeof(int));
        if (is_save_X)
        {
            // read X
            X.resize(G, GridVec(N));
            for (int k = 0; k < G; k++)
            {
                for (int n = 0; n < N; n++)
                {
                    if (!in.read((char*)&X[k][n], sizeof(T)))
                        return 2;
                }
            }
            if (verbose)
                cerr << "load X done" << endl;
        }
        // read S
        S.resize(G);
        for (int k = 0; k < G; k++)
        {
            size_t sz;
            if (!in.read((char*)&sz, sizeof(sz)) || sz < 1)
                return 2;
            S[k].resize(sz);
            for (size_t i = 0; i < sz; i++)
            {
                if (!in.read((char*)&S[k][i], sizeof(T)))
                    return 2;
            }
        }
        if (verbose)
            cerr << "load S done" << endl;
        // read A
        A.resize(G, vector<int>(N));
        for (int k = 0; k < G; k++)
        {
            for (int n = 0; n < N; n++)
            {
                if (!in.read((char*)&A[k][n], sizeof(int)))
                    return 2;
            }
        }
        if (verbose)
            cerr << "load A done" << endl;
        if (is_save_D)
        {
            // read D
            D.resize(G, vector<int>(N));
            for (int k = 0; k < G; k++)
            {
                for (int n = 0; n < N; n++)
                {
                    if (!in.read((char*)&D[k][n], sizeof(int)))
                        return 2;
                }
            }
            if (verbose)
                cerr << "load D done" << endl;
        }
        // read C
        C.resize(G);
        for (int k = 0; k < G; k++)
        {
            size_t sz;
            if (!in.read((char*)&sz, sizeof(sz)) || sz < 1)
                return 2;
            C[k].resize(sz);
            for (size_t i = 0; i < sz; i++)
            {
                if (!in.read((char*)&C[k][i], sizeof(int)))
                    return 2;
            }
        }
        if (verbose)
            cerr << "load C done" << endl;
        // read W
        W.resize(G);
        for (int k = 0; k < G; k++)
        {
            size_t sz;
            if (!in.read((char*)&sz, sizeof(sz)) || sz < 1)
                return 2;
            W[k].resize(sz);
            for (size_t i = 0; i < sz; i++)
            {
                size_t sz1;
                if (!in.read((char*)&sz1, sizeof(sz1)) || sz1 < 1)
                    return 2;
                W[k][i].resize(sz1);
                for (size_t j = 0; j < sz1; j++)
                {
                    if (!in.read((char*)&W[k][i][j], sizeof(int)))
                        return 2;
                }
            }
        }
        if (verbose)
            cerr << "load W done" << endl;

        return 0;
    }

    void query(const std::string& vcfquery, const std::string& samples, const std::string& region,
               int viewk = 0)
    {
        BcfReader vcf(vcfquery, samples, region);
        BcfRecord var(vcf.header);
        vector<bool> gt;
        vector<GridVec> Z(vcf.nsamples * 2, GridVec(G));
        size_t i{0};
        int k{0}, m{0};
        while (vcf.getNextVariant(var))
        {
            var.getGenotypes(gt);
            if (!var.isNoneMissing() || !var.allPhased())
                continue;
            for (i = 0; i < gt.size(); i++)
                Z[i][k] = (Z[i][k] << 1) | (gt[i] != 0);
            m++;
            if (m % B == 0)
            {
                for (i = 0; i < gt.size(); i++)
                    Z[i][k] = reverseBits(Z[i][k]); // reverset bits
                k++;                                // update next grid
            }
        }
        assert(m == M);
        if (G == k + 1)
        {
            for (i = 0; i < gt.size(); i++)
            {
                Z[i][k] <<= G * B - M;
                Z[i][k] = reverseBits(Z[i][k]); // reverset bits
            }
        }
        else if (G == k)
        {
            if (verbose)
                cerr << "no need padding\n";
        }
        else
        {
            throw std::runtime_error("something wrong\n");
        }
        cerr << "insert Z[1] now " << endl;
        view_zg(Z[1], viewk);
    }

    void report_v1(const GridVec& zg, IntMapU& haplens, IntMapU& hapends, int L = 5)
    {
        size_t i, s;
        int e{0}, f{0}, g{0}, n{0};
        for (int k = 0; k < G; k++)
        {
            auto kzs = std::lower_bound(S[k].begin(), S[k].end(), zg[k]);
            // s could be S[k].size() if kzs == S[k].end()
            s = std::fmin(std::distance(S[k].begin(), kzs), S[k].size() - 1);
            if (k > 0)
            {
                auto fzk = std::lower_bound(W[k][s].begin(), W[k][s].end(), f);
                auto gzk = std::lower_bound(W[k][s].begin(), W[k][s].end(), g);
                auto f1 =
                    C[k][s] + std::fmin(std::distance(W[k][s].begin(), fzk), W[k][s].size() - 1);
                auto g1 =
                    C[k][s] + std::fmin(std::distance(W[k][s].begin(), gzk), W[k][s].size() - 1);
                if (f1 >= g1 || k == G - 1)
                {
                    if (k == G - 1)
                    {
                        f = f1;
                        g = g1;
                    }
                    // report previous matches first
                    for (i = f; i <= g; i++)
                    {
                        n = A[k][i];
                        if (haplens.count(n) == 0)
                        {
                            haplens[n] = k - e;
                            hapends[n] = e;
                        }
                        else if (haplens[n] > k - e)
                        {
                            haplens[n] = k - e;
                            hapends[n] = e;
                        }
                    }
                }
                if (f1 < g1)
                {
                    f = f1;
                    g = g1;
                    continue;
                }
                // find new e, f, g
                i = 1;
                while (k >= i && X[k - i + 1][A[k][f1]] == zg[k - i + 1] &&
                       X[k - i][A[k][f1]] < zg[k - i])
                {
                    f1++;
                    if (X[k - i][A[k][f1]] == zg[k - i])
                        i++;
                }
                g = std::fmin(f1 + L, N - 1);
                f = f1;
                e = k - i;
                // TODO: use divergence D to do the follows efficiently
                // while (L--)
                // {
                //     i = k - e;
                //     while (i--)
                //     {
                //         if (X[k - i][A[k][g]] != zg[k - i])
                //             break;
                //     }
                //     if (X[k - i][A[k][g]] == zg[k - i])
                //         g++;
                // }
            }
            else
            {
                f = C[k][s];
                g = C[k][s] + W[k][s].size() - 1;
                e = 0;
            }
            // cerr << k << "," << f << "," << g << "," << e << endl;
        }
    }

    void report_setmaximal(IntMapU& haplens, IntMapU& hapends, const GridVec& zg, int L = 32,
                           int Step = 2)
    {
        int k, s, n, i, e, f, g, e1, f1, g1;
        for (k = 0; k < G; k++)
        {
            auto kzs = std::lower_bound(S[k].begin(), S[k].end(), zg[k]);
            // s could be S[k].size() if kzs == S[k].end()
            s = std::fmin(std::distance(S[k].begin(), kzs), S[k].size() - 1);
            if (k == 0)
            {
                f = C[k][s];
                g = C[k][s] + W[k][s].size() - 1;
                e = 0;
            }
            else
            {
                auto fzk = std::lower_bound(W[k][s].begin(), W[k][s].end(), f);
                auto gzk = std::lower_bound(W[k][s].begin(), W[k][s].end(), g);
                f1 = C[k][s] + std::fmin(std::distance(W[k][s].begin(), fzk), W[k][s].size() - 1);
                g1 = C[k][s] + std::fmin(std::distance(W[k][s].begin(), gzk), W[k][s].size() - 1);
                if (f1 < g1)
                {
                    e1 = e; // continue forwarding
                    f = f1;
                    g = g1;
                    e = e1;
                    if (k == G)
                    {
                        for (i = f; i < g; i++)
                        {
                            n = A[k][i];
                            if (haplens.count(n) == 0)
                            {
                                haplens[n] = k - e;
                                hapends[n] = e;
                            }
                            else if (k - e >= haplens[n])
                            {
                                haplens[n] = k - e;
                                hapends[n] = e;
                            }
                        }
                    }
                }
                else if (f1 == g1)
                {
                    // report matches from e to k for [f, g]
                    for (i = f; i < g; i++)
                    {
                        n = A[k][i];
                        if (haplens.count(n) == 0)
                        {
                            haplens[n] = k - e;
                            hapends[n] = e;
                        }
                        else if (k - e >= haplens[n])
                        {
                            haplens[n] = k - e;
                            hapends[n] = e;
                        }
                    }
                    // finding new e1, f1, g1
                    e1 = D[k][f1];
                    g1 = f1 + 1;
                    while (e1 >= 0 && zg[e1] == X[e1][A[k][f1]])
                        e1--;
                    while (g1 < N && D[k][g1] <= e1)
                        g1++;
                    g = g1;
                    f = f1;
                    e = e1;
                }
                else
                {
                    cerr << "something may be wrong with divergence matrix D\n";
                }
            }
        }
    }

    void report_neighourings(IntMapU& haplens, IntMapU& hapends, IntMapU& hapnindicies,
                             const GridVec& zg, int L = 32)
    {
        int k, s, klen, i, j, n, l, Gi, ki, idi, ni{-1};
        for (idi = 0; idi < nindices; idi++)
        {
            auto Gv = seq_by(idi, G - 1, nindices);
            Gi = Gv.size();
            int zak_prev, zak_curr, valid_grid_start = 0;
            bool first_valid_grid_start = true;
            for (ki = 0; ki < Gi; ki++)
            {
                k = Gv[ki];
                ni++;
                auto kzs = std::lower_bound(S[ni].begin(), S[ni].end(), zg[k]);
                s = std::fmin(std::distance(S[ni].begin(), kzs), S[ni].size() - 1);
                zak_curr = C[ni][s];

                if (zg[k] == S[ni][s])
                {
                    if (first_valid_grid_start)
                        valid_grid_start = ki;
                    first_valid_grid_start = false;
                }
                else
                {
                    if (verbose)
                        cerr << "skip: " << ki << endl;
                    // if zg[k] symbol not exists, skip this grid and start over.
                    first_valid_grid_start = true;
                    continue;
                }

                if (ki > valid_grid_start)
                {
                    if (zak_prev >= zak_curr)
                    {
                        auto kzi = std::lower_bound(W[ni][s].begin(), W[ni][s].end(), zak_prev);
                        zak_curr +=
                            std::fmin(std::distance(W[ni][s].begin(), kzi), W[ni][s].size() - 1);
                    }

                    // // // re-confirm new position with longest matches
                    // i = 1;
                    // while (ki >= i && X[Gv[ki - i + 1]][A[ni][zak_curr]] == zg[Gv[ki - i + 1]]
                    // &&
                    //        X[Gv[ki - i]][A[ni][zak_curr]] < zg[Gv[ki - i]])
                    // {
                    //     zak_curr++;
                    //     if (X[Gv[ki - i]][A[ni][zak_curr]] == zg[Gv[ki - i]])
                    //         i++;
                    // }


                    for (l = 0; l < L; l++)
                    {
                        n = A[ni][std::fmax(zak_curr - l - 1, 0)];
                        if (X[Gv[ki]][n] != zg[Gv[ki]])
                            break;
                        klen = 0;
                        for (j = ki; j >= 0; j--)
                        {
                            if (X[Gv[j]][n] == zg[Gv[j]])
                            {
                                klen++;
                            }
                            else
                            {
                                if (haplens.count(n) == 0)
                                {
                                    haplens[n] = klen;
                                    hapends[n] = ki;
                                }
                                else if (klen >= haplens[n])
                                {
                                    haplens[n] = klen;
                                    hapends[n] = ki;
                                }
                                if (hapnindicies.count(n))
                                    hapnindicies[n] =
                                        idi >= hapnindicies[n] ? (idi + 1) : hapnindicies[n];
                                break;
                            }
                        }
                        if (zak_curr - l - 1 == 0)
                            break;
                    }

                    for (l = 0; l < L; l++)
                    {
                        n = A[ni][std::fmin(zak_curr + l, N - 1)];
                        if (X[Gv[ki]][n] != zg[Gv[ki]])
                            break;
                        klen = 0;
                        for (j = ki; j >= 0; j--)
                        {
                            if (X[Gv[j]][n] == zg[Gv[j]])
                            {
                                klen++;
                            }
                            else
                            {
                                if (haplens.count(n) == 0)
                                {
                                    haplens[n] = klen;
                                    hapends[n] = ki;
                                }
                                else if (klen >= haplens[n])
                                {
                                    haplens[n] = klen;
                                    hapends[n] = ki;
                                }
                                if (hapnindicies.count(n))
                                    hapnindicies[n] =
                                        idi >= hapnindicies[n] ? (idi + 1) : hapnindicies[n];
                                break;
                            }
                        }
                        if (zak_curr + l == N - 1)
                            break;
                    }
                }
                zak_prev = zak_curr;
            }
        }
    }

    void report_neighourings_v1(IntMapU& haplens, IntMapU& hapends, const GridVec& zg, int L = 32,
                                int Step = 2)
    {
        int k, s, klen, klen_long, i, n, l, j;
        int zak_prev, zak_curr, valid_grid_start;
        bool first_valid_grid_start = true;
        for (k = 0; k < G; k++)
        {
            auto kzs = std::lower_bound(S[k].begin(), S[k].end(), zg[k]);
            s = std::fmin(std::distance(S[k].begin(), kzs), S[k].size() - 1);
            if (zg[k] == S[k][s])
            {
                if (first_valid_grid_start)
                    valid_grid_start = k;
                first_valid_grid_start = false;
                zak_curr = C[k][s];
            }
            else
            {
                if (verbose)
                    cerr << "skip: " << k << endl;
                // if zg[k] symbol not exists, skip this grid and start over.
                first_valid_grid_start = true;
                // hapi_curr = hapi_prev;
                continue;
            }
            if (k > valid_grid_start)
            {
                if (zak_prev >= zak_curr)
                {
                    auto kzi = std::lower_bound(W[k][s].begin(), W[k][s].end(), zak_prev);
                    zak_curr += std::fmin(std::distance(W[k][s].begin(), kzi), W[k][s].size() - 1);
                }
                i = 1;
                while (k >= i && X[k - i + 1][A[k][zak_curr]] == zg[k - i + 1] &&
                       X[k - i][A[k][zak_curr]] < zg[k - i])
                {
                    zak_curr++;
                    if (X[k - i][A[k][zak_curr]] == zg[k - i])
                        i++;
                }
                zak_prev = zak_curr;

                // scan up L haps
                for (l = 0; l < L; l++)
                {
                    n = A[k][std::fmax(zak_curr - l - 1, 0)];
                    klen = 0;
                    for (j = k; j > 0; j--)
                    {
                        if (X[j][n] == zg[j])
                        {
                            klen++;
                        }
                        else
                        {
                            if (haplens.count(n) == 0)
                            {
                                haplens[n] = klen;
                                hapends[n] = k;
                            }
                            else if (klen >= haplens[n])
                            {
                                haplens[n] = klen;
                                hapends[n] = k;
                            }
                            break;
                        }
                    }
                    if (zak_curr + l == N - 1)
                        break;
                }

                klen_long = i;
                // scan down L haps
                for (l = 0; l < L; l++)
                {
                    n = A[k][std::fmin(zak_curr + l, N - 1)];
                    klen = i;
                    // klen = 0;
                    for (j = k - i; j > 0; j--)
                    {
                        if (X[j][n] == zg[j])
                        {
                            klen++;
                        }
                        else
                        {
                            if (haplens.count(n) == 0)
                            {
                                haplens[n] = klen;
                                hapends[n] = k;
                            }
                            else if (klen >= haplens[n])
                            {
                                haplens[n] = klen;
                                hapends[n] = k;
                            }
                            break;
                        }
                    }
                    if (klen > klen_long)
                    {
                        zak_prev = zak_curr + l; // update
                        klen_long = klen;
                    }
                    if (zak_curr + l == N - 1)
                        break;
                }
            }
            else
            {
                zak_prev = zak_curr;
            }
        }
    }

    void report_neighourings_v2(IntMapU& haplens, IntMapU& hapends, const GridVec& zg, int L = 32,
                                int Step = 2)
    {
        int k, s, klen, i, n, l, j;
        int zak_prev, zak_curr, valid_grid_start;
        bool first_valid_grid_start = true;
        for (k = 0; k < G; k++)
        {
            auto kzs = std::lower_bound(S[k].begin(), S[k].end(), zg[k]);
            s = std::fmin(std::distance(S[k].begin(), kzs), S[k].size() - 1);
            if (zg[k] == S[k][s])
            {
                if (first_valid_grid_start)
                    valid_grid_start = k;
                first_valid_grid_start = false;
                zak_curr = C[k][s];
            }
            else
            {
                if (verbose)
                    cerr << "skip: " << k << endl;
                // if zg[k] symbol not exists, skip this grid and start over.
                first_valid_grid_start = true;
                // hapi_curr = hapi_prev;
                continue;
            }
            if (k > valid_grid_start)
            {
                if (zak_prev >= zak_curr)
                {
                    auto kzi = std::lower_bound(W[k][s].begin(), W[k][s].end(), zak_prev);
                    zak_curr += std::fmin(std::distance(W[k][s].begin(), kzi), W[k][s].size() - 1);
                }

                // // re-confirm new position with longest matches
                i = 1;
                while (k >= i && X[k - i + 1][A[k][zak_curr]] == zg[k - i + 1] &&
                       X[k - i][A[k][zak_curr]] < zg[k - i])
                {
                    zak_curr++;
                    if (X[k - i][A[k][zak_curr]] == zg[k - i])
                        i++;
                }

                // scan up L haps
                for (l = 0; l < L; l++)
                {
                    n = A[k][std::fmax(zak_curr - l - 1, 0)];
                    klen = 0;
                    for (j = k; j > 0; j--)
                    {
                        if (X[j][n] == zg[j])
                        {
                            klen++;
                        }
                        else
                        {
                            if (haplens.count(n) == 0)
                            {
                                haplens[n] = klen;
                                hapends[n] = k;
                            }
                            else if (klen >= haplens[n])
                            {
                                haplens[n] = klen;
                                hapends[n] = k;
                            }
                            break;
                        }
                    }
                    if (zak_curr - l - 1 == 0)
                        break;
                }

                // scan down L haps
                for (l = 0; l < L; l++)
                {
                    n = A[k][std::fmin(zak_curr + l, N - 1)];
                    klen = 0;
                    for (j = k; j > 0; j--)
                    {
                        if (X[j][n] == zg[j])
                        {
                            klen++;
                        }
                        else
                        {
                            if (haplens.count(n) == 0)
                            {
                                haplens[n] = klen;
                                hapends[n] = k;
                            }
                            else if (klen >= haplens[n])
                            {
                                haplens[n] = klen;
                                hapends[n] = k;
                            }
                            break;
                        }
                    }
                    if (zak_curr + l == N - 1)
                        break;
                }
            }
            zak_prev = zak_curr;
        }
    }

    vector<int> insert_and_match(const GridVec& zg)
    {
        vector<int> za(G);
        size_t i, s;
        int k;
        for (k = 0; k < G; k++)
        {
            auto kzs = std::lower_bound(S[k].begin(), S[k].end(), zg[k]);
            s = std::fmin(std::distance(S[k].begin(), kzs), S[k].size() - 1);
            za[k] = C[k][s];
            if (k > 0)
            {
                if (za[k - 1] >= za[k])
                {
                    auto kzi = std::lower_bound(W[k][s].begin(), W[k][s].end(), za[k - 1]);
                    za[k] += std::fmin(std::distance(W[k][s].begin(), kzi), W[k][s].size() - 1);
                }
                i = 1;
                while (k >= i && X[k - i + 1][A[k][za[k]]] == zg[k - i + 1] &&
                       X[k - i][A[k][za[k]]] < zg[k - i])
                {
                    za[k]++;
                    if (X[k - i][A[k][za[k]]] == zg[k - i])
                        i++;
                }
                cerr << k << "," << i << endl;
            }
        }
        return za;
    }

    vector<int> insert(const GridVec& zg, bool start_with_last = false)
    {
        vector<int> za(G);
        int k, s;
        for (k = 0; k < G; k++)
        {
            auto kzs = std::lower_bound(S[k].begin(), S[k].end(), zg[k]);
            s = std::fmin(std::distance(S[k].begin(), kzs), S[k].size() - 1);
            if (start_with_last && k == 0)
                za[k] = C[k][s] + W[k][s].size() - 1;
            else
                za[k] = C[k][s] + 1;
            if (k > 0)
            {
                if (za[k - 1] >= za[k])
                {
                    auto kzi = std::lower_bound(W[k][s].begin(), W[k][s].end(), za[k - 1]);
                    za[k] += std::fmin(std::distance(W[k][s].begin(), kzi), W[k][s].size() - 1);
                }
            }
        }
        return za;
    }

    void view_panel(int k, bool bit = true)
    {
        for (size_t i = 0; i < X[0].size(); i++)
        {
            for (int j = 0; j <= k + 3; j++)
            {
                if (bit)
                {
                    auto rb = reverseBits(X[j][A[k][i]]);
                    cout << std::bitset<sizeof(T) * 8>(rb) << " ";
                }
                else
                {
                    cout << X[j][A[k][i]] << " ";
                }
            }
            cout << endl;
        }
        for (size_t i = 0; i < D[0].size(); i++)
        {
            for (size_t j = 0; j < D.size(); j++)
            {
                cerr << D[j][i] << " ";
            }
            cerr << endl;
        }
    }

    void view_zg(const GridVec& zg, int k, bool bit = true, int L = 0)
    {
        auto za1 = insert(zg, true);
        auto za2 = insert_and_match(zg);
        for (size_t i = 0; i < X[0].size(); i++)
        {
            if (((i < za1[k] - L) || (i > za1[k] + L)) && L != 0)
                continue;
            for (int j = 0; j <= k + 2; j++)
            {
                if (bit)
                {
                    auto rb = reverseBits(X[j][A[k][i]]);
                    cout << std::bitset<sizeof(T) * 8>(rb) << " ";
                }
                else
                {
                    cout << X[j][A[k][i]] << " ";
                }
            }
            cout << endl;
            if (i == za1[k])
            {
                // print out original Z bits
                cout << "========= fisrt zg is inserting here ========  k=" << k
                     << ", za1[k]=" << za1[k] << endl;
                for (int j = 0; j <= k + 2; j++)
                {
                    if (bit)
                    {
                        auto rb = reverseBits(zg[j]);
                        cout << std::bitset<sizeof(T) * 8>(rb) << " ";
                    }
                    else
                    {
                        cout << zg[j] << " ";
                    }
                }
                cout << endl;
                cout << "========= fisrt zg is inserting here ========  k=" << k
                     << ", za1[k]=" << za1[k] << endl;
            }
            if (i == za2[k])
            {
                // print out original Z bits
                cout << "========= last zg is inserting here ========  k=" << k
                     << ", za2[k]=" << za2[k] << endl;
                for (int j = 0; j <= k + 2; j++)
                {
                    if (bit)
                    {
                        auto rb = reverseBits(zg[j]);
                        cout << std::bitset<sizeof(T) * 8>(rb) << " ";
                    }
                    else
                    {
                        cout << zg[j] << " ";
                    }
                }
                cout << endl;
                cout << "========= last zg is inserting here ========  k=" << k
                     << ", za2[k]=" << za2[k] << endl;
            }
        }
    }
};


#endif // MSPBWT_H_