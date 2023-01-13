#include "mspbwt/mspbwt.h"
#include "utils/timer.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
void mspbwt_build(const std::string& binfile, const std::string& vcfpanel,
                  const std::string& samples, const std::string& region, double maf = 0)
{
    msPBWT<uint64_t> msp;
    msp.build(vcfpanel, samples, region, maf);
    msp.save(binfile);
}

//' @export
// [[Rcpp::export]]
SEXP mspbwt_load(const std::string& binfile)
{
    /* creating a pointer to msPBWT */
    msPBWT<uint64_t>* msp = new msPBWT<uint64_t>();
    msp->load(binfile);

    /* wrap the pointer as an external pointer */
    /* this automatically protected the external pointer from R garbage
     collection until xp goes out of scope. */
    XPtr<msPBWT<uint64_t>> xp(msp, true);

    /* return it back to R, since xp goes out of scope after the return
     the external pointer is no more protected by xp, but it gets
     protected by being on the R side */
    return (xp);
}


//' @export
// [[Rcpp::export]]
List mspbwt_report(SEXP xp_, const IntegerVector& z, int pbwtL, int pbwtS)
{
    Timer tm;
    tm.clock();

    Rcpp::XPtr<msPBWT<uint64_t>> xp(xp_);
    vector<int> zc = as<vector<int>>(z);
    auto zg = xp->encodezg(zc);
    IntMapU haplens, hapstarts;
    xp->report_neighourings(haplens, hapstarts, zg, pbwtL, pbwtS);
    int n = haplens.size();
    vector<int> haps(n), lens(n), starts(n);
    n = 0;
    for (auto const& h : haplens)
    {
        haps[n] = h.first;
        lens[n] = h.second;
        starts[n] = hapstarts[h.first];
        n++;
    }

    Rcout << "elapsed time of mspbwt insert: " << tm.abstime() << " milliseconds" << endl;
    List out = List::create(Named("haps", haps), Named("lens", lens), Named("ends", starts));

    return out;
}
