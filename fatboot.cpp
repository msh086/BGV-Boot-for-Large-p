/* Copyright (C) 2012-2020 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <NTL/BasicThreadPool.h>
#include <helib/FHE.h>
#include <helib/EncryptedArray.h>
#include <helib/matmul.h>
#include <helib/debugging.h>
#include <helib/fhe_stats.h>
#include <helib/sample.h>
#include <helib/ArgMap.h>

// #include "gtest/gtest.h"
// #include "test_common.h"
struct Parameters
{
    const long p; // plaintext base
    const long r; // exponent
    // p^r is the plaintext-space modulus
    const long c;        // number of columns in the key-switching matrices
    const long bits;     // # of bits in the modulus chain
    const long skHwt;    // Hamming weight of recryption secret key
    const long nthreads; // number of threads
    const long seed;     // random number seed
    const int useCache;  // 0: zzX cache, 1: DCRT cache
    const int force_bsgs;
    const int force_hoist;
    // int disable_intFactor // fhe_disable_intFactor (disabled by Victor)
    const int chen_han;
    const bool debug; // generate debugging output
    const int scale;  // scale parameter
    const NTL::Vec<long> global_gens;
    const NTL::Vec<long> global_ords;
    const NTL::Vec<long> global_mvec;
    const int c_m; // = 100;
    const long outer_rep;
    long inner_rep;

    Parameters(long p,
               long r,
               long c,
               long bits,
               long skHwt,
               long nthreads,
               long seed,
               long useCache,
               int c_m,
               int force_bsgs,
               int force_hoist,
               int chen_han,
               bool debug,
               int scale,
               const std::vector<long> &global_gens,
               const std::vector<long> &global_ords,
               const std::vector<long> &global_mvec,
               long outer_rep,
               long inner_rep) : p(p),
                                 r(r),
                                 c(c),
                                 bits(bits),
                                 skHwt(skHwt),
                                 nthreads(nthreads),
                                 seed(seed),
                                 useCache(useCache),
                                 force_bsgs(force_bsgs),
                                 force_hoist(force_hoist),
                                 chen_han(chen_han),
                                 debug(debug),
                                 scale(scale),
                                 global_gens(helib::convert<NTL::Vec<long>>(global_gens)),
                                 global_ords(helib::convert<NTL::Vec<long>>(global_ords)),
                                 global_mvec(helib::convert<NTL::Vec<long>>(global_mvec)),
                                 c_m(c_m),
                                 outer_rep(outer_rep),
                                 inner_rep(inner_rep)
    {
        if (global_gens.empty() || global_ords.empty() || global_mvec.empty())
            throw helib::LogicError("gens, ords, and mvec must be non-empty");
    };

    friend std::ostream &operator<<(std::ostream &os, const Parameters &params)
    {
        return os << "{"
                  << "p" << params.p << ","
                  << "r" << params.r << ","
                  << "c" << params.c << ","
                  << "bits" << params.bits << ","
                  << "skHwt" << params.skHwt << ","
                  << "nthreads" << params.nthreads << ","
                  << "seed" << params.seed << ","
                  << "useCache" << params.useCache << ","
                  << "force_bsgs" << params.force_bsgs << ","
                  << "force_hoist" << params.force_hoist << ","
                  << "chen_han" << params.chen_han << ","
                  << "debug" << params.debug << ","
                  << "scale" << params.scale << ","
                  << "global_gens" << params.global_gens << ","
                  << "global_ords" << params.global_ords << ","
                  << "global_mvec" << params.global_mvec << ","
                  << "c_m" << params.c_m << "}";
    }
};

class GTestFatboot
{
    void preContextSetup()
    {

        if (!noPrint)
            helib::fhe_stats = true;

        if (!noPrint)
        {
            std::cout << "*** GTestFatboot";
            if (helib::isDryRun())
                std::cout << " (dry run)";
            std::cout << ": p=" << p << ", r=" << r << ", bits=" << bits
                      << ", t=" << skHwt << ", c=" << c << ", m=" << m
                      << " mvec=" << mvec << ", gens=" << helib::vecToStr(gens)
                      << ", ords=" << helib::vecToStr(ords) << std::endl;
            std::cout << "Computing key-independent tables..." << std::flush;
        }
        helib::setTimersOn();
        helib::setDryRun(
            false); // Need to get a "real context" to test bootstrapping
        time = -NTL::GetTime();
    }

    static void setGlobals(int force_bsgs, int force_hoist, int chen_han)
    {
        helib::fhe_test_force_bsgs = force_bsgs;
        helib::fhe_test_force_hoist = force_hoist;
        helib::fhe_force_chen_han = chen_han;
    };

    void cleanupBootstrappingGlobals()
    {
        helib::fhe_test_force_bsgs = old_fhe_test_force_bsgs;
        helib::fhe_test_force_hoist = old_fhe_test_force_hoist;
        helib::fhe_force_chen_han = old_fhe_force_chen_han;
    }

    static void setSeedIfNeeded(const long seed)
    {
        if (seed)
            SetSeed(NTL::ZZ(seed));
    };

    static void checkPM(const long p, const long m)
    {
        helib::assertTrue(NTL::GCD(p, m) == 1, "GCD(p, m) == 1");
    }

    const int old_fhe_test_force_bsgs;
    const int old_fhe_test_force_hoist;
    const int old_fhe_force_chen_han;
    const long p;
    const long r;
    const long c;
    const long bits;
    const long skHwt;
    const long nthreads;
    const long seed;
    const long useCache;
    const int force_bsgs;
    const int force_hoist;
    const int chen_han;
    const bool debug;
    const int scale;
    const NTL::Vec<long> mvec;
    const std::vector<long> gens;
    const std::vector<long> ords;
    const int c_m; // = 100;
    const long outer_rep;
    const long inner_rep;
    const bool noPrint;
    const bool dry;
    // new params
    const long encapSkHwt;
    const long btsC;
    const long t;
    const double btsScale;
    const bool newBtsFlag;
    const bool isThick;

    const long m, phim;
    double time;
    helib::Context context;

public:
    GTestFatboot(Parameters param,
                 bool noPrint,
                 bool dry,
                 long encapSkHwt,
                 long btsC,
                 long t,
                 long btsScale,
                 long newBtsFlag,
                 bool newKSFlag,
                 bool isThick) : old_fhe_test_force_bsgs(helib::fhe_test_force_bsgs),
                                 old_fhe_test_force_hoist(helib::fhe_test_force_hoist),
                                 old_fhe_force_chen_han(helib::fhe_force_chen_han),
                                 p((setGlobals(param.force_bsgs,
                                               param.force_hoist,
                                               param.chen_han),
                                    param.p)),
                                 r(param.r),
                                 c(param.c),
                                 bits(param.bits),
                                 skHwt(param.skHwt),
                                 nthreads((NTL::SetNumThreads(param.nthreads), param.nthreads)),
                                 seed((setSeedIfNeeded(param.seed), param.seed)),
                                 useCache(param.useCache),
                                 force_bsgs(param.force_bsgs),
                                 force_hoist(param.force_hoist),
                                 chen_han(param.chen_han),
                                 debug(param.debug),
                                 scale(param.scale),
                                 mvec(param.global_mvec),
                                 gens(helib::convert<std::vector<long>>(param.global_gens)),
                                 ords(helib::convert<std::vector<long>>(param.global_ords)),
                                 c_m(param.c_m),
                                 outer_rep(param.outer_rep),
                                 inner_rep(param.inner_rep),
                                 noPrint(noPrint),
                                 dry(dry),
                                 encapSkHwt(encapSkHwt),
                                 btsC(btsC),
                                 t(t),
                                 btsScale(btsScale),
                                 newBtsFlag(newBtsFlag),
                                 isThick(isThick),
                                 m(helib::computeProd(mvec)),
                                 phim((checkPM(p, m), helib::phi_N(m))),
                                 time(0),
                                 context((preContextSetup(),
                                          helib::ContextBuilder<helib::BGV>()
                                              .m(m)
                                              .p(p)
                                              .r(r)
                                              .gens(gens)
                                              .ords(ords)
                                              .scale(scale ? scale : 10 /*default is 10.*/)
                                              // new params
                                              .encapSkHwt(encapSkHwt)
                                              .btsScale(btsScale)
                                              .btsC(btsC)
                                              .t(t)
                                              .newKS(newKSFlag)
                                              .newBts(newBtsFlag)
                                              // remove the hack
                                              .bits(bits)
                                              .c(c)
                                              .bootstrappable(true)
                                              .skHwt(skHwt)
                                              .mvec(mvec)
                                              .buildCache(useCache)
                                              .setThick(isThick)
                                              //    .buildModChain(false)
                                              .build()))
    {
    }

    void TearDown()
    {
        if (noPrint)
        {
            helib::printAllTimers();
        }
        if (helib::fhe_stats)
            helib::print_stats(std::cout);

        cleanupBootstrappingGlobals();
        helib::cleanupDebugGlobals();
    }

    void run()
    {
        // context.buildModChain(bits,
        //                       c,
        //                       /*willBeBootstrappable=*/true,
        //                       /*t=*/skHwt);

        if (!noPrint)
        {
            std::cout << "security=" << context.securityLevel() << std::endl;
            std::cout << "# small primes = " << context.getSmallPrimes().card()
                      << std::endl;
            std::cout << "# ctxt primes = " << context.getCtxtPrimes().card()
                      << std::endl;
            std::cout << "# bits in ctxt primes = "
                      << long(context.logOfProduct(context.getCtxtPrimes()) / log(2.0) +
                              0.5)
                      << std::endl;
            std::cout << "# special primes = " << context.getSpecialPrimes().card()
                      << std::endl;
            std::cout << "# bits in special primes = "
                      << long(context.logOfProduct(context.getSpecialPrimes()) /
                                  log(2.0) +
                              0.5)
                      << std::endl;
            std::cout << "scale=" << context.getScale() << std::endl;
        }

        // context.enableBootStrapping(mvec, useCache);
        time += NTL::GetTime();

        if (!noPrint)
        {
            std::cout << " done in " << time << " seconds" << std::endl;
            std::cout << "  e=" << context.getRcData().e
                      << ", e'=" << context.getRcData().ePrime
                      << ", t=" << context.getRcData().skHwt << "\n"
                      << "  ";
            context.printout();
        }
        helib::setDryRun(
            dry); // Now we can set the dry-run flag if desired

        long p2r = context.getAlMod().getPPowR();
        helib::BootBench meanBenchmarker;

        for (long numkey = 0; numkey < outer_rep; numkey++)
        { // test with 3 keys
            time = -NTL::GetTime();
            if (!noPrint)
                std::cout << "Generating keys, " << std::flush;
            helib::SecKey secretKey(context);
            helib::PubKey &publicKey = secretKey;
            secretKey.GenSecKey(); // A +-1/0 secret key
            helib::addSome1DMatrices(
                secretKey); // compute key-switching matrices that we need
            helib::addFrbMatrices(secretKey);
            if (!noPrint)
                std::cout << "computing key-dependent tables..." << std::flush;
            secretKey.genRecryptData();
            time += NTL::GetTime();
            if (!noPrint)
                std::cout << " done in " << time << " seconds\n";

            NTL::ZZX ptxt_poly1, ptxt_poly;
            std::vector<NTL::ZZX> val1;
            if (isThick)
            {
                NTL::zz_p::init(p2r);
                NTL::zz_pX poly_p = NTL::random_zz_pX(context.getPhiM());
                helib::zzX poly_p1 = helib::balanced_zzX(poly_p);
                ptxt_poly = helib::convert<NTL::ZZX>(poly_p1);
                helib::PolyRed(ptxt_poly1, ptxt_poly, p2r, true);
            }
            else
            {
                long d = context.getOrdP();
                long phim = context.getPhiM();
                long nslots = phim / d;

                NTL::zz_p::init(p2r);
                NTL::Vec<NTL::zz_p> val0(NTL::INIT_SIZE, nslots);

                bool debug_ordered = false;
                if (debug_ordered)
                    for (long i = 0; i < nslots; i++)
                        val0[i] = i + 1;
                else
                    for (auto &x : val0)
                        random(x);

                val1.resize(nslots);
                for (long i = 0; i < nslots; i++)
                    val1[i] = NTL::conv<NTL::ZZX>(NTL::conv<NTL::ZZ>(rep(val0[i])));
            }
            // this is the format produced by decryption

            helib::setupDebugGlobals(&secretKey, context.shareEA());

            NTL::ZZX poly2;
            std::vector<NTL::ZZX> val2;
            helib::Ctxt c1(publicKey);

            if (isThick)
                secretKey.Encrypt(c1, ptxt_poly, p2r);
            else
                context.shareEA()->encrypt(c1, publicKey, val1);

            helib::resetAllTimers();
            for (long num = 0; num < inner_rep; num++)
            { // multiple tests with same key
                if (isThick)
                {
                    auto benchmarker = publicKey.reCrypt(c1);
                    // print the current results
                    printf("time for linear1 = %f, linear2 = %f, extract = %f, total = %f\n",
                            benchmarker.time_linear_1, benchmarker.time_linear_2,
                            benchmarker.time_extract, benchmarker.time_total);
                    printf("bits for linear1 = %f, linear2 = %f, extract = %f, min cap = %f, after cap = %f, after prod = %f\n",
                            benchmarker.bits_down_linear_1, benchmarker.bits_down_linear_2,
                            benchmarker.bits_down_extract, context.getMinCap(), benchmarker.bits_final, benchmarker.bits_after_inner_prod);
                    printf("cap / time = %f\n",
                           (benchmarker.bits_final - context.getMinCap()) / benchmarker.time_total);
                    // add to the mean benchmarker
                    meanBenchmarker += benchmarker;
                    secretKey.Decrypt(poly2, c1);
                    for (long i = 0; i < context.getPhiM(); i++)
                    {
                        if (ptxt_poly1[i] != poly2[i])
                            std::cout << i << " is bad\n";
                    }
                    helib::assertEq(ptxt_poly1, poly2, "fat results not match");
                }
                else
                {
                    helib::BootBench benchmarker;
                    if (t >= 0 && newBtsFlag)
                        benchmarker = publicKey.thinReCryptRefine(c1);
                    else
                        benchmarker = publicKey.thinReCrypt(c1);
                    printf("time for linear1 = %f, linear2 = %f, extract = %f, total = %f\n",
                            benchmarker.time_linear_1, benchmarker.time_linear_2,
                            benchmarker.time_extract, benchmarker.time_total);
                    printf("bits for linear1 = %f, linear2 = %f, extract = %f, min cap = %f, after cap = %f, after prod = %f\n",
                            benchmarker.bits_down_linear_1, benchmarker.bits_down_linear_2,
                            benchmarker.bits_down_extract, context.getMinCap(), benchmarker.bits_final, benchmarker.bits_after_inner_prod);
                    printf("cap / time = %f\n",
                           (benchmarker.bits_final - context.getMinCap() - benchmarker.bits_down_linear_1) / benchmarker.time_total);
                    meanBenchmarker += benchmarker;
                    context.shareEA()->decrypt(c1, secretKey, val2);
                    for (long i = 0; i < val1.size(); i++)
                    {
                        if (val1[i] != val2[i])
                            std::cout << i << " is bad\n";
                    }
                    helib::assertEq(val1, val2, "thin results not match");
                }
            }
        }
        meanBenchmarker.Mult(1.0 / (inner_rep * outer_rep));            
        printf("\n### summary over %ld test\n", inner_rep * outer_rep);        
        printf("time for linear1 = %f, linear2 = %f, extract = %f, total = %f\n",
                meanBenchmarker.time_linear_1, meanBenchmarker.time_linear_2,
                meanBenchmarker.time_extract, meanBenchmarker.time_total);
        printf("bits for linear1 = %f, linear2 = %f, extract = %f, min cap = %f, after cap = %f, after prod = %f\n",
                meanBenchmarker.bits_down_linear_1, meanBenchmarker.bits_down_linear_2,
                meanBenchmarker.bits_down_extract, context.getMinCap(), meanBenchmarker.bits_final, meanBenchmarker.bits_after_inner_prod);
        if (isThick)
            printf("cap / time = %f\n",
                           (meanBenchmarker.bits_final - context.getMinCap()) / meanBenchmarker.time_total);
        else
            printf("cap / time = %f\n",
                (meanBenchmarker.bits_final - context.getMinCap() - meanBenchmarker.bits_down_linear_1) / meanBenchmarker.time_total);

        printf("\n\n### bts finished, everything ok ###\n\n");
        if (!noPrint)
            helib::printAllTimers();
#if (defined(__unix__) || defined(__unix) || defined(unix))
        struct rusage rusage;
        getrusage(RUSAGE_SELF, &rusage);
        if (!noPrint)
            std::cout << "  rusage.ru_maxrss=" << rusage.ru_maxrss << std::endl;
#endif
        if (helib::fhe_stats)
            helib::print_stats(std::cout);
    }
};

int main(int argc, char *argv[])
{
    helib::ArgMap amap;

    long i_arg = 0;
    amap.arg("i", i_arg, "index of the chosen parameter set");

    long h_arg = 0;
    amap.arg("h", h_arg, "hwt of encapsulated key");

    long t_arg = 0;
    amap.arg("t", t_arg, "parameter t used in new bts");

    long newBts_flag = 0;
    amap.arg("newbts", newBts_flag, "if new bts is used");

    long newKS_flag = 1;
    amap.arg("newks", newKS_flag, "if new ks is used");

    long thick_flag = 0;
    amap.arg("thick", thick_flag, "if thick bts is used");

    long repeat_arg = 5;
    amap.arg("repeat", repeat_arg, "number of tests");

    amap.parse(argc, argv);
    // some candidate primes:
    // 2^14 - 3
    // 2^13 - 1
    // p,r,c, bits,skHwt,nthreads,seed,useCache,c_m,force_bsgs,force_hoist,chen_han,debug,scale,global_gens,global_ords,global_mvec,outer_rep, inner_rep
    bool force_chen_han = false;
    // clang-format off
    Parameters toy_params[] = {
        // NOTE: toy param sets, only for debug
        Parameters(
                // p, r, c, bits, h
                    17, 4, 3, 800, 0, 
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {262, 146, 181}, {4, 6, 7}, {5, 9, 29}, 1, 5), // d = 4
        Parameters(
                // p, r, c, bits, h
                127, 2, 3, 800, 0, 
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {686, 838}, {6, 17}, {9, 137}, 1, 5), // d = 4
        Parameters(
                // p, r, c, bits, h
                65537, 1, 3, 800, 0, 
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {923, 475}, {2, 23}, {3, 461}, 1, 5), // d = 3
    };
    Parameters real_params[] = {
        // NOTE: parameters used in GIKV22
        // TODO: update bits
        // NOTE: we don't have the code of GIKV23... 
        //  if we use encapsulated sparse key, this is no way to compare with them
        // XXX: always set the mvec entries as primes, otherwise bound (5) may fail
        // Parameters( // d = 40, C = 149, good dim, m = 45551
        //         // log(qks*R) ~= 72, encapHwt = 22(12), seclvl = 129.1(82.3)
        //         // main seclvl = 78, total bits = (1593, 84.01793423200083)
        //         // p, r, c, bits, h
        //         17, 4, 3, 1100, 0, 
        //         // nthreads,seed,useCache,c_m,
        //             1, 0, 1, 100, 
        //         // force_bsgs,force_hoist,chen_han,debug,scale
        //             0, 0, force_chen_han, 0, 0, 
        //         // global_gens,global_ords,global_mvec,outer_rep,inner_rep
        //         {19394, 37270}, {100, 10}, {101, 11, 41}, 1, 5),
        // Parameters( // d = 14, C = 161, bad dim, m = 32551
        //         // log(qks*R) ~= 66, encapHwt = 24(10), seclvl = 134.7(70.9)
        //         // main seclvl = 66, total bits = (1593, 69.76804235959568)
        //         // p, r, c, bits, h
        //         127, 2, 3, 1100, 0, 
        //         // nthreads,seed,useCache,c_m,
        //             1, 0, 1, 100, 
        //         // force_bsgs,force_hoist,chen_han,debug,scale
        //             0, 0, force_chen_han, 0, 0, 
        //         // global_gens,global_ords,global_mvec,outer_rep,inner_rep
        //         {7571, 28768}, {42, -54}, {43, 757}, 1, 5),
        // NOTE: new params for p = 17 and 127
/*A*/   Parameters( // d = 24, C = 105, good dim, m = 38309
                // t = 2, log(qks*R) ~= 71, encapHwt = 14(24), seclvl = 89.8(133.1)
                // native, log(qks*R) ~= 65(64), encapHwt = 12(24), seclvl = 81.3(136.0)
                // main seclvl = 82.5, total bits <= 1500
                // p, r, c, bits, h
                17, 4, 3, 1070, 0, // 1462 bits
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {26421, 23781}, {28, 55}, {29, 1321}, 1, 5),
/*B*/   Parameters( // d = 45, C = 113, good dim, m = 56647
                // t = 1, log(qks*R) ~= 67, encapHwt = 12(22), seclvl = 85.7(135.4)
                // native, log(qks*R) ~= 61(70), encapHwt = 12(22), seclvl = 87.3(134.1)
                // TODO: larger t's
                // main seclvl = 82.6, total bits <= 2250
                // p, r, c, bits, h
                127, 2, 3, 1630, 0, // 2253 bits
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {12249, 49026}, {36, 34}, {37, 1531}, 1, 5),
        // NOTE: parameters with larger p
/*C*/   Parameters( // d = 28, C = 159, bad dim, m = 55427
                // t = 1, log(qks*R) ~= 71, encapHwt = 12(22), seclvl = 85.5(133.4)
                // native, log(qks*R) ~= 65(65), encapHwt = 12(22), seclvl = 87.4(135.5)
                // main seclvl = 82.8, total bits = 2200
                // p, r, c, bits, h
                257, 2, 3, 1610, 0,  // 2176 bits
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {52850, 2581}, {42, -46}, {43, 1289}, 1, 5),
/*D*/   Parameters( // d = 14, C = 129, good dim, m = 45193
                // t = 1, log(qks*R) ~= 73, encapHwt = 12(24), seclvl = 81.9(136.7)
                // t = -1, log(qks*R) ~= 63, encapHwt = 12(22), seclvl = 83.8(131.8)
                // native, log(qks*R) ~= 67(67), encapHwt = 12(22), seclvl = 83.3(130.2)
                // main seclvl = 82.3, total bits = 1800
                // p, r, c, bits, h
                8191, 1, 3, 1320, 0, // 1803 bits for t = 1, 1787 bits for t = -1
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {26276, 32595}, {42, 75}, {43, 1051}, 1, 5),
/*E*/   Parameters( // d = 18, C = 141, good dim, m = 50731
                // t = 1, log(qks*R) ~= 82, encapHwt = 12(24), seclvl = 82.3(136.8)
                // t = -1, log(qks*R) ~= 68, encapHwt = 12(22), seclvl = 84.8(132.8)
                // native, log(qks*R) ~= 76(75), encapHwt = 12(22), seclvl = 83.0(130.6)
                // main seclvl = 82.3, total bits = 2050
                // p, r, c, bits, h
                65537, 1, 3, 1500, 0, // 2036 bits for t = 1, 2018 for t = -1
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {48117, 5239}, {96, 29}, {97, 523}, 1, 5),

        };
    // clang-format on
    bool noPrint = false;
    bool dry = false; // don't set to true...
    long encapSkHwt = h_arg;
    long btsC = 3;
    long t = t_arg; // -1 for non-power-of-p aux, 0 for auto-decided power-of-p aux, >0 for aux = p^t
    double btsScale = 8;
    bool newBtsFlag = newBts_flag;
    bool newKSFlag = newKS_flag;
    bool isThick = thick_flag;
    // for ordinary key with 128 bit security, use h = 492 (no, achieving 128 is too expensive)
    // use h = 120 for the ordinary key
    std::cout << "chosen idx = " << i_arg << "\n";
    std::cout << "encapHwt = " << encapSkHwt << "\n";
    std::cout << "t for new bts is " << t << "\n";
    std::cout << "new bts ? " << newBtsFlag << "\n";
    std::cout << "new ks ? " << newKSFlag << "\n";
    std::cout << "thick bts ?" << isThick << "\n";
    std::cout << "repeat times = " << repeat_arg << "\n";
    real_params[i_arg].inner_rep = repeat_arg;
    GTestFatboot test(real_params[i_arg], noPrint, dry, encapSkHwt, btsC, t, btsScale, newBtsFlag, newKSFlag, isThick);
    test.run();
    test.TearDown();
    return 0;
}

// LEGACY TEST DEFAULT PARAMETERS:
// long p=2;
// long r=1;
// long c=3;
// long bits=600;
// long t=64;
// long nthreads=1;
// long seed=0;
// long useCache=1;
// int c_m = 100;
