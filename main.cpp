//
// run a full protocol
//
#include <iostream>
#include <chrono>
#include <fstream>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

#include "properties.h"
#include "utils/Logger.h"
#include "utils/aux.h"

//helib
#include <helib/binaryArith.h>
#include <helib/binaryCompare.h>


int main() {
    /*      Bootstrapping System    */
    //  logger
    auto t0_main = std::chrono::high_resolution_clock::now();
    Logger logger;
    logger.log(log_trace, "Starting Protocol");
    //  Keys Server
    KeysServer keysServer;
//    helib::Context context = keysServer.getContext();
    const helib::EncryptedArray &ea = keysServer.getEA();

    // Choose two random n-bit integers
    long pa = NTL::RandomBits_long(bitSize);
    long pb = NTL::RandomBits_long(bitSize + 1);
    long pMax = std::max(pa, pb);
    long pMin = std::min(pa, pb);
    bool pMu = pa > pb;
    bool pNi = pa < pb;

    // Encrypt the individual bits
    NTL::Vec<helib::Ctxt> eMax, eMin, enca, encb;

    const helib::SecKey& encryptionKey = keysServer.getSecKey();
    helib::Ctxt mu(encryptionKey), ni(encryptionKey);
    resize(enca, long(bitSize), mu);
    resize(encb, long(bitSize) + 1, ni);
    for (long i = 0; i <= bitSize; i++) {
        if (i < bitSize)
            encryptionKey.Encrypt(enca[i], NTL::ZZX((pa >> i) & 1));
        encryptionKey.Encrypt(encb[i], NTL::ZZX((pb >> i) & 1));
        if (helib_bootstrap) { // put them at a lower level
            if (i < bitSize)
                enca[i].bringToSet(keysServer.getCtxtPrimes(5));
            encb[i].bringToSet(keysServer.getCtxtPrimes(5));
        }
    }
#ifdef HELIB_DEBUG
    decryptAndPrint((std::cout << " before comparison: "),
                  encb[0],
                  secKey,
                  ea,
                  0);
#endif

    std::vector<long> slotsMin, slotsMax, slotsMu, slotsNi;
    // comparison only
    compareTwoNumbers(mu,
                      ni,
                      helib::CtPtrs_VecCt(enca),
                      helib::CtPtrs_VecCt(encb),
                      false);
    ea.decrypt(mu, encryptionKey, slotsMu);
    ea.decrypt(ni, encryptionKey, slotsNi);
    auto muniPair = std::make_pair(slotsMu[0], slotsNi[0]);
    auto pmupniPair = std::make_pair((long) pMu, (long) pNi);
    std::cout << "Comparison (without min max) error: a=" << pa << ", b=" << pb
         << ", mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << std::endl;
    std::cout << std::get<0>(muniPair) << ", " << std::get<1>(muniPair) << std::endl;
    std::cout << std::get<0>(pmupniPair) << ", " << std::get<1>(pmupniPair) << std::endl;
    if (VERBOSE) {
        std::cout << "Comparison (without min max) succeeded: ";
        std::cout << '(' << pa << ',' << pb << ")=> mu=" << slotsMu[0]
                  << ", ni=" << slotsNi[0] << std::endl;
    }

    {
        helib::CtPtrs_VecCt wMin(eMin),
                wMax(eMax); // A wrappers around output vectors
        // comparison with max and min
        compareTwoNumbers(wMax,
                          wMin,
                          mu,
                          ni,
                          helib::CtPtrs_VecCt(enca),
                          helib::CtPtrs_VecCt(encb),
                          false);
        decryptBinaryNums(slotsMax, wMax, encryptionKey, ea);
        decryptBinaryNums(slotsMin, wMin, encryptionKey, ea);
    } // get rid of the wrapper
    ea.decrypt(mu, encryptionKey, slotsMu);
    ea.decrypt(ni, encryptionKey, slotsNi);

    auto muniTuple = std::make_tuple(pMax, pMin, pMu, pNi);
    auto pmupniTuple = std::make_tuple(slotsMax[0], slotsMin[0], slotsMu[0], slotsNi[0]);
    std::cout
            << "Comparison (with min max) error: a=" << pa << ", b=" << pb
            << ", but min=" << slotsMin[0] << ", max=" << slotsMax[0]
            << ", mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << std::endl;
    std::cout << std::get<0>(muniTuple) << ", " << std::get<1>(muniTuple) << ", " << std::get<2>(muniTuple) << ", " << std::get<3>(muniTuple)
         << std::endl;
    std::cout << std::get<0>(pmupniTuple) << ", " << std::get<1>(pmupniTuple) << ", " << std::get<2>(pmupniTuple) << ", "
         << std::get<3>(pmupniTuple) << std::endl;
    if (VERBOSE) {
        std::cout << "Comparison (with min max) succeeded: ";
        std::cout << '(' << pa << ',' << pb << ")=>(" << slotsMin[0] << ','
                  << slotsMax[0] << "), mu=" << slotsMu[0] << ", ni=" << slotsNi[0]
                  << std::endl;
    }

#ifdef HELIB_DEBUG
    const helib::Ctxt* minLvlCtxt = nullptr;
  long minLvl = 1000;
  for (const helib::Ctxt& c : eMax) {
    long lvl = c.logOfPrimeSet();
    if (lvl < minLvl) {
      minLvlCtxt = &c;
      minLvl = lvl;
    }
  }
  decryptAndPrint((std::cout << " after comparison: "),
                  *minLvlCtxt,
                  secKey,
                  ea,
                  0);
  std::cout << std::endl;
#endif

    if (VERBOSE) helib::printAllTimers(std::cout);


    //-----------------------------------------------------------

    printDuration(t0_main, "Main");
    logger.print_log(log_error, false);


    return 0;
}







////// GLOBAL VARIABLES good example
//int g_x{}; // global variable g_x
//void doSomething() {
//    // global variables can be seen and used everywhere in the file
//    g_x = 3;
//    std::cout << g_x << '\n';
//    int g_x = 2;
//    std::cout << g_x << '\n';
//}
//int main() {
//    std::cout << g_x << '\n';
//    doSomething();
//    std::cout << g_x << '\n';
//    g_x = 5;
//    std::cout << g_x << '\n';

//    int mynum;
//    std::cout << "enter number->" << std::endl;
//    std::cin >> mynum;
//    printNameVal(mynum);

//}
