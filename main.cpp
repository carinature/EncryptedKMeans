//
// run a full protocol
//
#include <iostream>
#include <chrono>
#include <fstream>

using std::cout;
using std::endl;

#include "properties.h"
#include "utils/Logger.h"
#include "utils/aux.h"

//helib
#include <helib/binaryArith.h>
#include <helib/binaryCompare.h>
#include <bitset>



int main() {
    //  loging
    auto t0_main = std::chrono::high_resolution_clock::now();
    Logger logger;
    logger.log(log_trace, "Starting Protocol");
    //  Keys Server
    KeysServer keysServer;
    logger.log(log_trace, printDuration(t0_main, "KeysServer Initialization"));

    const helib::EncryptedArray &ea = keysServer.getEA();    //    helib::Context context = keysServer.getContext();

    // Choose two random n-bit integers
    long pa = NTL::RandomBits_long(bitSize);
    long pb = NTL::RandomBits_long(bitSize + 1);
    long pMax = std::max(pa, pb);
    long pMin = std::min(pa, pb);
    bool pMu = pa > pb;
    bool pNi = pa < pb;

    printNameVal(pa);
    printNameVal(pb);
    printNameVal(pMax);
    printNameVal(pMin);
    printNameVal(pMu);
    printNameVal(pNi);

    std::bitset<LONG_BIT> ba(pa);
    std::bitset<LONG_BIT> bb(pb);
    std::bitset<LONG_BIT> bMax(pMax);
    std::bitset<LONG_BIT> bMin(pMin);
    std::bitset<LONG_BIT> bMu(pMu);
    std::bitset<LONG_BIT> bNi(pNi);
    printNameVal(ba);
    printNameVal(bb);
    printNameVal(bMax);
    printNameVal(bMin);
    printNameVal(bMu);
    printNameVal(bNi);
    std::bitset<LONG_BIT> ba0((pa >> 0) & 1);
    std::bitset<LONG_BIT> ba1((pa >> 1) & 1);
    std::bitset<LONG_BIT> ba2((pa >> 2) & 1);
    std::bitset<LONG_BIT> ba3((pa >> 3) & 1);
    printNameVal(ba0);
    printNameVal(ba1);
    printNameVal(ba2);
    printNameVal(ba3);

    NTL::Vec<helib::Ctxt> eMax, eMin, enca, encb;
    //// This part should be done in Client/Point
    const helib::SecKey &encryptionKey = keysServer.getSecKey();
    helib::Ctxt mu(encryptionKey), ni(encryptionKey);
    resize(enca, long(bitSize), mu);
    resize(encb, long(bitSize) + 1, ni);

    // Encrypt the individual bits
    for (long i = 0; i <= bitSize; i++) {
        //todo in future remove this `if`, currently serves as an example for init/enc of 2 different-sized numbers
        if (i < bitSize)
            encryptionKey.Encrypt(enca[i], NTL::ZZX((pa >> i) & 1));
        encryptionKey.Encrypt(encb[i], NTL::ZZX((pb >> i) & 1));
        //        if (helib_bootstrap) { // put them at a lower level
        //            if (i < bitSize)
        //                enca[i].bringToSet(keysServer.getCtxtPrimes(5));
        //            encb[i].bringToSet(keysServer.getCtxtPrimes(5));
        //        }
    }

#ifdef HELIB_DEBUG
    decryptAndPrint((cout << " before comparison: "),
                  encb[0],
                  secKey,
                  ea,
                  0);
#endif

    std::vector<long> slotsMin, slotsMax, slotsMu, slotsNi;
    // comparison only
    compareTwoNumbers(mu, ni,
                      helib::CtPtrs_VecCt(enca),
                      helib::CtPtrs_VecCt(encb),
                      false);
    ea.decrypt(mu, encryptionKey, slotsMu);
    ea.decrypt(ni, encryptionKey, slotsNi);
    auto muniPair = std::make_pair(slotsMu[0], slotsNi[0]);
    auto pmupniPair = std::make_pair((long) pMu, (long) pNi);
    cout << "Comparison (without min max) : a=" << pa << ", b=" << pb
              << ", mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << endl;
    cout << std::get<0>(muniPair) << ", " << std::get<1>(muniPair) << endl;
    cout << std::get<0>(pmupniPair) << ", " << std::get<1>(pmupniPair) << endl;
    if (VERBOSE)
        cout << "Comparison (without min max) : " << '(' << pa << ',' << pb << ")"
                  << "=> mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << endl;

    {
        helib::CtPtrs_VecCt wMin(eMin),
                wMax(eMax); // A wrappers around output vectors
        // comparison with max and min
        compareTwoNumbers(wMax, wMin,
                          mu, ni,
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
    cout << "Comparison (with min max) a=" << pa << ", b=" << pb
              << ", min=" << slotsMin[0] << ", max=" << slotsMax[0]
              << ", mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << endl;
    cout << std::get<0>(muniTuple) << ", " << std::get<1>(muniTuple) << ", "
              << std::get<2>(muniTuple) << ", " << std::get<3>(muniTuple) << endl;
    cout << std::get<0>(pmupniTuple) << ", " << std::get<1>(pmupniTuple) << ", "
              << std::get<2>(pmupniTuple) << ", " << std::get<3>(pmupniTuple) << endl;
    if (VERBOSE)
        cout << "Comparison (with min max) : " << '(' << pa << ',' << pb << ")=>(" << slotsMin[0] << ','
                  << slotsMax[0] << "), mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << endl;


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
  decryptAndPrint((cout << " after comparison: "),
                  *minLvlCtxt,
                  secKey,
                  ea,
                  0);
  cout << endl;
#endif

    if (VERBOSE) helib::printAllTimers(cout);


    //-----------------------------------------------------------

    printDuration(t0_main, "Main");
    logger.print_log(log_error, false);


    return 0;
}



