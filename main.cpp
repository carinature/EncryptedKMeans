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
#include "src/Client.h"

//helib
#include <helib/binaryArith.h>
#include <helib/binaryCompare.h>
#include <bitset>


int main() {
    //  lLgging
    auto t0_main = std::chrono::high_resolution_clock::now();
    Logger logger;
    logger.log("Starting Protocol", log_trace);
    //  Keys Server
    KeysServer keysServer;
    logger.log(printDuration(t0_main, "KeysServer Initialization"));

    const helib::EncryptedArray &ea = keysServer.getEA();    //    helib::Context context = keysServer.getContext();

    cout << "start client" << endl;

    const long arr[] = {1L, 2l};
//    Client client1(keysServer, arr);
//    Client client2(keysServer, arr);
    Client client1(keysServer);
    Client client2(keysServer);
    client1.encryptPoint(arr);
    client2.encryptPoint(arr);

    cout << "fin client" << endl;

    const helib::PubKey &publicKey = keysServer.getSecKey();
    long n = keysServer.getContext().getNSlots();


    //// This part should be done in Client/Point
    NTL::Vec<helib::Ctxt> eMax, eMin, enca, encb;
    const helib::SecKey &encryptionKey = keysServer.getSecKey();
    helib::Ctxt mu(encryptionKey), ni(encryptionKey);
    resize(enca, long(BIT_SIZE), mu);
    resize(encb, long(BIT_SIZE) + 1, ni);

#ifdef alt
    // Choose two random n-bit integers
    long pa = NTL::RandomBits_long(BIT_SIZE);
    long pb = NTL::RandomBits_long(BIT_SIZE + 1);

    long pMax = std::max(pa, pb);
    long pMin = std::min(pa, pb);
    bool pMu = pa > pb;
    bool pNi = pa < pb;

    //    auto pv = helib::PtrVector<helib::Ctxt>(pa);

    // Encrypt the individual bits
    for (long i = 0; i <= BIT_SIZE; i++) {
        //todo in future remove this `if`,
        // currently serves as an example for init/enc of 2 different-sized numbers
        if (i < BIT_SIZE)
            encryptionKey.Encrypt(enca[i], NTL::ZZX((pa >> i) & 1));
        encryptionKey.Encrypt(encb[i], NTL::ZZX((pb >> i) & 1));
        //        if (helib_bootstrap) { // put them at a lower level
        //            if (i < BIT_SIZE)
        //                enca[i].bringToSet(keysServer.getCtxtPrimes(5));
        //            encb[i].bringToSet(keysServer.getCtxtPrimes(5));
        //        }
    }

//
//    printNameVal(n);
//    printNameVal(n);
//    printNameVal(n);
//    printNameVal(n);
//    printNameVal(n);
//    printNameVal(n);
//    std::vector<long> v0(n);
//    for (long i = 0; i < n; i++)
//        v0[i] = i;
//    helib::PtxtArray ptxtArray(keysServer.getContext(), v0);
//    helib::Ctxt c0(publicKey);
//    ptxtArray.encrypt(c0);
//
//    compareTwoNumbers(mu, ni, helib::CtPtrs_VecCt(c0), helib::CtPtrs_VecCt(c0), false, nullptr);
#else
    // Each bit of the binary number is encoded into a single ciphertext. Thus
    // for a 16 bit binary number, we will represent this as an array of 16
    // unique ciphertexts.
    // i.e. b0 = [0] [0] [0] ... [0] [0] [0]        ciphertext for bit 0
    //      b1 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 1
    //      b2 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 2
    // These 3 ciphertexts represent the 3-bit binary number 110b = 6

    // Note: several numbers can be encoded across the slots of each ciphertext
    // which would result in several parallel slot-wise operations.
    // For simplicity we place the same data into each slot of each ciphertext,
    // printing out only the back of each vector.
    // NB: fifteenOrLess4Four max is 15 bits. Later in the code we pop the MSB.
    long bitSize = 16;
    long outSize = 2 * bitSize;
    long a_data = NTL::RandomBits_long(bitSize);
    long b_data = NTL::RandomBits_long(bitSize);
    long c_data = NTL::RandomBits_long(bitSize);

    // Use a scratch ciphertext to populate vectors.
    const helib::PubKey &public_key = encryptionKey;
//    const helib::PubKey &public_key = keysServer.getSecKey();
    helib::Ctxt scratch(public_key);
    std::vector<helib::Ctxt> encrypted_a(bitSize, scratch);
    std::vector<helib::Ctxt> encrypted_b(bitSize, scratch);
    std::vector<helib::Ctxt> encrypted_c(bitSize, scratch);
//    cout << "odin" << endl;
    // Encrypt the data in binary representation.
    for (long i = 0; i < bitSize; ++i) {
//        cout << "dva" << endl;
        std::vector<long> a_vec(ea.size());
        std::vector<long> b_vec(ea.size());
        std::vector<long> c_vec(ea.size());
        // Extract the i'th bit of a,b,c.
        for (auto &slot : a_vec)
            slot = (a_data >> i) & 1;
        for (auto &slot : b_vec)
            slot = (b_data >> i) & 1;
        for (auto &slot : c_vec)
            slot = (c_data >> i) & 1;
        ea.encrypt(encrypted_a[i], public_key, a_vec);
        ea.encrypt(encrypted_b[i], public_key, b_vec);
        ea.encrypt(encrypted_c[i], public_key, c_vec);
//        cout << "tri" << endl;

        std::vector<long> slotsMin, slotsMax, slotsMu, slotsNi;
        // comparison only
        compareTwoNumbers(mu, ni,
                          helib::CtPtrs_vectorCt(encrypted_a),
                          helib::CtPtrs_vectorCt(encrypted_a),
                          false);
//        cout << "chetiri" << endl;
        ea.decrypt(mu, encryptionKey, slotsMu);
        ea.decrypt(ni, encryptionKey, slotsNi);
//        for (auto s : slotsMu) printNameVal(s);
//        for (auto s : slotsNi) printNameVal(s);

    }
    // Although in general binary numbers are represented here as
    // std::vector<helib::Ctxt> the binaryArith APIs for HElib use the PtrVector
    // wrappers instead, e.g. helib::CtPtrs_vectorCt. These are nothing more than
    // thin wrapper classes to standardise access to different vector types, such
    // as NTL::Vec and std::vector. They do not take ownership of the underlying
    // object but merely provide access to it.
    //
    // helib::CtPtrMat_vectorCt is a wrapper for
    // std::vector<std::vector<helib::Ctxt>>, used for representing a list of
    // encrypted binary numbers.
#endif

#ifdef HELIB_DEBUG
    decryptAndPrint((cout << " before comparison: "),
                  encb[0],
                  encryptionKey,
                  ea,
                  0);
#endif

#if alt
    //// this should stay in main but use `Client` class compare
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
        cout << "Comparison (with min max) : " << '(' << pa << ',' << pb << ")=>(" << slotsMin[0]
             << ','
             << slotsMax[0] << "), mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << endl;
#endif

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
                  encryptionKey,
                  ea,
                  0);
  cout << endl;
#endif

    if (VERBOSE) helib::printAllTimers(cout);


    //-----------------------------------------------------------

    printDuration(t0_main, "Main");
    logger.print_log(log_trace);//, false);
    //    clientLogger.print_log(log_error, false);

    return 0;
}



