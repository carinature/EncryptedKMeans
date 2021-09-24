//
// The KeysServer functionallity.
// Distibution of keys for Clients and the DataServer
//

#include "KeysServer.h"
#include "properties.h"
#include "utils/Logger.h"
//#include "utils/aux.h"


#include <iostream>
#include <helib/Context.h>
#include <helib/intraSlot.h> //for buildUnpackSlotEncoding
#include <helib/debugging.h> //for setupDebugGlobals

using std::cout;
using std::endl;

Logger keysServerLogger(log_debug,"keysServerLogger");//todo change to log_trace

// Validates the prm value, throwing if invalid
// [prm] Corresponds to the number of entry in mValues table
long KeysServer::validatePrm(long prm) {
    if (prm < 0 || prm >= 5) throw std::invalid_argument("prm must be in the interval [0, 4]");
    return prm;
};

long KeysServer::correctBitSize(long minimum, long oldBitSize) {
    long newBitSize;
    if (oldBitSize <= 0) newBitSize = minimum;
    else if (oldBitSize > 32) newBitSize = 32;
    else newBitSize = oldBitSize;
    return newBitSize;
};

NTL::Vec<long> KeysServer::calculateMvec(const long *vals) {
    NTL::Vec<long> mvec;
    append(mvec, vals[4]);
    if (vals[5] > 1) append(mvec, vals[5]);
    if (vals[6] > 1) append(mvec, vals[6]);
    return mvec;
};

std::vector<long> KeysServer::calculateGens(const long *vals) {
    std::vector<long> gens;
    gens.push_back(vals[7]);
    if (vals[8] > 1) gens.push_back(vals[8]);
    if (vals[9] > 1) gens.push_back(vals[9]);
    return gens;
};

std::vector<long> KeysServer::calculateOrds(const long *vals) {
    std::vector<long> ords;
    ords.push_back(vals[10]);
    if (abs(vals[11]) > 1) ords.push_back(vals[11]);
    if (abs(vals[12]) > 1) ords.push_back(vals[12]);
    return ords;
};

long KeysServer::calculateLevels(bool bootstrap, long bitSize) {
    return bootstrap
           ? 900
           : 30 * (7 + NTL::NumBits(bitSize + 2)); // that should be enough
};

helib::Context &KeysServer::prepareContext(helib::Context &contxt) {
    if (VERBOSE) {
        cout << "input bitSize=" << bitSize << endl;
        if (nthreads > 1) cout << "  using " << NTL::AvailableThreads() << " threads\n";
        cout << "computing key-independent tables..." << std::flush;
    }
    contxt.buildModChain(L, c, /*willBeBootstrappable=*/bootstrap);
    if (bootstrap) contxt.enableBootStrapping(mvec);

    buildUnpackSlotEncoding(KeysServer::unpackSlotEncoding, contxt.getEA());

    if (VERBOSE) {
        cout << " done.\n";
        contxt.printout();
    }

    return contxt;
}

void KeysServer::prepareSecKey(helib::SecKey &key) const {
    if (VERBOSE) {
        cout << "\ncomputing key-dependent tables..." << std::flush;
    }
    key.GenSecKey();
    addSome1DMatrices(key); // compute key-switching matrices
    addFrbMatrices(key);
    if (bootstrap) key.genRecryptData();
    if (VERBOSE) cout << " done\n";
};

std::vector<helib::zzX> KeysServer::unpackSlotEncoding; //todo move? already defined in class - check what a 2nd def does

//! KeysServer c'tor: default values
// NOTE: The parameters used in this example code are for demonstration only.
// They were chosen to provide the best performance of execution while
// providing the context to demonstrate how to use the HElib's FHE encryption.
// (The parameters may not provide the security level that might be required by real use/application scenarios.)
KeysServer::KeysServer(long prm, long bitSize, bool bootstrap, long seed, long nthreads) :
        prm(validatePrm(prm)),
        bitSize(correctBitSize(5, bitSize)),
        bootstrap(bootstrap),
        seed(seed),
        nthreads(nthreads),
        vals(mValues[prm]),
        p(vals[0]),
        m(vals[2]),
        mvec(calculateMvec(vals)),
        gens(calculateGens(vals)),
        ords(calculateOrds(vals)),
        c(vals[14]),
        L(calculateLevels(bootstrap, bitSize)),
        context(helib::ContextBuilder<helib::BGV>()
                        .m(m)
                        .p(p)
                        .r(1)
                        .gens(gens)
                        .ords(ords)
                        .buildModChain(false)
                        .build()),
        secKey(prepareContext(context)),
        // In HElib, the SecKey class is actually a subclass if the PubKey class.  So
        // one way to initialize a public key object is like this:

        // TECHNICAL NOTE: Note the "&" in the declaration of publicKey. Since the
        // SecKey class is a subclass of PubKey, this particular PubKey object is
        // ultimately a SecKey object, and through the magic of C++ polymorphism,
        // encryptions done via publicKey will actually use the secret key, which has
        // certain advantages.  If one left out the "&", then encryptions done via
        // publicKey will NOT use the secret key.
        pubKey(secKey) {

    if (seed) NTL::SetSeed(NTL::ZZ(seed));
    if (nthreads > 1) NTL::SetNumThreads(nthreads);

    prepareSecKey(secKey);

    helib::activeContext = &context; // make things a little easier sometimes

    helib::setupDebugGlobals(&secKey, context.shareEA());

}


//long KeysServer::decrypt(const helib::Ctxt& cBit) {
//    //    long pBit = 0;
//    NTL::ZZX pp;
//    secKey.Decrypt(pp, cBit);
//    return IsOne(pp);
//}
//
//long KeysServer::decrypt(const std::vector<helib::Ctxt> cNum) {
//    long pNum = 0;
//    NTL::ZZX pp;
//    for (int bit = 0; bit < bitSize; ++bit) {
//        secKey.Decrypt(pp, cNum[bit]);
//        printNameVal(pp);
//        if (IsOne(pp)) pNum += std::pow(2, bit);
//        printNameVal(pNum);
//    }
//    return pNum;
//}
//
//std::vector<long> KeysServer::decrypt(const Point& p) {
//    cout << "KeysServer::decryptCoordinates" << endl;
//    keysServerLogger.log("decryptCoordiantes", log_debug);
//    std::vector<long> dCoordinates(DIM);
//    if (!p[0][0].isEmpty()) return dCoordinates;
//    for (int dim = 0; dim < DIM; ++dim) dCoordinates[dim] = decrypt(p[dim]);
//    return dCoordinates;
//}

//std::vector<long> KeysServer::decryptCtxt(helib::Ctxt ctxt) {
//#if VERBOSE
//    cout << "decryptCoordiantes" << endl;
//#endif
//    keysServerLogger.log("decryptCoordiantes", log_debug);
//    std::vector<long> decrypted_result(); // todo redundant, decide which of the two to keep
//    for (int i = 0; i < DIM; ++i) {
//        helib::decryptBinaryNums(decrypted_result,
//                                 helib::CtPtrs_vectorCt(ctxt),
//                                 secKey, getEA());
//        dCoordinates[i] = decrypted_result[i].back();
//    }
//#if VERBOSE
//    for (auto dec_coordinate :decrypted_result) printNameVal(dec_coordinate.back());
//        cout << endl << endl;
//#endif
//    return dCoordinates;
//}










