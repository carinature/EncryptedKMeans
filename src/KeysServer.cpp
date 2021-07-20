//
// The KeysServer functionallity.
// Distibution of keys for Clients and the DataServer
//

#include "KeysServer.h"
#include "../properties.h"
#include "../utils/Logger.h"
#include "../utils/aux.h"


#include <iostream>
#include <helib/Context.h>
#include <helib/intraSlot.h> //for buildUnpackSlotEncoding
#include <helib/debugging.h> //for setupDebugGlobals

// fixme remove both
#include <helib/binaryArith.h>
#include <helib/binaryCompare.h>

using namespace std; //fixme remove

// Corresponds to the number of entry in mValues table
//long prm = 0;

//p,      phi(m),     m,          d,        m1, m2, m3,         g1, g2, g3,       ord1, ord2, ord3,         B, c
//2,      48,         105,        12,       3,  35,  0,         71, 76, 0,        2,  2,  0,                25, 2}


// Validates the prm value, throwing if invalid
long KeysServer::validatePrm(long prm) {
    if (prm < 0 || prm >= 5) throw std::invalid_argument("prm must be in the interval [0, 4]");
    return prm;
};

long KeysServer::correctBitSize(long minimum, long oldBitSize) {
    long newBitSize;
    if (oldBitSize <= 0)
        newBitSize = minimum;
    else if (oldBitSize > 32)
        newBitSize = 32;
    else
        newBitSize = oldBitSize;
    return newBitSize;
};

NTL::Vec<long> KeysServer::calculateMvec(const long *vals) {
    NTL::Vec<long> mvec;
    append(mvec, vals[4]);
    if (vals[5] > 1)
        append(mvec, vals[5]);
    if (vals[6] > 1)
        append(mvec, vals[6]);
    return mvec;
};

std::vector<long> KeysServer::calculateGens(const long *vals) {
    std::vector<long> gens;
    gens.push_back(vals[7]);
    if (vals[8] > 1)
        gens.push_back(vals[8]);
    if (vals[9] > 1)
        gens.push_back(vals[9]);
    return gens;
};

std::vector<long> KeysServer::calculateOrds(const long *vals) {
    std::vector<long> ords;
    ords.push_back(vals[10]);
    if (abs(vals[11]) > 1)
        ords.push_back(vals[11]);
    if (abs(vals[12]) > 1)
        ords.push_back(vals[12]);
    return ords;
};

long KeysServer::calculateLevels(bool bootstrap, long bitSize) {
    return bootstrap
           ? 900
           : 30 * (7 + NTL::NumBits(bitSize + 2)); // that should be enough
};

helib::Context &KeysServer::prepareContext(helib::Context &context) {
    if (VERBOSE) {
        std::cout << "input bitSize=" << bitSize << std::endl;
        if (nthreads > 1)
            std::cout << "  using " << NTL::AvailableThreads() << " threads\n";
        std::cout << "computing key-independent tables..." << std::flush;
    }
    context.buildModChain(L, c, /*willBeBootstrappable=*/bootstrap);
    if (bootstrap) {
        context.enableBootStrapping(mvec);
    }
    buildUnpackSlotEncoding(KeysServer::unpackSlotEncoding, context.getEA());

    if (VERBOSE) {
        std::cout << " done.\n";
        context.printout();
    }

    return context;
}

void KeysServer::prepareSecKey(helib::SecKey &secKey) const {
    if (VERBOSE) {
        std::cout << "\ncomputing key-dependent tables..." << std::flush;
    }
    secKey.GenSecKey();
    addSome1DMatrices(secKey); // compute key-switching matrices
    addFrbMatrices(secKey);
    if (bootstrap)
        secKey.genRecryptData();
    if (VERBOSE)
        std::cout << " done\n";
};

std::vector<helib::zzX> KeysServer::unpackSlotEncoding; //todo move? already defined in class - check what a 2nd def does

//::Values(
//Parameters(1, 5, false, 0, 1)            // SLOW
// Parameters(0, 5, false, 0, 1)            // FAST
//),
KeysServer::KeysServer() :
        prm(validatePrm(0)),
        bitSize(correctBitSize(5, 5)),
        bootstrap(false),
        seed(0),
        nthreads(0),
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
        secKey(prepareContext(context)) {
    //#endif

    cout << "KeysServer" << endl;
    auto t1 = std::chrono::high_resolution_clock::now();


    if (seed)
        NTL::SetSeed(NTL::ZZ(seed));
    if (nthreads > 1)
        NTL::SetNumThreads(nthreads);

    prepareSecKey(secKey);

    helib::activeContext = &context; // make things a little easier sometimes

    helib::setupDebugGlobals(&secKey, context.shareEA());

    //todo should init and return(?) "global" vars like context , pubKey etc.
    Logger logger;
    //    logger.log(log_trace, "KeysServer init finished");
    logger.log(log_trace, printDuration(t1, "KeysServer()"));

}









