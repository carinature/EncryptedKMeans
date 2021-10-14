//
// The KeysServer functionallity.
// Distibution of keys for Clients and the DataServer
//

#include "utils/aux.h"
#include "KeysServer.h"
#include "Point.h"

Logger keysServerLogger(log_debug, "keysServerLogger");//todo change to log_trace

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
        cout << "input BIT_SIZE=" << bitSize << endl;
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
//KeysServer::KeysServer(long prm, long bitSize, bool bootstrap, long seed, long nthreads)


long KeysServer::decryptCtxt(const helib::Ctxt &cBit) {
    if (cBit.isEmpty()) return -1; //todo should return 0?
    NTL::ZZX pp;
    secKey.Decrypt(pp, cBit);
    return IsOne(pp);
}

//long KeysServer::decryptNum(std::vector<helib::Ctxt> cNum, bool isProduct) {
long KeysServer::decryptNum(const std::vector<helib::Ctxt> &cNum) {
    if (!cNum.size()) return -1; //todo should return 0?
    long pNum = 0;
    NTL::ZZX pp;
    //    int out_size = (1 + isProduct) * BIT_SIZE;
    //    for (int bit = 0; bit < out_size; ++bit) {
    for (int bit = 0; bit < cNum.size(); ++bit) {
        secKey.Decrypt(pp, cNum[bit]);
        //        printNameVal(pp);
        if (IsOne(pp)) pNum += std::pow(2, bit);
        //        printNameVal(pNum);
    }
    return pNum;
}

long KeysServer::decryptSize(const std::vector<helib::Ctxt> &cSize) {
    long size = 0;
    NTL::ZZX pp;
    for (const Ctxt &ctxt: cSize) {
        secKey.Decrypt(pp, ctxt);
        //        printNameVal(pp);
        size += IsOne(pp);
    }
    return size;
}


const Point KeysServer::scratchPoint() const {
    cout << " scratchPoint" << endl;
    return Point(getPublicKey());//, nullptr);
}

const Point KeysServer::tinyRandomPoint() const {
    long arr[DIM];
    for (short dim = 0; dim < DIM; ++dim) arr[dim] = rand() % 1 * epsilon;
    // note that despite the rand illusion, currently this always returns 0
    // which is perfectly for us, but the "real" solution will be ok too
    // to create the "real" solution we'll change the `rand()` function (and `long` to `double`)
    return Point(getPublicKey(), arr);
}



