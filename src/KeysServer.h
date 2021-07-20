//
// The KeysServer functionallity.
// Distribution of keys for Clients and the DataServer.
//  Big parts of the skfj c'tor and init uses GTestBinaryCompare code snipets
//


#ifndef ENCKMEAN_KEYSSERVER_H
#define ENCKMEAN_KEYSSERVER_H


#include <vector>
#include <helib/zzX.h>
#include <helib/Context.h>
#include <helib/keys.h>

class KeysServer {

protected:
    static std::vector<helib::zzX> unpackSlotEncoding;
    constexpr static long mValues[][15] = {
            // { p, phi(m),   m,   d, m1, m2, m3,    g1,   g2,   g3, ord1,ord2,ord3, B,c}
            {2, 48,    105,   12, 3,  35,  0,   71,    76,    0,     2,  2,  0,   25, 2},
            {2, 600,   1023,  10, 11, 93,  0,   838,   584,   0,     10, 6,  0,   25, 2},
            {2, 2304,  4641,  24, 7,  3,   221, 3979,  3095,  3760,  6,  2,  -8,  25, 3},
            {2, 15004, 15709, 22, 23, 683, 0,   4099,  13663, 0,     22, 31, 0,   25, 3},
            // clang-format off
            {2, 27000, 32767, 15, 31, 7,   151, 11628, 28087, 25824, 30, 6,  -10, 28, 4}
            // clang-format on
    };

    const long prm;
    const long bitSize;
    const bool bootstrap;
    const long seed;
    const long nthreads;
    const long *vals;
    const long p;
    const long m;
    const NTL::Vec<long> mvec;
    const std::vector<long> gens;
    const std::vector<long> ords;
    const long c;
    const long L;
    helib::Context context;
protected:
    helib::SecKey secKey;
public:
    KeysServer();
    //    explicit KeysServer(const long prm=0);

    // Getters
    [[nodiscard]] const helib::Context &getContext() const {
        return context;
    }

    [[nodiscard]] const helib::EncryptedArray &getEA() const {
        return context.getEA();
    }

    [[nodiscard]] helib::IndexSet getCtxtPrimes(long nprimes = 5) const {
        return context.getCtxtPrimes(nprimes);
    }

    [[nodiscard]] const helib::SecKey &getSecKey() const {
        return secKey;
    }


private:
    // used for c'tor and "sys"/KeysServer bootstrapping
    static long validatePrm(long prm);

    static long correctBitSize(long minimum, long oldBitSize);

    static NTL::Vec<long> calculateMvec(const long *vals);

    static std::vector<long> calculateGens(const long *vals);

    static std::vector<long> calculateOrds(const long *vals);

    static long calculateLevels(bool bootstrap, long bitSize);

    helib::Context &prepareContext(helib::Context &context);

    void prepareSecKey(helib::SecKey &secKey) const;

};


#endif //ENCKMEAN_KEYSSERVER_H