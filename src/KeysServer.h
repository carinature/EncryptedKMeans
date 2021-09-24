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

/**
 * @class KeysServer
 * @brief A wrapper for the Keys Server (CA)
 * Creates the context and keys for the encryption services.
 * Manages protocols for creating Client and Data-Centers (DataServer) shared keys.
 * */
class KeysServer {
public:    static std::vector<helib::zzX> unpackSlotEncoding;

protected:
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

    const long prm; // parameter size (0-tiny,...,4-huge)
    const long bitSize; // itSize of input integers (<=32)
    const bool bootstrap; // comparison with bootstrapping (??)
    const long seed; // PRG seed
    const long nthreads; // number of threads
    const long *vals;
    const long p;
    const long m;
    const NTL::Vec<long> mvec;
    const std::vector<long> gens;
    const std::vector<long> ords;
    const long c;
    const long L;
    helib::Context context;
    helib::SecKey secKey; //private?

public:
    // TECHNICAL NOTE: Note the "&" in the declaration of publicKey. Since the
    // SecKey class is a subclass of PubKey, this particular PubKey object is
    // ultimately a SecKey object, and through the magic of C++ polymorphism,
    // encryptions done via publicKey will actually use the secret key, which has
    // certain advantages.  If one left out the "&", then encryptions done via
    // publicKey will NOT use the secret key.
     helib::PubKey &pubKey;
    helib::PubKey & getPublicKey() {
        return pubKey;
    }

    //::Values(     prm,    bitSize,    bootstrap,  seed,   nthreads)
    // Parameters(  1,      5,          false,      0,      1)            // SLOW
    // Parameters(  0,      5,          false,      0,      1)            // FAST
    explicit KeysServer(long prm = 0, // parameter size (0-tiny,...,4-huge)
                        long bitSize = 5, // itSize of input integers (<=32)
                        bool bootstrap = false, // comparison with bootstrapping (??)
                        long seed = 0, // PRG seed
                        long nthreads = 1); // number of threads

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

    helib::Ctxt enccryptCtxt(long l){
        NTL::ZZX pl(l);
        helib::Ctxt cl(pubKey);
        secKey.Encrypt(cl, pl);
        return cl;
    }

    long decryptCtxt(helib::Ctxt);

    long decryptNum(std::vector<helib::Ctxt> cNum, bool isProduct=false);
    /*{
        long pNum = 0;
            NTL::ZZX pp;
            for (int bit = 0; bit < BIT_SIZE; ++bit) {
                secKey.Decrypt(pp, cNum[bit]);
                            printNameVal(pp);
                if (IsOne(pp)) pNum+=std::pow(2, bit);
                            printNameVal(pNum);
            }
        return pNum;
    }*/

protected:


    //    [[nodiscard]] const helib::SecKey &getSecKey() const {
    [[nodiscard]]helib::SecKey getSecKey() const {
        return secKey;
    }
    friend class TestPoint;
    //fixme -https://stackoverflow.com/questions/3903180/make-a-friend-class-have-only-special-access-to-1-function-of-another-class
    friend class Client;


private:
    // Following are methods from helib
    // All used for c'tor and "sys"/KeysServer bootstrapping
    static long validatePrm(long prm);

    static long correctBitSize(long minimum, long oldBitSize);

    static NTL::Vec<long> calculateMvec(const long *vals);

    static std::vector<long> calculateGens(const long *vals);

    static std::vector<long> calculateOrds(const long *vals);

    static long calculateLevels(bool bootstrap, long bitSize);

    helib::Context &prepareContext(helib::Context &contxt);

    void prepareSecKey(helib::SecKey &key) const;

};


#endif //ENCKMEAN_KEYSSERVER_H
