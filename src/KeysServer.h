//
// The KeysServer functionallity.
// Distribution of keys for Clients and the DataServer.
//  Big parts of the skfj c'tor and init uses GTestBinaryCompare code snipets
//

#ifndef ENCKMEAN_KEYSSERVER_H
#define ENCKMEAN_KEYSSERVER_H

#include "utils/aux.h"

/**
 * @class KeysServer
 * @brief A wrapper for the Keys Server (CA)
 * Creates the context and keys for the encryption services.
 * Manages protocols for creating Client and Data-Centers (DataServer) shared keys.
 * */
class KeysServer {
public:
    static std::vector<helib::zzX> unpackSlotEncoding;
    helib::PubKey &public_key;
    //    helib::PubKey &pubKeyRef;
    //    publicKey.writeTo(str));
    //    std::shared_ptr<helib::PubKey> deserialized_pkp;// =            std::make_shared<helib::PubKey>(helib::PubKey::readFrom(str, context));
protected:
    friend class TestPoint;

    friend class TestDataServer;

    //fixme -https://stackoverflow.com/questions/3903180/make-a-friend-class-have-only-special-access-to-1-function-of-another-class
    friend class Client;


    constexpr static long mValues[][15] = { // todo note it's 15 param for the 5 first line and 14 for the rest (because you took them from 2 different files)
            // { p, phi(m),   m,   d, m1, m2, m3,    g1,   g2,   g3, ord1,ord2,ord3, B,c}
            {2,   48,    105,   12, 3,   35,   0,   71,    76,    0,     2,   2,   0,   25, 2},
            {2,   600,   1023,  10, 11,  93,   0,   838,   584,   0,     10,  6,   0,   25, 2},
            {2,   2304,  4641,  24, 7,   3,    221, 3979,  3095,  3760,  6,   2,   -8,  25, 3},
            {2,   15004, 15709, 22, 23,  683,  0,   4099,  13663, 0,     22,  31,  0,   25, 3},
            // clang-format off
            {2,   27000, 32767, 15, 31,  7,    151, 11628, 28087, 25824, 30,  6,   -10, 28, 4},
            // clang-format on
            {7,   36,    57,    3,  3,   19,   0,   20,    40,    0,     2,   -6,  0,   100}, // m=3*(19) :-( m/phim(m)=1.58 C=14 D=3 E=0

            {17,  48,    105,   12, 3,   35,   0,   71,    76,    0,     2,   2,   0,   100}, // m=3*(5)*{7} m/phim(m)=2.18 C=14 D=2 E=2
            {17,  576,   1365,  12, 7,   3,    65,  976,   911,   463,   6,   2,   4,   100}, // m=3*(5)*7*{13} m/phim(m)=2.36  C=22  D=3
            {17,  18000, 21917, 30, 101, 217,  0,   5860,  5455,  0,     100, 6,   0,   100}, // m=(7)*{31}*101 m/phim(m)=1.21  C=134 D=2
            {17,  30000, 34441, 30, 101, 341,  0,   2729,  31715, 0,     100, 10,  0,   100}, // m=(11)*{31}*101 m/phim(m)=1.14 C=138 D=2
            {17,  40000, 45551, 40, 101, 451,  0,   19394, 7677,  0,     100, 10,  0,   200}, // m=(11)*{41}*101 m/phim(m)=1.13 C=148 D=2
            {17,  46656, 52429, 36, 109, 481,  0,   46658, 5778,  0,     108, 12,  0,   100}, // m=(13)*{37}*109 m/phim(m)=1.12 C=154 D=2
            {17,  54208, 59363, 44, 23,  2581, 0,   25811, 5199,  0,     22,  56,  0,   100}, // m=23*(29)*{89} m/phim(m)=1.09  C=120 D=2
            {17,  70000, 78881, 10, 101, 781,  0,   67167, 58581, 0,     100, 70,  0,   100}, // m=(11)*{71}*101 m/phim(m)=1.12 C=178 D=2

            {127, 576,   1365,  12, 7,   3,    65,  976,   911,   463,   6,   2,   4,   100}, // m=3*(5)*7*{13} m/phim(m)=2.36   C=22  D=3
            {127, 1200,  1925,  20, 11,  175,  0,   1751,  199,   0,     10,  6,   0,   100}, //  m=(5^2)*{7}*11 m/phim(m)=1.6   C=34 D=2
            {127, 2160,  2821,  30, 13,  217,  0,   652,   222,   0,     12,  6,   0,   100}, // m=(7)*13*{31} m/phim(m)=1.3     C=46 D=2
            {127, 18816, 24295, 28, 43,  565,  0,   16386, 16427, 0,     42,  16,  0,   100}, // m=(5)*43*{113} m/phim(m)=1.29   C=84  D=2
            {127, 26112, 30277, 24, 17,  1781, 0,   14249, 10694, 0,     16,  68,  0,   100}, // m=(13)*17*{137} m/phim(m)=1.15  C=106 D=2
            {127, 31752, 32551, 14, 43,  757,  0,   7571,  28768, 0,     42,  54,  0,   100}, // m=43*(757) :-( m/phim(m)=1.02   C=161 D=3
            {127, 46656, 51319, 36, 37,  1387, 0,   48546, 24976, 0,     36,  -36, 0,   200}, //m=(19)*37*{73}:-( m/phim(m)=1.09 C=141 D=3
            {127, 49392, 61103, 28, 43,  1421, 0,   1422,  14234, 0,     42,  42,  0,   200}, // m=(7^2)*{29}*43 m/phim(m)=1.23  C=110 D=2
            {127, 54400, 61787, 40, 41,  1507, 0,   30141, 46782, 0,     40,  34,  0,   100}, // m=(11)*41*{137} m/phim(m)=1.13  C=112 D=2
            {127, 72000, 77531, 30, 61,  1271, 0,   7627,  34344, 0,     60,  40,  0,   100}  // m=(31)*{41}*61 m/phim(m)=1.07   C=128 D=2
    };

    const long prm; // parameter size (0-tiny,...,4-huge) //todo this says which row is chosen from mValues
    const long bitSize; // itSize of input integers (<=32)
    const bool bootstrap; // comparison with bootstrapping (??)
    const long seed; // PRG seed
    const long nthreads; // number of threads
    const long *vals;   //todo this is initialized w/ the row #prm chosen from mValues
    const long p;
    const long m;
    const NTL::Vec<long> mvec;
    const std::vector<long> gens;
    const std::vector<long> ords;
    const long c;
    const long L;
    helib::Context context;
    helib::SecKey secKey; //private? //reference?
    helib::SecKey &secKeyRef; //private? //reference?

public:

    explicit KeysServer(
            long prm = 0, // parameter size (0-tiny,...,4-huge) //  CT bigger is slower...
            long bitSize = BIT_SIZE, // bitSize of input integers (<=32)
            bool bootstrap = true, // comparison with bootstrapping
            // (KT-26.oct.21) definitely make bootstrap true - for cmp w/ min/max (and for huge number of points(?))
            long seed = 0, // PRG seed
            long nthreads = N_Threads // number of threads
    )
            :
            prm(validatePrm(prm)),//todo this says which row is chosen from mValues
            bitSize(correctBitSize(5, bitSize)),
            bootstrap(bootstrap),
            // (KT-26.oct.21) definitely make bootstrap true - for cmp w/ min/max (and for huge number of points(?))
            seed(seed),
            nthreads(nthreads),
            vals(mValues[prm]),   //todo this is initialized w/ the row #prm chosen from mValues
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
            //            fla(prepareSecKey(secKey)),
            secKeyRef(secKey),
            // In HElib, the SecKey class is actually a subclass if the PubKey class.  So
            // one way to initialize a public key object is like this:

            // TECHNICAL NOTE: Note the "&" in the declaration of publicKey. Since the
            // SecKey class is a subclass of PubKey, this particular PubKey object is
            // ultimately a SecKey object, and through the magic of C++ polymorphism,
            // encryptions done via publicKey will actually use the secret key, which has
            // certain advantages. If one left out the "&", then encryptions done via
            // publicKey will NOT use the secret key.
            public_key(secKey) {

        if (seed) NTL::SetSeed(NTL::ZZ(seed));
        if (nthreads > 1) NTL::SetNumThreads(nthreads);

        prepareSecKey(secKey);

        helib::activeContext = &context; // make things a little easier sometimes

        helib::setupDebugGlobals(&secKey, context.shareEA());

    }

    // TECHNICAL NOTE: Note the "&" in the declaration of publicKey. Since the
    // SecKey class is a subclass of PubKey, this particular PubKey object is
    // ultimately a SecKey object, and through the magic of C++ polymorphism,
    // encryptions done via publicKey will actually use the secret key, which has
    // certain advantages.  If one left out the "&", then encryptions done via
    // publicKey will NOT use the secret key.
    helib::PubKey &getPublicKey() const {
        return public_key;
        //        return deserialized_pkp;
    }

    /* * *  for DBG    * * */
    helib::Ctxt encryptCtxt(bool b) const {
        NTL::ZZX pl(b);
        helib::Ctxt cl(public_key);
        public_key.Encrypt(cl, pl);
        return cl;
    }

    EncryptedNum encryptNum(long l) const {
        helib::PubKey & public_key = getPublicKey();
        EncryptedNum cl(BIT_SIZE, helib::Ctxt(public_key));
        for (long bit = 0; bit < BIT_SIZE; ++bit)
            public_key.Encrypt(cl[bit],
                                     NTL::to_ZZX((l >> bit) & 1));
        return cl;
    }

    long decryptCtxt(const helib::Ctxt &ctxt) const;

    long decryptNum(const EncryptedNum &cNum) const;

    long decryptSize(const std::vector<helib::Ctxt> &size) const;
    /* * *  end for DBG    * * */


protected:
    helib::SecKey &getSecKey() const { // return CONST SecKey?
        return secKeyRef;
    }

    [[nodiscard]] const helib::Context &getContext() const {
        return context;
    }

    [[nodiscard]] const helib::EncryptedArray &getEA() const {
        return context.getEA();
    }

    [[nodiscard]] helib::IndexSet getCtxtPrimes(long nprimes = 5) const {
        return context.getCtxtPrimes(nprimes);
    }


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

public:
    /*  services    */ //todo consider moving back to DataServer
    /**
     * @brief return a scratch point. for cases of ...(?)
     * */
    const Point scratchPoint() const;

    /**
     * @brief  returns a tiny point, as close to zero point as possible,
     * within EPSILON margin.
     * used, by data server, for classification of null (zero) points -
     * so as to not assign then to one of the cells
     * */
    const Point tinyRandomPoint() const;

    const Point getQuotientPoint(const Point &point, const std::vector<Ctxt> &sizeBitVector,
                                 const short repsNum) const;

    const EncryptedNum getQuotient(const EncryptedNum &encryptedNum, const long num) const;
};


#endif //ENCKMEAN_KEYSSERVER_H
