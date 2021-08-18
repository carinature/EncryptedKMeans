
#ifndef ENCKMEAN_CLIENT_H
#define ENCKMEAN_CLIENT_H

#include "utils/aux.h" // for including KeysServer.h


class Point {
    std::vector<std::vector<helib::Ctxt> > cCoordinates;
    const helib::PubKey &public_key;// = encryptionKey;

public:
    Point(const helib::EncryptedArray &ea,
          const helib::PubKey &public_key,
          const long coordinates[] = nullptr) :
            public_key(public_key) {
        std::vector<long> a_vec(ea.size());
        for (int dim = 0; dim < DIM; ++dim) {
            cCoordinates.emplace_back(bitSize, helib::Ctxt(public_key));
            for (long bit = 0; bit < bitSize; ++bit) {
                // Extract the i'th bit of coordinates[dim].
                for (auto &slot : a_vec) slot = (coordinates[dim] >> bit) & 1;
                if (coordinates) ea.encrypt(cCoordinates[dim][bit], public_key, a_vec);
            }
        }
    }
};


/**
 * @class Client
 * @brief the Client can encrypt info using the key from the Keys Server. decryption is currently only for dbg purposes.
 * */
class Client {
    const helib::SecKey encryptionKey;
    const helib::EncryptedArray ea;
    const helib::Ctxt scratch;        // Use a scratch ciphertext to populate vectors.

protected:
    //    NTL::Vec<helib::Ctxt> cCoordinatesNTL;
    std::vector<std::vector<helib::Ctxt> > cCoordinatesStd;
    const helib::PubKey &public_key;// = encryptionKey;
    //#if DBG
    //private:
    //    long *pCoordinates;
    //#endif

public:
    /**
     * Constructor for \class{Client},
     * @param keysServer binds to the \class{KeysServer} responsible for the distributing the appropriate key
     * @brief simulates a mini-protocol between the keys server and the client
     * */
    explicit Client(KeysServer &keysServer);

    /**
     *
     * */
    std::vector<std::vector<helib::Ctxt>> encryptPoint(const long coordinates[] = nullptr);

    /**
     *
     * */
    std::vector<long> decryptCoordinates();

    /**
     *
     * */
    std::vector<Point> points;

    [[maybe_unused]] helib::PubKey &getPublicKey() {
        return (helib::PubKey &) encryptionKey;
    }

};


#endif //ENCKMEAN_CLIENT_H
