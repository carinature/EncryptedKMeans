
#ifndef ENCKMEAN_CLIENT_H
#define ENCKMEAN_CLIENT_H

#include "utils/aux.h" // for including KeysServer.h


class Point {
    std::vector<std::vector<helib::Ctxt> > cCoordinates;
    const helib::PubKey &public_key;// = encryptionKey;
    friend class Client;
    friend class KeysServer;

public:
    Point(const helib::EncryptedArray &ea,
          const helib::PubKey &public_key,
          const long coordinates[] = nullptr) :
            public_key(public_key) {
        cout << " Point"<< endl;
        std::vector<long> a_vec(ea.size());
        cout << "ehad" << endl;
        if (coordinates)
            for (int dim = 0; dim < DIM; ++dim) {
            cout << "shtaim" << endl;
            cCoordinates.emplace_back(bitSize, helib::Ctxt(public_key));
            cout << "shalosh" << endl;
            for (long bit = 0; bit < bitSize; ++bit) {
                // Extract the i'th bit of coordinates[dim].
                cout << "arba" << endl;
                for (auto &slot : a_vec) slot = (coordinates[dim] >> bit) & 1;
                cout << "hamesh" << endl;
                if (coordinates) ea.encrypt(cCoordinates[dim][bit], public_key, a_vec);
            }
        }
    }

private:
    std::vector<long> decrypt(const helib::EncryptedArray &ea,
                              const helib::PubKey &public_key){
        std::vector<std::vector<long>> decrypted_result(DIM);
        std::vector<long> dCoordinates(DIM); // todo redundant, decide which of the two to keep
        for (int i = 0; i < DIM; ++i) {
//            helib::SecKey sk;
//            sk.Decrypt()
//            helib::decryptBinaryNums(decrypted_result[i],
//                                     helib::CtPtrs_vectorCt(cCoordinates[i]),
//                                     public_key, ea);
//            dCoordinates[i] = decrypted_result[i].back();
        }
#if VERBOSE
        for (auto dec_coordinate :decrypted_result) printNameVal(dec_coordinate.back());
        cout << endl << endl;
#endif
        return dCoordinates;
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
