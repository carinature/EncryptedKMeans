
#ifndef ENCKMEAN_CLIENT_H
#define ENCKMEAN_CLIENT_H

#include "utils/aux.h" // for including KeysServer.h


class Point {
    std::vector<std::vector<helib::Ctxt> > cCoordinates;
    const helib::PubKey &public_key;// = encryptionKey;
    friend class Client;

    friend class KeysServer;

public:

    // Each bit of the binary number is encoded into a single ciphertext. Thus
    // for a 16 bit binary number, we will represent this as an array of 16
    // unique ciphertexts.
    // i.e. b0 = [0] [0] [0] ... [0] [0] [0]        ciphertext for bit 0
    //      b1 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 1
    //      b2 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 2
    // These 3 ciphertexts represent the 3-bit binary number 110b = 6
    // todo why all slots the same? redundant? (KT)
    // Note: several numbers can be encoded across the slots of each ciphertext
    // which would result in several parallel slot-wise operations.
    // For simplicity we place the same data into each slot of each ciphertext,
    // printing out only the back of each vector.
    // NB: fifteenOrLess4Four max is 15 bits. Later in the code we pop the MSB.
    //        long bitSize = 16;
    //        long outSize = 2 * bitSize;
    explicit Point(const helib::PubKey &public_key, const long coordinates[] = nullptr) :
            public_key(public_key),
            cCoordinates(DIM, std::vector(bitSize, helib::Ctxt(public_key))) {
        cout << " Point Init" << endl;
        if (coordinates)
            for (int dim = 0; dim < DIM; ++dim)
                for (long bit = 0; bit < bitSize; ++bit) // Extract the i'th bit of coordinates[dim]
                    this->public_key.Encrypt(cCoordinates[dim][bit],
                                             NTL::to_ZZX((coordinates[dim] >> bit) & 1));
    }

private:
    std::vector<long> decrypt(const helib::EncryptedArray &ea,
                              const helib::PubKey &public_key) {
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
