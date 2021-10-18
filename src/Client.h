
#ifndef ENCKMEAN_CLIENT_H
#define ENCKMEAN_CLIENT_H

#include "Point.h"

/**
 * @class Client
 * @brief the Client can encrypt info using the key from the Keys Server. decryption is currently only for dbg purposes.
 * */
class Client {
    const helib::SecKey encryptionKey; //reference?
    const helib::EncryptedArray ea;
    const helib::Ctxt scratch;        // Use a scratch ciphertext to populate vectors.

    /*  todo for DBG. remove */
    std::vector<std::vector<long> > pCoordinatesDBG;
    const helib::PubKey *pubKeyPtrDBG;

protected:
    //    NTL::Vec<helib::Ctxt> cCoordinatesNTL;
    //    std::vector<std::vector<helib::Ctxt> > cCoordinatesStd;
    const helib::PubKey &public_key;// = encryptionKey;
//    helib::PubKey pubKey;// = encryptionKey;
    //#if DBG
    //private:
    //    long *pCoordinatesDBG;
    //#endif


public:
    /**
     * Constructor for \class{Client},
     * @param keysServer binds to the \class{KeysServer} responsible for the distributing the appropriate key
     * @brief simulates a mini-protocol between the keys server and the client
     * */
    explicit Client(const KeysServer &keysServer);

    /**
     *
     * */
    //todo consider returning arr instad of vector (DIM is constant throughout the program)
    std::vector<std::vector<helib::Ctxt>> encryptPoint(const long coordinates[] = nullptr);

    Client &addEncryptedPoint(Point &point);

    /**
     *
     * */
    std::vector<long> decryptCoordinate(int i = 0);

    /**
     *
     * */
    std::vector<Point> points;

    const std::vector<Point> &getPoints() const{
        return points;
    }


    const helib::PubKey &getPublicKey() const{
        return (helib::PubKey &) encryptionKey;
    }

};


#endif //ENCKMEAN_CLIENT_H
