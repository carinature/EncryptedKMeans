
#ifndef ENCKMEAN_CLIENT_H
#define ENCKMEAN_CLIENT_H

#include <utils/aux.h>

#include "Point.h"

/**
 * @class ClientDevice
 * @brief the ClientDevice can encrypt info using the key from the Keys Server. decryption is currently only for dbg purposes.
 * */
class ClientDevice {
    const helib::SecKey encryptionKey; //reference?
    const helib::EncryptedArray ea;
    const helib::Ctxt scratch;        // Use a scratch ciphertext to populate vectors.

    std::vector<std::vector<long> > pCoordinatesDBG;
    const helib::PubKey *pubKeyPtrDBG;

protected:
    helib::PubKey &public_key;// = encryptionKey;

public:
    /**
     * Constructor for \class{ClientDevice},
     * @brief simulates a mini-protocol between the keys server and the client
     * @param keysServer binds to the \class{KeysServer} responsible for the distributing the appropriate key
     * */
    explicit ClientDevice(const KeysServer &keysServer);

    /**
     *
     * */
    //todo consider returning arr instad of vector (DIM is constant throughout the program)
    std::vector<std::vector<helib::Ctxt> > encryptPoint(const long coordinates[] = nullptr);

    /**
     * @brief add an encrypted point to client's list of data points
     * */
    ClientDevice &addEncryptedPoint(Point &point);

    //  todo remove? doesn't fit buisness logic
    std::vector<long> decryptCoordinate(int i = 0);

    /**
     *
     * */
    std::vector<Point> points;

    /**
     * @return a list of clients data points points
     * */
    const std::vector<Point> &getPoints() const {
        return points;
    }

    /**
     * @return clients public key, with which data can be encrypted, but not decrypted
     * */
    const helib::PubKey &getPublicKey() const {
        return (helib::PubKey &) encryptionKey;
    }

};

#endif //ENCKMEAN_CLIENT_H
