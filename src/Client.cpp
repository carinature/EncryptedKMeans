

#include "Client.h"

using std::cout;
using std::endl;

Client::Client(KeysServer &keysServer, const long *coordinates)
        :
#if DBG
        pCoordinates(coordinates),
#endif
        encryptionKey(keysServer.getSecKey()),
        ea(keysServer.getEA()) {
    Logger logger;
    logger.log("Client()");
    cout << "Client()" << endl;
    NTL::Vec<helib::Ctxt> eMax, eMin, enca, encb;
    helib::Ctxt mu(encryptionKey), ni(encryptionKey);
    resize(enca, long(bitSize), mu);
    resize(cCoordinates, long(bitSize), mu);
    logger.log("odin");
    cout << "odin" << endl;
    for (int i = 0; i < dim; ++i)
        for (long b = 0; b < bitSize; b++) {
            logger.log("dva");
            cout << "dva" << endl;
            encryptionKey.Encrypt(enca[b], NTL::ZZX((coordinates[i] >> b) & 1));
            cout << "tree" << endl;
            encryptionKey.Encrypt(cCoordinates[b], NTL::ZZX((coordinates[i] >> b) & 1));
            cout << "chetiri" << endl;
        }
}

