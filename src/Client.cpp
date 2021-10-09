
#include "Client.h"


//helib
#include <helib/binaryArith.h>
//#include <helib/binaryCompare.h>

using std::cout;
using std::endl;
static Logger clientLogger(log_debug, "clientLogger");

Client::Client(const KeysServer &keysServer) :
        encryptionKey(keysServer.getSecKey()),
        //        public_key(encryptionKey),
        public_key(keysServer.getPublicKey()),
        pubKeyPtrDBG(&public_key),
        ea(keysServer.getEA()),
        scratch(encryptionKey)
//        cCoordinatesStd(DIM)
{
    clientLogger.log("Initializing Client Protocol Finished");
    cout << "Initializing Client Protocol Finished" << endl;
}

std::vector<std::vector<helib::Ctxt>> Client::encryptPoint(const long coordinates[]) {
#if VERBOSE
    printNameVal(ea.size());
        cout << "encryptPoint for coordinates: " << endl;
        for (int i = 0; i < DIM; ++i) printNameVal(coordinates[i]);
#endif
    //    std::vector<long> a_vec(ea.size());
//    pCoordinatesDBG.push_back(std::vector<long>(DIM));
//    for (int dim = 0; dim < DIM; ++dim) pCoordinatesDBG.back()[dim] = coordinates[dim];
    points.emplace_back(public_key, coordinates); //be careful when changing to `emplace_back`
    return points.back().cCoordinates;
}

Client &Client::addEncryptedPoint(Point &point) {
    points.emplace_back(point);
    return *this;
}

std::vector<long> Client::decryptCoordinate(int i) {
    cout << "Client::decryptCoordinate" << endl;
    clientLogger.log("decryptCoordiantes", log_debug);
    std::vector<long> dCoordinates(DIM);
    if (points[i][0][0].isEmpty()) return dCoordinates;
    for (int dim = 0; dim < DIM; ++dim) {
        NTL::ZZX pp;
        for (int bit = 0; bit < BIT_SIZE; ++bit) {
            encryptionKey.Decrypt(pp, points[i].cCoordinates[dim][bit]);
            if (IsOne(pp)) dCoordinates[dim] += std::pow(2, bit);
        }
    }
    return dCoordinates;
}

