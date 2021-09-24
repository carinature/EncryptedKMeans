
#include "Client.h"


//helib
#include <helib/binaryArith.h>
//#include <helib/binaryCompare.h>

using std::cout;
using std::endl;
static Logger clientLogger(log_debug, "clientLogger");

Client::Client(KeysServer &keysServer) :
        encryptionKey(keysServer.getSecKey()),
        public_key(encryptionKey),
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
    std::vector<long> a_vec(ea.size());
    points.emplace_back(Point(public_key, coordinates));
    return points[0].cCoordinates;
}

Client &Client::addEncryptedPoint(Point &point) {
    points.emplace_back(point);
    return *this;
}

std::vector<long> Client::decryptCoordinates(int i) {
    cout << "Client::decryptCoordinates" << endl;
    clientLogger.log("decryptCoordiantes", log_debug);
    std::vector<long> dCoordinates(DIM);
    if (points[i][0][0].isEmpty()) return dCoordinates;
    for (int dim = 0; dim < DIM; ++dim) {
        NTL::ZZX pp;
        for (int bit = 0; bit < BIT_SIZE; ++bit) {
            encryptionKey.Decrypt(pp, points[i].cCoordinates[dim][bit]);
//            printNameVal(pp);
            if (IsOne(pp)) dCoordinates[dim] += std::pow(2, bit);
//            printNameVal(dCoordinates[dim]);
        }
    }
    return dCoordinates;
}

