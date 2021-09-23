
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
//    for (int dim = 0; dim < DIM; ++dim) {
//        cCoordinatesStd.emplace_back(bitSize, scratch);
//        if (coordinates)
//            for (long bit = 0; bit < bitSize; ++bit) {
//                // Extract the i'th bit of coordinates[dim].
//                for (auto &slot : a_vec) slot = (coordinates[dim] >> bit) & 1;
//                ea.encrypt(cCoordinatesStd[dim][bit], encryptionKey, a_vec);
//            }
//    }
    points.emplace_back(Point(public_key, coordinates));
    return points[0].cCoordinates;
}


std::vector<long> Client::decryptCoordinates() {
    cout << "Client::decryptCoordinates" << endl;
    clientLogger.log("decryptCoordiantes", log_debug);
    std::vector<long> dCoordinates(DIM);
//    if (!cCoordinatesStd[0][0].isEmpty())
    for (int dim = 0; dim < DIM; ++dim) {
        NTL::ZZX pp;
        for (int bit = 0; bit < bitSize; ++bit) {
            encryptionKey.Decrypt(pp, points[0].cCoordinates[dim][bit]);
//            printNameVal(pp);
            if (IsOne(pp)) dCoordinates[dim]+=std::pow(2, bit);
//            printNameVal(dCoordinates[dim]);
        }
    }
    return dCoordinates;
}
