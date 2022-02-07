
#include "ClientDevice.h"

using std::cout;
using std::endl;
static Logger loggerClient(log_debug, "loggerClient");

ClientDevice::ClientDevice(const KeysServer &keysServer) :
        encryptionKey(keysServer.getSecKey()),
        //        public_key(encryptionKey),
        public_key(keysServer.getPublicKey()),
        pubKeyPtrDBG(&public_key),
        ea(keysServer.getEA()),
        scratch(encryptionKey)
//        cCoordinatesStd(DIM)
{
    loggerClient.log("Initializing ClientDevice Protocol Finished");
    cout << "Initializing ClientDevice Protocol Finished" << endl;
}

std::vector<std::vector<helib::Ctxt> > ClientDevice::encryptPoint(const long coordinates[]) {
#if VERBOSE
    printNameVal(ea.size());
        cout << "encryptPoint for coordinates: " << endl;
        for (int i = 0; i < DIM; ++i) printNameVal(coordinates[i]);
#endif
    //    std::vector<long> a_vec(ea.size());
//    pCoordinatesDBG.push_back(std::vector<long>(DIM));
//    for (short dim = 0; dim < DIM; ++dim) pCoordinatesDBG.back()[dim] = coordinates[dim];
    points.emplace_back(public_key, coordinates); //be careful when changing to `emplace_back`
    return points.back().cCoordinates;
}

ClientDevice &ClientDevice::addEncryptedPoint(Point &point) {
    points.emplace_back(point);
    return *this;
}

std::vector<long> ClientDevice::decryptCoordinate(int i) {
    cout << "ClientDevice::decryptCoordinate" << endl;
    loggerClient.log("decryptCoordiantes", log_debug);
    std::vector<long> dCoordinates(DIM);
    if (points[i][0][0].isEmpty()) return dCoordinates;
    for (short dim = 0; dim < DIM; ++dim) {
        NTL::ZZX pp;
        for (int bit = 0; bit < BIT_SIZE; ++bit) {
            encryptionKey.Decrypt(pp, points[i].cCoordinates[dim][bit]);
            if (IsOne(pp)) dCoordinates[dim] += std::pow(2, bit);
        }
    }
    return dCoordinates;
}

