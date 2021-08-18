
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
std::vector<std::vector<helib::Ctxt>> Client::encryptPoint(const long coordinates[]) {
#if VERBOSE
    printNameVal(ea.size());
        cout << "encryptPoint for coordinates: " << endl;
        for (int i = 0; i < DIM; ++i) printNameVal(coordinates[i]);
#endif
    std::vector<long> a_vec(ea.size());
    for (int dim = 0; dim < DIM; ++dim) {
        cCoordinatesStd.emplace_back(bitSize, scratch);
        if (coordinates)
            for (long bit = 0; bit < bitSize; ++bit) {
                // Extract the i'th bit of coordinates[dim].
                for (auto &slot : a_vec) slot = (coordinates[dim] >> bit) & 1;
                ea.encrypt(cCoordinatesStd[dim][bit], encryptionKey, a_vec);
                points.emplace_back(Point(ea, public_key, coordinates));
            }
    }
    return cCoordinatesStd;
}


std::vector<long> Client::decryptCoordinates() {
#if VERBOSE
    cout << "decryptCoordiantes" << endl;
#endif
    clientLogger.log("decryptCoordiantes", log_debug);
    std::vector<std::vector<long>> decrypted_result(DIM);
    std::vector<long> dCoordinates(DIM); // todo redundant, decide which of the two to keep
    for (int i = 0; i < DIM; ++i) {
        helib::decryptBinaryNums(decrypted_result[i],
                                 helib::CtPtrs_vectorCt(cCoordinatesStd[i]),
                                 encryptionKey, ea);
        dCoordinates[i] = decrypted_result[i].back();
    }
#if VERBOSE
    for (auto dec_coordinate :decrypted_result) printNameVal(dec_coordinate.back());
        cout << endl << endl;
#endif
    return dCoordinates;
}
