
#include "Client.h"

//helib
#include <helib/binaryArith.h>
#include <helib/binaryCompare.h>

using std::cout;
using std::endl;
//Logger clientLogger;

//Client::Client(KeysServer &keysServer, const long coordinates) :
//        encryptionKey(keysServer.getSecKey()),
//        ea(keysServer.getEA()) {
//    clientLogger.log("Client()");
//    cout << "Client()" << endl;
//    NTL::Vec<helib::Ctxt> eMax, eMin, enca, encb;
//    helib::Ctxt mu(encryptionKey), ni(encryptionKey);
//    resize(enca, long(bitSize), mu);
//    resize(cCoordinatesNTL, long(bitSize), mu);
//    resize(cCoordinatesStd, long(bitSize), mu);
//}
//
//std::vector<helib::Ctxt> Client::encryptPoint(const long coordinates[]) {
//    // Each bit of the binary number is encoded into a single ciphertext. Thus
//    // for a 16 bit binary number, we will represent this as an array of 16
//    // unique ciphertexts.
//    // i.e. b0 = [0] [0] [0] ... [0] [0] [0]        ciphertext for bit 0
//    //      b1 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 1
//    //      b2 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 2
//    // These 3 ciphertexts represent the 3-bit binary number 110b = 6
//
//    // Note: several numbers can be encoded across the slots of each ciphertext
//    // which would result in several parallel slot-wise operations.
//    // For simplicity we place the same data into each slot of each ciphertext,
//    // printing out only the back of each vector.
//    // NB: fifteenOrLess4Four max is 15 bits. Later in the code we pop the MSB.
//    long bitSize = 16;
//    long outSize = 2 * bitSize;
//    //    long a_data = NTL::RandomBits_long(bitSize);
//    //    long b_data = NTL::RandomBits_long(bitSize);
//    //    long c_data = NTL::RandomBits_long(bitSize);
//
//    // Use a scratch ciphertext to populate vectors.
//    const helib::PubKey &public_key = encryptionKey;
//    //    const helib::PubKey &public_key = keysServer.getSecKey();
//    helib::Ctxt scratch(public_key);
//    std::vector<helib::Ctxt> encrypted_a(bitSize, scratch);
//    //    cout << "odin" << endl;
//    // Encrypt the data in binary representation.
//    int DIM = 1;
//    for (int c = 0; c < DIM; ++c)
//        for (long i = 0; i < bitSize; ++i) {
//            std::vector<long> a_vec(ea.size());
//            printNameVal(ea.size());
//            // Extract the i'th bit of a,b,c.
//            for (auto &slot : a_vec) slot = (coordinates[c] >> i) & 1;
//            ea.encrypt(encrypted_a[i], public_key, a_vec);
//        }
//    for (int i = 0; i < DIM; ++i) {
//        std::vector<long> slots;
//        ea.decrypt(encrypted_a[i], encryptionKey, slots);
//        for (auto s : slots) printNameVal(s);
//    }
//    return encrypted_a;
//}

//helib::Ctxt Client::compare(Client &client) {
//    cout << "compare(Client &client) " << endl;
//    helib::Ctxt isBiggerFlag(encryptionKey), isSmallerFlag(encryptionKey);
//    cout << "odin" << endl;
//    compareTwoNumbers(isBiggerFlag, isSmallerFlag,
//                      helib::CtPtrs_VecCt(this->cCoordinatesNTL),
//                      helib::CtPtrs_VecCt(client.cCoordinatesNTL),
//                      false, nullptr);
//    cout << "onaes" << endl;
//    return isBiggerFlag;
//}
//
//void Client::decryptCoordinates() {
//    cout << "decryptCoordiantes" << endl;
//    //    decryptAndPrint((cout << " before comparison: "),
//    //                    this->cCoordinatesNTL[0], encryptionKey, ea, 0);
//    std::vector<long> slots;
//    ea.decrypt(this->cCoordinatesNTL[0], encryptionKey, slots);
//    for (auto slot:slots) printNameVal(slot);
//    cout << "-------" << endl;
//    ea.decrypt(this->cCoordinatesNTL[1], encryptionKey, slots);
//    for (auto slot:slots) printNameVal(slot);
//
//    cout << endl << endl;
//}
