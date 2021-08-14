

#ifndef ENCKMEAN_CLIENT_H
#define ENCKMEAN_CLIENT_H


#include "utils/aux.h"

//fixme remove //todo move to h_file
using std::cout;
using std::endl;
static Logger clientLogger;

#include <helib/binaryCompare.h>
#include <helib/debugging.h>

class Client {
    const helib::SecKey encryptionKey;
    const helib::EncryptedArray ea;
    //    // Use a scratch ciphertext to populate vectors.
    const helib::Ctxt scratch;//(public_key);

protected:
    NTL::Vec<helib::Ctxt> cCoordinatesNTL;
    std::vector<std::vector<helib::Ctxt> > cCoordinatesStd;

    const helib::PubKey &public_key;// = encryptionKey;

private:
#if DBG
    long *pCoordinates; // for dbg
#endif
    //
public:
    //    explicit Client(KeysServer &keysServer, const long coordinates[] = nullptr);
    explicit Client(KeysServer &keysServer, const long coordinates = 0) :
            encryptionKey(keysServer.getSecKey()),
            public_key(encryptionKey),
            ea(keysServer.getEA()),
            scratch(encryptionKey)
    //            cCoordinatesStd(DIM)
    {
        clientLogger.log("Client()");
        cout << "Client()" << endl;
    }

    [[maybe_unused]] helib::PubKey getPublicKey() {
        return encryptionKey;
    }

    std::vector<std::vector<helib::Ctxt>> encryptPoint(const long coordinates[]) {
        // Each bit of the binary number is encoded into a single ciphertext. Thus
        // for a 16 bit binary number, we will represent this as an array of 16
        // unique ciphertexts.
        // i.e. b0 = [0] [0] [0] ... [0] [0] [0]        ciphertext for bit 0
        //      b1 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 1
        //      b2 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 2
        // These 3 ciphertexts represent the 3-bit binary number 110b = 6

        // Note: several numbers can be encoded across the slots of each ciphertext
        // which would result in several parallel slot-wise operations.
        // For simplicity we place the same data into each slot of each ciphertext,
        // printing out only the back of each vector.
        // NB: fifteenOrLess4Four max is 15 bits. Later in the code we pop the MSB.
        //        long bitSize = 16;
        long outSize = 2 * bitSize;
        // Use a scratch ciphertext to populate vectors.
        const helib::PubKey &public_key = encryptionKey;
        for (int dimi = 0; dimi < DIM; ++dimi) {
            std::vector<long> a_vec(ea.size());
            printNameVal(coordinates[dimi]);
            for (long biti = 0; biti < bitSize; ++biti) {
                printNameVal(ea.size());
                // Extract the biti'th bit of a,b,dimi.
                for (auto &slot : a_vec) slot = (coordinates[dimi] >> biti) & 1;
                //                std::vector<helib::Ctxt> ctxt(bitSize, scratch);
                cCoordinatesStd.emplace_back(bitSize, scratch);
                ea.encrypt(cCoordinatesStd[dimi][biti], encryptionKey, a_vec);
            }
        }
        cout << endl << "----------decrupt--------" << endl;
        for (int i = 0; i < DIM; ++i) {
            cout << endl;
            printNameVal(i);
            std::vector<long> slots;
            for (int j = 0; j < bitSize; ++j)
                ea.decrypt(cCoordinatesStd[i][j], encryptionKey, slots);
            for (auto s : slots) printNameVal(s);
            std::vector<long> decrypted_result;
            helib::CtPtrs_vectorCt resultWrapper(cCoordinatesStd[i]);
            helib::decryptBinaryNums(decrypted_result, resultWrapper, encryptionKey, ea);
            std::cout << "a*b+c = " << decrypted_result.back() << std::endl;

        }
        return cCoordinatesStd;
    }

    helib::Ctxt compare(Client &client) {
        cout << "compare(Client &client) " << endl;
        helib::Ctxt isBiggerFlag(encryptionKey), isSmallerFlag(encryptionKey);
        cout << "odin" << endl;
        compareTwoNumbers(isBiggerFlag, isSmallerFlag,
                          helib::CtPtrs_VecCt(this->cCoordinatesNTL),
                          helib::CtPtrs_VecCt(client.cCoordinatesNTL),
                          false, nullptr);
        cout << "onaes" << endl;
        return isBiggerFlag;
    }


    void decryptCoordinates() {
        cout << "decryptCoordiantes" << endl;
        clientLogger.log("decryptCoordiantes", log_debug);
        //        std::vector<long> slots;
        //        ea.decrypt(this->cCoordinatesStd[0], encryptionKey, slots);
        //        for (auto slot:slots) printNameVal(slot);
        //        cout << "-------" << endl;
        //        ea.decrypt(this->cCoordinatesStd[1], encryptionKey, slots);
        //        for (auto slot:slots) printNameVal(slot);

        cout << endl << endl;
    }


};


#endif //ENCKMEAN_CLIENT_H
