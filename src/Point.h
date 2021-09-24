
#ifndef ENCRYPTEDKMEANS_POINT_H
#define ENCRYPTEDKMEANS_POINT_H

#include <iostream>
#include <helib/helib.h>

#include "utils/aux.h" // for including KeysServer.h


class Point {
    const helib::PubKey &public_key;// = encryptionKey;
    friend class Client;

    friend class KeysServer;

public:
    std::vector<std::vector<helib::Ctxt> > cCoordinates;

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
    explicit Point(const helib::PubKey &public_key, const long coordinates[] = nullptr) :
            public_key(public_key),
            cCoordinates(DIM, std::vector(bitSize, helib::Ctxt(public_key))) {
        cout << " Point Init" << endl;
        if (coordinates)
            for (int dim = 0; dim < DIM; ++dim)
                for (long bit = 0; bit < bitSize; ++bit) // Extract the i'th bit of coordinates[dim]
                    this->public_key.Encrypt(cCoordinates[dim][bit],
                                             NTL::to_ZZX((coordinates[dim] >> bit) & 1));
    }


    std::vector<helib::Ctxt> operator[](short int i) const {
        return cCoordinates[i];
    }


    [[maybe_unused]] bool isEmpty() {
        return cCoordinates[0][0].isEmpty();
    }

    Point operator+(Point &point) {
        //                if (point.isEmpty()) return this;
        //                if (this.isEmpty()) return point;
        const long arr[] = {0, 0};
        Point sum(this->public_key, arr);
        ////        for (int dim = 0; dim < DIM; ++dim) {
        ////            sum[dim] = cCoordinates[dim] + point.cCoordinates[dim];
        ////        }
        //        sum.cCoordinates[0] = cCoordinates[0] + point.cCoordinates[0];
        cout << "whahad" << endl;
        std::vector<Ctxt> temp(bitSize, helib::Ctxt(point.public_key));
//        helib::CtPtrs_vectorCt result_wrapper(                sum.cCoordinates[0]);
                std::vector<helib::Ctxt> encrypted_product;
        helib::CtPtrs_vectorCt result_wrapper(encrypted_product);
        cout << "tnin" << endl;
        // Perform the multiplication first and put it in encrypted_product.
        helib::addTwoNumbers(
                result_wrapper,
                                helib::CtPtrs_vectorCt(point.cCoordinates[0]),
                                helib::CtPtrs_vectorCt(this->cCoordinates[0])
//                outSize    // sizeLimit=0 means use as many bits as needed.
//                &(KeysServer::unpackSlotEncoding) // Information needed for bootstrapping.
        );
        cout << "talate" << endl;
        sum.cCoordinates[0] = encrypted_product;
//        sum.cCoordinates[0] = result_wrapper.v;
        return sum;
        //        return Point(public_key);
    }

};


#endif //ENCRYPTEDKMEANS_POINT_H
