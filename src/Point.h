
#ifndef ENCRYPTEDKMEANS_POINT_H
#define ENCRYPTEDKMEANS_POINT_H

#include <iostream>
#include <helib/helib.h>
#include <helib/binaryCompare.h>
#include <helib/binaryArith.h>

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
    //        long BIT_SIZE = 16;
    //        long OUT_SIZE = 2 * BIT_SIZE;
    explicit Point(const helib::PubKey &public_key, const long coordinates[] = nullptr) :
            public_key(public_key),
            cCoordinates(DIM, std::vector(BIT_SIZE, helib::Ctxt(public_key))) {
        cout << " Point Init" << endl;
        if (coordinates)
            for (int dim = 0; dim < DIM; ++dim)
                for (long bit = 0;
                     bit < BIT_SIZE; ++bit) // Extract the i'th bit of coordinates[dim]
                    this->public_key.Encrypt(cCoordinates[dim][bit],
                                             NTL::to_ZZX((coordinates[dim] >> bit) & 1));
    }


    const std::vector<helib::Ctxt> &operator[](short int i) const {
        return cCoordinates[i];
    }


    bool isEmpty() const {
        return cCoordinates[0][0].isEmpty();
    }

    Point operator+(Point &point) {
        if (point.isEmpty()) return *this; //todo consider
        if (this->isEmpty()) return point;
        const long arr[] = {0, 0};
        Point sum(this->public_key, arr);
        for (int dim = 0; dim < DIM; ++dim) {
            helib::CtPtrs_vectorCt result_wrapper(sum.cCoordinates[dim]);
            helib::addTwoNumbers(
                    result_wrapper,
                    helib::CtPtrs_vectorCt(point.cCoordinates[dim]),
                    helib::CtPtrs_vectorCt(this->cCoordinates[dim]),
                    OUT_SIZE,   // sizeLimit=0 means use as many bits as needed.
                    &(KeysServer::unpackSlotEncoding) // Information needed for bootstrapping.
            );
        }
        return sum;
    }

    //todo notice that IT IS 3-for-2 (whoohhoo)
    //// Calculates the sum of many numbers using the 3-for-2 method
    static Point addManyPoints(std::vector<Point> points) {
        //        if (points.empty()) return static_cast<Point>(nullptr);
        const long arr[] = {0, 0};
        //        std::vector<std::vector<helib::Ctxt>>> c
        Point sum(points.back().public_key);//, arr);
        std::vector<std::vector<helib::Ctxt>> summands; // = {encrypted_a, encrypted_b, encrypted_c};
        std::vector<helib::Ctxt> encrypted_result;
        for (short dim = 0; dim < DIM; ++dim) {
//            if (80 < points.size()) {
//                cout << "lskdjflksdjflksdjflksfjslkfj" << endl;
//                for (Point &point : points) {
//                    summands.push_back(point.cCoordinates[0]);
//                    summands.push_back(point.cCoordinates[0]);
//                }
//            } else
                for (Point &point : points) summands.push_back(point.cCoordinates[dim]);
            helib::CtPtrMat_vectorCt summands_wrapper(summands);
            //            helib::CtPtrs_vectorCt result_wrapper(sum.cCoordinates[dim]);
            helib::CtPtrs_vectorCt result_wrapper(encrypted_result);
            // Calculates the sum of many numbers using the 3-for-2 method
            addManyNumbers(
                    result_wrapper,
                    summands_wrapper
                    //                    ,
                    //                    points.size() * BIT_SIZE * 2, // sizeLimit=0 means use as many bits as needed.
                    //                    &(KeysServer::unpackSlotEncoding) // Information needed for bootstrapping.
            );
            sum.cCoordinates[dim] = encrypted_result;
            encrypted_result.clear();
            summands.clear();
        }
        return sum;
    }

    /*
     * multiply  an encrypted point by an encrypted bit
     * */
    Point operator*(Ctxt &bit) const {
        if (bit.isEmpty()) return *this; //todo consider
        if (this->isEmpty()) return Point(this->public_key);
        //        const long arr[] = {0, 0};
        Point product(*this);
        for (int dim = 0; dim < DIM; ++dim) {
            helib::CtPtrs_vectorCt result_wrapper(product.cCoordinates[dim]);
            //                vecCopy(product.cCoordinates[dim], this->cCoordinates[dim], BIT_SIZE);
            binaryMask(result_wrapper, bit);
            //                for (long i = 0; i < resSize; i++)
            //                    productCoors[i]->multiplyBy(*(lhs[0]));
            //                return;
        }
        return product;

        /*  todo notice this two functions from @file binaryArith.cpp :
                //! Apply mask across the vector of bits slot-wise.
                void binaryMask(CtPtrs& bits, const Ctxt& mask)        {
                    for (long i = 0; i < bits.size(); ++i)
                        bits[i]->multiplyBy(mask);
                }
                //! Implementation of output = cond ? trueValue : falseValue
                void binaryCond(CtPtrs& output,
                                const Ctxt& cond,
                                const CtPtrs& trueValue,
                                const CtPtrs& falseValue)
                {...}
         */
    }

    Point operator*(Point &point) {
        if (point.isEmpty()) return *this; //todo consider
        if (this->isEmpty()) return point;
        const long arr[] = {0, 0};
        Point product(this->public_key, arr);
        for (int dim = 0; dim < DIM; ++dim) {
            helib::CtPtrs_vectorCt result_wrapper(product.cCoordinates[dim]);
            helib::multTwoNumbers(
                    result_wrapper,
                    helib::CtPtrs_vectorCt(point.cCoordinates[dim]),
                    helib::CtPtrs_vectorCt(this->cCoordinates[dim]),
                    false,
                    OUT_SIZE,   // sizeLimit=0 means use as many bits as needed.
                    &(KeysServer::unpackSlotEncoding) // Information needed for bootstrapping.
            );
        }
        return product;
    }

    helib::Ctxt isBiggerThan(Point &p, short currentDim = DIM - 1) {
        Ctxt mu(public_key), ni(public_key);
        compareTwoNumbers(mu,
                          ni,
                          helib::CtPtrs_vectorCt(this->cCoordinates[currentDim]),
                          helib::CtPtrs_vectorCt(p.cCoordinates[currentDim]),
                          false,
                          &KeysServer::unpackSlotEncoding
        );
        return mu;
        /* todo notice there is also // comparison with max and min
         *  maybe useful later
         *      compareTwoNumbers(wMax,
                      wMin,
                      mu,
                      ni,
                      helib::CtPtrs_VecCt(enca),
                      helib::CtPtrs_VecCt(encb),
                      false,
                      &unpackSlotEncoding);
                      */
    }


};


#endif //ENCRYPTEDKMEANS_POINT_H
