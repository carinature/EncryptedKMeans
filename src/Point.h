
#ifndef ENCRYPTEDKMEANS_POINT_H
#define ENCRYPTEDKMEANS_POINT_H

#include <iostream>
#include <helib/helib.h>
#include <helib/binaryCompare.h>
#include <helib/binaryArith.h>

#include "utils/aux.h" // for including KeysServer.h
#include "KeysServer.h"

static long counter = 0; //fixme move to the cpp file

class Point {
    friend class Client;

    friend class KeysServer;

    /*  for DBG */
    std::vector<long> pCoordinatesDBG;
    bool isCopyDBG = false;
    bool isEmptyDBG = true;
    const helib::PubKey *pubKeyPtrDBG;
    long cmpCounter, addCounter, multCounter; //todo

public:
    const helib::PubKey &public_key;// = encryptionKey;
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
            cmpCounter(0),  //todo maybe better to init to 0, depending future impl & use
            addCounter(0),  //todo maybe better to init to 0, depending future impl & use
            multCounter(0), //todo maybe better to init to 0, depending future impl & use
    //            public_key(public_key){
            id(counter++),
            public_key(public_key),
            pubKeyPtrDBG(&public_key),
            cCoordinates(DIM, std::vector(BIT_SIZE, helib::Ctxt(this->public_key))) {
        //        this->cCoordinates.reserve(DIM);
        pCoordinatesDBG.reserve(DIM);
        //        cout << " Point Init" << endl;
        if (coordinates) {
            isEmptyDBG = false;
            for (short dim = 0; dim < DIM; ++dim) {
                pCoordinatesDBG.push_back(coordinates[dim]);
                // Extract the i'th bit of coordinates[dim]
                for (long bit = 0; bit < BIT_SIZE; ++bit)
                    this->public_key.Encrypt(cCoordinates[dim][bit],
                                             NTL::to_ZZX((coordinates[dim] >> bit) & 1));
            }
        }
    }

    bool isEmpty() const {
        return cCoordinates[0][0].isEmpty();
    }

    Point(const Point &point) :
            cmpCounter(point.cmpCounter),  //todo maybe better to init to 0, depending future impl & use
            addCounter(point.addCounter),  //todo maybe better to init to 0, depending future impl & use
            multCounter(point.multCounter), //todo maybe better to init to 0, depending future impl & use
            id(point.id), //todo make sure this approach won't cuase some unpredictable behaviour later (e.g in cases of EXPLICIT copy (instead of implicit, in which this way makes sense))
            public_key(point.public_key),
            pubKeyPtrDBG(&(point.public_key)),
            //            cCoordinates(DIM, std::vector(BIT_SIZE, helib::Ctxt(public_key))),
            cCoordinates(point.cCoordinates),
            pCoordinatesDBG(point.pCoordinatesDBG),
            isEmptyDBG(point.isEmptyDBG),
            isCopyDBG(true) {
        //        cout << " Point copy copy" << endl; //this print is important for later. efficiency...
        if (!point.isEmpty())
            for (short dim = 0; dim < DIM; ++dim) cCoordinates[dim] = point.cCoordinates[dim];

    }

    //  todo need?
    //    Point operator=(Point &point) {
    //        cout << " Point assign" << endl;
    //        public_key. = point.public_key;
    //        if (point.isEmpty()) return
    //            for (short dim = 0; dim < DIM; ++dim)
    //                cCoordinates[dim] = point.cCoordinates[dim];
    //    }


    const std::vector<helib::Ctxt> &operator[](short int i) const {
        //        if (isEmpty()) return std::vector<helib::Ctxt>(helib::Ctxt(public_key));
        return cCoordinates[i];
    }

    std::vector<helib::Ctxt> operator[](short int i) {
        //        if (isEmpty()) return std::vector<helib::Ctxt>(helib::Ctxt(public_key));
        return cCoordinates[i];
    }

    Point operator+(Point &point) {
        if (point.isEmpty()) return *this; //todo consider
        if (this->isEmpty()) return point;
        const long arr[] = {0, 0};
        Point sum(this->public_key, arr);
        for (short dim = 0; dim < DIM; ++dim) {
            helib::CtPtrs_vectorCt result_wrapper(sum.cCoordinates[dim]);
            /*
 * @brief Adds two numbers in binary representation where each ciphertext of the
 * input vector contains a bit.
 * @param sum result of the addition operation.
 * @param lhs left hand side of the addition.
 * @param rhs right hand side of the addition.
 * @param sizeLimit number of bits to compute on, taken from the least
 * significant end.
 * @param unpackSlotEncoding vector of constants for unpacking, as used in
 * bootstrapping.
 **/
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
    static Point addManyPoints(std::vector<Point> &points) {
        //        if (points.empty()) return static_cast<Point>(nullptr);
        //        const long arr[] = {0, 0};
        //        std::vector<std::vector<helib::Ctxt>>> c
        Point sum(points.back().public_key);//, arr);
        std::vector<std::vector<helib::Ctxt>> summands; // = {encrypted_a, encrypted_b, encrypted_c};
        for (short dim = 0; dim < DIM; ++dim) {
            for (Point &point : points) summands.push_back(point.cCoordinates[dim]);
            helib::CtPtrMat_vectorCt summands_wrapper(summands);
            helib::CtPtrs_vectorCt result_wrapper(sum.cCoordinates[dim]);
            //            std::vector<helib::Ctxt> encrypted_result;
            //            helib::CtPtrs_vectorCt result_wrapper(encrypted_result);
            /*
 * @brief Sum an arbitrary amount of numbers in binary representation.
 * @param sum result of the summation.
 * @param numbers values of which to sum.
 * @param sizeLimit number of bits to compute on, taken from the least
 * significant end.
 * @param unpackSlotEncoding vector of constants for unpacking, as used in
 * bootstrapping.
 *
 * Calculates the sum of many numbers using the 3-for-2 method.
 **/
            // Calculates the sum of many numbers using the 3-for-2 method
            addManyNumbers(
                    result_wrapper,
                    summands_wrapper,
                    2 * points.size() * BIT_SIZE, // sizeLimit=0 means use as many bits as needed.
                    &(KeysServer::unpackSlotEncoding) // Information needed for bootstrapping.
            );
            //            sum.cCoordinates[dim] = encrypted_result;
            //            encrypted_result.clear();
            summands.clear();
        }
        return sum;
    }

    /**
     * @brief Multiplies an encrypted point by an encrypted bit
     * */
    Point operator*(Ctxt &bit) const {
        if (bit.isEmpty()) return *this; //todo consider
        if (this->isEmpty()) return Point(this->public_key);
        //        const long arr[] = {0, 0};
        Point product(*this);
        for (short dim = 0; dim < DIM; ++dim) {
            //todo use operator[]
            helib::CtPtrs_vectorCt result_wrapper(product.cCoordinates[dim]);
            binaryMask(result_wrapper, bit);
            //                for (long i = 0; i < resSize; i++)
            //                    productCoors[i]->multiplyBy(*(lhs[0]));
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
        for (short dim = 0; dim < DIM; ++dim) {
            //todo use operator[]
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

    //    friend std::ostream& operator<<(std::ostream& os, const Point& dt);
    //    std::ostream& operator<<(std::ostream& os, const Point& p){
    //        os << "( ";
    //        for (short dim = 0; dim < DIM - 1; ++dim)
    //            os << keysServer.decryptNum(p[dim]) << ",";
    //        os << keysServer.decryptNum(p[short(DIM) - 1]) << " ) " << endl;
    //
    //        //    os << dt.mo << '/' << dt.da << '/' << dt.yr;
    //        return os;
    //    }


    /**
     * @brief compares 2 points by comparing the values of 2 coordinates, in a specified dimention.
     * @param point - a point with encrypted coordinate values.
     * @param currentDim - the index of the coordinates to be compared.
     * @returns a tuple that answers - ((p1[d]>p2[d]), (p2[d]>p1[d])). value are encrypted.
     * @return std::vector<helib::Ctxt>
     * */
    const std::vector<helib::Ctxt> isBiggerThan(const Point &point, short currentDim = DIM - 1) const {
        Ctxt mu(public_key), ni(public_key);
        if (!(isEmpty() || point.isEmpty())) {
            if (point.id == id) {
                // good solution but results in representatienves being picked twice
                //  - once for their own group and once for the group above fixme
                public_key.Encrypt(mu, NTL::to_ZZX((1)));
                public_key.Encrypt(ni, NTL::to_ZZX((1))); // make sure
            } else {
                std::vector<helib::Ctxt> c = (*this)[currentDim];
                std::vector<helib::Ctxt> coor = point[currentDim];
                compareTwoNumbers(mu,
                                  ni,
                                  helib::CtPtrs_vectorCt(c),
                                  helib::CtPtrs_vectorCt(coor),
                                  false,
                                  &KeysServer::unpackSlotEncoding
                );
            }
        }
        return std::vector<helib::Ctxt>{mu, ni};
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

    bool operator==(const Point &point) const {
        return id == point.id;
    }

public:
    //! @var long id
    //! used in createCmpDict for comparison
    const long id; // = 0;  //fixme
};


struct cmpPoints {
    long id;
    bool operator()(const Point &a, const Point &b) const {
        return a.id > b.id;
    }
    bool operator==(const Point &p) const {
        return id == p.id;
    }
};

struct hashPoints {
    std::size_t operator()(const Point &point) const {
        return std::hash<long>()(point.id);
    }
//    bool operator==(const Point &p) const {
//        return id == p.id;
//    }
};

namespace std {

    template <>
    struct hash<const Point>
    {
        std::size_t operator()(const Point& point) const
        {
            using std::size_t;
            using std::hash;

            return hash<long>()(point.id);
//            return point.id;
        }
    };

}


#endif //ENCRYPTEDKMEANS_POINT_H
