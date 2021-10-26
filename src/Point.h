
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

public:

    //! @var long id
    //! used in createCmpDict for comparison
    const long id;

    //     helib::PubKey &public_key;// = encryptionKey;
    const helib::PubKey &public_key;// = encryptionKey;
    std::vector<std::vector<helib::Ctxt> > cCoordinates;

    /*
     Each bit of the binary number is encoded into a single ciphertext. Thus
     for a 16 bit binary number, we will represent this as an array of 16
     unique ciphertexts.
     i.e. b0 = [0] [0] [0] ... [0] [0] [0]        ciphertext for bit 0
          b1 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 1
          b2 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 2
     These 3 ciphertexts represent the 3-bit binary number 110b = 6
     todo why all slots the same? redundant? (KT)
     Note: several numbers can be encoded across the slots of each ciphertext
     which would result in several parallel slot-wise operations.
     For simplicity we place the same data into each slot of each ciphertext,
     printing out only the back of each vector.
     NB: fifteenOrLess4Four max is 15 bits. Later in the code we pop the MSB.
            long BIT_SIZE = 16;
            long OUT_SIZE = 2 * BIT_SIZE;
     */
    explicit Point(const helib::PubKey &public_key, const long coordinates[] = nullptr) :
            cmpCounter(0),  //todo maybe better to init to 0, depending future impl & use
            addCounter(0),  //todo maybe better to init to 0, depending future impl & use
            multCounter(0), //todo maybe better to init to 0, depending future impl & use
            id(counter++),
            public_key(public_key),
            pubKeyPtrDBG(&public_key),
            cCoordinates(DIM, std::vector(BIT_SIZE, helib::Ctxt(public_key))) {
        //        this->public_key = public_key.;
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

        //        todo notice this helibs function:  (CT 24.oct.2021)
        //        /**
        //         * @brief Returns a number as a vector of bits with LSB on the left.
        //         * @param num Number to be converted.
        //         * @param bitSize Number of bits of the input and output.
        //         * @return Bit vector representation of num.
        //         * @note `bitSize` must be non-negative.
        //         **/
        //        std::vector<long> longToBitVector(long num, long bitSize);
    }

    explicit Point(const std::vector<std::vector<helib::Ctxt> > &cCoordinates) :
            cmpCounter(0),  //todo maybe better to init to 0, depending future impl & use
            addCounter(0),  //todo maybe better to init to 0, depending future impl & use
            multCounter(0), //todo maybe better to init to 0, depending future impl & use
            id(counter++),
            public_key(cCoordinates[0][0].getPubKey()),
            pubKeyPtrDBG(&public_key),
            pCoordinatesDBG(DIM),
            cCoordinates(cCoordinates)
    //  cCoordinates(DIM, std::vector(BIT_SIZE, helib::Ctxt(public_key)))
    {
        //        public_key.keyExists()
        //        this->pCoordinatesDBG.reserve(DIM);
        //        this->cCoordinates.reserve(DIM);
        //        for (int dim = 0; dim < DIM; ++dim) {
        ////            this->cCoordinates.push_back(cCoordinates[dim]);
        //            vecCopy(this->cCoordinates[dim], cCoordinates[dim]);
        //            //            this->cCoordinates[dim] = cCoordinates[dim];
        //        }
        //        cout << "Point c'tor from encrypted data" << endl;

    }

    bool isEmpty() const {
        return cCoordinates[0][0].isEmpty();
    }

    Point(const Point &point) :
    //todo maybe better to init to 0, depending future impl & use
            cmpCounter(point.cmpCounter),
            addCounter(point.addCounter),
            multCounter(point.multCounter),
            //todo make sure this approach won't cuase some unpredictable behaviour later
            // (e.g in cases of EXPLICIT copy (instead of implicit, in which this way makes sense))
            id(point.id),
            public_key(point.public_key),
            pubKeyPtrDBG(&(point.public_key)),
            cCoordinates(point.cCoordinates),
            pCoordinatesDBG(point.pCoordinatesDBG),
            isEmptyDBG(point.isEmptyDBG),
            isCopyDBG(true) {
        //        cout << " Point copy copy" << endl; //this print is important for later. efficiency...
        if (!point.isEmpty())
            for (short dim = 0; dim < DIM; ++dim) {
                vecCopy(cCoordinates[dim], point.cCoordinates[dim]); //helibs version of vec copy
                //                cCoordinates[dim] = point.cCoordinates[dim];
            }

    }

    Point &operator=(const Point &point) {
        //            cout << " Point assign" << endl;
        if (&point == this || point.isEmpty()) return *this;
        for (short dim = 0; dim < DIM; ++dim) {
            //   cCoordinates[dim] = point.cCoordinates[dim];
            vecCopy(cCoordinates[dim], point.cCoordinates[dim]); //helibs version of vec copy
            pCoordinatesDBG[dim] = point.pCoordinatesDBG[dim];
        }
        return *this;
    }

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
        Point sum(this->public_key);
        for (short dim = 0; dim < DIM; ++dim) {
            helib::CtPtrs_vectorCt result_wrapper(sum.cCoordinates[dim]);
            //             * @brief Adds two numbers in binary representation where each ciphertext of the
            //             * input vector contains a bit.
            //             * @param sum result of the addition operation.
            //             * @param lhs left hand side of the addition.
            //             * @param rhs right hand side of the addition.
            //             * @param sizeLimit number of bits to compute on, taken from the least
            //             * significant end.
            //             * @param unpackSlotEncoding vector of constants for unpacking, as used in
            //             * bootstrapping.
            //             *
            helib::addTwoNumbers(
                    result_wrapper,
                    helib::CtPtrs_vectorCt(point.cCoordinates[dim]),
                    helib::CtPtrs_vectorCt(this->cCoordinates[dim]),
                    OUT_SIZE,   // sizeLimit=0 means use as many bits as needed.
                    &(KeysServer::unpackSlotEncoding) // Information needed for bootstrapping.
            );
            sum.pCoordinatesDBG[dim] += point.pCoordinatesDBG[dim];
        }
        return sum;
    }

    //// Calculates the sum of many numbers using the 3-for-2 method
    static Point addManyPoints(const std::vector<Point> &points) {
        //        if (points.empty()) return static_cast<Point>(nullptr);
        //        Point sum(points.back().public_key);
        std::vector<std::vector<helib::Ctxt> > summandsVec[DIM];
        std::vector<std::vector<helib::Ctxt> > encrypted_results(DIM);

        for (short dim = 0; dim < DIM; ++dim) {
            summandsVec[dim].reserve(points.size());
            for (const Point &point : points) {
                std::vector<helib::Ctxt> copyVec;
                summandsVec[dim].push_back(point.cCoordinates[dim]);
                //                sum.pCoordinatesDBG[dim] += point.pCoordinatesDBG[dim];
            }
            helib::CtPtrMat_vectorCt summands_wrapper(summandsVec[dim]);

            helib::CtPtrs_vectorCt result_wrapper(encrypted_results[dim]);
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
                    BIT_SIZE * points.size() *
                    BIT_SIZE, // sizeLimit=0 means use as many bits as needed.
                    &(KeysServer::unpackSlotEncoding) // Information needed for bootstrapping.
            );
            //            sum.cCoordinates[dim] = encrypted_result;
            //            vecCopy(sum.cCoordinates[dim], encrypted_results[dim]);
        }
        Point sum(encrypted_results);

        for (short dim = 0; dim < DIM; ++dim)
            for (const Point &point : points)
                sum.pCoordinatesDBG[dim] += point.pCoordinatesDBG[dim];

        return sum;
    }

    /**
     * @brief Multiplies an encrypted point by an encrypted bit
     * */
    const Point operator*(const Ctxt &bit) const {
        //        if (bit.isEmpty()) return *this; //todo consider
        //        if (this->isEmpty()) return Point(this->public_key);
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

    Point &operator*=(const Ctxt &bit) {
        //        if (bit.isEmpty()) return *this; //todo consider
        //        if (this->isEmpty()) return Point(this->public_key);
        Point product(*this);
        for (short dim = 0; dim < DIM; ++dim) {
            //todo use operator[]
            helib::CtPtrs_vectorCt result_wrapper(this->cCoordinates[dim]);
            binaryMask(result_wrapper, bit);
            //                for (long i = 0; i < resSize; i++)
            //                    productCoors[i]->multiplyBy(*(lhs[0]));

        }
        return *this;

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

    /**
     * @brief compares 2 points by comparing the values of 2 coordinates, in a specified dimention.
     * @param point - a point with encrypted coordinate values.
     * @param currentDim - the index of the coordinates to be compared.
     * @returns a tuple that answers - ((p1[d]>p2[d]), (p2[d]>p1[d])). value are encrypted.
     * @return std::vector<helib::Ctxt>
     * */
    [[nodiscard]] std::vector<helib::Ctxt>
    isBiggerThan(const Point &point, short currentDim = DIM - short(1)) const {
        Ctxt mu(public_key), ni(public_key);
        if (!(isEmpty() || point.isEmpty())) {
            if (point.id == id) {
                // good solution but results in representatienves being picked twice
                //  - once for their own group and once for the group above fixme
                public_key.Encrypt(mu, NTL::to_ZZX((true)));
                public_key.Encrypt(ni, NTL::to_ZZX((true))); // make sure
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

                todo and this one too:
                    //                    /**
                    // * @brief Compute a bitwise NOT of `input`.
                    // * @param output Result of bit-flipping `input`.
                    // * @param input Binary number to be bit-flipped.
                    // * @note The size of `output` and `input` must be the same.
                    //
                    //        void bitwiseNot(CtPtrs& output, const CtPtrs& input);
                    */
    }

    bool operator==(const Point &point) const {
        return id == point.id;
    }

    /**
     * @brief return the encrypted (square of the) distance
     * @param point from which we measure our distance
     * @return encrypted (square of the) distance from @param point
     * @return std::vector<Ctxt>
     * */
    EncryptedNum
    calculateDistanceFromPoint(Point &point, const KeysServer &keysServer) {
        auto t0 = CLOCK::now();

        // c1 = p1.coor[dim]      
        // c2 = p2.coor[dim]      

        /*
                std::vector<std::vector<helib::Ctxt> > sqaredDiffs(DIM);
                for (int dim = 0; dim < DIM; ++dim) {
                    // subtract: c1 - c2
                    std::vector<helib::Ctxt> sub_vector(BIT_SIZE, helib::Ctxt(public_key));
                    helib::CtPtrs_vectorCt sub_wrapper(sub_vector);
                    helib::subtractBinary(sub_wrapper,
                                          helib::CtPtrs_vectorCt(this->cCoordinates[dim]),
                                          helib::CtPtrs_vectorCt(point.cCoordinates[dim]));
                    printNameVal(keysServer.decryptNum(sub_vector));

                    // fixme this creates a problem when the diff between the original values is negative
                    //   (meaning when c1 is smaller then c2) -
                    //   the result of `subtractBinary` is a (positive) 2's compliment.
                    //   in other words it would be
                    //   sub_result = [ NUMBERS_RANGE + (c1-c2) ] mod NUMBERS_RANGE
                    // todo
                    //  two options to fix:
                    //  1. use helib cmp w/ min/max ptrs
                    //  2. calculate (c1^2)+(c2^2)-2(c1*c2) = (c1-c2)^2  <<---------------- follows this imp
                    //  3. consult adidanny
                    // square: ( c1 - c2 )^2
                    helib::CtPtrs_vectorCt sqr_wrapper(sqaredDiffs[dim]);
                    helib::multTwoNumbers(sqr_wrapper,
                                          sub_wrapper,
                                          sub_wrapper);
                    printNameVal(keysServer.decryptNum(sqaredDiffs[dim]));
                }
                // sum: SUM[ ( c1-c2 )^2 | for all dim ]
                helib::CtPtrMat_vectorCt summands_wrapper(sqaredDiffs);
                std::vector<helib::Ctxt> result_vector;//(BIT_SIZE, helib::Ctxt(public_key));
                helib::CtPtrs_vectorCt output_wrapper(result_vector);
                helib::addManyNumbers(output_wrapper, summands_wrapper);
                printNameVal(keysServer.decryptNum(result_vector));
                */

        std::vector<std::vector<helib::Ctxt> >
                sqaredDiffs(DIM,
                            std::vector<helib::Ctxt>(2 * BIT_SIZE,
                                                     helib::Ctxt(point.public_key)));

        for (int dim = 0; dim < DIM; ++dim) {
            // subtract: c1 - c2
            helib::CtPtrs_vectorCt p1c(this->cCoordinates[dim]);
            helib::CtPtrs_vectorCt p2c(point.cCoordinates[dim]);

            //  (c1) ^ 2
            std::vector<helib::Ctxt> c1sqr;
            helib::CtPtrs_vectorCt sqr_wrapper1(c1sqr);
            helib::multTwoNumbers(sqr_wrapper1,
                                  p1c,
                                  p1c,
                                  false,
                                  2 * BIT_SIZE
            );
            printNameVal(keysServer.decryptNum(c1sqr));

            //  (c2) ^ 2
            std::vector<helib::Ctxt> c2sqr;
            helib::CtPtrs_vectorCt sqr_wrapper2(c2sqr);
            helib::multTwoNumbers(sqr_wrapper2,
                                  p2c,
                                  p2c,
                                  false,
                                  2 * BIT_SIZE
            );
            printNameVal(keysServer.decryptNum(c2sqr));

            //  (c1)^2 + (c2)^2
            std::vector<helib::Ctxt> sumSqrs;
            helib::CtPtrs_vectorCt sum_wrapper(sumSqrs);
            helib::addTwoNumbers(sum_wrapper,
                                 sqr_wrapper1,
                                 sqr_wrapper2,
                                 2 * BIT_SIZE
            );
            printNameVal(keysServer.decryptNum(sumSqrs));

            //  (c1) * (c2)
            std::vector<helib::Ctxt> multC1C2;
            helib::CtPtrs_vectorCt mult_wrapper(multC1C2);
            helib::multTwoNumbers(mult_wrapper,
                                  p1c,
                                  p2c,
                                  false,
                                  2 * BIT_SIZE
            );
            printNameVal(keysServer.decryptNum(multC1C2));

            //  2 * [(c1) * (c2)]
            std::vector<helib::Ctxt> dblMultC1C2;
            helib::CtPtrs_vectorCt dbl_wrapper(dblMultC1C2);
            helib::addTwoNumbers(dbl_wrapper,
                                 mult_wrapper,
                                 mult_wrapper,
                                 2 * BIT_SIZE
            );
            printNameVal(keysServer.decryptNum(dblMultC1C2));

            // sum: (c1^2) + (c2^2) - 2(c2*c2)
            helib::CtPtrs_vectorCt output(sqaredDiffs[dim]);
            helib::subtractBinary(output,
                                  sum_wrapper,
                                  dbl_wrapper);
            printNameVal(keysServer.decryptNum(sqaredDiffs[dim]));

        }

        // sum: SUM[ ( c1-c2 )^2 | for all dim ]
        helib::CtPtrMat_vectorCt summands_wrapper(sqaredDiffs);
        std::vector<helib::Ctxt> result_vector;
        helib::CtPtrs_vectorCt output_wrapper(result_vector);
        helib::addManyNumbers(output_wrapper,
                              summands_wrapper);
        printNameVal(keysServer.decryptNum(result_vector));

        printDuration(t0, "dist");

        return result_vector;
    }


    std::vector<std::tuple<Point, Point, EncryptedNum> >
    collectDistancesFromMeans(
            const std::vector<Point> &means,
            const KeysServer &keysServer
    ) {
        //  Collect Means
        std::vector<std::tuple<Point, Point, EncryptedNum> > distances;
        distances.reserve(means.size());
        for (const Point &mean:means) {
            //            this.
        }
        return distances;
    }


    /*  for DBG */
    std::vector<long> pCoordinatesDBG;

    bool isCopyDBG = false;
    bool isEmptyDBG = true;
    const helib::PubKey *pubKeyPtrDBG;
    long cmpCounter, addCounter, multCounter; //todo
};



/** hash functionality for std::unordered_map of Point  */
namespace std {

    template<>
    struct hash<const Point> {
        std::size_t operator()(const Point &point) const {
            using std::size_t;
            using std::hash;

            return hash<long>()(point.id);
            //            return point.id;
        }
    };

}

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


#endif //ENCRYPTEDKMEANS_POINT_H
