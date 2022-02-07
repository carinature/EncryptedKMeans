
#ifndef ENCRYPTEDKMEANS_POINT_H
#define ENCRYPTEDKMEANS_POINT_H

#include <helib/binaryCompare.h>
#include <helib/binaryArith.h>
#include <NTL/ZZX.h>

#include "utils/aux.h"
#include "KeysServer.h"


class Point {
    friend class ClientDevice;

    friend class KeysServer;

public:

    //! @var long id
    //! used in createCmpDict for comparison
    const long id;
    Point *originalPointAddress;
    EncryptedNum cid;
    const helib::PubKey &public_key;// = encryptionKey;
    std::vector<EncryptedNum> cCoordinates;

    /*
     Each bit of the binary number is encoded into a single ciphertext. Thus
     for a 16 bit binary number, we will represent this as an array of 16
     unique ciphertexts.
     i.e. b0 = [0] [0] [0] ... [0] [0] [0]        ciphertext for bit 0
          b1 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 1
          b2 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 2
     These 3 ciphertexts represent the 3-bit binary number 110b = 6
     Note: several numbers can be encoded across the slots of each ciphertext
     which would result in several parallel slot-wise operations.
     For simplicity we place the same data into each slot of each ciphertext,
     printing out only the back of each vector.
     NB: fifteenOrLess4Four max is 15 bits. Later in the code we pop the MSB.
            long BIT_SIZE = 16;
            long OUT_SIZE = 2 * BIT_SIZE;
     */
    /**
     * @brief Construct a point from plaintext coordinates
     * @param public_key The key used to encrypt the data (can't be used for decryption)
     * @param coordinates Plaintext data (currently coordinates is an array - will be changed to vector?)
     * @note an "empty" point is created when no coordinates are passed
     * */
    explicit Point(const helib::PubKey &public_key, const long coordinates[] = nullptr);

    /**
     * @brief Construct a point from ciphertext coordinates
     * @param cCoordinates Encrypted data
     * @note an "empty" point is created when no coordinates are passed
     * */
    explicit Point(const std::vector<EncryptedNum> &cCoordinates);

    bool isEmpty() const {
        return cCoordinates[0][0].isEmpty();
    }

    /**
     * @brief Copy Constructor
     * @note value of point's id is copied
     * @note value of point's isCopyDBG is set to `true`
     * */
    Point(const Point &point);

    Point &operator=(const Point &point);

    /**
     * @param i The coordinate for which we request the data
     * @return The const reference to the data in coordinate i
     * */
    const EncryptedNum &operator[](short int i) const;

    /**
     * @param i The coordinate for which we request the data
     * @return The data in coordinate i
     * */
    EncryptedNum operator[](short int i);

    /**
     * @brief for each i returns this[i]+point[i]
     * */
    Point operator+(Point &point);

    Point operator+(Point point);

    //// Calculates the sum of many numbers using the 3-for-2 method
    static Point addManyPoints(const std::vector<Point> &points, const KeysServer &keysServer);

    /**
     * @brief Multiplies an encrypted point by an encrypted bit
     * */
    Point operator*(const Ctxt &bit) const;

    Point &operator*=(const Ctxt &bit);

    // todo maybe unused. delete?
    Point operator*(Point &point) {
        if (point.isEmpty()) return *this;
        if (this->isEmpty()) return point;
        const long arr[] = {0, 0};
        Point product(this->public_key, arr);
        for (short dim = 0; dim < DIM; ++dim) {
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
     * @return EncryptedNum
     * */
    std::vector<CBit>
    isBiggerThan(const Point &point, short int currentDim = DIM - 1) const;

    bool operator==(const Point &point) const {
        //        return cid == point.cid;
        return id == point.id;
    }

    /**
     * @brief return the encrypted (square of the) distance
     * @param point from which we measure our distance
     * @return encrypted (square of the) distance from point
     * @return EncryptedNum
     * */
    EncryptedNum
    distanceFrom(
            const Point &point,
            const KeysServer &keysServer
    ) const;

    /**
     * @brief find closest point, from a list, and minimal distance from it
     * @param points list of points from which we measure our distance
     * @return minimal distance and corresponding closest point
     * @return std::pair<Point, EncryptedNum>
     * */
    std::pair<Point, EncryptedNum>
    findMinDistFromMeans(
            const std::vector<Point> &points,
            const KeysServer &keysServer
    ) const;


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
        //        return this == p;
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

//int counter = 0;
//
//class PointCompact {
//    helib::PubKey public_key;
//    std::vector<helib::Ctxt> coordinates;
//    long id;
//    helib::Ctxt cid;
//    /*  for DBG */
//    std::vector<long> pCoordinates;
//    /*  end for DBG */
//
///**
// * Ctor
// * */
//    PointCompact(const helib::PubKey &public_key, const std::vector<long> coordinates) :
//            public_key(public_key),
//            id(counter++),
//            cid(public_key),
//            coordinates(DIM, helib::Ctxt(public_key)),
//            pCoordinates(coordinates) {
//        for
//    }
//};

#endif //ENCRYPTEDKMEANS_POINT_H
