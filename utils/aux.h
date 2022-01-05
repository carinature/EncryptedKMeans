
#ifndef ENCKMEAN_AUX_H
#define ENCKMEAN_AUX_H

/** @file aux
 * */

#include <iostream>

#include <helib/FHE.h>
#include <helib/helib.h>
#include <helib/intraSlot.h>

#include <vector>
#include <helib/zzX.h>
#include <helib/Context.h>
#include <helib/keys.h>
#include <helib/debugging.h>

#include "properties.h"
#include "Logger.h"

using std::cout;
using std::endl;

using helib::Ctxt;
//using EncryptedBit = helib::Ctxt;
using EncryptedNum = std::vector<helib::Ctxt>;

#include <random>
//static std::mt19937 mt;
static std::random_device rd;
static std::mt19937 mt(rd()); //    for extra randomness
static std::uniform_int_distribution<long> randomLongInRange(0, NUMBERS_RANGE);
#define giveMeRandomLong() randomLongInRange(mt)

class KeysServer;

class Point;

class Client;

//  print both the value and it's name. comfy for dgb  // best. macro. EVA!
#define printNameVal(val)   cout << # val << ": " << (val) << endl
// TODO save this somewhere (list of useful tricks)

/* * for DBG * */
void printPoint(const Point &p, const KeysServer &keysServer);

std::vector<long> decryptPoint(const Point &p, const KeysServer &keysServer);

void printPoints(const std::vector<Point> &points, const KeysServer &keysServer);

void printNonEmptyPoints(const std::vector<Point> &points, const KeysServer &keysServer);

/** Writing points to a specified.
 *      The points are first decrypted and returned into a double form (from long) **/
void
decAndWriteToFile(const std::vector<Point> & points, const std::string & filename, const KeysServer &keysServer ); // <----------

/* * fin * */

std::chrono::time_point<std::chrono::system_clock> NowTime();

std::string
printDuration(const std::chrono::time_point<std::chrono::system_clock> &t1,
              const std::string &funcName = "");

/**
 * @brief generate a database of Clients and their data Points (randomized)
 * @param keysServer - a reference to the CA
 * @returns a list of data Clients
 * @return std::vector<Client>
 * */
std::vector<Client> generateDataClients(const KeysServer &keysServer);

//struct PointTuple {
//    const Point &point;
//    const CBit isIn;
//
//    PointTuple(const Point &point, const CBit &isIn) : point(point), isIn(isIn) {};
//};
//using PointTuple = std::tuple<Point, CBit>;
using PointTuple = std::pair<Point, CBit>;

/**
 * and aux struct to put some order in the
 * */
struct Slice {
    std::vector<Point> reps;
    std::vector<Point> points;
    std::vector<helib::Ctxt> counter;
    std::vector<PointTuple> pointTuples;

    Slice() {
        //        cout << "init cell" << endl;
        reps.reserve(DIM); //   should be one rep per dimension
        points.reserve(NUMBER_OF_POINTS); //   should be one rep per dimension
        counter.reserve(NUMBER_OF_POINTS); //   should be one rep per dimension
        pointTuples.reserve(NUMBER_OF_POINTS); //   should be one rep per dimension
    }

    Slice &addRep(const Point &point) {
        reps.push_back(point);
        return *this;
    }

    Slice &addReps(const std::vector<Point> &repPoints) {
        for (const Point &rep:repPoints) reps.push_back(rep);
        return *this;
    }

    Slice &addPoint(const Point &point, const helib::Ctxt &isIncluded);

    void printSlice(const KeysServer &keysServer) const;

    void clear() {
        reps.clear();
        points.clear();
        counter.clear();
        pointTuples.clear();
    }

    ~Slice() {
        this->clear();
    }

};

/**
 * Implementing the alg in the article:
 * "Setup-Free Secure Search on Encrypted Data: Faster and Post-Processing Free"
 */
std::vector<Ctxt> prefix(const std::vector<Ctxt> &v, long k);
std::vector<Ctxt> suffix(const std::vector<Ctxt> &v, long len);

Ctxt isEqual(const std::vector<Ctxt> &a,
             const std::vector<Ctxt> &b,
             long w=BIT_SIZE
             );

Ctxt isGreater(const std::vector<Ctxt> &a, const std::vector<Ctxt> &b, long w);

#endif //ENCKMEAN_AUX_H
