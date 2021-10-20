
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

//#include "src/KeysServer.h"

using std::cout;
using std::endl;

using helib::Ctxt;
using EncryptedBit = helib::Ctxt;
using EncryptedNum = std::vector<helib::Ctxt>;

class KeysServer;
class Point;
class Client;

//  print both the value and it's name. comfy for dgb  // best. macro. EVA! //TODO save this somewhere (list of useful tricks)
#define printNameVal(val)   cout << # val << ": " << (val) << endl

/* * for DBG * */
void printPoint(const Point &p,  const KeysServer &keysServer);

void printPoints(const std::vector<Point> &points, const KeysServer &keysServer);

void printNonEmptyPoints(const std::vector<Point> &points, const KeysServer &keysServer);
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
std::vector<Client> generateDataClients(const KeysServer &keysServer) ;

/**
 * and aux struct to put some order in the
 * */
struct Slice{
    std::vector<Point> reps;
    std::vector<Point> includedPoints;
    std::vector<helib::Ctxt> included;

    Slice(){
//        cout << "init cell" << endl;
        reps.reserve(DIM); //   should be one rep per dimension
        includedPoints.reserve(NUMBER_OF_POINTS); //   should be one rep per dimension
        included.reserve(NUMBER_OF_POINTS); //   should be one rep per dimension
    }

    Slice & addRep(const Point &point){
        reps.push_back(point);
        return *this;
    }
    Slice & addPoint(const Point &point, const helib::Ctxt &isIncluded){
        includedPoints.push_back(point);
        included.push_back(isIncluded);
        return *this;
    }

    void printSlice(const KeysServer&keysServer)const ;
};

#endif //ENCKMEAN_AUX_H
