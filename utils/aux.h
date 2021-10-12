
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

//  print both the value and it's name. comfy for dgb  // best. macro. EVA! //TODO save this somewhere (list of useful tricks)
#define printNameVal(val)   cout << # val << ": " << (val) << endl

std::chrono::time_point<std::chrono::system_clock> NowTime();

std::string
printDuration(const std::chrono::time_point<std::chrono::system_clock> &t1,
              const std::string &funcName = "");

/*
 * for DBG
 * */
void printPoint(const Point &p, KeysServer &keysServer);

void printPoints(const std::vector<Point> &points, KeysServer &keysServer);

#endif //ENCKMEAN_AUX_H
