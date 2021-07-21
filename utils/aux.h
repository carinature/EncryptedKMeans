

#ifndef ENCKMEAN_AUX_H
#define ENCKMEAN_AUX_H


#include <iostream>

#include <helib/FHE.h>
#include <helib/helib.h>
#include <helib/binaryArith.h>
#include <helib/intraSlot.h>

//#include "../properties.h"

#include "../src/KeysServer.h"

//  print both the value and it's name. comfy for dgb  // best. macro. EVA! //TODO save this somewhere (list of useful tricks)
#define printNameVal(val)   cout << # val << ": " << (val) << endl

std::chrono::time_point<std::chrono::system_clock> NowTime();

std::string
printDuration(const std::chrono::time_point<std::chrono::system_clock> &t1, const std::string &funcName = "");


#endif //ENCKMEAN_AUX_H
