

#ifndef ENCKMEAN_AUX_H
#define ENCKMEAN_AUX_H


#include <iostream>

#include <helib/FHE.h>
#include <helib/helib.h>
#include <helib/binaryArith.h>
#include <helib/intraSlot.h>

//#include "../properties.h"

//  print both the value and it's name. comfy for dgb
#define printNameVal(val)   std::cout << # val << ": " << (val) << std::endl

int check_DBG(); // fixme remove

std::chrono::time_point<std::chrono::system_clock> NowTime();
void printDuration( const std::chrono::time_point<std::chrono::system_clock> & t1, const std::string & funcName = "" );


#endif //ENCKMEAN_AUX_H
