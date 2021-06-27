//
// Created by karina on 19/06/2021.
//

#ifndef ENCKMEAN_AUX_H
#define ENCKMEAN_AUX_H



#include <iostream>

#include <helib/FHE.h>
#include <helib/helib.h>
#include <helib/binaryArith.h>
#include <helib/intraSlot.h>

int check_DBG();
#endif //ENCKMEAN_AUX_H

std::chrono::time_point<std::chrono::system_clock> NowTime();
void printDuration( const std::chrono::time_point<std::chrono::system_clock> & t1, const std::string & funcName = "" );

