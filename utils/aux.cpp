

#include "aux.h"
#include "../properties.h"
#include "Logger.h"

#include <sstream>      // std::stringstream

int check_DBG() {
    std::cout << "DBG is "; //ON" <<endl;

#ifdef DBG
    std::cout << "defined ";
#if DBG
    std::cout << "TRUE" << std::endl;
    return 1;
#else
    std::cout<< "False" <<std::endl;
    return 0;
#endif

#else
    std::cout << "NOT defined" << std::endl;
#endif

    return -1;
}


std::chrono::time_point<std::chrono::system_clock> NowTime() {
    return std::chrono::high_resolution_clock::now();
}

std::string printDuration(const std::chrono::time_point<std::chrono::system_clock> &t1, const std::string &funcName) {
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::string str = "\'" + funcName + "\' Finished in " + std::to_string(duration) + " seconds.";
    std::cout << str << std::endl;
    fcout << str << std::endl; //todo remove?
    return str;
}



