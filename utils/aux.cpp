

#include "aux.h"
#include "properties.h"
#include "Logger.h"

#include <sstream>      // std::stringstream

using std::cout;
using std::endl;

std::chrono::time_point<std::chrono::system_clock> NowTime() {
    return std::chrono::high_resolution_clock::now();
}

std::string printDuration(const std::chrono::time_point<std::chrono::system_clock> &t1, const std::string &funcName) {
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::string str = "\'" + funcName + "\' Finished in " + std::to_string(duration) + " seconds.";
    cout << str << endl;
    fcout << str << endl; //todo remove?
    return str;
}



