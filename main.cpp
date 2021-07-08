//
// run a full protocol
//
#include <iostream>
#include <chrono>
#include <fstream>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

#include "properties.h"
#include "utils/Logger.h"
#include "utils/aux.h"

#define printNameVal(val)   std::cout << # val << ": " << (val) << std::endl

int main() {
    Logger logger;
    logger.log(log_info, "Starting Protocol");
    auto t1 = std::chrono::high_resolution_clock::now();
    check_DBG();

    //// JSON
    std::cout << "jsonConfig " << jsonConfig << std::endl;
    std::cout << "range_lim_comment " << range_lim_comment << std::endl;
    printNameVal(range_lim_comment);

    printDuration(t1, "Main");
    logger.print_log(log_error, false);
    return 0;
}



////// GLOBAL VARIABLES good example
//int g_x{}; // global variable g_x
//void doSomething() {
//    // global variables can be seen and used everywhere in the file
//    g_x = 3;
//    std::cout << g_x << '\n';
//    int g_x = 2;
//    std::cout << g_x << '\n';
//}
//int main() {
//    std::cout << g_x << '\n';
//    doSomething();
//    std::cout << g_x << '\n';
//    g_x = 5;
//    std::cout << g_x << '\n';

//    int mynum;
//    std::cout << "enter number->" << std::endl;
//    std::cin >> mynum;
//    printNameVal(mynum);

//}
