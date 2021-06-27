//
// run a full protocol
//
#include "properties.h"
#include "utils/aux.h"
#include "utils/Logger.h"

int main(){
    Logger logger;
    logger.log(log_info, "Starting Protocol");
    auto t1 = std::chrono::high_resolution_clock::now();
    check_DBG();

    printDuration( t1, "Main" );
    return 0;
}