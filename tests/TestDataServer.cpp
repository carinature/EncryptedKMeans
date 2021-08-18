
#include "TestDataServer.h"

#include "utils/Logger.h"
#include "src/DataServer.h"

static Logger loggerTestDataServer(log_debug, "loggerTestDataServer");

void TestDataServer::testConstructor() {
    loggerTestDataServer.log("testConstructor");
    KeysServer keysServer;
    DataServer DataServer(keysServer);
}

//void TestDataServer::testDecryptCoordinates() {
//
//}
