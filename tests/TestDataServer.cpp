
#include "TestDataServer.h"

#include "utils/Logger.h"
#include "src/DataServer.h"

static Logger loggerTestDataServer(log_debug, "loggerTestDataServer");

void TestDataServer::testConstructor() {
    loggerTestDataServer.log("testConstructor");
    KeysServer keysServer;
    DataServer dataServer(keysServer);
}

void TestDataServer::testScratchPoint() {
    loggerTestDataServer.log("testEncryptScratchPoint");
    KeysServer keysServer;
    DataServer dataServer(keysServer);
    Point scratchPoint = dataServer.scratchPoint();
//    cCoordinates.emplace_back(bitSize, helib::Ctxt(public_key));
}

//void TestDataServer::testDecryptCoordinates() {
//
//}
