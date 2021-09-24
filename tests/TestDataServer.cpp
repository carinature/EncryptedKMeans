
#include "TestDataServer.h"

#include "utils/Logger.h"
#include "src/DataServer.h"

static Logger loggerTestDataServer(log_debug, "loggerTestDataServer");

void TestDataServer::testConstructor() {
    loggerTestDataServer.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestDataServer::testScratchPoint() {
    loggerTestDataServer.log("testEncryptScratchPoint");
    cout << " ------ testEncryptScratchPoint ------ " << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);
    Point scratchPoint = dataServer.scratchPoint();
    //    cCoordinates.emplace_back(BIT_SIZE, helib::Ctxt(public_key));
    cout << " ------ testEncryptScratchPoint finished ------ " << endl << endl;
}

void TestDataServer::testCompareClients() {

}

//void TestDataServer::testDecryptCoordinates() {
//
//}
