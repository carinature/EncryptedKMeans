
#include "TestKeysServer.h"

#include "utils/Logger.h"
#include "src/Client.h"

void TestKeysServer::testConstructor() {

    //    loggerTestClient.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl ;
    KeysServer server;
//    Client client(server);
    //    loggerTestClient.print_log();
    cout << " ------ testConstructor finished ------ " << endl << endl;
}
