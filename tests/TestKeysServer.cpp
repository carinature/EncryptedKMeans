
#include "TestKeysServer.h"

#include "utils/Logger.h"
#include "src/Client.h"

void TestKeysServer::testConstructor() {

    //    loggerTestClient.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl;
    KeysServer server;
    //    loggerTestClient.print_log();
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestKeysServer::testEncryptCtxt() {
    cout << " ------ testConstructor ------ " << endl;
    KeysServer server;
    server.encryptCtxt(0);
    server.encryptCtxt(1);
    //    loggerTestClient.print_log();
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestKeysServer::testDecryptCtxt() {

    cout << " ------ testDecryptCtxt ------ " << endl;
    KeysServer server;
    long bit0 = 0;
    long bit1 = 1;
    helib::Ctxt cbit0 = server.encryptCtxt(bit0);
    helib::Ctxt cbit1 = server.encryptCtxt(bit1);
    long dbit0 = server.decryptCtxt(cbit0);
    long dbit1 = server.decryptCtxt(cbit1);
    assert(bit0 == dbit0);
    assert(bit1 == dbit1);
    //    loggerTestClient.print_log();
    cout << " ------ testDecryptCtxt finished ------ " << endl << endl;
}

//todo
void TestKeysServer::testEncryptNum() {
    cout << " ------ testEncryptNum ------ " << endl;
    KeysServer server;
    //    Client client(server);
    //    loggerTestClient.print_log();
    cout << " ------ testEncryptNum finished ------ " << endl << endl;
}

//todo
void TestKeysServer::testDecryptNum() {

    cout << " ------ testDecryptNum ------ " << endl;
    KeysServer server;
    //    Client client(server);
    //    loggerTestClient.print_log();
    cout << " ------ testDecryptNum finished ------ " << endl << endl;
}
