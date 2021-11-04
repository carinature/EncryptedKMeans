
#include "TestKeysServer.h"

#include "utils/Logger.h"
#include "src/Client.h"

void TestKeysServer::testConstructor() {
    //    loggerTestClient.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl;
    KeysServer keysServer;
    //    loggerTestClient.print_log();
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestKeysServer::testEncryptCtxt() {
    cout << " ------ testConstructor ------ " << endl;
    KeysServer keysServer;
    keysServer.encryptCtxt(0);
    keysServer.encryptCtxt(1);
    //    loggerTestClient.print_log();
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestKeysServer::testDecryptCtxt() {
    cout << " ------ testDecryptCtxt ------ " << endl;
    KeysServer keysServer;
    helib::PubKey &publicKey = keysServer.getPublicKey();

    long bit0 = 0;
    long bit1 = 1;

    helib::Ctxt cbit0(publicKey);
    helib::Ctxt cbit1(publicKey);

    publicKey.Encrypt(cbit0, NTL::ZZX(bit0));
    publicKey.Encrypt(cbit1, NTL::ZZX(bit1));

    long dbit0 = keysServer.decryptCtxt(cbit0);
    long dbit1 = keysServer.decryptCtxt(cbit1);

    assert(bit0 == dbit0);
    assert(bit1 == dbit1);
    //    loggerTestClient.print_log();
    cout << " ------ testDecryptCtxt finished ------ " << endl << endl;
}

#include <random>

void TestKeysServer::testEncryptNum() {
    cout << " ------ testEncryptNum ------ " << endl;
    KeysServer keysServer;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<long> dist(1, NUMBERS_RANGE);
    const long l = dist(mt);
    keysServer.encryptNum(l);
    //    Client client(keysServer);
    //    loggerTestClient.print_log();
    cout << " ------ testEncryptNum finished ------ " << endl << endl;
}

//todo
void TestKeysServer::testDecryptNum() {
    cout << " ------ testDecryptNum ------ " << endl;
    KeysServer keysServer;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<long> dist(1, NUMBERS_RANGE);
    const long l = dist(mt);
    EncryptedNum cl = keysServer.encryptNum(l);
    const long pl = keysServer.decryptNum(cl);
    assert(l==pl);
    //    Client client(keysServer);
    //    loggerTestClient.print_log();
    cout << " ------ testDecryptNum finished ------ " << endl << endl;
}

void TestKeysServer::testScratchPoint() {
    cout << " ------ testEncryptScratchPoint ------ " << endl;
//    loggerTestDataServer.log("testEncryptScratchPoint");
    KeysServer keysServer;
//    DataServer dataServer(keysServer);
    Point scratchPoint = keysServer.scratchPoint();
    //    cCoordinates.emplace_back(BIT_SIZE, helib::Ctxt(public_key));
    printNameVal(scratchPoint.isEmpty());
    cout << " ------ testEncryptScratchPoint finished ------ " << endl << endl;
}

void TestKeysServer::testTinyRandomPoint() {
    cout << " ------ testTinyRandomPoint ------ " << endl;
    KeysServer keysServer;
    for (int i = 0; i < 10; ++i) {
        Point tiny(keysServer.tinyRandomPoint());
        printPoint(tiny, keysServer);
    }
    //    loggerTestClient.print_log();
    cout << " ------ testTinyRandomPoint finished ------ " << endl << endl;
}
