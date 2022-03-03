
#include "TestKeysServer.h"

#include "utils/Logger.h"
#include "src/ClientDevice.h"

void TestKeysServer::testConstructor() {
    cout << " ------ testConstructor ------ " << endl;

    KeysServer keysServer;

    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestKeysServer::testEncryptCtxt() {
    cout << " ------ testConstructor ------ " << endl;

    KeysServer keysServer;

    keysServer.encryptCtxt(0);
    keysServer.encryptCtxt(1);

    long num = GiveMeRand;
    keysServer.encryptCtxtCompact(num);
    keysServer.encryptCtxtCompact(num);

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

    for (int test_i = 0; test_i < 10; ++test_i) {
        long randomLong = GiveMeRand;
        helib::Ptxt<helib::BGV> ptxt(publicKey.getContext(), NTL::ZZX(randomLong));
        Ctxt ctxt(keysServer.encryptCtxtCompact(randomLong));
        long decryptedNum = keysServer.decryptCtxtCompact(ctxt);
        printNameVal(randomLong);
        printNameVal(decryptedNum);
        assert(randomLong==decryptedNum);
    }

    cout << " ------ testDecryptCtxt finished ------ " << endl << endl;
}

void TestKeysServer::testEncryptNum() {
    cout << " ------ testEncryptNum ------ " << endl;

    KeysServer keysServer;

    const long l = randomLongInRange(mt);

    keysServer.encryptNum(l);

    cout << " ------ testEncryptNum finished ------ " << endl << endl;
}

void TestKeysServer::testDecryptNum() {
    cout << " ------ testDecryptNum ------ " << endl;

    KeysServer keysServer;

    const long l = randomLongInRange(mt);
    EncryptedNum cl = keysServer.encryptNum(l);

    const long pl = keysServer.decryptNum(cl);

    assert(l == pl);

    cout << " ------ testDecryptNum finished ------ " << endl << endl;
}

void TestKeysServer::testScratchPoint() {
    cout << " ------ testEncryptScratchPoint ------ " << endl;

    KeysServer keysServer;

    Point scratchPoint = keysServer.scratchPoint();

    printNameVal(scratchPoint.isEmpty());

    cout << " ------ testEncryptScratchPoint finished ------ " << endl << endl;
}

void TestKeysServer::testTinyRandomPoint() {
    cout << " ------ testTinyRandomPoint ------ " << endl;

    KeysServer keysServer;

    for (int i = 0; i < 10; ++i) {
        Point tiny(keysServer.tinyRandomPoint());
        cout << printPoint(tiny, keysServer);
    }

    cout << " ------ testTinyRandomPoint finished ------ " << endl << endl;
}
