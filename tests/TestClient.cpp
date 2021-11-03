
#include "TestClient.h"

#include "utils/Logger.h"
#include "src/Client.h"

//static Logger loggerTestClient(log_debug, "loggerTestClient");

void TestClient::testConstructor() {
    //    loggerTestClient.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl;
    KeysServer keysServer;
    Client client(keysServer);
    //    loggerTestClient.print_log();
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestClient::testEncryptCoordinates() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testEncryptCoordinates ------ " << endl;
    KeysServer keysServer;
    long arr[DIM];
    for (short dim = 0; dim < DIM; ++dim) arr[dim] = rand() % NUMBERS_RANGE;
    Client client(keysServer);
    client.encryptPoint(arr);
    cout << " ------ testEncryptCoordinates finished ------ " << endl << endl;
}


void TestClient::testDecryptCoordinates() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testDecryptCoordinates ------ " << endl;
    KeysServer keysServer;
    long arr[DIM];
    for (short dim = 0; dim < DIM; ++dim) arr[dim] = rand() % NUMBERS_RANGE;
    Client client(keysServer);
    client.encryptPoint(arr);
    std::vector<long> decryptCoordinates = client.decryptCoordinate(0);
    for (int i = 0; i < DIM; ++i) {
        //        printNameVal(arr[i]);
        //        printNameVal(decryptCoordinate[i]);
        assert(decryptCoordinates[i] == arr[i]);
    }
    //    loggerTestClient.print_log();
    cout << " ------ testDecryptCoordinates finished ------ " << endl << endl;
}


void TestClient::testEncryptScratchPoint() {
    cout << " ------ testEncryptScratchPoint ------ " << endl << endl;
    KeysServer keysServer;// = KeysServer();
    Client client(keysServer);
    client.encryptPoint();
    std::vector<long> decryptCoordinates = client.decryptCoordinate(0);
    for (auto ds:decryptCoordinates) assert(0 == ds);
    cout << " ------ testEncryptScratchPoint finished ------ " << endl << endl;
}


void TestClient::testCompare() {
    //    loggerTestClient.log("testComparePoints");
    cout << " ------ testCompare ------ " << endl << endl;
    KeysServer keysServer;
    long arr1[DIM], arr2[DIM];
    for (short dim = 0; dim < DIM; ++dim) {
        arr1[dim] = rand() % NUMBERS_RANGE;
        arr2[dim] = rand() % NUMBERS_RANGE;
    }

    Client client1(keysServer);
    client1.encryptPoint(arr1);
    Point point1(client1.getPoints().back());

    Client client2(keysServer);
    client2.encryptPoint(arr2);
    Point point2(client2.getPoints().back());

    for (int i = 0; i < 100; ++i)
        for (short dim = 0; dim < DIM; ++dim) {
            helib::Ctxt res = point1.isBiggerThan(point2, dim)[0];
            assert((arr1[dim] > arr2[dim]) == keysServer.decryptCtxt(res));

            helib::Ctxt res2 = point2.isBiggerThan(point1, dim)[0];
            assert((arr2[dim] > arr1[dim]) == keysServer.decryptCtxt(res2));
        }

    cout << " ------ testCompare finished ------ " << endl << endl;
}

void TestClient::testAddition() {
    //    loggerTestClient.log("testAddition");
}

void TestClient::testMultiplication() {
    //    loggerTestClient.log("testMultiplication");
}