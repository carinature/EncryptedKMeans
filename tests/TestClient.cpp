
#include "TestClient.h"

#include "utils/Logger.h"
#include "src/Client.h"

//static Logger loggerTestClient(log_debug, "loggerTestClient");

void TestClient::testConstructor() {
    //    loggerTestClient.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl;
    KeysServer server;
    Client client(server);
    //    loggerTestClient.print_log();
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestClient::testEncryptCoordinates() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testEncryptCoordinates ------ " << endl;
    KeysServer server;
    long arr[DIM];
    for (short dim = 0; dim < DIM; ++dim) arr[dim] = rand() % NUMBERS_RANGE;
    Client client(server);
    client.encryptPoint(arr);
    cout << " ------ testEncryptCoordinates finished ------ " << endl << endl;
}


void TestClient::testDecryptCoordinates() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testDecryptCoordinates ------ " << endl;
    KeysServer server;
    long arr[DIM];
    for (short dim = 0; dim < DIM; ++dim) arr[dim] = rand() % NUMBERS_RANGE;
    Client client(server);
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
    KeysServer server = KeysServer();
    Client client(server);
    client.encryptPoint();
    std::vector<long> decryptCoordinates = client.decryptCoordinate(0);
    for (auto ds:decryptCoordinates) assert(0 == ds);
    cout << " ------ testEncryptScratchPoint finished ------ " << endl << endl;
}


void TestClient::testCompare() {
    //    loggerTestClient.log("testComparePoints");
    cout << " ------ testCompare ------ " << endl << endl;
    KeysServer server;
    long arr1[DIM], arr2[DIM];
    for (short dim = 0; dim < DIM; ++dim) {
        arr1[dim] = rand() % NUMBERS_RANGE;
        arr2[dim] = rand() % NUMBERS_RANGE;
    }

    Client client1(server);
    client1.encryptPoint(arr1);
    Point point1(client1.points.back());

    Client client2(server);
    client2.encryptPoint(arr2);
    Point point2(client2.points.back());

    for (int i = 0; i < 100; ++i)
        for (short dim = 0; dim < DIM; ++dim) {
            helib::Ctxt res = point1.isBiggerThan(point2, dim)[0];
            assert((arr1[dim] > arr2[dim]) == server.decryptCtxt(res));

            helib::Ctxt res2 = point2.isBiggerThan(point1, dim)[0];
            assert((arr2[dim] > arr1[dim]) == server.decryptCtxt(res2));
        }

    //    loggerTestClient.print_log();
    cout << " ------ testCompare finished ------ " << endl << endl;
}

void TestClient::testAddition() {
    //    loggerTestClient.log("testAddition");
}

void TestClient::testMultiplication() {
    //    loggerTestClient.log("testMultiplication");
}