
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
    for (int dim = 0; dim < DIM; ++dim) arr[dim] = rand() % NUMBERS_RANGE;
    Client client(server);
    client.encryptPoint(arr);
    cout << " ------ testEncryptCoordinates finished ------ " << endl << endl;
}


void TestClient::testDecryptCoordinates() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testDecryptCoordinates ------ " << endl;
    KeysServer server;
    long arr[DIM];
    for (int dim = 0; dim < DIM; ++dim) arr[dim] = rand() % NUMBERS_RANGE;
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
    cout << " ------ testEncryptScratchPoint ------ " << endl << endl;
    KeysServer server = KeysServer();
    Client client1(server);
    Client client2(server);
    long arr[DIM], arr1[DIM], arr2[DIM];
    for (int dim = 0; dim < DIM; ++dim) {
        arr[dim] = rand() % NUMBERS_RANGE;
        arr1[dim] = rand() % NUMBERS_RANGE;
        arr2[dim] = rand() % NUMBERS_RANGE;
    }
    client1.encryptPoint(arr);
    client1.encryptPoint(arr1);
    client2.encryptPoint(arr2);
    //            client1.compare(client2);
    //    helib::decryptBinaryNums();
    //    loggerTestClient.print_log();
    cout << " ------ testEncryptScratchPoint finished ------ " << endl << endl;
}

void TestClient::testAddition() {
    //    loggerTestClient.log("testAddiiton");
}

void TestClient::testMultiplication() {
    //    loggerTestClient.log("testMultiplication");
}