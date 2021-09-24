
#include "TestClient.h"

#include "utils/Logger.h"
#include "src/Client.h"

//static Logger loggerTestClient(log_debug, "loggerTestClient");

void TestClient::testConstructor() {
    //    loggerTestClient.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl;
    KeysServer server;// = KeysServer();
    Client client(server);
    //    loggerTestClient.print_log();
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestClient::testEncryptCoordinates() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testEncryptCoordinates ------ " << endl;
    KeysServer server;// = KeysServer();
    long arr[DIM];
    for (int dim = 0; dim < DIM; ++dim) arr[dim] = rand() % bitSizeRange;
    Client client(server);
    client.encryptPoint(arr);
    cout << " ------ testEncryptCoordinates finished ------ " << endl << endl;
}


void TestClient::testDecryptCoordinates() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testDecryptCoordinates ------ " << endl;
    KeysServer server;// = KeysServer();
    long arr[DIM];
    for (int dim = 0; dim < DIM; ++dim) arr[dim] = rand() % bitSizeRange;
    Client client(server);
    client.encryptPoint(arr);
    std::vector<long> decryptCoordinates = client.decryptCoordinates(0);
    for (int i = 0; i < DIM; ++i) {
        printNameVal(arr[i]);
        printNameVal(decryptCoordinates[i]);
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
    std::vector<long> decryptCoordinates = client.decryptCoordinates(0);
    for (auto ds:decryptCoordinates) assert(0 == ds);
    cout << " ------ testEncryptScratchPoint finished ------ " << endl << endl;
}


void TestClient::testCompare() {
    //    loggerTestClient.log("testCompareClients");
    KeysServer server = KeysServer();
    Client client1(server);
    Client client2(server);
    const long arr1[] = {1L, 2L};
    const long arr2[] = {2L, 2L};
    client1.encryptPoint(arr1);
    client2.encryptPoint(arr2);
    //        client1.compare(client2);
    //    helib::decryptBinaryNums();
    //    loggerTestClient.print_log();
}

void TestClient::testAddition() {
    //    loggerTestClient.log("testAddiiton");
}

void TestClient::testMultiplication() {
    //    loggerTestClient.log("testMultiplication");
}