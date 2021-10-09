
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
    cout << " ------ testCompare ------ " << endl << endl;
    KeysServer server; // = KeysServer();
    Client client1(server);
    Client client2(server);
    long arr[DIM], arr1[DIM], arr2[DIM];
    std::vector<Client> clients;//(2, Client(server));
    //    for (Client &c:clients) {
    clients.reserve(2);
    for (int i = 0; i < 2; ++i) {
        Client c(server);
        c.encryptPoint(arr);
        clients.push_back(c);
    }
    for (int dim = 0; dim < DIM; ++dim) {
        arr[dim] = rand() % NUMBERS_RANGE;
        arr1[dim] = rand() % NUMBERS_RANGE;
        arr2[dim] = rand() % NUMBERS_RANGE;
    }
    client1.encryptPoint(arr);
    client1.encryptPoint(arr1);
    client2.encryptPoint(arr2);
    std::vector<long*> arrs;//(2, long[DIM]);
//    for (Client &c:clients) {
    arrs.reserve(2);
    for (int i = 0; i < 2; ++i) {
//        arrs.push_back(long[DIM]);
        long arr0[DIM];
        for (int dim = 0; dim < DIM; ++dim) {
            arr0[dim] = rand() % NUMBERS_RANGE;
        }
        arrs.push_back(arr0);
    }
    std::vector<Point> points;
    points.reserve(2);
    for (Client &c:clients) {
        for (Point &p:c.points) {
            points.push_back(Point(p));
        }
    }
    cout << "one ";
    Point point1(client1.points[0]);
    cout << "two ";
    Point point2(client2.points[0]);
    cout << "three ";
    Ctxt isBigger = point1.isBiggerThan(point2);
    cout << "four ";
    bool well = server.decryptCtxt(isBigger);
    printNameVal(well);

    //            client1.compare(client2);
    //    helib::decryptBinaryNums();
    //    loggerTestClient.print_log();
    cout << " ------ testCompare finished ------ " << endl << endl;
}

void TestClient::testAddition() {
    //    loggerTestClient.log("testAddiiton");
}

void TestClient::testMultiplication() {
    //    loggerTestClient.log("testMultiplication");
}