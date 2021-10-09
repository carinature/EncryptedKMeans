
#include "TestDataServer.h"

//#include "utils/Logger.h"
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
    printNameVal(scratchPoint.isEmpty());
    cout << " ------ testEncryptScratchPoint finished ------ " << endl << endl;
}

void TestDataServer::testRetrievePoints() {
    //    loggerTestClient.log("testComparePoints");
    cout << " ------ testRetrievePoints ------ " << endl << endl;
    KeysServer server;
    std::vector<Client> clients = DataServer::generateDataClients(server);
    std::vector<Point> points = DataServer::retrievePoints(clients);

    cout << "===========" << endl;
    //    printNameVal(uniquePointsNum);
    printNameVal(points.size());
    for (Point &p:points) {
        cout << "( ";
        for (int dim = 0; dim < DIM - 1; ++dim)
            cout << server.decryptNum(p[dim]) << ",";
        cout << server.decryptNum(p[DIM - 1]) << " ) " << endl;
    }
    // todo worth reading again, may improve efficiency:
    //  https://stackoverflow.com/questions/25108854/initializing-the-size-of-a-c-vector
    cout << " ------ testRetrievePoints finished ------ " << endl << endl;
}



void TestDataServer::testComparePoints() {
    //    loggerTestClient.log("testComparePoints");
    cout << " ------ testComparePoints ------ " << endl << endl;

    KeysServer server;
    long arr[DIM], arr2[DIM];
    for (int dim = 0; dim < DIM; ++dim) {
        arr[dim] = rand() % NUMBERS_RANGE;
        arr2[dim] = rand() % NUMBERS_RANGE;
    }
    Point point(server.getPublicKey(), arr);
    Point point2(server.getPublicKey(), arr2);
    for (short dim = 0; dim < DIM; ++dim) {
        for (int i = 0; i < 100; ++i) {
            helib::Ctxt res = point.isBiggerThan(point2, dim);
            assert((arr[dim] > arr2[dim]) == server.decryptCtxt(res));

            helib::Ctxt res2 = point2.isBiggerThan(point, dim);
            assert((arr2[dim] > arr[dim]) == server.decryptCtxt(res2));
        }
    }
    cout << " ------ testComparePoints finished ------ " << endl << endl;
}

//void TestDataServer::testDecryptCoordinates() {
//    cout << " ------ testDecryptCoordinates ------ " << endl << endl;
//    cout << " ------ testDecryptCoordinates finished ------ " << endl << endl;
//
//}
