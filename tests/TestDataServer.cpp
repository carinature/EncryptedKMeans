
#include "TestDataServer.h"

#include "utils/Logger.h"
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
    cout << " ------ testEncryptScratchPoint finished ------ " << endl << endl;
}

void TestDataServer::testRetrievePoints() {
    //    loggerTestClient.log("testComparePoints");
    cout << " ------ testRetrievePoints ------ " << endl << endl;
    KeysServer server;// = KeysServer(); todo
    int pointsNum = 3 + rand() % 10;
    long arrs[pointsNum][DIM];
    for (auto &arr: arrs)
        for (int dim = 0; dim < DIM; ++dim)
            arr[dim] = rand() % NUMBERS_RANGE;

    int clientsNum = pointsNum;
    // todo worth reading again, to maybe improve on efficiency:
    //  https://stackoverflow.com/questions/25108854/initializing-the-size-of-a-c-vector
    std::vector<Client> clients;
    clients.reserve(clientsNum);

    /*
     *      client [0] stays empty
     *      client [1] has 1 point - {point0}
     *      client [2] has 2 points - {point0, point1}
     *      ...
     *      client [n] has n points -  {point0, point1, ... , pointN}
     */
    for (int c = 0; c < clientsNum; ++c) {
        clients.emplace_back(Client(server));
        // fixme - without this line, at THIS point, the next lines cause segfault (bcs publicKey).
        Client client(server);         //  WHY??? help needed
        for (int a = 0; a < c; ++a) clients.back().encryptPoint(arrs[a]);
    }

    std::vector<Point> points;
    points.reserve(pointsNum * pointsNum / 2); // preallocate memory
    for (Client &client:clients)
        for (Point &point:client.points)
            points.push_back(point);
    //    points.shrink_to_fit();

    cout << "===========" << endl;
    printNameVal(pointsNum);
    printNameVal(points.size());
    Point sum2 = Point::addManyPoints(points);
    for (int dim = 0; dim < DIM - 1; ++dim) printNameVal(server.decryptNum(sum2[dim]));

    //    for (Point &point:points) {
    //        cout << "( ";
    //        //        for (int dim = 0; dim < DIM - 1; ++dim)
    //        //            cout << server.decryptNum(point.cCoordinates[dim]) << ",";
    ////        cout << server.decryptNum(point[0]) << " ) " << endl;
    ////        server.decryptNum(point[0]);
    //    }


    cout << " ------ testRetrievePoints finished ------ " << endl << endl;
}

void TestDataServer::testComparePoints() {
    //    loggerTestClient.log("testComparePoints");
    cout << " ------ testComparePoints ------ " << endl << endl;

    KeysServer server = KeysServer();
    Client client1(server);
    Client client2(server);
    const long arr1[] = {1L, 2L};
    const long arr2[] = {2L, 2L};
    client1.encryptPoint(arr1);
    client2.encryptPoint(arr2);


    //            client1.compare(client2);
    //    helib::decryptBinaryNums();
    //    loggerTestClient.print_log();
    cout << " ------ testComparePoints finished ------ " << endl << endl;
}

//void TestDataServer::testDecryptCoordinates() {
//    cout << " ------ testDecryptCoordinates ------ " << endl << endl;
//    cout << " ------ testDecryptCoordinates finished ------ " << endl << endl;
//
//}
