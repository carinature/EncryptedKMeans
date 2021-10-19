
#include <random>

#include "TestDataServer.h"

//#include "utils/Logger.h"
#include "utils/aux.h"
#include "src/DataServer.h"

//static Logger loggerTestDataServer(log_debug, "loggerTestDataServer");

void TestDataServer::testConstructor() {
    //    loggerTestDataServer.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestDataServer::testComparePoints() {
    //    loggerTestClient.log("testComparePoints");
    cout << " ------ testComparePoints ------ " << endl << endl;

    KeysServer server;
    long arr[DIM], arr2[DIM];
    for (short dim = 0; dim < DIM; ++dim) {
        arr[dim] = rand() % NUMBERS_RANGE;
        arr2[dim] = rand() % NUMBERS_RANGE;
    }
    Point point(server.getPublicKey(), arr);
    Point point2(server.getPublicKey(), arr2);

    for (short dim = 0; dim < DIM; ++dim) {
        //        for (int i = 0; i < number_of_points; ++i) {
        helib::Ctxt res = point.isBiggerThan(point2, dim)[0];
        assert((arr[dim] > arr2[dim]) == server.decryptCtxt(res));

        helib::Ctxt res2 = point2.isBiggerThan(point, dim)[0];
        assert((arr2[dim] > arr[dim]) == server.decryptCtxt(res2));
        //        }
    }
    cout << " ------ testComparePoints finished ------ " << endl << endl;
}

void TestDataServer::testAddition() {
    cout << " ------ testAddition ------ " << endl << endl;
    KeysServer server;
    std::vector<Client> clients = generateDataClients(server);
    std::vector<Point> points = DataServer::retrievePoints(clients);
    printPoints(points, server);
    Point sum = Point::addManyPoints(points);
    for (short dim = 0; dim < DIM; ++dim) {
        printNameVal(server.decryptNum(sum[dim]));
    }
    cout << " ------ testAddition finished ------ " << endl << endl;
}

void TestDataServer::testMultiplication() {
    cout << " ------ testMultiplication ------ " << endl << endl;
    cout << " ------ testMultiplication finished ------ " << endl << endl;
}

void TestDataServer::testGenerateDataClients() {
    cout << " ------ testGenerateDataClients ------ " << endl << endl;
    KeysServer keysServer;

    std::vector<Client> clients = generateDataClients(keysServer);
    //    std::vector<Point> points = DataServer::retrievePoints(clients);
    //
    //    cout << " --- Points  ---" << endl;
    //    printPoints(points, keysServer);
    //    cout << " --- --- --- --- ---" << endl;
    cout << " ------ testGenerateDataClients finished ------ " << endl << endl;
}

void TestDataServer::testRetrievePoints() {
    cout << " ------ testRetrievePoints ------ " << endl << endl;
    KeysServer keysServer;

    std::vector<Client> clients = generateDataClients(keysServer);
    std::vector<Point> points = DataServer::retrievePoints(clients);

    cout << " --- Points  ---" << endl;
    printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    cout << " ------ testRetrievePoints finished ------ " << endl << endl;
}

void TestDataServer::testminirand() {
    // Fill vector with numbers 0,1,2,3...,kMaxValue
    static const int kMaxValue = 7;
    std::vector<int> v(kMaxValue + 1);
    for (size_t i = 0; i < v.size(); ++i)
        v[i] = i;

    // Random seed generator
    std::random_device rd;

    // Psuedo random number generator
    std::mt19937 prng(rd());

    // Shuffle the 0,1,2,3...,kMaxValue integer sequence
    std::shuffle(v.begin(), v.end(), prng);

    // Print random sequence
    for (int x : v)
        std::cout << x << ' ';
    std::cout << std::endl;

    cout << " ==== and another one ===  " << endl;
    auto rng = std::default_random_engine{rd()};
    std::vector<int> v2(kMaxValue + 1);
    for (size_t i = 0; i < v2.size(); ++i)
        v2[i] = i;
    std::shuffle(std::begin(v2), std::end(v2), rng);
    // Print random sequence
    for (int x : v2)
        std::cout << x << ' ';
    std::cout << std::endl;

    /*
        cout << " ==== and with points ===  " << endl;
        KeysServer keysServer;
        std::vector<Client> clients = generateDataClients(keysServer);
        std::vector<Point> points = DataServer::retrievePoints(clients);
        std::shuffle(std::begin(points), std::end(points), rng);
        // Print random sequence
        for (Point &point: points) printPoint(point, keysServer);
        std::cout << std::endl;
    //    error: no matching function for call to ‘swap(Point&, Point&)’
    //    swap(*__a, *__b);
    */


}

void TestDataServer::testPickRandomPoints() {
    cout << " ------ testPickRandomPoints ------ " << endl << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<Client> clients = generateDataClients(keysServer);
    std::vector<Point> points = DataServer::retrievePoints(clients);
    cout << " --- All Points  ---" << endl;
    printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    std::vector<std::vector<Point>> randomPoints = dataServer.pickRandomPoints(points,
                                                                               1 / epsilon);
    cout << " --- Random Points  ---" << endl;
    for (auto vec :randomPoints) printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    cout << " ------ testPickRandomPoints finished ------ " << endl << endl;
}

void TestDataServer::testCreateCmpDict() {
    cout << " ------ testCreateCmpDict ------ " << endl << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<Client> clients = generateDataClients(keysServer);
    std::vector<Point> points = DataServer::retrievePoints(clients);
    cout << " --- All Points  ---" << endl;
    printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    std::vector<std::vector<Point>> randomPoints = dataServer.pickRandomPoints(points,
                                                                               1 / epsilon);
    cout << " --- Random Points  ---" << endl;
    for (auto vec :randomPoints) printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    std::vector<
            std::unordered_map<
                    const Point,
                    std::unordered_map<
                            const Point,
                            helib::Ctxt> > >
            cmp = dataServer.createCmpDict(points, randomPoints);

    cout << "The Dictionary: " << endl;
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "    ======   ";
        printNameVal(DIM);// << " ======" << endl;
        for (auto const&[point, map] : cmp[dim]) {
            printPoint(point, keysServer);
            cout << endl;
            for (auto const&[point2, val]: map) {
                printPoint(point2, keysServer);
                printNameVal(keysServer.decryptCtxt(val));
            }
            printNameVal(map.size());
            cout << " --- --- ---" << endl;
        }
        printNameVal(cmp[dim].size());
        cout << " === === ===" << endl;
    }

    cout << " ------ testCreateCmpDict finished ------ " << endl << endl;
}

void TestDataServer::testSplit() {
    cout << " ------ testSplit ------ " << endl << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<Client> clients = generateDataClients(keysServer);
    std::vector<Point> points = DataServer::retrievePoints(clients);

    cout << " --- All Points  ---" << endl;
    printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    const Point &tinyRandomPoint = keysServer.tinyRandomPoint();
    std::vector<std::vector<Point>>
            randomPoints = dataServer.pickRandomPoints(points, 1 / epsilon);
    cout << " --- Random Points  ---" << endl;
    for (auto vec :randomPoints) printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    std::vector<
            std::unordered_map<
                    const Point,
                    std::unordered_map<
                            const Point,
                            helib::Ctxt> > >
            cmpDict = dataServer.createCmpDict(points, randomPoints);

    /*
     * cout << "The Dictionary: " << endl;
        for (int dim = 0; dim < DIM; ++dim) {
            cout << "    ======   ";
            printNameVal(DIM);// << " ======" << endl;
            for (auto const&[point, map] : cmp[dim]) {
                printPoint(point, keysServer);
                cout << endl;
                for (auto const&[point2, val]: map) {
                    printPoint(point2, keysServer);
                    printNameVal(keysServer.decryptCtxt(val));
                }
                printNameVal(map.size());
                cout << " --- --- ---" << endl;
            }
            printNameVal(cmp[dim].size());
            cout << " === === ===" << endl;
        }
        */

    std::map<int, //DIM
            std::vector< //current slices for approp dimension
                    Slice
            >
    >
            groups = dataServer.splitIntoEpsNet(points, randomPoints, cmpDict, keysServer);
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "   ---   For dim " << dim << "  --- " << endl;
        for (Slice &cell: groups[dim]) {
            cell.printSlice(keysServer);
        }
        cout << "   ---     --- " << endl;
        cout << endl;
    }

    cout << " ------ testSplit finished ------ " << endl << endl;


}
