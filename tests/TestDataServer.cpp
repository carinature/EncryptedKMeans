
#include <src/coreset/run1meancore.h>
#include "TestDataServer.h"

#include "src/DataServer.h"

void TestDataServer::testConstructor() {
    //    loggerTestDataServer.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestDataServer::testComparePoints() {
    cout << " ------ testComparePoints ------ " << endl << endl;

    KeysServer keysServer;
    long arr[DIM], arr2[DIM];
    for (short dim = 0; dim < DIM; ++dim) {
        arr[dim] = randomLongInRange(mt);
        arr2[dim] = randomLongInRange(mt);
    }
    Point point(keysServer.getPublicKey(), arr);
    Point point2(keysServer.getPublicKey(), arr2);

    for (short dim = 0; dim < DIM; ++dim) {
        //        for (int i = 0; i < NUMBER_OF_POINTS; ++i) {
        helib::Ctxt res = point.isBiggerThan(point2, dim)[0];
        assert((arr[dim] > arr2[dim]) == keysServer.decryptCtxt(res));

        helib::Ctxt res2 = point2.isBiggerThan(point, dim)[0];
        assert((arr2[dim] > arr[dim]) == keysServer.decryptCtxt(res2));
        //        }
    }
    cout << " ------ testComparePoints finished ------ " << endl << endl;
}

void TestDataServer::testRetrievePoints() {
    cout << " ------ testRetrievePoints ------ " << endl << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<ClientDevice> clients = generateDataClients(keysServer);
    std::vector<Point> points = dataServer.retrievePoints(clients);

    cout << " --- Points  ---" << endl;
    cout << printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    decAndWriteToFile(points, POINTS_FILE, keysServer);
    cout << " ------ testRetrievePoints finished ------ " << endl << endl;
}

void TestDataServer::testRetrievePoints_Threads() {
    cout << " ------ testRetrievePoints_Threads ------ " << endl << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<ClientDevice> clients = generateDataClients(keysServer);
    std::vector<Point> points = dataServer.retrievePoints(clients);

    cout << " --- Points  ---" << endl;
    cout << printPoints(points, keysServer);
    printNameVal(points.size());
    cout << " --- --- --- --- ---" << endl;

    std::vector<Point> points_withThreads = dataServer.retrievePoints_WithThreads(clients);
    std::vector<Point> points_withThreads_2 = dataServer.retrievePoints_WithThreads(clients, NUMBER_OF_CLIENTS);


    cout << " --- Points 1 ---" << endl;
    cout << printPoints(points_withThreads, keysServer);
    printNameVal(points_withThreads.size());
    cout << " --- Points 2 ---" << endl;
    cout << printPoints(points_withThreads_2, keysServer);
    printNameVal(points_withThreads_2.size());
    cout << " --- --- --- --- ---" << endl;

    cout << " ------ testRetrievePoints_Threads finished ------ " << endl << endl;
}

void TestDataServer::testPickRandomPoints() {
    cout << " ------ testPickRandomPoints ------ " << endl << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<ClientDevice>
            clients = generateDataClients(keysServer);
    std::vector<Point>
            points = dataServer.retrievePoints(clients);

    cout << " --- All Points  ---" << endl;
    cout << printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    std::vector<std::vector<Point>>
            randomPoints = dataServer.pickRandomPoints(points); //, 1 / EPSILON);

    cout << " --- Random Points  ---" << endl;
    for (auto vec :randomPoints) cout << printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    cout << " ------ testPickRandomPoints finished ------ " << endl << endl;
}

void TestDataServer::testCreateCmpDict() {
    cout << " ------ testCreateCmpDict ------ " << endl << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<ClientDevice>
            clients = generateDataClients(keysServer);
    std::vector<Point>
            points = dataServer.retrievePoints(clients);

    cout << " --- All Points  ---" << endl;
    cout << printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    std::vector<std::vector<Point>>
            randomPoints = dataServer.pickRandomPoints(points);//, 1 / EPSILON);

    cout << " --- Random Points  ---" << endl;
    for (auto vec :randomPoints) cout << printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    std::vector<
            std::unordered_map<
                    const Point,
                    std::unordered_map<
                            const Point,
                            helib::Ctxt> > >
            cmp = dataServer.createCmpDict(points, randomPoints);

    cout << "The Dictionary: " << endl;
    for (short dim = 0; dim < DIM; ++dim) {
        cout << "    ======   ";
        printNameVal(dim);// << " ======" << endl;
        for (auto const&[point, map] : cmp[dim]) {
            cout << printPoint(point, keysServer);
            cout << endl;

            long p1c = keysServer.decryptNum(point[dim]);

            for (auto const&[point2, val]: map) {

                long p2c = keysServer.decryptNum(point2[dim]);

                long pVal = keysServer.decryptCtxt(val);
                assert(pVal == (p1c > p2c) || pVal == (p1c == p2c));

                cout << printPoint(point2, keysServer);
                printNameVal(pVal);

            }
            printNameVal(map.size());
            cout << " --- --- ---" << endl;
        }
        printNameVal(cmp[dim].size());
        cout << " === === ===" << endl;
    }

    cout << " ------ testCreateCmpDict finished ------ " << endl << endl;
}

void TestDataServer::testCreateCmpDict_Threads() {
    cout << " ------ testCreateCmpDict_Threads ------ " << endl << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<ClientDevice>
            clients = generateDataClients(keysServer);
    std::vector<Point>
            points = dataServer.retrievePoints(clients);

    cout << " --- All Points  ---" << endl;
    cout << printPoints(points, keysServer);
    printNameVal(points.size());
    cout << " --- --- --- --- ---" << endl;

    dataServer.retrievePoints_WithThreads(clients);

    cout << " --- Points  ---" << endl;
    cout << printPoints(dataServer.retrievedPoints, keysServer);
    printNameVal(dataServer.retrievedPoints.size());
    cout << " --- --- --- --- ---" << endl;

    std::vector<std::vector<Point>>
            randomPoints = dataServer.pickRandomPoints(points);//, 1 / EPSILON);

    cout << " --- Random Points  ---" << endl;
    for (const auto &vec :randomPoints) cout << printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    std::vector<
            std::unordered_map<
                    const Point,
                    std::unordered_map<
                            const Point,
                            helib::Ctxt> > >
            cmp = dataServer.createCmpDict(points, randomPoints);
    dataServer.createCmpDict_WithThreads(points, randomPoints, NUMBER_OF_THREADS);

    cout << "The Dictionary: " << endl;
    for (short dim = 0; dim < DIM; ++dim) {
        cout << "    ======   ";
        printNameVal(dim);// << " ======" << endl;
        for (auto const&[point, map] : cmp[dim]) {
            cout << printPoint(point, keysServer);
            cout << endl;

            long p1c = keysServer.decryptNum(point[dim]);

            for (auto const&[point2, val]: map) {

                long p2c = keysServer.decryptNum(point2[dim]);

                long pVal = keysServer.decryptCtxt(val);
                assert(pVal == (p1c > p2c) || pVal == (p1c == p2c));

                cout << printPoint(point2, keysServer);
                printNameVal(pVal);

            }
            printNameVal(map.size());
            cout << " --- --- ---" << endl;
        }
        printNameVal(cmp[dim].size());
        cout << " === === ===" << endl;
    }

    cout << "The Dictionary - With Threads: " << endl;
    for (short dim = 0; dim < DIM; ++dim) {
        cout << "    ======   ";
        printNameVal(dim);// << " ======" << endl;
        for (auto const&[point, map] : dataServer.cmpDict[dim]) {
            cout << printPoint(point, keysServer);
            cout << endl;

            long p1c = keysServer.decryptNum(point[dim]);

            for (auto const&[point2, val]: map) {

                long p2c = keysServer.decryptNum(point2[dim]);

                long pVal = keysServer.decryptCtxt(val);
                assert(pVal == (p1c > p2c) || pVal == (p1c == p2c));

                cout << printPoint(point2, keysServer);
                printNameVal(pVal);

            }
            printNameVal(map.size());
            cout << " --- --- ---" << endl;
        }
        printNameVal(dataServer.cmpDict[dim].size());
        cout << " === === ===" << endl;
    }

    cout << " ------ testCreateCmpDict_Threads finished ------ " << endl << endl;
}

void TestDataServer::testSplitIntoEpsNet() {
    cout << " ------ testSplitIntoEpsNet ------ " << endl;// << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<ClientDevice>
            clients = generateDataClients(keysServer);
    std::vector<Point>
            points = dataServer.retrievePoints(clients);
    std::vector<std::vector<Point>>
            randomPoints = dataServer.pickRandomPoints(points);//, 1 / EPSILON);

    std::vector<
            std::unordered_map<
                    const Point,
                    std::unordered_map<
                            const Point,
                            helib::Ctxt> > >
            cmpDict = dataServer.createCmpDict(points, randomPoints);

    cout << " --- All Points  ---" << endl;
    cout << printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    cout << " --- Random Points  ---" << endl;
    for (auto vec :randomPoints) cout << printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    cout << "The Dictionary: " << endl;
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "    ======   ";
        printNameVal(dim);// << " ======" << endl;
        for (auto const&[point, map] : cmpDict[dim]) {
            cout << printPoint(point, keysServer);
            cout << endl;
            for (auto const&[point2, val]: map) {
                cout << printPoint(point2, keysServer);
                printNameVal(keysServer.decryptCtxt(val));
            }
            printNameVal(map.size());
            cout << " --- --- ---" << endl;
        }
        printNameVal(cmpDict[dim].size());
        cout << " === === ===" << endl;
    }

    std::map<int, std::vector<Slice> >
            slices = dataServer.splitIntoEpsNet(points, randomPoints, cmpDict);
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "   ---   For dim " << dim << "  --- " << endl;
        for (Slice &cell: slices[dim]) cell.printSlice(keysServer);
        cout << "   ---     --- " << endl;
        cout << endl;
    }

    cout << " ------ testSplitIntoEpsNet finished ------ " << endl << endl;

}

void TestDataServer::testSplitIntoEpsNet_WithThreads() {
    cout << " ------ testSplitIntoEpsNet_WithThreads ------ " << endl;// << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<ClientDevice> clients = generateDataClients(keysServer);

    std::vector<Point> points = dataServer.retrievePoints(clients);
    dataServer.retrievePoints_WithThreads(clients);

    std::vector<std::vector<Point>> randomPoints = dataServer.pickRandomPoints(points);

    std::vector<
            std::unordered_map<
                    const Point,
                    std::unordered_map<
                            const Point,
                            helib::Ctxt> > >
            cmpDict = dataServer.createCmpDict(points, randomPoints);
    dataServer.createCmpDict_WithThreads(points, randomPoints, NUMBER_OF_THREADS);

    cout << " --- All Points  ---" << endl;
    cout << printPoints(points, keysServer);
    printNameVal(points.size());
    cout << " --- --- --- --- ---" << endl;

    cout << " --- Random Points  ---" << endl;
    for (auto vec :randomPoints) cout << printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    cout << "The Dictionary: " << endl;
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "    ======   ";
        printNameVal(dim);// << " ======" << endl;
        for (auto const&[point, map] : cmpDict[dim]) {
            cout << printPoint(point, keysServer);
            cout << endl;
            for (auto const&[point2, val]: map) {
                cout << printPoint(point2, keysServer);
                printNameVal(keysServer.decryptCtxt(val));
            }
            printNameVal(map.size());
            cout << " --- --- ---" << endl;
        }
        printNameVal(cmpDict[dim].size());
        cout << " === === ===" << endl;
    }

    std::map<int, std::vector<Slice> >
            slices = dataServer.splitIntoEpsNet(points, randomPoints, cmpDict);
    std::map<int, std::vector<Slice> >
            slices_Threads = dataServer.splitIntoEpsNet_WithThreads(points, randomPoints, cmpDict);

    for (int dim = 0; dim < DIM; ++dim) {
        cout << "   ---   For dim " << dim << "  --- " << endl;
        for (Slice &cell: slices[dim]) cell.printSlice(keysServer);
        cout << "   ---     --- " << endl;
        cout << endl;
    }
    cout << "==================\n";
    cout << "Slices With Treads\n";
    cout << "==================\n";
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "   ---   For dim " << dim << "  --- " << endl;
        for (Slice &cell: slices_Threads[dim]) cell.printSlice(keysServer);
        cout << "   ---     --- " << endl;
        cout << endl;
    }

    cout << " ------ testSplitIntoEpsNet_WithThreads finished ------ " << endl << endl;

}

void TestDataServer::testCalculateCellMeans() {
    cout << " ------ testCalculateCellMeans ------ " << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<ClientDevice> clients = generateDataClients(keysServer);
    std::vector<Point> points = dataServer.retrievePoints(clients);
    std::vector<std::vector<Point> > randomPoints = dataServer.pickRandomPoints(points);
    std::vector<std::unordered_map<const Point, std::unordered_map<const Point, helib::Ctxt> > >
            cmpDict = dataServer.createCmpDict(points, randomPoints);

    std::map<int, std::vector<Slice> >
            epsNet = dataServer.splitIntoEpsNet(points, randomPoints, cmpDict);

    std::vector<std::tuple<Point, Slice> >
            meanCellTuples = dataServer.calculateSlicesMeans(epsNet[DIM - 1]);

    cout << " === All Points  ===" << endl;
    cout << printPoints(points, keysServer);
    cout << " === === === === ===" << endl;
    cout << " === Random Points  ===" << endl;
    for (auto const &vec :randomPoints) cout << printPoints(vec, keysServer);
    cout << " === === === === ===" << endl;
    cout << " === Slices  ===" << endl;
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "   ---   For dim " << dim << "  --- " << endl;
        for (const Slice &slice: epsNet[dim]) slice.printSlice(keysServer);
        cout << endl;
    }
    cout << " === === === === ===" << endl;
    cout << " === === === === ===" << endl;
    cout << " === === === === ===" << endl;
    for (auto const &tup:meanCellTuples) {
        cout << "The Mean is: ";
        cout << printPoint(std::get<0>(tup), keysServer);
        cout << endl;
        std::get<1>(tup).printSlice(keysServer);
    }
    cout << endl;

    cout << " ------ testCalculateCellMeans finished ------ " << endl << endl;
}

void TestDataServer::testCalculateCellMeans_WithThreads() {
    cout << " ------ testCalculateCellMeans ------ " << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<ClientDevice> clients = generateDataClients(keysServer);
    //  collect points
    std::vector<Point> points = dataServer.retrievePoints(clients);
    dataServer.retrievePoints_WithThreads(clients);
    std::vector<Point> points_withThreads = dataServer.retrievedPoints;
    //  random points
    std::vector<std::vector<Point> > randomPoints = dataServer.pickRandomPoints(points);
    std::vector<std::vector<Point> > randomPoints_forThreads = dataServer.randomPointsList;
    //  compare dict
    std::vector<std::unordered_map<const Point, std::unordered_map<const Point, helib::Ctxt> > >
            cmpDict = dataServer.createCmpDict(points, randomPoints);
    CmpDict &cmpDict_withThreads = dataServer.createCmpDict_WithThreads(points, randomPoints, NUMBER_OF_THREADS);
    //  EPS net
    std::map<int, std::vector<Slice> >
            epsNet = dataServer.splitIntoEpsNet(points, randomPoints, cmpDict);
    std::map<int, std::vector<Slice> >
            epsNet_Threads = dataServer.splitIntoEpsNet_WithThreads(points, randomPoints, cmpDict);//points, randomPoints, cmpDict, keysServer);
    //  eps-net means
    std::vector<std::tuple<Point, Slice> >
            meanCellTuples = dataServer.calculateSlicesMeans(epsNet[DIM - 1]);
    std::vector<std::tuple<Point, Slice> >
            meanCellTuples_Threads = dataServer.calculateSlicesMeans_WithThreads(
            epsNet_Threads[DIM - 1]);

    cout << " === === === === ===" << endl;
    cout << " === All Points  ===" << endl;
    cout << " === === === === ===" << endl;
    cout << printPoints(points, keysServer);
    cout << " === === === === ===" << endl;
    cout << " === Random Points  ===" << endl;
    cout << " === === === === ===" << endl;
    for (auto const &vec :randomPoints) cout << printPoints(vec, keysServer);
    cout << " === === === === ===" << endl;
    cout << " === Slices  ===" << endl;
    cout << " === === === === ===" << endl;
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "   ---   For dim " << dim << "  --- " << endl;
        for (const Slice &slice: epsNet[dim]) slice.printSlice(keysServer);
        cout << endl;
    }
    cout << " === === === === ===" << endl;
    cout << " === Slices With Threads ===" << endl;
    cout << " === === === === ===" << endl;
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "   ---   For dim " << dim << "  --- " << endl;
        for (const Slice &slice: epsNet_Threads[dim]) slice.printSlice(keysServer);
        cout << endl;
    }
    cout << " === === === === ===" << endl;
    cout << " === === Means === ===" << endl;
    cout << " === === === === ===" << endl;
    for (auto const &tup:meanCellTuples) {
        cout << "The Mean is: ";
        cout << printPoint(std::get<0>(tup), keysServer);
        cout << endl;
        std::get<1>(tup).printSlice(keysServer);
    }
    cout << endl;
    cout << " === === === === === ===" << endl;
    cout << " === Means w/ Threads ===" << endl;
    cout << " === === === === === ===" << endl;
    for (auto const &tup:meanCellTuples_Threads) {
        cout << "The Mean is: ";
        cout << printPoint(std::get<0>(tup), keysServer);
        cout << endl;
        std::get<1>(tup).printSlice(keysServer);
    }
    cout << endl;

    cout << " ------ testCalculateCellMeans finished ------ " << endl << endl;
}

void TestDataServer::testGetMinimalDistances() {
    cout << " ------ testGetMinimalDistances ------ " << endl;
    const KeysServer keysServer;
    DataServer dataServer(keysServer);

    //  Creating Data
    int n = NUMBER_OF_POINTS, m = 1 / EPSILON;
    std::vector<Point> points, points2;
    points.reserve(n);
    points2.reserve(DIM);
    long tempArrs[n][DIM], tempArrs2[m][DIM];
    for (int i = 0; i < n; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs[i][dim] = randomLongInRange(mt);
        points.emplace_back(Point(keysServer.getPublicKey(), tempArrs[i]));
    }
    for (int i = 0; i < m; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs2[i][dim] = randomLongInRange(mt);
        points2.emplace_back(Point(keysServer.getPublicKey(), tempArrs2[i]));
    }
    std::vector<Point> dummyMeans(points2.begin(), points2.begin() + 1 / EPSILON);

    cout << endl << "Points: ";
    cout << printPoints(points, keysServer);
    cout << endl << "Dummies: ";
    cout << printPoints(dummyMeans, keysServer);
    cout << endl;

    //  Calculating Algorithm
    const std::vector<std::tuple<Point, Point, EncryptedNum> >
            minDistanceTuples =
            dataServer.collectMinimalDistancesAndClosestPoints(
                    points,
                    dummyMeans);

    //  Check results
    for (int i = 0; i < points.size(); ++i) {
        long currCoors[DIM];
        for (int dim = 0; dim < DIM; ++dim) currCoors[dim] = tempArrs[i][dim];
        long pMinDistance = DIM * std::pow(NUMBERS_RANGE, 2);
        int closestMinIndex = -1;
        for (int j = 0; j < dummyMeans.size(); ++j) {
            long tempsum = 0;
            for (int dim = 0; dim < DIM; ++dim)
                tempsum += pow(currCoors[dim] - tempArrs2[j][dim], 2);
            if (pMinDistance > tempsum) {
                pMinDistance = tempsum;
                closestMinIndex = j;
            }
        }

        const Point &closestMean = std::get<1>(minDistanceTuples[i]);
        const Point &pClosestMean = dummyMeans[closestMinIndex];
/*
        cout << " Point: ";
        cout << printPoint(std::get<0>(minDistanceTuples[i]), keysServer);
        cout << " pClosestMean: ";
        cout << printPoint(pClosestMean, keysServer);
        cout << " closestMean: ";
        cout << printPoint(closestMean, keysServer);
        cout << endl;
*/

        assert(keysServer.decryptNum(pClosestMean.cid) == keysServer.decryptNum(closestMean.cid));

        long minDistance = keysServer.decryptNum(std::get<2>(minDistanceTuples[i]));
        //        printNameVal(pMinDistance);
        //        printNameVal(minDistance);
        assert(pMinDistance == minDistance);
    }

    cout << " ------ testGetMinimalDistances finished ------ " << endl << endl;
}

void TestDataServer::testGetMinimalDistances_WithThreads() {
    cout << " ------ testGetMinimalDistances_WithThreads ------ " << endl;
    const KeysServer keysServer;
    DataServer dataServer(keysServer);

    //  Creating Data
    int n = NUMBER_OF_POINTS, m = 1 / EPSILON;
    std::vector<Point> points, points2;
    points.reserve(n);
    points2.reserve(DIM);
    long tempArrs[n][DIM], tempArrs2[m][DIM];
    for (int i = 0; i < m; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs2[i][dim] = randomLongInRange(mt);
        points2.emplace_back(Point(keysServer.getPublicKey(), tempArrs2[i]));
    }
    std::vector<Point> dummyMeans(points2.begin(), points2.begin() + DIM);


    //  Retrieve expected results
    std::vector<long> pMinDistances(n);
    std::vector<int> closestMinIndices(n);


    for (int i = 0; i < n; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs[i][dim] = randomLongInRange(mt);
        points.emplace_back(Point(keysServer.getPublicKey(), tempArrs[i]));
        //    }
        //    for (int i = 0; i < points.size(); ++i) {
        //        long currCoors[DIM];
        //        for (int dim = 0; dim < DIM; ++dim) currCoors[dim] = tempArrs[i][dim];
        long pMinDistance = DIM * std::pow(NUMBERS_RANGE, 2);
        int closestMinIndex = -1;
        for (int j = 0; j < dummyMeans.size(); ++j) {
            long tempsum = 0;
            for (int dim = 0; dim < DIM; ++dim)
                tempsum += pow(tempArrs[i][dim] - tempArrs2[j][dim], 2);
            if (pMinDistance > tempsum) {
                pMinDistance = tempsum;
                closestMinIndex = j;
            }
        }
        //        printNameVal(pMinDistance);
        pMinDistances[i] = pMinDistance;
        closestMinIndices[i] = closestMinIndex;
        /*                for (int j = 0; j < pMinDistances.size(); ++j) printNameVal(pMinDistances[j]);
                        printNameVal(pMinDistances[i]);
                        cout << "-------" << endl;*/

    }

    //  Calculating Algorithm
    const std::vector<std::tuple<Point, Point, EncryptedNum> >
            minDistanceTuples =
            dataServer.collectMinimalDistancesAndClosestPoints(
                    points,
                    dummyMeans);

    const std::vector<std::tuple<Point, Point, EncryptedNum> >
            minDistanceTuples_WithThreads =
            dataServer.collectMinimalDistancesAndClosestPoints_WithThreads(
                    points,
                    dummyMeans
            );
    /*
    cout << endl << "Points: ";
    cout << printPoints(points, keysServer);
    cout << endl << "Dummies: ";
    cout << printPoints(dummyMeans, keysServer);
    cout << endl;
    */

    //  Check results
    for (int i = 0; i < points.size(); ++i) {
        const Point &closestMean = std::get<1>(minDistanceTuples[i]);
        const Point &pClosestMean = dummyMeans[closestMinIndices[i]];
        const Point &pointWithThreads = std::get<0>(minDistanceTuples_WithThreads[i]);
        const Point &closestMean_WithThreads = std::get<1>(minDistanceTuples_WithThreads[i]);
        Point &approPoint = points[0];
        int appropI = -1;
        for (int p = 0; p < points.size(); ++p)
            if (points[p].id == pointWithThreads.id) {
                approPoint = points[p];
                appropI = p;
                break;
            }
        Point &pClosestMean_WithThreads = dummyMeans[closestMinIndices[appropI]];
        /*
           cout << " Point: ";
           cout << printPoint(std::get<0>(minDistanceTuples[i]), keysServer);
           cout << " pClosestMean: ";
           cout << printPoint(pClosestMean, keysServer);
           cout << " closestMean: ";
           cout << printPoint(closestMean, keysServer);
           cout << " Point: with threads:";
           cout << printPoint(pointWithThreads, keysServer);
           cout << " closestMean_withThreads: ";
           cout << printPoint(closestMean_WithThreads, keysServer);
           cout << " Point by id:";
           cout << printPoint(approPoint, keysServer);
           cout << " Point by index:";
           cout << printPoint(points[appropI], keysServer);
           printNameVal(appropI);
           cout << "pClosestMean_WithThreads: ";
           cout << printPoint(pClosestMean_WithThreads, keysServer);
           cout << endl;
           */
        long pClosestMeanCid = keysServer.decryptNum(pClosestMean.cid);
        long closestMeanCid = keysServer.decryptNum(closestMean.cid);
        long pClosestMean_WithThreads_Cid = keysServer.decryptNum(pClosestMean_WithThreads.cid);
        long closestMean_WithThreadsCid = keysServer.decryptNum(closestMean_WithThreads.cid);
        assert(pClosestMeanCid == closestMeanCid);
        assert(pClosestMean_WithThreads_Cid == closestMean_WithThreadsCid);

        long minDistance = keysServer.decryptNum(std::get<2>(minDistanceTuples[i]));
        long minDistance_WithThreads = keysServer.decryptNum(
                std::get<2>(minDistanceTuples_WithThreads[i]));
        /*
        printNameVal(pMinDistances[i]);
        printNameVal(pMinDistances[appropI]);
        printNameVal(minDistance);
        printNameVal(minDistance_WithThreads);
        */
        assert(pMinDistances[i] == minDistance);
        assert(pMinDistances[appropI] == minDistance_WithThreads);
    }

    cout << " ------ testGetMinimalDistances_WithThreads finished ------ " << endl << endl;
}

void TestDataServer::testCalculateThreshold() {
    cout << " ------ testCalculateThreshold ------ " << endl;
    const KeysServer keysServer;
    DataServer dataServer(keysServer);

    //  Creating Data
    int n = NUMBER_OF_POINTS, m = 1 / EPSILON;
    std::vector<Point> points, points2;
    points.reserve(n);
    points2.reserve(DIM);
    long tempArrs[n][DIM], tempArrs2[m][DIM];
    for (int i = 0; i < n; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs[i][dim] = randomLongInRange(mt);
        points.emplace_back(Point(keysServer.getPublicKey(), tempArrs[i]));
    }
    for (int i = 0; i < m; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs2[i][dim] = randomLongInRange(mt);
        points2.emplace_back(Point(keysServer.getPublicKey(), tempArrs2[i]));
    }
    std::vector<Point> dummyMeans(points2.begin(), points2.begin() + DIM);

    //  Calculating Algorithm
    const std::vector<std::tuple<Point, Point, EncryptedNum> >
            minDistanceTuples =
            dataServer.collectMinimalDistancesAndClosestPoints(
                    points,
                    dummyMeans);

    const EncryptedNum
            threshold =
            dataServer.calculateThreshold(
                    minDistanceTuples,
                    0);

    //  Print Results
    cout << endl << "Points: ";
    cout << printPoints(points, keysServer);
    cout << endl << "Dummies: ";
    cout << printPoints(dummyMeans, keysServer);
    cout << endl;

    for (int i = 0; i < points.size(); ++i) {
        cout << " Point: ";
        cout << printPoint(std::get<0>(minDistanceTuples[i]), keysServer);
        Point closestMean = std::get<1>(minDistanceTuples[i]);
        cout << " closestMean: ";
        cout << printPoint(closestMean, keysServer);
        cout << endl;
        long minDistance = keysServer.decryptNum(std::get<2>(minDistanceTuples[i]));
        printNameVal(minDistance);
    }
    printNameVal(keysServer.decryptNum(threshold));

    cout << " ------ testCalculateThreshold finished ------ " << endl << endl;
}

void TestDataServer::testChoosePointsByDistance() {
    cout << " ------ testChoosePointsByDistance ------ " << endl;
    const KeysServer keysServer;
    DataServer dataServer(keysServer);

    //  Creating Data
    int n = NUMBER_OF_POINTS, m = 1 / EPSILON;
    std::vector<Point> points, points2;
    points.reserve(n);
    points2.reserve(DIM);
    long tempArrs[n][DIM], tempArrs2[m][DIM];
    for (int i = 0; i < n; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs[i][dim] = randomLongInRange(mt);
        points.emplace_back(Point(keysServer.getPublicKey(), tempArrs[i]));
    }
    for (int i = 0; i < m; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs2[i][dim] = randomLongInRange(mt);
        points2.emplace_back(Point(keysServer.getPublicKey(), tempArrs2[i]));
    }
    std::vector<Point> dummyMeans(points2.begin(), points2.begin() + DIM);

    //  Calculating Algorithm
    const std::vector<std::tuple<Point, Point, EncryptedNum> >
            minDistanceTuples =
            dataServer.collectMinimalDistancesAndClosestPoints(
                    points,
                    dummyMeans);

    EncryptedNum
            threshold = dataServer.calculateThreshold(minDistanceTuples, 0);

    cout << endl << "Points: ";
    cout << printPoints(points, keysServer);
    cout << endl << "Dummies: ";
    cout << printPoints(dummyMeans, keysServer);
    cout << endl;

    for (int i = 0; i < points.size(); ++i) {
        cout << " Point: ";
        cout << printPoint(std::get<0>(minDistanceTuples[i]), keysServer);
        Point closestMean = std::get<1>(minDistanceTuples[i]);
        cout << " closestMean: ";
        cout << printPoint(closestMean, keysServer);
        cout << endl;
        long minDistance = keysServer.decryptNum(std::get<2>(minDistanceTuples[i]));
        printNameVal(minDistance);
    }
    printNameVal(keysServer.decryptNum(threshold));

    std::tuple<
            std::unordered_map<long, std::vector<std::pair<Point, CBit> > >,
            std::vector<std::pair<Point, CBit> >,
            std::vector<std::pair<Point, CBit> >
    > groups = dataServer.choosePointsByDistance(
            minDistanceTuples,
            dummyMeans,
            threshold
    );

    cout << " === Groups by Means === " << endl;
    for (auto const &[meanI, points] : std::get<0>(groups)) {
        printNameVal(meanI) << "\tClose Points: ";
        for (auto const &pair: points) {
            cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
            cout << printPoint(pair.first, keysServer);
        }
        cout << endl;
    }
    cout << " === === === === === " << endl << endl;

    cout << " === Closest Points === " << endl;
    for (auto const &pair : std::get<1>(groups)) {
        cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
        cout << printPoint(pair.first, keysServer);
    }
    cout << endl << " === === === === === " << endl << endl;

    cout << " === Farthest Points === " << endl;
    for (auto const &pair : std::get<2>(groups)) {
        cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
        cout << printPoint(pair.first, keysServer);
    }
    cout << endl << " === === === === === " << endl << endl;


    cout << " ------ testChoosePointsByDistance finished ------ " << endl << endl;
}

void TestDataServer::testChoosePointsByDistance_WithThreads() {
    cout << " ------ testChoosePointsByDistance_WithThreads ------ " << endl;
    const KeysServer keysServer;
    DataServer dataServer(keysServer);

    //  Creating Data
    int n = NUMBER_OF_POINTS, m = 1 / EPSILON;
    std::vector<Point> points, points2;
    points.reserve(n);
    points2.reserve(DIM);
    long tempArrs[n][DIM], tempArrs2[m][DIM];
    for (int i = 0; i < n; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs[i][dim] = randomLongInRange(mt);
        points.emplace_back(Point(keysServer.getPublicKey(), tempArrs[i]));
    }
    for (int i = 0; i < m; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs2[i][dim] = randomLongInRange(mt);
        points2.emplace_back(Point(keysServer.getPublicKey(), tempArrs2[i]));
    }
    std::vector<Point> dummyMeans(points2.begin(), points2.begin() + DIM);

    //  Calculating Algorithm
    const std::vector<std::tuple<Point, Point, EncryptedNum> >
            minDistanceTuples =
            dataServer.collectMinimalDistancesAndClosestPoints_WithThreads( //fixme ?
                    points,
                    dummyMeans
            );

    EncryptedNum
            threshold = dataServer.calculateThreshold(minDistanceTuples, 0);

    cout << endl << "Points: ";
    cout << printPoints(points, keysServer);
    cout << endl << "Dummies: ";
    cout << printPoints(dummyMeans, keysServer);
    cout << endl;

    for (int i = 0; i < points.size(); ++i) {
        cout << " Point: ";
        cout << printPoint(std::get<0>(minDistanceTuples[i]), keysServer);
        Point closestMean = std::get<1>(minDistanceTuples[i]);
        cout << " closestMean: ";
        cout << printPoint(closestMean, keysServer);
        cout << endl;
        long minDistance = keysServer.decryptNum(std::get<2>(minDistanceTuples[i]));
        printNameVal(minDistance);
    }
    printNameVal(keysServer.decryptNum(threshold));

    std::tuple<
            std::unordered_map<long, std::vector<std::pair<Point, CBit> > >,
            std::vector<std::pair<Point, CBit> >,
            std::vector<std::pair<Point, CBit> >
    > groups = dataServer.choosePointsByDistance(
            minDistanceTuples,
            dummyMeans,
            threshold
    );

/*
    std::tuple<
            std::unordered_map<long, std::vector<std::pair<Point, CBit> > >,
            std::vector<std::pair<Point, CBit> >
    > groups_withThreads1 = dataServer.choosePointsByDistance_WithThreads_slower(
            minDistanceTuples,
            dummyMeans,
            threshold
    );
*/


    std::tuple<
            std::unordered_map<long, std::vector<std::pair<Point, CBit> > >,
            std::vector<std::pair<Point, CBit> >
    > groups_withThreads2 = dataServer.choosePointsByDistance_WithThreads(
            minDistanceTuples,
            dummyMeans,
            threshold
    );

    cout << " === === === === === " << endl;
    cout << " === Groups by Means === " << endl;
    cout << " === === === === === " << endl;
    for (auto const &[meanI, points] : std::get<0>(groups)) {
        printNameVal(meanI) << "\tClose Points: ";
        for (auto const &pair: points) {
            cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
            cout << printPoint(pair.first, keysServer);
        }
        cout << endl;
    }

/*
    cout << " === === === === === === === === === " << endl;
    cout << " === Groups by Means - With Threads  - 111  === " << endl;
    cout << " === === === === === === === === === " << endl;
    for (auto const &[meanI, points] : std::get<0>(groups_withThreads1)) {
        printNameVal(meanI) << "\tClose Points: ";
        for (auto const &pair: points) {
            cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
            cout << printPoint(pair.first, keysServer);
        }
        cout << endl;
    }
*/

    cout << " === === === === === === === === === " << endl;
    cout << " === Groups by Means - With Threads  - 222  === " << endl;
    cout << " === === === === === === === === === " << endl;
    for (auto const &[meanI, points] : std::get<0>(groups_withThreads2)) {
        printNameVal(meanI) << "\tClose Points: ";
        for (auto const &pair: points) {
            cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
            cout << printPoint(pair.first, keysServer);
        }
        cout << endl;
    }
    cout << " === === === === === " << endl << endl;
    /*
        cout << " === === === === === " << endl;
        cout << " === Closest Points === " << endl;
        cout << " === === === === === " << endl;
        for (auto const &pair : std::get<1>(groups)) {
            cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
            cout << printPoint(pair.first, keysServer);
        }
        cout << endl << " === === === === === " << endl << endl;

        cout << " === === === === === === === === === " << endl;
        cout << " === Closest Points - With Threads === " << endl;
        cout << " === === === === === === === === === " << endl;
        for (auto const &pair : std::get<1>(groups_withThreads)) {
            cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
            cout << printPoint(pair.first, keysServer);
        }
        cout << endl << " === === === === === " << endl << endl;
    */
    cout << " === === === === === " << endl;
    cout << " === Farthest Points === " << endl;
    cout << " === === === === === " << endl;
    for (auto const &pair : std::get<2>(groups)) {
        cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
        cout << printPoint(pair.first, keysServer);
    }
    cout << endl << " === === === === === " << endl << endl;

/*
    cout << " === === === === === === === === === " << endl;
    cout << " === Farthest Points - With Threads  -  111  === " << endl;
    cout << " === === === === === === === === === " << endl;
    for (auto const &pair : std::get<1>(groups_withThreads1)) {
        cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
        cout << printPoint(pair.first, keysServer);
    }
    cout << endl << " === === === === === " << endl << endl;
*/

    cout << " === === === === === === === === === " << endl;
    cout << " === Farthest Points - With Threads  -  222  === " << endl;
    cout << " === === === === === === === === === " << endl;
    for (auto const &pair : std::get<1>(groups_withThreads2)) {
        cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
        cout << printPoint(pair.first, keysServer);
    }
    cout << endl << " === === === === === " << endl << endl;

    cout << " ------ testChoosePointsByDistance_WithThreads finished ------ " << endl << endl;

}

void TestDataServer::testChoosePointsByDistance_WithThreads_withYonis() {
    cout << " ------ testChoosePointsByDistance_WithThreads_withYonis ------ " << endl;
    auto t0_test = CLOCK::now();     //  for logging, profiling, DBG

    const KeysServer keysServer;
    DataServer dataServer(keysServer);

    //  Creating Data
    int n = NUMBER_OF_POINTS, m = 1 / EPSILON;
    std::vector<Point> points, points2;
    points.reserve(n);
    points2.reserve(DIM);
    long tempArrs[n][DIM], tempArrs2[m][DIM];
    for (int i = 0; i < n; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs[i][dim] = randomLongInRange(mt);
        points.emplace_back(Point(keysServer.getPublicKey(), tempArrs[i]));
    }
    for (int i = 0; i < m; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs2[i][dim] = randomLongInRange(mt);
        points2.emplace_back(Point(keysServer.getPublicKey(), tempArrs2[i]));
    }
    std::vector<Point> dummyMeans(points2.begin(), points2.begin() + DIM);

    //  Calculating Algorithm
    const std::vector<std::tuple<Point, Point, EncryptedNum> >
            minDistanceTuples =
            dataServer.collectMinimalDistancesAndClosestPoints_WithThreads(
                    points,
                    dummyMeans
            );

    EncryptedNum
            threshold = dataServer.calculateThreshold(minDistanceTuples, 0);

    for (int i = 0; i < points.size(); ++i) {
        cout << " Point: ";
        cout << printPoint(std::get<0>(minDistanceTuples[i]), keysServer);
        Point closestMean = std::get<1>(minDistanceTuples[i]);
        cout << " closestMean: ";
        cout << printPoint(closestMean, keysServer);
        cout << endl;
        long minDistance = keysServer.decryptNum(std::get<2>(minDistanceTuples[i]));
        printNameVal(minDistance);
    }
    printNameVal(keysServer.decryptNum(threshold));

    std::tuple<
            std::unordered_map<long, std::vector<std::pair<Point, CBit> > >,
            std::vector<std::pair<Point, CBit> >
    > groups = dataServer.choosePointsByDistance_WithThreads(
            minDistanceTuples,
            dummyMeans,
            threshold
    );

    cout << endl << "Points: ";
    cout << printPoints(points, keysServer);
    cout << endl << "Dummies: ";
    cout << printPoints(dummyMeans, keysServer);
    cout << endl;
    cout << endl;

    cout << " === Groups by Means === " << endl;
    for (auto const &[meanI, points] : std::get<0>(groups)) {
        cout << printPoint(dummyMeans[meanI], keysServer);
        printNameVal(meanI) << "Close Points: ";
        std::vector<Point> pointsGroup;
        pointsGroup.reserve(points.size());
        for (auto const &pair: points) {
            // todo do not remove yet
//            cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
//            cout << printPoint(pair.first, keysServer);
            pointsGroup.emplace_back(pair.first);
        }
        cout << endl;
        printNonEmptyPoints(pointsGroup, keysServer);
        cout << endl;
        decAndWriteToFile(pointsGroup, IO_DIR+to_string(meanI)+"_"+CHOSEN_FILE, keysServer);

    }
    cout << " === === === === === " << endl << endl;


    cout << " === === === === === === === === === " << endl;
    cout << " === Farthest Points - With Threads  === " << endl;
    cout << " === === === === === === === === === " << endl;
    std::vector<std::pair<Point, CBit> > farthest = std::get<1>(groups);
    std::vector<Point> leftover;
    leftover.reserve(points.size());
    for (auto const &[point, isIn] : farthest) {
        cout << "-" << keysServer.decryptCtxt(isIn) << "-";
        cout << printPoint(point, keysServer);
        leftover.emplace_back(point);
    }
    leftover.shrink_to_fit();
    cout << endl << " === === === === === " << endl << endl;


    printDuration(t0_test,
                  "testChoosePointsByDistance_WithThreads_withYonis - part 1");


    /**********   integration to coreset alg    ***********/
    //  for each mean-group in groups
    for (auto const &[meanI, points] : std::get<0>(groups)) {
        std::vector<std::vector<double>> pointsGroup;
        pointsGroup.reserve(points.size());
        std::vector<Point> forClean;
        // add all non-zero points
        for (auto const &pair: points) {
            const Point &currPoint(pair.first);
            std::vector<long> tempPoint = decryptPoint(currPoint, keysServer);
            std::vector<double> doublePoint(DIM);
            long checkNull = 0;
            for (const long &coordinate:decryptPoint(currPoint, keysServer)) {
                checkNull += coordinate;
                doublePoint.emplace_back(coordinate / CONVERSION_FACTOR);
            }
            if (checkNull) pointsGroup.emplace_back(doublePoint);
            forClean.emplace_back(currPoint);
        }
        pointsGroup.shrink_to_fit();
        cout << "For Group of Points:\t";
        printNonEmptyPoints(forClean, keysServer);
        //      call yoni's alg with mean-group
        cout<< "\nRunning 1-Mean Coreset Algorithm:"<<endl;

        auto t0_itr_rep = CLOCK::now();     //  for logging, profiling, DBG
        runCoreset(pointsGroup, pointsGroup.size(), DIM, EPSILON);  // <|-------------
        printDuration(t0_itr_rep, "runCoreset in iteration");

    }
    printDuration(t0_test, "testChoosePointsByDistance_WithThreads_withYonis");


    // CT for coreset.csv result for multiple iterations
    auto t = time(nullptr);
    auto tm = *localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y_%m_%d_%H_%M_%S");
    string timestamp = oss.str();
    //    oss.flush();
    oss.clear();
//    string filename = "io/" + timestamp + "_coreset.csv";
    string prefix = IO_DIR + timestamp + "_coreset.csv";

    decAndWriteToFile(points, IO_DIR+POINTS_FILE, keysServer);
//    decAndWriteToFile(randomPoints[DIM - 1], RANDS_FILE, keysServer);
    decAndWriteToFile(dummyMeans, IO_DIR+MEANS_FILE, keysServer);
    decAndWriteToFile(leftover, IO_DIR+LEFTOVER_FILE, keysServer);


    cout << " ------ testChoosePointsByDistance_WithThreads_withYonis finished ------ " << endl << endl;

}

