
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
        arr[dim] = dist(mt);
        arr2[dim] = dist(mt);
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

    std::vector<Client> clients = generateDataClients(keysServer);
    std::vector<Point> points = DataServer::retrievePoints(clients);

    cout << " --- Points  ---" << endl;
    printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    cout << " ------ testRetrievePoints finished ------ " << endl << endl;
}

void TestDataServer::testPickRandomPoints() {
    cout << " ------ testPickRandomPoints ------ " << endl << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<Client>
            clients = generateDataClients(keysServer);
    std::vector<Point>
            points = DataServer::retrievePoints(clients);

    cout << " --- All Points  ---" << endl;
    printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    std::vector<std::vector<Point>>
            randomPoints = dataServer.pickRandomPoints(points); //, 1 / EPSILON);

    cout << " --- Random Points  ---" << endl;
    for (auto vec :randomPoints) printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    cout << " ------ testPickRandomPoints finished ------ " << endl << endl;
}

void TestDataServer::testCreateCmpDict() {
    cout << " ------ testCreateCmpDict ------ " << endl << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<Client>
            clients = generateDataClients(keysServer);
    std::vector<Point>
            points = DataServer::retrievePoints(clients);

    cout << " --- All Points  ---" << endl;
    printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    std::vector<std::vector<Point>>
            randomPoints = dataServer.pickRandomPoints(points);//, 1 / EPSILON);

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
    for (short dim = 0; dim < DIM; ++dim) {
        cout << "    ======   ";
        printNameVal(dim);// << " ======" << endl;
        for (auto const&[point, map] : cmp[dim]) {
            printPoint(point, keysServer);
            cout << endl;

            long p1c = keysServer.decryptNum(point[dim]);

            for (auto const&[point2, val]: map) {

                long p2c = keysServer.decryptNum(point2[dim]);

                long pVal = keysServer.decryptCtxt(val);
                assert(pVal == (p1c > p2c) || pVal == (p1c == p2c));

                printPoint(point2, keysServer);
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

void TestDataServer::testSplitIntoEpsNet() {
    cout << " ------ testSplitIntoEpsNet ------ " << endl;// << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<Client>
            clients = generateDataClients(keysServer);
    std::vector<Point>
            points = DataServer::retrievePoints(clients);
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
    printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    cout << " --- Random Points  ---" << endl;
    for (auto vec :randomPoints) printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    cout << "The Dictionary: " << endl;
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "    ======   ";
        printNameVal(dim);// << " ======" << endl;
        for (auto const&[point, map] : cmpDict[dim]) {
            printPoint(point, keysServer);
            cout << endl;
            for (auto const&[point2, val]: map) {
                printPoint(point2, keysServer);
                printNameVal(keysServer.decryptCtxt(val));
            }
            printNameVal(map.size());
            cout << " --- --- ---" << endl;
        }
        printNameVal(cmpDict[dim].size());
        cout << " === === ===" << endl;
    }

    std::map<int, std::vector<Slice> >
            slices = dataServer.splitIntoEpsNet(points, randomPoints, cmpDict, keysServer);
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "   ---   For dim " << dim << "  --- " << endl;
        for (Slice &cell: slices[dim]) cell.printSlice(keysServer);
        cout << "   ---     --- " << endl;
        cout << endl;
    }

    cout << " ------ testSplitIntoEpsNet finished ------ " << endl << endl;

}

void TestDataServer::testCalculateCellMeans() {
    cout << " ------ testCalculateCellMeans ------ " << endl;
    KeysServer keysServer;
    DataServer dataServer(keysServer);

    std::vector<Client> clients = generateDataClients(keysServer);
    std::vector<Point> points = DataServer::retrievePoints(clients);
    std::vector<std::vector<Point> > randomPoints = dataServer.pickRandomPoints(points);
    std::vector<std::unordered_map<const Point, std::unordered_map<const Point, helib::Ctxt> > >
            cmpDict = dataServer.createCmpDict(points, randomPoints);

    std::map<int, std::vector<Slice> >
            epsNet = dataServer.splitIntoEpsNet(points, randomPoints, cmpDict, keysServer);

    std::vector<std::tuple<Point, Slice> >
            meanCellTuples = DataServer::calculateSlicesMeans(epsNet[DIM - 1], keysServer);

    cout << " === All Points  ===" << endl;
    printPoints(points, keysServer);
    cout << " === === === === ===" << endl;
    cout << " === Random Points  ===" << endl;
    for (auto const &vec :randomPoints) printPoints(vec, keysServer);
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
        printPoint(std::get<0>(tup), keysServer);
        cout << endl;
        std::get<1>(tup).printSlice(keysServer);
    }
    cout << endl;

    cout << " ------ testCalculateCellMeans finished ------ " << endl << endl;
}


void TestDataServer::testGetMinimalDistances() {
    cout << " ------ testGetMinimalDistances ------ " << endl;
    const KeysServer keysServer;
    const DataServer dataServer(keysServer);

    //  Creating Data
    int n = NUMBER_OF_POINTS, m = 1 / EPSILON;
    std::vector<Point> points, points2;
    points.reserve(n);
    points2.reserve(DIM);
    long tempArrs[n][DIM], tempArrs2[m][DIM];
    for (int i = 0; i < n; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs[i][dim] = dist(mt);
        points.emplace_back(Point(keysServer.getPublicKey(), tempArrs[i]));
    }
    for (int i = 0; i < m; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs2[i][dim] = dist(mt);
        points2.emplace_back(Point(keysServer.getPublicKey(), tempArrs2[i]));
    }
    std::vector<Point> dummyMeans(points2.begin(), points2.begin() + DIM);

    cout << endl << "Points: ";
    printPoints(points, keysServer);
    cout << endl << "Dummies: ";
    printPoints(dummyMeans, keysServer);
    cout << endl;

    //  Calculating Algorithm
    const std::vector<std::tuple<Point, Point, EncryptedNum> >
            minDistanceTuples =
            DataServer::collectMinimalDistancesAndClosestPoints(
                    points,
                    dummyMeans,
                    keysServer);

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

        //        cout << " Point: ";
        //        printPoint(std::get<0>(minDistanceTuples[i]), keysServer);
        Point closestMean = std::get<1>(minDistanceTuples[i]);
        Point pClosestMean = dummyMeans[closestMinIndex];
        //        cout << " pClosestMean: ";
        //        printPoint(pClosestMean, keysServer);
        //        cout << " closestMean: ";
        //        printPoint(closestMean, keysServer);
        //        cout << endl;
        assert(keysServer.decryptNum(pClosestMean.cid) == keysServer.decryptNum(closestMean.cid));

        long minDistance = keysServer.decryptNum(std::get<2>(minDistanceTuples[i]));
        //        printNameVal(pMinDistance);
        //        printNameVal(minDistance);
        assert(pMinDistance == minDistance);
    }

    cout << " ------ testGetMinimalDistances finished ------ " << endl << endl;
}

void TestDataServer::testCalculateThreshold() {

}
