

#include "TestPoint.h"

#include "src/Point.h"

void TestPoint::testConstructor() {
    cout << " ------ testConstructor ------ " << endl;
    KeysServer keysServer;

    Point point(keysServer.getPublicKey());
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestPoint::testEncryptCoordinates() {
    cout << " ------ testEncryptCoordinates ------ " << endl;
    KeysServer keysServer;

    long arr[DIM];
    for (long &a :arr) a = randomLongInRange(mt);

    Point point(keysServer.getPublicKey(), arr);
    cout << " ------ testEncryptCoordinates finished ------ " << endl << endl;
}

void TestPoint::testOperatorSubscript() {
    cout << " ------ testOperatorSubscript ------ " << endl;
    KeysServer keysServer;

    long arr[DIM];
    for (long &a :arr) a = randomLongInRange(mt);
    Point point(keysServer.getPublicKey(), arr);

    for (short dim = 0; dim < DIM; ++dim) {
        EncryptedNum cCoordiante(point[dim]);
        assert(cCoordiante == point.cCoordinates[dim]);
        long coordinate = keysServer.decryptNum(cCoordiante);
        assert(coordinate == arr[dim]);
    }
    cout << " ------ testOperatorSubscript finished ------ " << endl << endl;

}

void TestPoint::testIsEmpty() {
    //    loggerTestClient.log("testConstructor");
    cout << " ------ testIsEmpty ------ " << endl;
    KeysServer keysServer;
    Point emptyPoint(keysServer.getPublicKey());
    assert(true == emptyPoint.isEmpty());

    long arr[DIM];
    for (long &a :arr) a = randomLongInRange(mt);
    Point point(keysServer.getPublicKey(), arr);
    assert(false == point.isEmpty());

    cout << " ------ testIsEmpty finished ------ " << endl << endl;
}

void TestPoint::testAddition() {
    //    loggerTestClient.log("testAddition");
    cout << " ------ testAddition ------ " << endl;
    KeysServer keysServer;
    long arr[DIM], arr2[DIM], arrSum[DIM];
    for (short dim = 0; dim < DIM; ++dim) {
        arr[dim] = randomLongInRange(mt);
        arr2[dim] = randomLongInRange(mt);
        arrSum[dim] = arr[dim] + arr2[dim];
    }

    //  this option is good for readability
    //      `const helib::PubKey &public_key = keysServer.getPublicKey();`
    //      but the way it is now (asking the serv for a public_key for each point)
    //      makes more sense in a "real world" application
    Point point(keysServer.getPublicKey(), arr);
    const Point point2(keysServer.getPublicKey(), arr2);
    Point sum = point + point2;

    for (short dim = 0; dim < DIM; ++dim)
        assert(arrSum[dim] == keysServer.decryptNum(sum[dim]));

    cout << " ------ testAddition finished ------ " << endl << endl;
}

void TestPoint::testAddManyPoints() {
    //    loggerTestClient.log("testAddManyPoints");
    cout << " ------ testAddManyPoints ------ " << endl;
    KeysServer keysServer;
    long arr[DIM], arr2[DIM], arr3[DIM], arrSum[DIM];
    for (short dim = 0; dim < DIM; ++dim) {
        arr[dim] = randomLongInRange(mt);
        arr2[dim] = randomLongInRange(mt);
        arr3[dim] = randomLongInRange(mt);
        arrSum[dim] = arr[dim] + arr2[dim] + arr3[dim];
    }
    /*
    //    for (auto a: arr) printNameVal(a);
    //    for (auto a: arr2) printNameVal(a);
    //    for (auto a: arr3) printNameVal(a);
    //    for (auto a: arrSum) printNameVal(a);
*/
    std::vector<Point> points;
    points.reserve(4 * NUMBER_OF_POINTS);
    printNameVal(NUMBER_OF_POINTS);
    for (int i = 0; i < NUMBER_OF_POINTS; ++i) {
        points.emplace_back(Point(keysServer.getPublicKey(), arr));
        points.emplace_back(Point(keysServer.getPublicKey(), arr2));
        points.emplace_back(Point(keysServer.getPublicKey(), arr3));
    }
    Point sum(Point::addManyPoints(points, keysServer));

//    printNameVal(keysServer.decryptNum(sum[0]));
//    printNameVal(arrSum[0]);
//    printNameVal(keysServer.decryptNum(sum[1]));
//    printNameVal(arrSum[1]);

    for (short dim = 0; dim < DIM; ++dim) {
        cout << "----" << endl;
                printNameVal(arrSum[dim]);
                printNameVal(NUMBER_OF_POINTS * arrSum[dim]);
                printNameVal(keysServer.decryptNum(sum[dim]));
        assert(NUMBER_OF_POINTS * arrSum[dim] == keysServer.decryptNum(sum[dim]));
        //fixme this part of the test acts weird when running after _all_ the other tests.
        // incidentally, looks like the 2nd coordinate of sum is just a copy of the 1st.
        // could be problem w/ the decryption result?
    }
    cout << " ------ testAddManyPoints finished ------ " << endl << endl;
}

void TestPoint::testMultiplication() {
    cout << " ------ testMultiplication ------ " << endl;
    KeysServer keysServer;

    long arr[DIM], arr2[DIM], arrProd[DIM];
    for (short dim = 0; dim < DIM; ++dim) {
        arr[dim] = randomLongInRange(mt);
        arr2[dim] = randomLongInRange(mt);
        arrProd[dim] = arr[dim] * arr2[dim];
    }

    Point point(keysServer.getPublicKey(), arr);
    Point point2(keysServer.getPublicKey(), arr2);
    Point product = point * point2;

    for (short dim = 0; dim < DIM; ++dim)
        assert(arrProd[dim] == keysServer.decryptNum(product[dim]));
    cout << " ------ testMultiplication finished ------ " << endl << endl;
}

void TestPoint::testMultiplicationByBit() {
    cout << " ------ testMultiplicationByBit ------ " << endl;
    KeysServer keysServer;
    helib::PubKey &publicKey = keysServer.getPublicKey();

    long arr[DIM], arr0[] = {0, 0};
    for (long &a :arr) a = randomLongInRange(mt);

    Point point(keysServer.getPublicKey(), arr);
    helib::Ctxt bit0(publicKey);
    helib::Ctxt bit1(publicKey);

    publicKey.Encrypt(bit0, NTL::ZZX(false));
    publicKey.Encrypt(bit1, NTL::ZZX(true));

    Point product0 = point * bit0;
    Point product1 = point * bit1;
    for (short dim = 0; dim < DIM; ++dim) {
        assert(arr[dim] == keysServer.decryptNum(product1[dim]));
        assert(arr0[dim] == keysServer.decryptNum(product0[dim]));
    }
    cout << " ------ testMultiplicationByBit finished ------ " << endl << endl;
}

void TestPoint::testCompare() {
    cout << " ------ testCompare ------ " << endl;
    KeysServer keysServer;

    long arr[DIM], arr2[DIM];
    for (short dim = 0; dim < DIM; ++dim) {
        arr[dim] = randomLongInRange(mt);
        arr2[dim] = randomLongInRange(mt);
    }

    Point point(keysServer.getPublicKey(), arr);
    Point point2(keysServer.getPublicKey(), arr2);

    for (short dim = 0; dim < DIM; ++dim) {
        for (int i = 0; i < NUMBER_OF_POINTS; ++i) {
            helib::Ctxt res = point.isBiggerThan(point2, dim)[0];
            assert((arr[dim] > arr2[dim]) == keysServer.decryptCtxt(res));

            helib::Ctxt res2 = point2.isBiggerThan(point, dim)[0];
            assert((arr2[dim] > arr[dim]) == keysServer.decryptCtxt(res2));
        }
    }
    cout << " ------ testCompare finished ------ " << endl << endl;
}

void TestPoint::testCalculateDistanceFromPoint() {
    cout << " ------ testCalculateDistanceFromPoint ------ " << endl;
    KeysServer keysServer;
    int n = NUMBER_OF_POINTS;
    std::vector<Point> points;

    long arrs[n][DIM];
    for (int i = 0; i < n; ++i) {
        for (int dim = 0; dim < DIM; ++dim) arrs[i][dim] = randomLongInRange(mt);
        points.emplace_back(Point(keysServer.getPublicKey(), arrs[i]));
    }

    for (int i = 0; i < n - 1; ++i) {
        EncryptedNum distance = points[i].distanceFrom(points[i + 1], keysServer);
        long dDistance = keysServer.decryptNum(distance);

        long pDistSquared = 0;
        for (int dim = 0; dim < DIM; ++dim)
            pDistSquared += std::pow((arrs[i][dim] - arrs[i + 1][dim]), 2);

        //        printPoint(points[i], keysServer);
        //        printPoint(points[i + 1], keysServer);
        //        printNameVal(i) << "  ----------  " << endl;
        //        printNameVal(dDistance);
        //        printNameVal(pDistSquared);
        assert(pDistSquared == dDistance);

    }

    cout << " ------ testCalculateDistanceFromPoint finished ------ " << endl << endl;
}

void TestPoint::testFindMinimalDistancesFromMeans() {
    cout << " ------ testFindMinimalDistancesFromMeans ------ " << endl;
    KeysServer keysServer;
    int n = NUMBER_OF_POINTS;
    std::vector<Point> points;  //    points.reserve(n);
    std::vector<std::pair<Point, long> > distPairs; //    distPairs.reserve(n);

    long arrs[n][DIM], arr[DIM];
    for (int dim = 0; dim < DIM; ++dim) arr[dim] = randomLongInRange(mt);
    Point point(keysServer.getPublicKey(), arr);

    for (int i = 0; i < n; ++i) {
        for (int dim = 0; dim < DIM; ++dim) arrs[i][dim] = randomLongInRange(mt);
        points.push_back(Point(keysServer.getPublicKey(), arrs[i]));
    }

    std::pair<Point, EncryptedNum>
            minimalDistance = point.findMinDistFromMeans(points, keysServer);

    Point &minDistPoint = points[0];
    long pMinDist = DIM * pow(NUMBERS_RANGE, 2);
    long minId = -1;
    EncryptedNum minCid;//(BIT_SIZE, Ctxt(publicKey));

    for (int i = 0; i < n; ++i) {
        long pDist = 0;
        for (int dim = 0; dim < DIM; ++dim) pDist += std::pow(arrs[i][dim] - arr[dim], 2);

        if (pMinDist >= pDist) {
            minDistPoint = points[i];
            pMinDist = pDist;
            minId = points[i].id;
            minCid = points[i].cid;
        }
    }

    assert(pMinDist == keysServer.decryptNum(minimalDistance.second));
    //    printNameVal(minId);
    //    printNameVal(keysServer.decryptNum(minCid));
    //    printNameVal(keysServer.decryptNum(minimalDistance.first.cid));
    //    assert(minId == keysServer.decryptNum(minimalDistance.first.cid));
    assert(keysServer.decryptNum(minCid) == keysServer.decryptNum(minimalDistance.first.cid));
    //    assert(minDistPoint == minimalDistance.first);
    //    assert(minDistPoint.id == minimalDistance.first.id);
    //    assert(minDistPoint.cid == minimalDistance.first.cid);
    //    for (int bit = 0; bit < BIT_SIZE; ++bit) {
    //        printNameVal(minDistPoint.cid[bit] == minimalDistance.first.cid[bit]);
    //        printNameVal(minDistPoint.cid[bit].equalsTo(minimalDistance.first.cid[bit]));
    //    }

    //    printPoint(minDistPoint, keysServer);
    //    printNameVal(pMinDist);// << endl;
    //    printNameVal(minimalDistance.first.id);// << endl;
    //    printNameVal(minDistPoint.id);// << endl;
    //
    //    printNameVal(keysServer.decryptNum(minimalDistance.first.cid));// << endl;
    //    printNameVal(keysServer.decryptNum(minDistPoint.cid));// << endl;


    //  todo check how long it takes
    //  todo also consider creating a cmp dict for cid's
    const helib::PubKey &publicKey = point.public_key;
    Ctxt mu(publicKey), ni(publicKey);
    helib::CtPtrs_vectorCt origMinCid(minCid);
    helib::CtPtrs_vectorCt resMinCid(minimalDistance.first.cid);
    helib::compareTwoNumbers(mu, ni, resMinCid, origMinCid);
    //  mu==ni==0   means that  resMinCid==origMinCid
    assert((0 == keysServer.decryptCtxt(mu)) && (0 == keysServer.decryptCtxt(ni)));

    cout << " ------ testFindMinimalDistancesFromMeans finished ------ " << endl << endl;


}
