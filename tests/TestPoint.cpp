
#include <random>

//#include <src/Client.h>
#include "TestPoint.h"

#include "utils/Logger.h"
#include "src/Point.h"
//#include "src/Client.h"

void TestPoint::testConstructor() {
    //    loggerTestClient.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl;
    KeysServer keysServer;
    Point point(keysServer.getPublicKey());
    //    loggerTestClient.print_log();
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestPoint::testEncryptCoordinates() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testEncryptCoordinates ------ " << endl;
    KeysServer keysServer;
    long arr[DIM];
    for (long &c :arr) c = rand() % NUMBERS_RANGE;
    Point point(keysServer.getPublicKey(), arr);
    cout << " ------ testEncryptCoordinates finished ------ " << endl << endl;
}

void TestPoint::testOperatorSubscript() {
    //    loggerTestClient.log("testOperatorSubscript");
    cout << " ------ testOperatorSubscript ------ " << endl;
    KeysServer keysServer;
    long arr[DIM];
    for (long &c :arr) c = rand() % NUMBERS_RANGE;
    Point point(keysServer.getPublicKey(), arr);
    std::vector<long> coordinate;
    for (short dim = 0; dim < DIM; ++dim) {
        std::vector<Ctxt> cCoordiante = point[dim];
        helib::decryptBinaryNums(coordinate,
                                 helib::CtPtrs_vectorCt(cCoordiante),
                                 keysServer.getSecKey(),
                                 keysServer.getEA());
        assert(coordinate.back() == arr[dim]);
        assert(cCoordiante == point.cCoordinates[dim]);
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
    for (long &c :arr) c = rand() % NUMBERS_RANGE;
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
        arr[dim] = rand() % NUMBERS_RANGE;
        arr2[dim] = rand() % NUMBERS_RANGE;
        arrSum[dim] = arr[dim] + arr2[dim];
        //            for (auto a: arr) printNameVal(a);
    }
    //  this option is good for readability
    //      `const helib::PubKey &pubKey = keysServer.getPublicKey();`
    //      but the way it is now (asking the serv for a pubKey for each point)
    //      makes more sense in a "real world" application
    Point point(keysServer.getPublicKey(), arr);
    Point point2(keysServer.getPublicKey(), arr2);
    Point sum = point + point2;
    for (short dim = 0; dim < DIM; ++dim) {
        printNameVal(arrSum[dim]);
        printNameVal(keysServer.decryptNum(sum[dim]));
        assert(arrSum[dim] == keysServer.decryptNum(sum[dim]));
    }
    cout << " ------ testAddition finished ------ " << endl << endl;
}

void TestPoint::testAddManyPoints() {
    //    loggerTestClient.log("testAddManyPoints");
    cout << " ------ testAddManyPoints ------ " << endl;
    KeysServer keysServer;
    long arr[DIM], arr2[DIM], arr3[DIM], arrSum[DIM];
    for (short dim = 0; dim < DIM; ++dim) {
        arr[dim] = rand() % (NUMBERS_RANGE / 2); //fixme should be full range?
        arr2[dim] = rand() % (NUMBERS_RANGE / 2);
        arr3[dim] = rand() % (NUMBERS_RANGE / 2);
        arrSum[dim] = arr[dim] + arr2[dim] + arr3[dim];
    }
    //    for (auto a: arr) printNameVal(a);
    //    for (auto a: arr2) printNameVal(a);
    //    for (auto a: arr3) printNameVal(a);
    //    for (auto a: arrSum) printNameVal(a);
    /*
    std::vector<Point> points = {
            Point(keysServer.getPublicKey(), arr),
            Point(keysServer.getPublicKey(), arr2),
            Point(keysServer.getPublicKey(), arr3)
    };
    Point sum = Point::addManyPoints(points);
    for (short dim = 0; dim < DIM; ++dim) {
        printNameVal(arrSum[dim]);
        printNameVal(keysServer.decryptNum(sum[dim]));
        assert(arrSum[dim] == keysServer.decryptNum(sum[dim]));
    }
    */

    std::vector<Point> points;
    points.reserve(4 * NUMBER_OF_POINTS);
    printNameVal(NUMBER_OF_POINTS);
    for (int i = 0; i < NUMBER_OF_POINTS; ++i) {
        points.emplace_back(Point(keysServer.getPublicKey(), arr));
        points.emplace_back(Point(keysServer.getPublicKey(), arr2));
        points.emplace_back(Point(keysServer.getPublicKey(), arr3));
        points.emplace_back(Point(keysServer.getPublicKey(), arrSum));
    }
    Point sum(Point::addManyPoints(points));
    //        printNameVal(keysServer.decryptNum(sum[0]));
    //        printNameVal(arrSum[0]);
    //        printNameVal(keysServer.decryptNum(sum[1]));
    //        printNameVal(arrSum[1]);
    for (short dim = 0; dim < DIM; ++dim) {
        cout << "----" << endl;
        printNameVal(arrSum[dim]);
        printNameVal(2 * NUMBER_OF_POINTS * arrSum[dim]);
        printNameVal(keysServer.decryptNum(sum[dim]));
        cout << "----" << endl;
        assert(2 * NUMBER_OF_POINTS * arrSum[dim] == keysServer.decryptNum(sum[dim]));
        //fixme this part of the test acts weird when running after _all_ the other tests.
        // incidentally, looks like the 2nd coordinate of sum is just a copy of the 1st.
        // could be problem w/ the decryption result?
    }
    cout << " ------ testAddManyPoints finished ------ " << endl << endl;
}

void TestPoint::testMultiplication() {
    //    loggerTestClient.log("testMultiplication");
    cout << " ------ testMultiplication ------ " << endl;
    KeysServer keysServer;
    long arr[DIM], arr2[DIM], arrProd[DIM];
    for (short dim = 0; dim < DIM; ++dim) {
        arr[dim] = rand() % (NUMBERS_RANGE / 8);
        arr2[dim] = rand() % (NUMBERS_RANGE / 8);
        arrProd[dim] = arr[dim] * arr2[dim];
    }
    //  this option is good for readability
    //      `const helib::PubKey &pubKey = keysServer.getPublicKey();`
    //      but the way it is now (asking the serv for a pubKey for each point)
    //      makes more sense in a "real world" application
    Point point(keysServer.getPublicKey(), arr);
    Point point2(keysServer.getPublicKey(), arr2);
    Point product = point * point2;
    for (short dim = 0; dim < DIM; ++dim)
        assert(arrProd[dim] == keysServer.decryptNum(product[dim]));
    cout << " ------ testMultiplication finished ------ " << endl << endl;
}

void TestPoint::testMultiplicationByBit() {
    //    loggerTestClient.log("testMultiplicationByBit");
    cout << " ------ testMultiplicationByBit ------ " << endl;
    KeysServer keysServer;
    long arr[DIM], arr0[] = {0, 0};
    for (short dim = 0; dim < DIM; ++dim) arr[dim] = rand() % NUMBERS_RANGE;
    //  this option is good for readability
    //      `const helib::PubKey &pubKey = keysServer.getPublicKey();`
    //      but the way it is now (asking the serv for a pubKey for each point)
    //      makes more sense in a "real world" application
    Point point(keysServer.getPublicKey(), arr);
    //    Point point2(keysServer.getPublicKey(), arr2);
    helib::Ctxt bit0 = keysServer.encryptCtxt(false);
    helib::Ctxt bit1 = keysServer.encryptCtxt(true);
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
        arr[dim] = rand() % NUMBERS_RANGE;
        arr2[dim] = rand() % NUMBERS_RANGE;
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

void TestPoint::testCalculateDistanceFromAnotherPoint() {
    cout << " ------ testCalculateDistanceFromAnotherPoint ------ " << endl;
    KeysServer keysServer;
    std::random_device rd;
    std::mt19937 mt(rd());
    //    std::uniform_real_distribution<double> dist(0, NUMBERS_RANGE);
    std::uniform_int_distribution<long> dist(0, NUMBERS_RANGE);
    std::vector<Point> points;
    int n = 4;
    points.reserve(n);
    long tempArrs[n][DIM];
    for (int i = 0; i < n; ++i) {
        for (int dim = 0; dim < DIM; ++dim) tempArrs[i][dim] = dist(mt);
        points.push_back(Point(keysServer.getPublicKey(), tempArrs[i]));
    }
    printPoints(points, keysServer);
    cout << endl;
    for (int i = 0; i < n - 1; ++i) {
        printNameVal(i)<<"  ----------  "<<endl;
        EncryptedNum distSquared =
                points[i].calculateDistanceFromPoint(points[i + 1],
                                                     keysServer);
        long dDistSquared = keysServer.decryptNum(distSquared);
        printPoint(points[i], keysServer);
        printPoint(points[i+1], keysServer);
        printNameVal(dDistSquared);

        long pDistSquared = 0;
        for (int dim = 0; dim < DIM; ++dim) {
            long diff = tempArrs[i][dim] - tempArrs[i + 1][dim];
            printNameVal(tempArrs[i][dim]);
            printNameVal(tempArrs[i + 1][dim]);
            printNameVal(diff);
            pDistSquared += std::pow(diff, 2);
            printNameVal(pDistSquared);
        }
        printNameVal(pDistSquared);
        assert(pDistSquared == dDistSquared);

    }

    cout << " ------ testCalculateDistanceFromAnotherPoint finished ------ " << endl << endl;
}

