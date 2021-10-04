//
// Created by karina on 29/08/2021.
//

#include <src/Client.h>
#include "TestPoint.h"

#include "utils/Logger.h"
#include "src/Point.h"
//#include "src/Client.h"

void TestPoint::testConstructor() {
    //    loggerTestClient.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl;
    KeysServer server;
    Point point(server.getPublicKey());
    //    loggerTestClient.print_log();
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestPoint::testEncryptCoordinates() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testEncryptCoordinates ------ " << endl;
    KeysServer server;
    long arr[DIM];
    for (long &c :arr) c = rand() % NUMBERS_RANGE;
    Point point(server.getPublicKey(), arr);
    cout << " ------ testEncryptCoordinates finished ------ " << endl << endl;
}

void TestPoint::testOperatorSubscript() {
    //    loggerTestClient.log("testOperatorSubscript");
    cout << " ------ testOperatorSubscript ------ " << endl;
    KeysServer server;
    long arr[DIM];
    for (long &c :arr) c = rand() % NUMBERS_RANGE;
    Point point(server.getPublicKey(), arr);
    std::vector<long> coordinate;
    for (short int dim = 0; dim < DIM; ++dim) {
        std::vector<Ctxt> cCoordiante = point[dim];
        helib::decryptBinaryNums(coordinate,
                                 helib::CtPtrs_vectorCt(cCoordiante),
                                 server.getSecKey(),
                                 server.getEA());
        assert(coordinate.back() == arr[dim]);
        assert(cCoordiante == point.cCoordinates[dim]);
    }
    cout << " ------ testOperatorSubscript finished ------ " << endl << endl;

}

void TestPoint::testIsEmpty() {
    //    loggerTestClient.log("testConstructor");
    cout << " ------ testIsEmpty ------ " << endl;
    KeysServer server;
    Point emptyPoint(server.getPublicKey());
    assert(true == emptyPoint.isEmpty());

    long arr[DIM];
    for (long &c :arr) c = rand() % NUMBERS_RANGE;
    Point point(server.getPublicKey(), arr);
    assert(false == point.isEmpty());

    cout << " ------ testIsEmpty finished ------ " << endl << endl;
}

void TestPoint::testAddition() {
    //    loggerTestClient.log("testAddition");
    cout << " ------ testAddition ------ " << endl;
    KeysServer server;
    long arr[DIM], arr2[DIM], arrSum[DIM];
    for (int dim = 0; dim < DIM; ++dim) {
        arr[dim] = rand() % NUMBERS_RANGE;
        arr2[dim] = rand() % NUMBERS_RANGE;
        arrSum[dim] = arr[dim] + arr2[dim];
        //            for (auto a: arr) printNameVal(a);
    }
    //  this option is good for readability
    //      `const helib::PubKey &pubKey = server.getPublicKey();`
    //      but the way it is now (asking the serv for a pubKey for each point)
    //      makes more sense in a "real world" application
    Point point(server.getPublicKey(), arr);
    Point point2(server.getPublicKey(), arr2);
    Point sum = point + point2;
    for (int dim = 0; dim < DIM; ++dim)
        assert(arrSum[dim] == server.decryptNum(sum[dim]));
    cout << " ------ testAddition finished ------ " << endl << endl;
}

void TestPoint::testAddManyPoints() {
    //    loggerTestClient.log("testAddManyPoints");
    cout << " ------ testAddManyPoints ------ " << endl;
    KeysServer server;
    long arr[DIM], arr2[DIM], arr3[DIM], arrSum[DIM];
    for (int dim = 0; dim < DIM; ++dim) {
        arr[dim] = rand() % (NUMBERS_RANGE / 2); //fixme should be full range?
        arr2[dim] = rand() % (NUMBERS_RANGE / 2);
        arr3[dim] = arr[dim] + arr2[dim];
        arrSum[dim] = arr[dim] + arr2[dim] + arr3[dim];
    }
    //    for (auto a: arr) printNameVal(a);
    //    for (auto a: arr2) printNameVal(a);
    //    for (auto a: arr3) printNameVal(a);
    //    for (auto a: arrSum) printNameVal(a);
    Point sum = Point::addManyPoints({
                                             Point(server.getPublicKey(), arr),
                                             Point(server.getPublicKey(), arr2),
                                             Point(server.getPublicKey(), arr3)
                                     });
    for (int dim = 0; dim < DIM; ++dim) {
                printNameVal(arrSum[dim]);
                printNameVal(server.decryptNum(sum[dim]));
        assert(arrSum[dim] == server.decryptNum(sum[dim]));
    }

    std::vector<Point> points;
    for (int i = 0; i < 40; ++i) {
        points.emplace_back(Point(server.getPublicKey(), arr));
        points.emplace_back(Point(server.getPublicKey(), arr2));
        points.emplace_back(Point(server.getPublicKey(), arr3));
        points.emplace_back(Point(server.getPublicKey(), arrSum));
    }
    Point sum2 = Point::addManyPoints(points);
    printNameVal(server.decryptNum(sum2[0]));
    printNameVal(arrSum[0]);
    printNameVal(server.decryptNum(sum2[1]));
    printNameVal(arrSum[1]);
    for (int dim = 0; dim < DIM; ++dim) {
                printNameVal(arrSum[dim]);
                printNameVal(80*arrSum[dim]);
                printNameVal(server.decryptNum(sum2[dim]));
        assert(80*arrSum[dim] == server.decryptNum(sum2[dim]));
        //fixme this part of the test acts weird when running after _all_ the other tests.
        // incidentally, looks like the 2nd coordinate of sum is just a copy of the 1st.
        // could be problem w/ the decryption result?
    }
    cout << " ------ testAddManyPoints finished ------ " << endl << endl;
}

void TestPoint::testMultiplication() {
    //    loggerTestClient.log("testMultiplication");
    cout << " ------ testMultiplication ------ " << endl;
    KeysServer server;
    long arr[DIM], arr2[DIM], arrProd[DIM];
    for (int dim = 0; dim < DIM; ++dim) {
        arr[dim] = rand() % (NUMBERS_RANGE / 8);
        arr2[dim] = rand() % (NUMBERS_RANGE / 8);
        arrProd[dim] = arr[dim] * arr2[dim];
    }
    //  this option is good for readability
    //      `const helib::PubKey &pubKey = server.getPublicKey();`
    //      but the way it is now (asking the serv for a pubKey for each point)
    //      makes more sense in a "real world" application
    Point point(server.getPublicKey(), arr);
    Point point2(server.getPublicKey(), arr2);
    Point product = point * point2;
    for (int dim = 0; dim < DIM; ++dim)
        assert(arrProd[dim] == server.decryptNum(product[dim]));
    cout << " ------ testMultiplication finished ------ " << endl << endl;
}

void TestPoint::testMultiplicationByBit() {
    //    loggerTestClient.log("testMultiplicationByBit");
    cout << " ------ testMultiplicationByBit ------ " << endl;
    KeysServer server;
    long arr[DIM], arr0[] = {0, 0};
    for (int dim = 0; dim < DIM; ++dim) arr[dim] = rand() % NUMBERS_RANGE;
    //  this option is good for readability
    //      `const helib::PubKey &pubKey = server.getPublicKey();`
    //      but the way it is now (asking the serv for a pubKey for each point)
    //      makes more sense in a "real world" application
    Point point(server.getPublicKey(), arr);
    //    Point point2(server.getPublicKey(), arr2);
    helib::Ctxt bit0 = server.encryptCtxt(false);
    helib::Ctxt bit1 = server.encryptCtxt(true);
    Point product0 = point * bit0;
    Point product1 = point * bit1;
    for (int dim = 0; dim < DIM; ++dim) {
        assert(arr[dim] == server.decryptNum(product1[dim]));
        assert(arr0[dim] == server.decryptNum(product0[dim]));
    }
    cout << " ------ testMultiplicationByBit finished ------ " << endl << endl;
}

void TestPoint::testCompare() {    //    loggerTestClient.log("testCompare");
    cout << " ------ testCompare ------ " << endl;
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
    cout << " ------ testCompare finished ------ " << endl << endl;
}

