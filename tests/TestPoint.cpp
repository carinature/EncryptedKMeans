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
    for (long &c :arr) c = rand() % bitSizeRange;
    Point point(server.getPublicKey(), arr);
    cout << " ------ testEncryptCoordinates finished ------ " << endl << endl;
}

void TestPoint::testOperatorSubscript() {
    //    loggerTestClient.log("testOperatorSubscript");
    cout << " ------ testOperatorSubscript ------ " << endl;
    KeysServer server;
    long arr[DIM];
    for (long &c :arr) c = rand() % bitSizeRange;
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
    for (long &c :arr) c = rand() % bitSizeRange;
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
        arr[dim] = rand() % bitSizeRange;
        arr2[dim] = rand() % bitSizeRange;
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
    //    loggerTestClient.log("testAddition");
    cout << " ------ testAddManyPoints ------ " << endl;
    KeysServer server;
    long arr[DIM], arr2[DIM], arr3[DIM], arrSum[DIM];
    for (int dim = 0; dim < DIM; ++dim) {
        arr[dim] = rand() % (bitSizeRange/2); //fixme should be full range?
        arr2[dim] = rand() % (bitSizeRange/2);
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
//        printNameVal(arrSum[dim]);
//        printNameVal(server.decryptNum(sum[dim]));
        assert(arrSum[dim] == server.decryptNum(sum[dim]));
    }
    cout << " ------ testAddManyPoints finished ------ " << endl << endl;
}

void TestPoint::testMultiplication() {
    //    loggerTestClient.log("testMultiplication");
    cout << " ------ testMultiplication ------ " << endl;
    KeysServer server;
    long arr[DIM], arr2[DIM], arrProd[DIM];
    for (int dim = 0; dim < DIM; ++dim) {
        arr[dim] = rand() % (bitSizeRange / 8);
        arr2[dim] = rand() % (bitSizeRange / 8);
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
    cout << " ------ testAddition finished ------ " << endl << endl;
}

void TestPoint::testMultiplicationByBit() {
    //    loggerTestClient.log("testMultiplicationByBit");
    cout << " ------ testMultiplicationByBit ------ " << endl;
    KeysServer server;
    long arr[DIM], arr0[] = {0, 0};
    for (int dim = 0; dim < DIM; ++dim) {
        arr[dim] = rand() % bitSizeRange;
    }
    //  this option is good for readability
    //      `const helib::PubKey &pubKey = server.getPublicKey();`
    //      but the way it is now (asking the serv for a pubKey for each point)
    //      makes more sense in a "real world" application
    Point point(server.getPublicKey(), arr);
    //    Point point2(server.getPublicKey(), arr2);
    helib::Ctxt bit0 = server.encryptCtxt(0), bit1 = server.encryptCtxt(1);
    Point product = point * bit1;
    Point product0 = point * bit0;
    for (int dim = 0; dim < DIM; ++dim) {
        assert(arr[dim] == server.decryptNum(product[dim]));
        assert(arr0[dim] == server.decryptNum(product0[dim]));
    }
    cout << " ------ testMultiplicationByBit finished ------ " << endl << endl;
}

