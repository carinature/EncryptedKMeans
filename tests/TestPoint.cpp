//
// Created by karina on 29/08/2021.
//

#include "TestPoint.h"

#include "utils/Logger.h"
#include "src/Client.h"

void TestPoint::testConstructor() {
    //    loggerTestClient.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl;
    KeysServer server;// = KeysServer();
    Point point(server.getPublicKey());
    //    loggerTestClient.print_log();
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestPoint::testEncryptCoordinates() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testEncryptCoordinates ------ " << endl;
    KeysServer server;// = KeysServer();
    long arr[DIM];
    for (long &c :arr) c = rand() % bitSizeRange;
    Point point(server.getPublicKey(), arr);
    cout << " ------ testEncryptCoordinates finished ------ " << endl << endl;
}

void TestPoint::testOperatorSubscript() {
    //    loggerTestClient.log("testOperatorSubscript");
    cout << " ------ testOperatorSubscript ------ " << endl;
    KeysServer server;// = KeysServer();
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
        assert(cCoordiante==point.cCoordinates[dim]);
    }
    cout << " ------ testOperatorSubscript finished ------ " << endl << endl;

}

void TestPoint::testIsEmpty() {
    //    loggerTestClient.log("testConstructor");
    cout << " ------ testIsEmpty ------ " << endl;
    KeysServer server;// = KeysServer();
    Point emptyPoint(server.getPublicKey());
    assert(true == emptyPoint.isEmpty());

    long arr[DIM];
    for (long &c :arr) c = rand() % bitSizeRange;
    Point point(server.getPublicKey(), arr);
    assert(false == point.isEmpty());

    cout << " ------ testIsEmpty finished ------ " << endl << endl;
}

void TestPoint::testAddition() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testAddition ------ " << endl;
    KeysServer server;// = KeysServer();
    long arr[DIM];
    for (long &c :arr) c = rand() % bitSizeRange;
    for (auto a: arr) printNameVal(a);
    Point point(server.getPublicKey(), arr);
    Point point2(server.getPublicKey(), arr);
//    Point sum = point+point2;
    cout << " ------ testAddition finished ------ " << endl << endl;
//
//    // Let's encrypt something!
//    KeysServerCKKS serverCKKS;// = KeysServer();
//
//    const helib::Context &context = serverCKKS.getContext();
//    cout <<"odin"<<endl;
////    helib::PubKey publicKey = serverCKKS.getPublicKeyCKKS();
//cout <<"dva "<<endl;
//    long n = 6;//context.getNSlots();
//    double PI=3.14;
//    std::vector<double> v(n);
//    for (long i = 0; i < n; i++)
//        v[i] = sin(2.0 * PI * i / n);
//    cout <<"tri "<<endl;
//    helib::PtxtArray p(context, v);
//    cout <<"chetiri "<<endl;
////    Ctxt c(publicKey);
////    cout <<"pyats]"<<endl;
////    p.encrypt(c);
////    cout <<"shests "<<endl;

}

void TestPoint::testMultiplication() {}
