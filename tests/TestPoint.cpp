//
// Created by karina on 29/08/2021.
//

#include "TestPoint.h"

#include "utils/Logger.h"
#include "src/Client.h"

void TestPoint::testConstructor() {
    //    loggerTestClient.log("testConstructor");
    cout << " ------ testConstructor ------ " << endl ;
    KeysServer server;// = KeysServer();
    Point point(server.getPublicKey());
    //    loggerTestClient.print_log();
    cout << " ------ testConstructor finished ------ " << endl << endl;
}

void TestPoint::testEncryptCoordinates() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testEncryptCoordinates ------ " << endl ;
    KeysServer server;// = KeysServer();
    const long arr[] = {12, 34};
    Point point(server.getPublicKey(), arr);
    cout << " ------ testEncryptCoordinates finished ------ " << endl << endl;
}

void TestPoint::testAddition() {
    //    loggerTestClient.log("testDecryptCoordinates");
    cout << " ------ testEncryptCoordinates ------ " << endl ;
    KeysServer server;// = KeysServer();
    const long arr[] = {12, 34};
    Point point(server.getPublicKey(), arr);
    Point point2(server.getPublicKey(), arr);
    Point sum = point+point2;
    cout << " ------ testEncryptCoordinates finished ------ " << endl << endl;
}

void TestPoint::testMultiplication() {

}