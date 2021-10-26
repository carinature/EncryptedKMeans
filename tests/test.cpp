//
// Created by karina on 25/07/2021.
//

#include "TestKeysServer.h"
#include "TestPoint.h"
#include "TestClient.h"
#include "TestAux.h"
#include "TestDataServer.h"

#include "utils/aux.h"

int main() {

//    cout << " ============ Test KeysServer ============ " << endl;
//    TestKeysServer::testConstructor();
//    TestKeysServer::testEncryptCtxt();
//    TestKeysServer::testDecryptCtxt();
//    TestKeysServer::testEncryptNum();   //  todo
//    TestKeysServer::testDecryptNum();   //  todo
//    TestKeysServer::testScratchPoint();
//    TestKeysServer::testTinyRandomPoint();
//    cout << " ============ Test KeysServer Finished ============ " << endl << endl;
//
    cout << " ============ Test Point ============ " << endl;
//    TestPoint::testConstructor();
//    TestPoint::testEncryptCoordinates();
//    TestPoint::testOperatorSubscript();
//    TestPoint::testIsEmpty();
//    TestPoint::testAddition();
//    TestPoint::testAddManyPoints();
//    TestPoint::testMultiplication();
//    TestPoint::testMultiplicationByBit();
//    TestPoint::testCompare();
//    TestPoint::testCalculateDistanceFromPoint();
    TestPoint::testFindMinimalDistancesFromMeans();
//    cout << " ============ Test Point Finished ============ " << endl << endl;
//
//    cout << " ============ Test Client ============ " << endl;
//    TestClient::testConstructor();
//    TestClient::testEncryptCoordinates();
//    TestClient::testDecryptCoordinates();
//    TestClient::testEncryptScratchPoint();
//    TestClient::testCompare();
//    cout << " ============ Test Client Finished ============ " << endl << endl;
//
//    cout << " ============ Test Aux ============ " << endl;
//    TestAux::testGenerateDataClients();
//    cout << " ============ Test Client Finished ============ " << endl << endl;
//
    cout << " ============ Test DataServer ============ " << endl;
//    TestDataServer::testConstructor();
//    TestDataServer::testComparePoints();
//    TestDataServer::testAddition();
//    TestDataServer::testMultiplication();
//    TestDataServer::testGenerateDataClients();
//    TestDataServer::testRetrievePoints();
//    TestDataServer::testminirand();
//    TestDataServer::testPickRandomPoints();
//    TestDataServer::testCreateCmpDict();
//    TestDataServer::testSplitIntoEpsNet();
    TestDataServer::testCalculateCellMeans(); //FIXED! fixme weird bug when huge #points
    cout << " ============ Test DataServer Finished ============ " << endl << endl;


}