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
////    TestKeysServer::testConstructor();
////    TestKeysServer::testEncryptCtxt();
//    TestKeysServer::testDecryptCtxt();
////    TestKeysServer::testEncryptNum();
//    TestKeysServer::testDecryptNum();
////    TestKeysServer::testScratchPoint();
////    TestKeysServer::testTinyRandomPoint();
//    cout << " ============ Test KeysServer Finished ============ " << endl << endl;

//    cout << " ============ Test Point ============ " << endl;
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
//    TestPoint::testFindMinimalDistancesFromMeans();
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
//    TestAux::minitest();
//    TestAux::minitest2();
//    TestAux::testBGVPackedArithmetics_Original();
//    TestAux::testBGVPackedArithmetics__Comparison();
//    TestAux::testMultithreading();
//    cout << " ============ Test Client Finished ============ " << endl << endl;

    cout << " ============ Test DataServer ============ " << endl;
//    TestDataServer::testConstructor();
//    TestDataServer::testComparePoints();
//    TestDataServer::testRetrievePoints();
//    TestDataServer::testRetrievePoints_Threads();
//    TestDataServer::testPickRandomPoints();
//    TestDataServer::testCreateCmpDict();
//    TestDataServer::testCreateCmpDict_Threads();
//    TestDataServer::testSplitIntoEpsNet();
//    TestDataServer::testSplitIntoEpsNet_WithThreads();
//    TestDataServer::testCalculateCellMeans();
//    TestDataServer::testCalculateCellMeans_WithThreads();
//    TestDataServer::testGetMinimalDistances();
    TestDataServer::testGetMinimalDistances_WithThreads();
//    TestDataServer::testCalculateThreshold();
//    TestDataServer::testChoosePointsByDistance_WithThreads();
    cout << " ============ Test DataServer Finished ============ " << endl << endl;

}