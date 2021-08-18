//
// Created by karina on 25/07/2021.
//

#include "TestClient.h"
#include "TestDataServer.h"

int main(){
    TestClient::testConstructor();
    TestClient::testEncryptCoordinates();
    TestClient::testDecryptCoordinates();
    TestClient::testEncryptScratchPoint();

    TestDataServer::testConstructor();
    TestDataServer::testScratchPoint();
    TestDataServer::testCompareClients();


}