//
// Created by karina on 25/07/2021.
//

#include "utils/Logger.h"
#include "TestClient.h"
#include "src/Client.h"

Logger logger;

void TestClient::testConstructor() {
    logger.log("testConstructor");
    KeysServer server;// = KeysServer();
//    const long arr[] = {1L, 2L};
    Client client(server);//, arr);
//    logger.print_log();
}

void TestClient::testDecryptCoordinates() {
    logger.log("testDecryptCoordinates");
    KeysServer server;// = KeysServer();
    const long arr[] = {12, 34};
    Client client(server);//, arr);
    client.encryptPoint(arr);
    client.decryptCoordinates();
//    logger.print_log();
}


void TestClient::testCompare() {
    logger.log("testCompare");
    KeysServer server = KeysServer();
    const long arr1[] = {1L, 2L};
    const long arr2[] = {2L, 2L};
    Client client1(server);//, arr1);
    Client client2(server);//, arr2);
    client1.encryptPoint(arr1);
    client2.encryptPoint(arr2);
    client1.compare(client2);
//    logger.print_log();
}

void TestClient::testAddiiton() {
    logger.log("testAddiiton");
}

void TestClient::testMultiplication() {
    logger.log("testMultiplication");
}

