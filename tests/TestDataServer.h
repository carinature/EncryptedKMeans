
#ifndef ENCRYPTEDKMEANS_TESTDATASERVER_H
#define ENCRYPTEDKMEANS_TESTDATASERVER_H


#include <vector>
#include <src/Client.h>

class TestDataServer {
public:
    static void testConstructor();

    static void testComparePoints();

    static void testAddition();

    static void testMultiplication();

    static void testgGenerateDataClients();
    static void testRetrievePoints();
    static void testPickRandomPoints();

    static void testminirand();
};


#endif //ENCRYPTEDKMEANS_TESTDATASERVER_H
