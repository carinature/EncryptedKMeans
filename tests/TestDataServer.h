
#ifndef ENCRYPTEDKMEANS_TESTDATASERVER_H
#define ENCRYPTEDKMEANS_TESTDATASERVER_H


#include <vector>
#include <src/Client.h>

class TestDataServer {
public:
    static void testConstructor();

    static void testComparePoints();

    static void testRetrievePoints();

    static void testPickRandomPoints();

    static void testCreateCmpDict();

    static void testSplitIntoEpsNet();

    static void testCalculateCellMeans();

    static void testCalculateThreshold();

    static void testGetMinimalDistances();

    static void testChoosePointsByDistance();
};


#endif //ENCRYPTEDKMEANS_TESTDATASERVER_H
