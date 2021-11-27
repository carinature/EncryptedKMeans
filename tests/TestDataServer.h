
#ifndef ENCRYPTEDKMEANS_TESTDATASERVER_H
#define ENCRYPTEDKMEANS_TESTDATASERVER_H


#include <vector>
#include <src/Client.h>

class TestDataServer {
public:
    static void testConstructor();

    static void testComparePoints();

    static void testRetrievePoints();

    static void testRetrievePoints_Threads();

    static void testPickRandomPoints();

    static void testCreateCmpDict();

    static void testCreateCmpDict_Threads();

    static void testSplitIntoEpsNet();

    static void testSplitIntoEpsNet_WithThreads();

    static void testCalculateCellMeans();

    static void testCalculateCellMeans_WithThreads();

    static void testCalculateThreshold();

    static void testGetMinimalDistances();
    static void testGetMinimalDistances_WithThreads();

    static void testChoosePointsByDistance();
};


#endif //ENCRYPTEDKMEANS_TESTDATASERVER_H
