#ifndef ENCRYPTEDKMEANS_TESTPOINT_H
#define ENCRYPTEDKMEANS_TESTPOINT_H

class TestPoint {
public:
    static void testConstructor();

    static void testEncryptCoordinates();

    static void testOperatorSubscript();

    static void testIsEmpty();

    static void testAddition();

    static void testAddManyPoints();

    static void testMultiplication();

    static void testMultiplicationByBit();

    static void testCompare();

    static void testCalculateDistanceFromPoint();

    static void testFindMinimalDistancesFromMeans();
};

#endif //ENCRYPTEDKMEANS_TESTPOINT_H
