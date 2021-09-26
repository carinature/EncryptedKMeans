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
};

#endif //ENCRYPTEDKMEANS_TESTPOINT_H
