
#ifndef ENCRYPTEDKMEANS_TESTDATASERVER_H
#define ENCRYPTEDKMEANS_TESTDATASERVER_H


#include <vector>
#include <src/Client.h>

class TestDataServer {
public:
    static void testConstructor();
    static void testScratchPoint();
    static void testComparePoints();
//    static void testAddiiton();
//    static void testMultiplication();
static void testRetrievePoints();

    static std::vector<Client> generateDataClients(const KeysServer &server);
};


#endif //ENCRYPTEDKMEANS_TESTDATASERVER_H
