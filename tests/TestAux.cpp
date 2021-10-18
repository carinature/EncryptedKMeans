
#include "utils/aux.h"

#include "TestAux.h"
#include "src/Client.h"

void TestAux::testGenerateDataClients() {
    cout << " ------ testGenerateDataClients ------ " << endl << endl;
    KeysServer keysServer;
    auto t0_main = CLOCK::now();

    std::vector<Client> clients = generateDataClients(keysServer);

    for (Client &client:clients)
        for (Point &point:client.points)
            printPoint(point, keysServer);

    printDuration(t0_main, "testGenerateDataClients");
    cout << " ------ testGenerateDataClients finished ------ " << endl << endl;

}
