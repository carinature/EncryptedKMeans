
#include "utils/aux.h"

#include "TestAux.h"
#include "src/Client.h"

void TestAux::testGenerateDataClients() {
    cout << " ------ testGenerateDataClients ------ " << endl << endl;
    KeysServer keysServer;
    auto t0_main = CLOCK::now();

    const std::vector<Client> clients = generateDataClients(keysServer);

    for (const Client &client:clients)
        for (const Point &point:client.getPoints())
            printPoint(point, keysServer);

    cout << endl<<printDuration(t0_main, "testGenerateDataClients");
    cout << " ------ testGenerateDataClients finished ------ " << endl << endl;

}
