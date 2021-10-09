

#ifndef ENCKMEAN_DATASERVER_H
#define ENCKMEAN_DATASERVER_H


#include "Client.h"

class DataServer {

    //    [[maybe_unused]] const helib::SecKey encryptionKey;
    //    [[maybe_unused]] const helib::PubKey encryptionKey;
    //    [[maybe_unused]] const helib::EncryptedArray ea;
    const helib::PubKey encryptionKey;
    //     const helib::EncryptedArray ea;
    const helib::Ctxt scratch;
protected:
    //    const helib::PubKey &public_key;// = encryptionKey;

public:
    /**
     * Constructor for \class{Client},
     * @param keysServer binds to the \class{KeysServer} responsible for the distributing the appropriate key
     * @brief simulates a mini-protocol between the keys server and the Data Server.
     * ks and ds will create a scratch that will allow to encrypt the results:
     * result bit in case of compare, and result num in case of add/multiplication.
     * */
    explicit DataServer(KeysServer &keysServer);

    //void BGV_binary_arithmetic();
    //void BGV_packed_arithmetic();
    //void KT_packed_arithmetic();
    //void add_numbers();
    //void mult_numbers();
    //void cmp_numbers();

    //    helib::Ctxt compare(Client &client1, Client &client2) {
    //#if VERBOSE
    //        cout << "compare(Client &client) " << endl;
    //#endif
    //        Ctxt &encryptionKey;
    //        helib::Ctxt isBiggerFlag(encryptionKey), isSmallerFlag(encryptionKey);
    //        cout << "odin" << endl;
    //        compareTwoNumbers(isBiggerFlag, isSmallerFlag,
    //                          helib::CtPtrs_VecCt(this->cCoordinatesNTL),
    //                          helib::CtPtrs_VecCt(client.cCoordinatesNTL),
    //                          false, nullptr);
    //        cout << "onaes" << endl;
    //        return isBiggerFlag;
    //    }
    Point scratchPoint();

    // todo move to aux ?
    static std::vector<Client> generateDataClients(const KeysServer &server) {
        int uniquePointsNum = 3 + rand() % 10, clientsNum = uniquePointsNum;
        //  init coordinate arrays
        long arrs[uniquePointsNum][DIM];
        for (auto &arr: arrs) for (int dim = 0; dim < DIM; ++dim) arr[dim] = rand() % NUMBERS_RANGE;
        /*  init clients vector
         *      client [0] stays empty
         *      client [1] has 1 point - {point0}
         *      client [2] has 2 points - {point0, point1}
         *      ...
         *      client [n] has n points -  {point0, point1, ... , pointN}
         */
        std::vector<Client> clients(clientsNum, Client(server));
        for (int i = 1; i < clients.size(); ++i)
            for (int j = 0; j < i; ++j)
                clients[i].encryptPoint(arrs[j]);
        return clients;
    }

    static std::vector<Point> retrievePoints(std::vector<Client> &clients) {
        //  retrieving points
        std::vector<Point> points;
        points.reserve(pow(clients.size(), 2) / 2); // preallocate memory
        for (Client &c:clients)
            for (Point &p:c.points)
                points.emplace_back(Point(p));
        //            points.insert(points.end(),c.points.begin(), c.points.end());
        points.shrink_to_fit();
        return points;
    }
};


#endif //ENCKMEAN_DATASERVER_H
