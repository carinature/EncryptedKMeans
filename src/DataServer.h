

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
    KeysServer &keysServer;
public:
    /**
     * Constructor for \class{Client},
     * @param keysServer binds to the \class{KeysServer} responsible for the distributing the appropriate key
     * @brief simulates a mini-protocol between the keys server and the Data Server.
     * ks and ds will create a scratch that will allow to encrypt the results:
     * result bit in case of compare, and result num in case of add/multiplication.
     * */
    explicit DataServer(KeysServer &keysServer);

    Point scratchPoint();

    /*  Generate random data.
     *  Returns a vector of clients with random points.
    *      client [0] stays empty
    *      client [1] has 1 point - {point0}
    *      client [2] has 2 points - {point0, point1}
    *      ...
    *      client [n] has n points -  {point0, point1, ... , pointN}
    */  // todo move to aux ?
    static std::vector<Client> generateDataClients(const KeysServer &server) {
        int uniquePointsNum = 3 + rand() % number_of_points,
                clientsNum = uniquePointsNum;
        /*        //  init coordinate arrays
                //        long arrs[uniquePointsNum][DIM];
                //        for (auto &arr: arrs)
                //            for (short dim = 0; dim < DIM; ++dim)
                //                arr[dim] = rand() % NUMBERS_RANGE;*/
        long arr[DIM];
        std::vector<Client> clients(clientsNum, Client(server));
        for (int i = 1; i < clients.size(); ++i)
            for (int j = 0; j < i; ++j) {
                //  pick random coordinates
                for (short dim = 0; dim < DIM; ++dim)
                    arr[dim] = rand() % NUMBERS_RANGE;
                //  encrypt coordinates and add to client
                clients[i].encryptPoint(arr);
                //                clients[i].encryptPoint(arrs[j]);}
            }
        return clients;
    }

    /* A simulated retrievel of data from clients.
     * Params:
     *  clients - a list of clients (chosen by the CA, to share a similar public key).
     * Returns a list of all the points.
     * */
    static std::vector<Point> retrievePoints(std::vector<Client> &clients) {
        std::vector<Point> points;
        if (clients.empty()) return points;
        points.reserve(pow(clients.size(), 2) / 2); // preallocate memory
        for (Client &c:clients)
            for (Point &p:c.points)
                points.emplace_back(Point(p));
        //            points.insert(points.end(),c.points.begin(), c.points.end());
        points.shrink_to_fit();
        return points;
    }

    /* Picks a group of random points to be used as both cell representative and max-limit.
     * Uses the Fisher–Yates shuffle for choosing k random points.
     * Params:
     *  points - a list of all points (in current group).
     *  numOfReps - the desired number of representatives, usually the number of desired data strips.
     * Returns a sub list of the original points, picked randomly using the .
     * */
    static std::vector<Point>
    pickRandomPoints(std::vector<Point> &points, int numOfStrips = int(1 / epsilon)) {
        // sanity check
        if (points.empty() || numOfStrips > points.size()) return std::vector<Point>();

        std::vector<Point> copy(points);
        auto begin = copy.begin();
        int k = numOfStrips, stripSize = numOfStrips;
        while (k--) {
            auto r = begin;
            // next line used to crushes the program with small(<20) number of points in file
            advance(r, random() % stripSize);
            swap(begin, r);
            ++begin;
            --stripSize;
        }
        //  todo see if there's a way to reduce this copy
        //        std::vector<Point> random(copy.begin(), copy.begin() + numOfStrips - 1);
        return std::vector<Point>(copy.begin(), copy.begin() + numOfStrips - 1);
    }

    static std::vector<
            std::tuple<
                    Point,
                    std::vector<Point>,
                    std::vector<Ctxt>
            >
    >
    split(std::vector<Point> &points, short dim, KeysServer &keysServer) {
        cout << "Split" << endl;
        auto t0_split = std::chrono::high_resolution_clock::now();
        std::vector<std::tuple<Point, std::vector<Point>, std::vector<Ctxt> > > groups;
        // sanity check
        if (points.empty() || 0 > dim || DIM < dim) return groups;

        std::vector<Point> randomPoints = DataServer::pickRandomPoints(points);

//        std::map<Point,
//                std::map<Point,
//                        std::vector<Bit>,
//                        cmpPoints
//                >,
//                cmpPoints
//        > cmp = createCmpDict(random, currStrip);

        std::vector<Ctxt> numPoints;//(BIT_SIZE);//, (points[0].public_key));
        numPoints.reserve(points.size());
        //        cmpDict

        for (Point &randomPoint: randomPoints) {
            auto t0_itr = std::chrono::high_resolution_clock::now();

            std::vector<Point> group;
            for (Point &point:points) {
                Ctxt ctxt(randomPoint.isBiggerThan(point, dim)[0]);
                group.push_back(point * ctxt);
                numPoints.emplace_back(ctxt);
            }

            // fixme need to implement hash function for Point to use map ?
            //      `in instantiation of member function 'std::less<Point>::operator()' requested here`
            groups.emplace_back(randomPoint, group, numPoints);
            printDuration(t0_itr, "split iteration");
        }

        // todo "tail" points - all the points that are bigger than all random reps


        //
        //        for (std::tuple tuple: groups) {
        //            cout << "curtuple" << endl;
        //            Point randomPoint = std::get<0>(tuple);
        //            std::vector<Point> group = std::get<1>(tuple);
        //            std::vector<Ctxt> groupSize = std::get<2>(tuple);
        //            cout << "for random point: ";
        //            printPoint(randomPoint, keysServer);
        //            cout << " these points will be included: \n";
        //            printPoints(group, keysServer);
        //            cout << "group size: " << keysServer.decryptSize(groupSize) << endl;
        //            //        for (Point &point:group) {
        //            //        }
        //        }

        printDuration(t0_split, "split ");
        return groups;
    }


//    std::map<Point,
//            std::map<Point,
//                    std::vector<Bit>,
//                    cmpPoints
//            >,
//            cmpPoints
//    >
//    createCmpDict(const std::vector<Point> &randomPoints,
//                  const std::vector<Point> &stripPoints,
//                  short dim) {
//        std::map<Point, std::map<Point, std::vector<Bit>, cmpPoints>, cmpPoints> cmpDict;
//        for (const Point &point : randomPoints) {
////            std::map<Point, std::vector<Bit>, cmpPoints> cmpDictMini; //, cmpDictMini2;
//            for (const Point &point2 : stripPoints) {
//                std::vector<helib::Ctxt> res = point.isBiggerThan(point2, dim);
//                // used to be cmpDictMini[point2].push_back(point > point2); and worked
//                cmpDict[point][point2].push_back(res[0]);   // point > point2
//                cmpDict[point2][point].push_back(res[1]);   // point2 > point
//            }
//        }
//        return cmpDict;
//    }

};


#endif //ENCKMEAN_DATASERVER_H
