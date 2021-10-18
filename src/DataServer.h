

#ifndef ENCKMEAN_DATASERVER_H
#define ENCKMEAN_DATASERVER_H


#include "Client.h"

//move to the cpp file
#include <algorithm> //for the random shuffle
#include <random>

static Logger loggerDataServer(log_debug, "loggerDataServer");

class DataServer {

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
    explicit DataServer(
            KeysServer &keysServer);

    /**
     * @brief A simulated retrievel of data from clients.
     * @param clients - a list of clients (chosen by the CA, to share a similar public key).
     * @returns a list of all the points.
    * @return std::vector<Point>
     * * */
    // TODO candidate for multithreading
    static std::vector<Point>
    retrievePoints(
            const std::vector<Client> &clients);

    /**
     * @brief request data from clients and conentrate into one list
     * @param points - all the points from all the clients in the data set
     * @param m - number of random representatives for each slice
     * @returns a list of #DIM lists - each containing m^d randomly chosen points
     * @return std::vector<Point>
     * */
    const std::vector<std::vector<Point>>
    pickRandomPoints(const std::vector<Point> &points,
                     int m,
                     const Point &tinyRandomPoint);

    // TODO candidate for multithreading
    /**
     * @brief create a comparison dict:
     *   for each 2 points a,b returns the answer to a[dim]>b[dim]
     *   | where a and b are from 2 different groups (one is supgroup of the other),
     *    and dim is the index of the desired coordinate.
     * @param points - a list of all points (in current group).
     * @param randomPoints - a sub group of all points (in current group).
     *  it is a vector of size #DIM, each node is a vector of size m^dim containing random reps
     * @param numOfReps - the desired number of representatives, usually the number of desired data slices.
     * @returns a vector of #DIM dictionaries.
     *   for each #dim the dictionary contains a pair of keys [point1,point2] and the encrypted value [point1[dim]>point[dim].
     * @return std::vector<Point>
     * */
    static
    std::vector<std::unordered_map<const Point, std::unordered_map<const Point, helib::Ctxt>>>
    createCmpDict(const std::vector<Point> &allPoints,
                  const std::vector<std::vector<Point>> &randomPoints,
                  const Point &tinyRandomPoint);


    // TODO candidate for multithreading
    /**
     * @brief Split into (1/eps) groups - each group is between 2 representative points.
     * @param points - a list of unordered points
     * @param randomPoints - a list of unordered points
     * @param dim - the index of coor by which the comparison is made. in other words - in what dimension the split is made.
     * @returns a list of pairs/tuples of a representative point and a list of the points in its Group (cell/slice).
     * */
    static
    std::map<int, //DIM
            std::vector< //current slices for approp dimension
                    Cell
            >
    >
    splitIntoEpsNet(
            const std::vector<Point> &points,
            const std::vector<std::vector<Point>> &randomPoints,
            const std::vector<
                    std::unordered_map<
                            const Point,
                            std::unordered_map<
                                    const Point,
                                    helib::Ctxt> > > &cmpDict,
            const KeysServer &keysServer // for dbg todo remove
    ) {
        auto t0_split = CLOCK::now();     //  for logging, profiling, DBG

        std::map<int, std::vector<Cell> > groups;
        if (points.empty()) return groups;        // sanity check

        // init initial 1st level of data
        Cell initCell;
        const helib::PubKey &publicKey = points[0].public_key;
        helib::Ctxt tempCtxt(publicKey);
        //        initCell.addRep(points[0]);
        for (auto const &point:points) initCell.addPoint(point, tempCtxt);
        groups[-1].push_back(initCell);

        for (int dim = 0; dim < DIM; ++dim) {
            auto t0_itr_dim = CLOCK::now();     //  for logging, profiling, DBG

            for (const Cell &cellFromPrevDim:groups[dim - 1]) {

                for (const Point &R: randomPoints[dim]) {
                    auto t0_itr_rep = CLOCK::now();     //  for logging, profiling, DBG

                    Cell newCell;
                    const std::vector<Point> &prevDimReps = cellFromPrevDim.reps;
                    for (auto const &rep:prevDimReps) newCell.addRep(rep);

                    newCell.addRep(R);
                    const std::unordered_map<const Point, std::unordered_map<const Point, helib::Ctxt>> &cmpDim = cmpDict[dim];
                    const auto &atR = cmpDim.at(R);
                    const Ctxt &atRR = atR.at(R);
                    helib::Ctxt isRepInPrevCell(atRR);
                    if (0 < dim) isRepInPrevCell = cmpDict[dim - 1].at(prevDimReps[dim - 1]).at(R);

                    for (const Point &p:cellFromPrevDim.includedPoints) {

                        // p < R
                        helib::Ctxt isBelowCurrentRep(cmpDict[dim].at(R).at(p));
                        //                printNameVal(keysServer.decryptCtxt(isBelowCurrentRep));
                        helib::Ctxt isAboveAboveSmallerReps(isBelowCurrentRep); //todo other init
                        //                printNameVal(keysServer.decryptCtxt(isAboveAboveSmallerReps));

                        for (const Point &r: randomPoints[dim]) {
                            //                    cout << "           --- other r: ";
                            //                    printPoint(r, keysServer);
                            if (r.id == R.id) continue;
                            //todo here's a good place to check for "tail" p - is curr p bigger than all the random points
                            //todo handle cases of p that is equal to representative

                            //  r > R (in which case we don't care about cmpDict results of p and r)
                            Bit otherRepIsAboveCurrentRep = cmpDict[dim].at(r).at(R);
                            // results in: Bit otherRepIsAboveCurrentRep = (r > R)
                            //                    printNameVal(keysServer.decryptCtxt(otherRepIsAboveCurrentRep));

                            //  R > r    AND     p > r
                            Bit pIsAboveOtherRep(
                                    cmpDict[dim].at(R).at(r)); //pIsAboveOtherRep = (R > r)
                            pIsAboveOtherRep.multiplyBy(cmpDict[dim].at(p).at(r));
                            // results in: Bit pIsAboveOtherRep = (R > r) * (p > r)
                            //                    printNameVal(keysServer.decryptCtxt(pIsAboveOtherRep));

                            //   R < r  OR  [ R > r    AND     p > r ]
                            Bit pIsBelowCurrentRepAndAboveOtherRep = pIsAboveOtherRep;
                            pIsBelowCurrentRepAndAboveOtherRep.multiplyBy(
                                    otherRepIsAboveCurrentRep);
                            pIsBelowCurrentRepAndAboveOtherRep.negate();
                            pIsBelowCurrentRepAndAboveOtherRep += pIsAboveOtherRep;
                            pIsBelowCurrentRepAndAboveOtherRep += otherRepIsAboveCurrentRep;
                            // results in: Bit pIsBelowCurrentRepAndAboveOtherRep =
                            // pIsAboveOtherRep + otherRepIsAboveCurrentRep
                            // - pIsAboveOtherRep * otherRepIsAboveCurrentRep
                            //                    printNameVal(keysServer.decryptCtxt(pIsBelowCurrentRepAndAboveOtherRep));

                            // this will hold the Product[ (R > r) && (p > r) | foreach r in randomPoints ]
                            isAboveAboveSmallerReps *= pIsBelowCurrentRepAndAboveOtherRep;
                            //                    printNameVal(keysServer.decryptCtxt(isAboveAboveSmallerReps));
                        }
                        Bit isInGroup = isBelowCurrentRep;
                        isInGroup *= isAboveAboveSmallerReps;

                        //fixme todo eod change! may need remove
                        isInGroup *= isRepInPrevCell;
                        //end fixme todo eod change! may need remove

                        Point pointIsInCell = p * isInGroup;
                        /*

                                                //                    long pIsInGroup = keysServer.decryptCtxt(isInGroup);
                                                //                    if (pIsInGroup) {
                                                //                        cout << "       Point is in group of       ";
                                                //                        printPoint(pointIsInCell, keysServer);
                                                //                    }
                                                //                cout << "       Point is in group?       "
                                                //                     << pIsInGroup ? "yes    " : "no   ";
                                                //                printPoint(pointIsInCell, keysServer);
                                                //                cout << "\n";
                        */

                        newCell.addPoint(pointIsInCell, isInGroup);
                    }

                    groups[dim].emplace_back(newCell);
                    //                groups.emplace_back(R, group, groupSize);

                    loggerDataServer.log(printDuration(
                            t0_itr_rep,
                            "split iteration (for single random point)"));
                }
            }
            loggerDataServer.log(printDuration(
                    t0_itr_dim,
                    "split iteration (for #" + std::to_string(dim) + " dimention)"));
        }
        // todo "tail" points - all the points that are bigger than all random reps
        loggerDataServer.log(printDuration(t0_split, "split total"));

        return groups;
    }




    // TODO candidate for multithreading
    /**
     * @brief Split into (1/eps) groups - each group is between 2 representative points.
     * @param points - a list of unordered points
     * @param randomPoints - a list of unordered points
     * @param dim - the index of coor by which the comparison is made. in other words - in what dimension the split is made.
     * @returns a list of pairs/tuples of a representative point and a list of the points in its Group (cell/slice).
     * */
    static
    std::vector<
            std::tuple< //should be replaced with std::pair in production
                    Point,
                    std::vector<Point>,
                    std::vector<Ctxt>
            >
    >
    split(
            const std::vector<Point> &points,
            const std::vector<std::vector<Point>> &randomPoints,
            const std::vector<
                    std::unordered_map<
                            const Point,
                            std::unordered_map<
                                    const Point,
                                    helib::Ctxt> > > &cmpDict) {

        auto t0_split = CLOCK::now();     //  for logging, profiling, DBG

        //        std::vector<std::tuple<Point, std::vector<Point>, std::vector<Ctxt> > > groups;
        std::vector<
                std::tuple<
                        Point,
                        std::vector<Point>,
                        std::vector<Ctxt> > > groups;
        if (points.empty()) return std::vector<std::tuple<Point, std::vector<Point>, std::vector<Ctxt> > >();        // sanity check

        //        //        cout << "Random Points in this group: ";
        //        printPoints(points, keysServer);
        //        for (int dim = 0; dim < DIM; ++dim) {
        for (int dim = 0; dim < 1; ++dim) { //fixme

            for (const Point &R: randomPoints[dim]) {
                auto t0_itr = CLOCK::now();     //  for logging, profiling, DBG
                //            cout << "\n   ############### current R: ";
                //            printPoint(R, keysServer);

                std::vector<Point> group;
                group.reserve(points.size() * 2 * epsilon);
                std::vector<Ctxt> groupSize;
                groupSize.reserve(points.size() * 2 * epsilon);

                for (const Point &p:points) {
                    //                cout << "       ============= current p: ";
                    //                printPoint(p, keysServer);

                    // p < R
                    helib::Ctxt isBelowCurrentRep(cmpDict[dim].at(R).at(p));
                    //                printNameVal(keysServer.decryptCtxt(isBelowCurrentRep));
                    helib::Ctxt isAboveAboveSmallerReps(isBelowCurrentRep); //todo other init
                    //                printNameVal(keysServer.decryptCtxt(isAboveAboveSmallerReps));

                    for (const Point &r: randomPoints[dim]) {
                        //                    cout << "           --- other r: ";
                        //                    printPoint(r, keysServer);
                        if (r.id == R.id) continue;
                        //todo here's a good place to check for "tail" p - is curr p bigger than all the random points
                        //todo handle cases of p that is equal to representative

                        //  r > R (in which case we don't care about cmpDict results of p and r)
                        Bit otherRepIsAboveCurrentRep = cmpDict[dim].at(r).at(R);
                        // results in: Bit otherRepIsAboveCurrentRep = (r > R)
                        //                    printNameVal(keysServer.decryptCtxt(otherRepIsAboveCurrentRep));

                        //  R > r    AND     p > r
                        Bit pIsAboveOtherRep(cmpDict[dim].at(R).at(r)); //pIsAboveOtherRep = (R > r)
                        pIsAboveOtherRep.multiplyBy(cmpDict[dim].at(p).at(r));
                        // results in: Bit pIsAboveOtherRep = (R > r) * (p > r)
                        //                    printNameVal(keysServer.decryptCtxt(pIsAboveOtherRep));

                        //   R < r  OR  [ R > r    AND     p > r ]
                        Bit pIsBelowCurrentRepAndAboveOtherRep = pIsAboveOtherRep;
                        pIsBelowCurrentRepAndAboveOtherRep.multiplyBy(otherRepIsAboveCurrentRep);
                        pIsBelowCurrentRepAndAboveOtherRep.negate();
                        pIsBelowCurrentRepAndAboveOtherRep += pIsAboveOtherRep;
                        pIsBelowCurrentRepAndAboveOtherRep += otherRepIsAboveCurrentRep;
                        // results in: Bit pIsBelowCurrentRepAndAboveOtherRep =
                        // pIsAboveOtherRep + otherRepIsAboveCurrentRep
                        // - pIsAboveOtherRep * otherRepIsAboveCurrentRep
                        //                    printNameVal(keysServer.decryptCtxt(pIsBelowCurrentRepAndAboveOtherRep));

                        // this will hold the Product[ (R > r) && (p > r) | foreach r in randomPoints ]
                        isAboveAboveSmallerReps *= pIsBelowCurrentRepAndAboveOtherRep;
                        //                    printNameVal(keysServer.decryptCtxt(isAboveAboveSmallerReps));
                    }
                    Bit isInGroup = isBelowCurrentRep;
                    isInGroup *= isAboveAboveSmallerReps;

                    Point pointIsInCell = p * isInGroup;

                    //                    long pIsInGroup = keysServer.decryptCtxt(isInGroup);
                    //                    if (pIsInGroup) {
                    //                        cout << "       Point is in group of       ";
                    //                        printPoint(pointIsInCell, keysServer);
                    //                    }
                    //                cout << "       Point is in group?       "
                    //                     << pIsInGroup ? "yes    " : "no   ";
                    //                printPoint(pointIsInCell, keysServer);
                    //                cout << "\n";

                    group.push_back(pointIsInCell);
                    groupSize.emplace_back(isInGroup);
                }

                groups.emplace_back(R, group, groupSize);

                loggerDataServer.log(
                        printDuration(t0_itr, "split iteration (for single random point)"));
            }

        }
        // todo "tail" points - all the points that are bigger than all random reps
        loggerDataServer.log(printDuration(t0_split, "split total"));

        return groups;
    }

};


#endif //ENCKMEAN_DATASERVER_H

/*
    */
/*  @brief Generate random data.
     *  Returns a vector of clients with random points.
    *      client [0] stays empty
    *      client [1] has 1 point - {point0}
    *      client [2] has 2 points - {point0, point1}
    *      ...
    *      client [n] has n points -  {point0, point1, ... , pointN}
    *//*
  // todo move to aux ?
    static std::vector<Client> generateDataClients(const KeysServer &keysServer) {
        */
/*        //        Note that rand() is considered harmful, and is discouraged in C++14
                int uniquePointsNum = 3 + random() % number_of_points,
                        clientsNum = uniquePointsNum;
                printNameVal(number_of_points);
                printNameVal(uniquePointsNum);
                *//*
*/
/*        //  init coordinate arrays
                //        long arrs[uniquePointsNum][DIM];
                //        for (auto &arr: arrs)
                //            for (short dim = 0; dim < DIM; ++dim)
                //                arr[dim] = rand() % NUMBERS_RANGE;*//*
*/
/*
        long arr[DIM];
        std::vector<Client> clients(clientsNum, Client(keysServer));
        for (int i = 1; i < clients.size(); ++i)
            for (int j = 0; j < i; ++j) {
                //  pick random coordinates
                for (short dim = 0; dim < DIM; ++dim)
                    arr[dim] = random() % NUMBERS_RANGE;
                //  encrypt coordinates and add to client
                clients[i].encryptPoint(arr);
                //                clients[i].encryptPoint(arrs[j]);}
            }*//*


        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0, NUMBERS_RANGE);
        long tempArr[DIM];
        std::vector<Client> clients(number_of_clients, Client(keysServer));
        for (Client &client:clients) {
            for (int dim = 0; dim < DIM; ++dim) tempArr[dim] = dist(mt);
            client.encryptPoint(tempArr);
        }
        return clients;
    }
*/

