

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
     * @brief Picks a group of random points to be used as both cell representative and max-limit.
     * Uses the Fisher–Yates shuffle for choosing m random points. (fixme not anymore. make sure)
     * @param points - a list of all points (in current group).
     * @param numOfReps - the desired number of representatives, usually the number of desired data slices.
     * @returns a sub list of the original points, picked randomly using the .
     * @return std::vector<Point>
     * */
    const std::vector<std::vector<Point>>
    pickRandomPoints(
            const std::vector<Point> &points,
            int m = int(1 / epsilon));

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
    std::vector<
            std::unordered_map<
                    const Point,
                    std::unordered_map<
                            const Point,
                            helib::Ctxt
                    >
            >
    >
//    std::vector<
//            std::unordered_map<
//                    const Point,
//                    std::unordered_map<
//                            const Point,
//                            helib::Ctxt,
//                            hashPoints
//                    >,
//                    hashPoints
//            >
//    >
    createCmpDict(
            const std::vector<std::vector<Point>> &randomPoints,
            const std::vector<Point> &allPoints) {
        auto t0_cmpDict = CLOCK::now();//for logging

//        std::vector<std::unordered_map<const Point, std::unordered_map<const Point, helib::Ctxt, hashPoints>, hashPoints> >
        std::vector<std::unordered_map<const Point, std::unordered_map<const Point, helib::Ctxt> > >
                cmpsDict(DIM);

        for (short dim = 0; dim < DIM; ++dim) {

//            std::unordered_map<const Point, std::unordered_map<const Point, helib::Ctxt, hashPoints>, hashPoints>
            std::unordered_map<const Point, std::unordered_map<const Point, helib::Ctxt> >
                    &cmpDictForCurrDim = cmpsDict[dim];
//            cmpDictForCurrDim.reserve(allPoints.size());

            for (const Point &point : randomPoints[dim]) {

                //                cmpDictForCurrDim[point].reserve(pow(sliceSize, dim));
                for (const Point &point2 : allPoints) {
                    //                if (!cmpsDict[point].empty() && !cmpsDict[point][point2].isEmpty()))
                    //                if (point.id == point2.id) continue; //this is checked inside isBigger function
                    // todo in the future, for efficiency
                    //  - need to check if exist
                    //  - find a way to utilize the 2nd result of the `isBiggerThan()`
                    //      since you get it for free in helibs cmp
                    //  - maybe useful to use helibs `min/max` somehow?
                    const std::vector<helib::Ctxt> res = point.isBiggerThan(point2, dim);
                    cmpDictForCurrDim[point].emplace(point2, res[0]);   //  point > point2
                    cmpDictForCurrDim[point2].emplace(point, res[1]);   //  point < point2
                }
            }

        }

        loggerDataServer.log(printDuration(t0_cmpDict, "CmpDict Creation"));
        return cmpsDict;

    }
/*


    // TODO candidate for multithreading
    */
/**
     * @brief Split into (1/eps) groups - each group is between 2 representative points.
     * @param points - a list of unordered points
     * @param dim - the index of coor by which the comparison is made. in other words - in what dimension the split is made.
     * @param keysServer - the appointed CA
     * @returns a list of pairs/tuples of a representative point and a list of the points in its Group (cell/slice).
     * *//*

    std::vector<
            std::tuple< //should be replaced with std::pair in production
                    Point,
                    std::vector<Point>,
                    std::vector<Ctxt>
            >
    >
    split(
            std::vector<Point> &points,
            short dim) {

        auto t0_split = CLOCK::now();
        std::vector<std::tuple<Point, std::vector<Point>, std::vector<Ctxt> > > groups;
        if (points.empty() || 0 > dim || DIM < dim) return groups;        // sanity check
        //        if (points.empty() || 0 > dim || DIM < dim) return std::vector<std::tuple<Point, std::vector<Point>, std::vector<Ctxt> > >();        // sanity check

        const std::vector<std::vector<Point>> randomPoints = pickRandomPoints(points, 1 / epsilon);

        //  create points-comparing dict - for every 2 points (p1,p2) answers p1[dim]>p2[dim]
        std::vector<
                std::unordered_map<
                        const Point,
                        std::unordered_map<const Point, helib::Ctxt,
                                hashPoints
                        >,
                        hashPoints
                >
        >
                cmp = createCmpDict(randomPoints, points);
        //        cout << "Random Points in this group: ";
        printPoints(points, keysServer);
        */
/*

                for (const Point &R: randomPoints) {
                    auto t0_itr = CLOCK::now();
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
                        helib::Ctxt isBelowCurrentRep(cmp[R].at(p));
                        //                printNameVal(keysServer.decryptCtxt(isBelowCurrentRep));
                        helib::Ctxt isAboveAboveSmallerReps(isBelowCurrentRep); //todo other init
                        //                printNameVal(keysServer.decryptCtxt(isAboveAboveSmallerReps));

                        for (const Point &r: randomPoints) {
                            //                    cout << "           --- other r: ";
                            //                    printPoint(r, keysServer);
                            if (r.id == R.id) continue;
                            //todo here's a good place to check for "tail" p - is curr p bigger than all the random points
                            //todo handle cases of p that is equal to representative

                            //  r > R (in which case we don't care about cmp results of p and r)
                            Bit otherRepIsAboveCurrentRep = cmp[r].at(R);
                            // results in: Bit otherRepIsAboveCurrentRep = (r > R)
                            //                    printNameVal(keysServer.decryptCtxt(otherRepIsAboveCurrentRep));

                            //  R > r    AND     p > r
                            Bit pIsAboveOtherRep(cmp[R].at(r)); //pIsAboveOtherRep = (R > r)
                            pIsAboveOtherRep.multiplyBy(cmp[p].at(r));
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

                        long pIsInGroup = keysServer.decryptCtxt(isInGroup);
                        if (pIsInGroup) {
                            cout << "       Point is in group of       ";
                            printPoint(pointIsInCell, keysServer);
                        }
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
        *//*


        // todo "tail" points - all the points that are bigger than all random reps
        loggerDataServer.log(printDuration(t0_split, "split total"));

        return groups;
    }
*/

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

