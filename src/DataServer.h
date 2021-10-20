

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
    const Point tinyRandomPoint;
public:
    /**
     * Constructor for \class{Client},
     * @param keysServer binds to the \class{KeysServer} responsible for the distributing the appropriate key
     * @brief simulates a mini-protocol between the keys server and the Data Server.
     * ks and ds will create a scratch that will allow to encrypt the results:
     * result bit in case of compare, and result num in case of add/multiplication.
     * */
    explicit DataServer(KeysServer &keysServer) :
            keysServer(keysServer),
            tinyRandomPoint(keysServer.tinyRandomPoint()),
            encryptionKey(keysServer.getPublicKey()),
            //        public_key(keysServer.getPublicKey()),
            //        public_key(encryptionKey),
            //        ea(keysServer.getEA()),
            scratch(keysServer.getPublicKey()) {
        //        dataServerLogger.log("DataServer()");
        cout << "DataServer()" << endl;
    }

    // TODO candidate for multithreading
    /**
     * @brief A simulated retrievel of data from clients.
     * @param clients - a list of clients (chosen by the CA, to share a similar public key).
     * @returns a list of all the points.
    * @return std::vector<Point>
     * * */
    static
    std::vector<Point>
    retrievePoints(
            const std::vector<Client> &clients);

    /**
     * @brief request data from clients and conentrate into one list
     * @param points - all the points from all the clients in the data set
     * @param m - number of random representatives for each slice
     * @returns a list of #DIM lists - each containing m^d randomly chosen points
     * @return std::vector<Point>
     * */
    std::vector<std::vector<Point>>
    pickRandomPoints(
            const std::vector<Point> &points,
            int m = 1/EPSILON   //  -1
            //            const Point & tinyRandomPoint
    );

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
    //    static
    std::vector<std::unordered_map<const Point, std::unordered_map<const Point, helib::Ctxt>>>
    createCmpDict(
            const std::vector<Point> &allPoints,
            const std::vector<std::vector<Point>> &randomPoints
            //            const Point & tinyRandomPoint
    );


    // TODO candidate for multithreading
    /**
     * @brief Split into (1/eps) groups - each group is between 2 representative points.
     * @param points - a list of unordered points
     * @param randomPoints - a list of unordered points
     * @param dim - the index of coor by which the comparison is made. in other words - in what dimension the split is made.
     * @returns a list of pairs/tuples of a representative point and a list of the points in its Group (slice/slice).
     * */
    //    static
    std::map<int, //DIM
            std::vector< //current slices for approp dimension
                    Slice
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

        std::map<int, std::vector<Slice> > groups;
        if (points.empty()) return groups;        // sanity check

        // init initial 1st level of data
        Slice startingSlice;
        const helib::PubKey &publicKey = points[0].public_key;
        helib::Ctxt tempCtxt(publicKey);
        //        startingSlice.addRep(points[0]);
        for (auto const &point:points) startingSlice.addPoint(point, tempCtxt);
        groups[-1].push_back(startingSlice);

        for (int dim = 0; dim < DIM; ++dim) {
            auto t0_itr_dim = CLOCK::now();     //  for logging, profiling, DBG

            for (const Slice &sliceFromPrevDim:groups[dim - 1]) {

                for (const Point &R: randomPoints[dim]) {
                    auto t0_itr_rep = CLOCK::now();     //  for logging, profiling, DBG

                    Slice newSlice;
                    //  const std::vector<Point> &prevDimReps = sliceFromPrevDim.reps;
                    for (auto const &rep:sliceFromPrevDim.reps) newSlice.addRep(rep);

                    const std::unordered_map<const Point, std::unordered_map<const Point, helib::Ctxt>> &cmpDim = cmpDict[dim];
                    helib::Ctxt isRepInPrevSlice(cmpDim.at(R).at(R));
                    if (0 < dim && !sliceFromPrevDim.reps.empty())
                        // todo why cmp at prev dim and not current?
                        isRepInPrevSlice *= cmpDict[dim - 1].at(sliceFromPrevDim.reps[dim - 1]).at(
                                R);
                    newSlice.addRep(R);

                    //  tailSlice.addRep(R);

                    for (const Point &p:sliceFromPrevDim.includedPoints) {
                        /*                        if (p.id == R.id) continue;
                                                // todo what about cases of p (non random point) that is equal to representative?
                                                // make sure with adi and dan current solution makes sense
                                                Bit inGroup = cmpDim.at(p).at(R);
                                                //  inGroup *= isRepInPrevSlice; //maybe-fixme eod change! may need remove. check in the morning
                                                Point pointIsInSlice = p * inGroup;
                                                newSlice.addPoint(p * inGroup,  inGroup);*/

                        // p < R
                        helib::Ctxt isBelowCurrentRep(cmpDict[dim].at(R).at(p));
                        //    printNameVal(keysServer.decryptCtxt(isBelowCurrentRep));
                        helib::Ctxt isAboveAboveSmallerReps(isBelowCurrentRep); //todo other init
                        //    printNameVal(keysServer.decryptCtxt(isAboveAboveSmallerReps));

                        for (const Point &r: randomPoints[dim]) {
                            //      cout << "           --- other r: ";
                            //      printPoint(r, keysServer);
                            if (r.id == R.id) continue;
                            // make sure with adi and dan current solution makes sense

                            //  r > R (in which case we don't care about cmpDict results of p and r)
                            Bit otherRepIsAboveCurrentRep = cmpDict[dim].at(r).at(R);
                            // results in: Bit otherRepIsAboveCurrentRep = (r > R)
                            //    printNameVal(keysServer.decryptCtxt(otherRepIsAboveCurrentRep));

                            //  R > r    AND     p > r
                            Bit pIsAboveOtherSmallerRep(
                                    cmpDict[dim].at(R).at(r)); //pIsAboveOtherSmallerRep = (R > r)
                            pIsAboveOtherSmallerRep *= (cmpDict[dim].at(p).at(r));
                            // results in: Bit pIsAboveOtherSmallerRep = (R > r) * (p > r)
                            //    printNameVal(keysServer.decryptCtxt(pIsAboveOtherSmallerRep));

                            //   R < r  OR  [ R > r    AND     p > r ]
                            Bit pIsBelowCurrentRepAndAboveOtherRep = pIsAboveOtherSmallerRep;
                            pIsBelowCurrentRepAndAboveOtherRep *= otherRepIsAboveCurrentRep;
                            pIsBelowCurrentRepAndAboveOtherRep.negate();
                            pIsBelowCurrentRepAndAboveOtherRep += pIsAboveOtherSmallerRep;
                            pIsBelowCurrentRepAndAboveOtherRep += otherRepIsAboveCurrentRep;
                            // results in: Bit pIsBelowCurrentRepAndAboveOtherRep =
                            // pIsAboveOtherSmallerRep + otherRepIsAboveCurrentRep
                            // - pIsAboveOtherSmallerRep * otherRepIsAboveCurrentRep
                            //    printNameVal(keysServer.decryptCtxt(pIsBelowCurrentRepAndAboveOtherRep));

                            // this will hold the Product[ (R > r) && (p > r) | foreach r in randomPoints ]
                            isAboveAboveSmallerReps *= pIsBelowCurrentRepAndAboveOtherRep;
                            //    printNameVal(keysServer.decryptCtxt(isAboveAboveSmallerReps));
                        }
                        Bit isInGroup = isBelowCurrentRep;
                        isInGroup *= isAboveAboveSmallerReps;

                        isInGroup *= isRepInPrevSlice; //maybe-fixme eod change! may need remove. check in the morning

                        Point pointIsInSlice = p * isInGroup;
                        /*
                        //                    long pIsInGroup = keysServer.decryptCtxt(isInGroup);
                        //                    if (pIsInGroup) {
                        //                        cout << "       Point is in group of       ";
                        //                        printPoint(pointIsInSlice, keysServer);
                        //                    }
                        //                cout << "       Point is in group?       "
                        //                     << pIsInGroup ? "yes    " : "no   ";
                        //                printPoint(pointIsInSlice, keysServer);
                        //                cout << "\n";
                        */
                        newSlice.addPoint(pointIsInSlice, isInGroup);
                    }

                    groups[dim].emplace_back(newSlice);

                    //                    loggerDataServer.log(printDuration(
                    //                            t0_itr_rep,
                    //                            "split iteration (for single random point)"));
                }

                // todo handle tail points - points bigger than all the random points at current slice
                // init separate slice for tail points 
//                Slice tailSlice;
//                for (const Point &p:sliceFromPrevDim.includedPoints) {
//                    // init separate counter for tail points
//                    Bit pIsAboveAllReps = cmpDict[dim].at(p).at(tinyRandomPoint);
//                    for (const Point &R: randomPoints[dim])
//                        pIsAboveAllReps *= cmpDict[dim].at(p).at(R);
//                    tailSlice.addPoint(p * pIsAboveAllReps, pIsAboveAllReps);
//                }
//                groups[dim].emplace_back(tailSlice);

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
     * @brief calculate cell-means
     * @param cells a list of Cells (each cell is a list of encrypted points and encrypted size)
     * @param keysServer
     * @return a list of cells and their corresponding means
     * @returns std::vector<std::tuple<Point, Slice> >
     * */
    static std::vector<std::tuple<Point, Slice> >
    caculateCellMeans(
            const std::vector<Slice> &cells,
            const KeysServer &keysServer
    ) {
        helib::PubKey pubKey = cells[0].includedPoints[0].public_key;

        std::vector<std::tuple<Point, Slice> > cellsMeans;//(cells.size());

        for (const Slice &slice: cells) {
            const std::vector<Point> slicePoints = slice.includedPoints;
            const std::vector<Ctxt> sliceSize = slice.included;
            Point sum = Point::addManyPoints(slicePoints);
//            Point sum(cells[0].includedPoints[0].public_key);
//            for (Point &point: slicePoints) sum = sum + point;
            const Point mean = keysServer.getQuotientPoint(sum, sliceSize);
            Ctxt size(sliceSize[0].getPubKey());
            for (const Ctxt &ctxt: sliceSize) size += ctxt;
            printNameVal(keysServer.decryptCtxt(size));
            printNameVal(keysServer.decryptSize(sliceSize));

            cout << "in calcMeans, the currnt mean: " ;
            printPoint(mean, keysServer);
            cellsMeans.emplace_back(mean, slice);
        }
        return cellsMeans;
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
                int uniquePointsNum = 3 + random() % NUMBER_OF_POINTS,
                        clientsNum = uniquePointsNum;
                printNameVal(NUMBER_OF_POINTS);
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
        std::vector<Client> clients(NUMBER_OF_CLIENTS, Client(keysServer));
        for (Client &client:clients) {
            for (int dim = 0; dim < DIM; ++dim) tempArr[dim] = dist(mt);
            client.encryptPoint(tempArr);
        }
        return clients;
    }
*/

