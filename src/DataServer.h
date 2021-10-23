

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
            int m = 1 / EPSILON   //  -1
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

        const helib::PubKey &publicKey = points[0].public_key;
        std::map<int, std::vector<Slice> > slices;
        if (points.empty()) return slices;        // sanity check

        // initialize base level of data
        Slice startingSlice;
        for (auto const &point:points)
            startingSlice.addPoint(point, cmpDict[0].at(point).at(tinyRandomPoint));
        slices[-1].push_back(startingSlice);

        /**     for DBG  (todo remove)    **/
        long PisRepInPrevSlice;// = keysServer.decryptCtxt(isRepInPrevSlice);
        long PisInGroup;// = keysServer.decryptCtxt(isInGroup);
        long PpIsBelowCurrentRep;// = keysServer.decryptCtxt(pIsBelowCurrentRep);
        long PpIsAboveAllSmallerReps;// = keysServer.decryptCtxt(pIsAboveAllSmallerReps);
        long PotherRepIsAboveCurrentRep;// = keysServer.decryptCtxt(otherRepIsAboveCurrentRep);
        long PpIsAboveOtherSmallerRep;// = keysServer.decryptCtxt(pIsAboveOtherSmallerRep);
        long PpIsBelowCurrentRepAndAboveOtherRep;// = keysServer.decryptCtxt(pIsBelowCurrentRepAndAboveOtherRep);

        for (int dim = 0; dim < DIM; ++dim) {
            auto t0_itr_dim = CLOCK::now();     //  for logging, profiling, DBG

            for (const Slice &baseSlice:slices[dim - 1]) {
                /*
                // for each slice from the previous iteration
                // slices[-1][m^0] = { p | p form all_points }
                // slices[0] [m^1] = { p < Ri | p from slices[-1] | Ri from random_points[0] }
                // slices[1] [m^2] = { p < Ri & p < Rj | p from slices[0] | Ri from random_points[0]
                //                                                          | Rj from random_points[1] }
                //  ...
                // slices[DIM-1][m^DIM]  = { p < Ri & p < Rj & .... & p < Rj | p from slices[1]
                //                                                          | Rj from random_points[0]
                //                                                          | Rj from random_points[1]
                //                                                          |   ...
                //                                                          | Rj from random_points[DIM-1] }
                */

                for (const Point &R: randomPoints[dim]) {
                    auto t0_itr_rep = CLOCK::now();     //  for logging, profiling, DBG

                    cout << endl << endl;
                    printNameVal(dim) << "Base Slice (prev) reps: ";
                    for (auto const &baseRep: baseSlice.reps)
                        printPoint(baseRep, keysServer);
                    cout << endl << " ========== R: ";
                    printPoint(R, keysServer);
                    cout << " ========== " << endl;

                    Slice newSlice;
                    newSlice.addReps(baseSlice.reps);
                    newSlice.addRep(R);

                    CBit isRepInPrevSlice(cmpDict[dim].at(R).at(
                            R)); //todo or cmpDict[dim].at(R).at(tinyRandPoint)

                    if (0 < dim && !baseSlice.reps.empty())
                        // does this rep belong to the slice
                        isRepInPrevSlice *= cmpDict[dim - 1].at(baseSlice.reps[dim - 1]).at(R);
                    PisRepInPrevSlice = keysServer.decryptCtxt(isRepInPrevSlice);
                    // todo why cmp at prev dim and not current?

                    //  tailSlice.addRep(R);

//                    for (const Point &p:baseSlice.points) {
                    for (const PointTuple &pointTuple:baseSlice.pointTuples) {

                        //                        const Point &p = pointTuple.point;
                        const Point &p = pointTuple.first;

//                        CBit isPointInPrevSlice(pointTuple.isIn);
                        CBit isPointInPrevSlice(pointTuple.second);

                        /*
                         if (p.id == R.id) continue;
                        // todo what about cases of p (non random point) that is equal to representative?
                        // make sure with adi and dan current solution makes sense
                        */

                        cout << endl;

//                        if (7 == p.id) cout << "*******";// << endl;
                        cout << " ~~~~~~ current point: { id=" << p.id
                             << " [" << p.pCoordinatesDBG[0] << "," << p.pCoordinatesDBG[1] << "] ";
                        printPoint(p, keysServer);
                        cout << "}";

                        CBit isInGroup = isRepInPrevSlice;
                         isInGroup *= isPointInPrevSlice; //fixme new   <--------------

                        // p < R
                        CBit pIsBelowCurrentRep(cmpDict[dim].at(R).at(p));
                        PpIsBelowCurrentRep = keysServer.decryptCtxt(pIsBelowCurrentRep);

                        CBit pIsAboveAllSmallerReps(pIsBelowCurrentRep); //todo other init

                        for (const Point &r: randomPoints[dim]) {
                            cout << "           --- other r: ";
                            printPoint(r, keysServer);
                            if ((r.id == R.id) || (p.id == r.id)) continue;
                            // make sure with adi and dan current solution makes sense

                            //  r > R (in which case we don't care about cmpDict results of p and r)
                            CBit otherRepIsAboveCurrentRep = cmpDict[dim].at(r).at(R);
                            // results in: CBit otherRepIsAboveCurrentRep = (r > R)

                            //  R > r    AND     p > r
                            CBit pIsAboveOtherSmallerRep(cmpDict[dim].at(R).at(r));
                            // results in: CBit pIsAboveOtherSmallerRep = (R > r)
                            pIsAboveOtherSmallerRep *= (cmpDict[dim].at(p).at(r));
                            // results in: CBit pIsAboveOtherSmallerRep = (R > r) * (p > r)

                            //   [ R > r    AND     p > r ]   OR   r > R
                            //  should be   [ (R > r) * (p > r) ] + (r > R) - [ (R > r) * (p > r) ] * (r > R)
                            //  actually is [ (R > r) * (p > r) ] + (r > R) - [ (R > r) * (p > r) ] * (r > R)
                            CBit pIsBelowCurrentRepAndAboveOtherRep = pIsAboveOtherSmallerRep;
                            // results in: CBit pIsBelowCurrentRepAndAboveOtherRep
                            // = (R > r) * (p > r)
                            pIsBelowCurrentRepAndAboveOtherRep *= otherRepIsAboveCurrentRep;
                            // results in: CBit pIsBelowCurrentRepAndAboveOtherRep
                            // = [(R > r) * (p > r)] * (r > R)
                            pIsBelowCurrentRepAndAboveOtherRep.negate();
                            // results in: CBit pIsBelowCurrentRepAndAboveOtherRep
                            // =  - [(R > r) * (p > r)] * (r > R)
                            pIsBelowCurrentRepAndAboveOtherRep += pIsAboveOtherSmallerRep;
                            // results in: CBit pIsBelowCurrentRepAndAboveOtherRep
                            // = [(R > r) * (p > r)] - [(R > r) * (p > r)] * (r > R)
                            pIsBelowCurrentRepAndAboveOtherRep += otherRepIsAboveCurrentRep;
                            PpIsBelowCurrentRepAndAboveOtherRep = keysServer.decryptCtxt(
                                    pIsBelowCurrentRepAndAboveOtherRep);

                            // results in: CBit pIsBelowCurrentRepAndAboveOtherRep
                            // = (r > R) + [(R > r) * (p > r)]- [(R > r) * (p > r)] * (r > R)
                            // results in: CBit pIsBelowCurrentRepAndAboveOtherRep =
                            // pIsAboveOtherSmallerRep + otherRepIsAboveCurrentRep
                            // - pIsAboveOtherSmallerRep * otherRepIsAboveCurrentRep
                            PotherRepIsAboveCurrentRep = keysServer.decryptCtxt(
                                    otherRepIsAboveCurrentRep);
                            PpIsAboveOtherSmallerRep = keysServer.decryptCtxt(
                                    pIsAboveOtherSmallerRep);

                            // this will hold the Product[ (R > r) && (p > r) | foreach r in randomPoints ]
                            pIsAboveAllSmallerReps *= pIsBelowCurrentRepAndAboveOtherRep;
                            PpIsAboveAllSmallerReps = keysServer.decryptCtxt(
                                    pIsAboveAllSmallerReps);
                        }

                        isInGroup *= pIsBelowCurrentRep;
                        isInGroup *= pIsAboveAllSmallerReps;
                        PisInGroup = keysServer.decryptCtxt(isInGroup);

                        Point pointIsInSlice = p * isInGroup;
//                                                Point pointIsInSlice(p);
//                                                pointIsInSlice*= isInGroup;

                        const std::vector<long>
                                &decP = decryptPoint(p, keysServer);
                        const std::vector<long>
                                &decPInSlice = decryptPoint(pointIsInSlice, keysServer);

                        /**     for DBG     (todo remove)   **/
                        if (PisInGroup) {
                            /*
                             cout << endl;
                            cout << endl << "  ******  ";
                            printPoint(pointIsInSlice, keysServer);
                            cout << "      Point is in group of       ";
                            printNameVal(dim);
                            cout << "Base Slice (prev) reps: ";
                            for (auto const &baseRep: baseSlice.reps)
                                printPoint(baseRep, keysServer);
                            //                            cout << endl;
                            cout << " ========== current R: ";
                            printPoint(R, keysServer);
                            cout << "  ******  " << endl;

                            long psum = 0;
                            for (auto item:decPInSlice)
                                psum += item;
                            if (!psum)
                                cout << "$$$$" << endl;*/
                            cout << "\t\t$$$$$$";
                        }

                        newSlice.addPoint(pointIsInSlice, isInGroup);
                    }

                    slices[dim].emplace_back(newSlice);

                    //                    loggerDataServer.log(printDuration(
                    //                            t0_itr_rep,
                    //                            "split iteration (for single random point)"));
                }
                /*
                // todo handle tail points - points bigger than all the random points at current slice
                // init separate slice for tail points
                Slice tailSlice;
                for (const Point &p:baseSlice.points) {
                    // init separate counter for tail points
                    CBit pIsAboveAllReps = cmpDict[dim].at(p).at(tinyRandomPoint);
                    for (const Point &R: randomPoints[dim])
                        pIsAboveAllReps *= cmpDict[dim].at(p).at(R);
                    tailSlice.addPoint(p * pIsAboveAllReps, pIsAboveAllReps);
                }
                slices[dim].emplace_back(tailSlice);
                 */
            }
            cout << endl;
            loggerDataServer.log(printDuration(
                    t0_itr_dim,
                    "split iteration (for #" + std::to_string(dim) + " dimention)"));
        }
        loggerDataServer.log(printDuration(t0_split, "split total"));

        /*       todo  checkout this function, from @file Ctxt.h, could be used for all iterarion of 2nd rep at once?
                // set out=prod_{i=0}^{n-1} v[j], takes depth log n and n-1 products
                // out could point to v[0], but having it pointing to any other v[i]
                // will make the result unpredictable.
                void totalProduct(Ctxt& out, const std::vector<Ctxt>& v);
                */
        return slices;
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
        helib::PubKey pubKey = cells[0].points[0].public_key;

        std::vector<std::tuple<Point, Slice> > cellsMeans;//(cells.size());

        for (const Slice &slice: cells) {
            const std::vector<Point> slicePoints = slice.points;
            const std::vector<Ctxt> sliceSize = slice.counter;
            Point sum = Point::addManyPoints(slicePoints);
            //            Point sum(cells[0].points[0].public_key);
            //            for (Point &point: slicePoints) sum = sum + point;
            const Point mean = keysServer.getQuotientPoint(sum, sliceSize);
            Ctxt size(sliceSize[0].getPubKey());
            for (const Ctxt &ctxt: sliceSize) size += ctxt;
            printNameVal(keysServer.decryptCtxt(size));
            printNameVal(keysServer.decryptSize(sliceSize));

            cout << "in calcMeans, the currnt mean: ";
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

