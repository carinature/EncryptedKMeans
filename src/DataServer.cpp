
#include "utils/aux.h"
#include "DataServer.h"

#include <algorithm> //for the random shuffle

using std::cout;
using std::endl;

static Logger loggerDataServer(log_debug, "loggerDataServer");

//DataServer::DataServer(KeysServer &keysServer) ;

/**
 * @brief request data from clients and conentrate into one list
 * @param clients
 * @return std::vector<Point>
 * @returns a list containing all the points from all the clients in the data set
 * */
std::vector<Point>
DataServer::retrievePoints(
        const std::vector<Client> &clients) {
    auto t0_retrievePoints = CLOCK::now();     //  for logging, profiling, DBG// logging

    std::vector<Point> points;
    if (clients.empty()) return points;
    points.reserve(pow(clients.size(), 2) / 2); // preallocate memory
    for (const Client &c:clients)
        for (const Point &p:c.getPoints())
            points.emplace_back(Point(p));
    //            points.insert(points.end(),c.getPoints().begin(), c.getPoints().end());
    points.shrink_to_fit();

    loggerDataServer.log(printDuration(t0_retrievePoints, "retrievePoints"));

    return points;

}


void
DataServer::retrievePoints_Thread(
        const std::vector<Client> &clients) {
    auto t0_retrievePoints = CLOCK::now();     //  for logging, profiling, DBG// logging
    if (clients.empty()) return;

    for (const Client &c:clients)
        for (const Point &p:c.getPoints()) {
            retrievedPointsLock.lock();
            retrievedPoints.emplace_back(Point(p));
            retrievedPointsLock.unlock();
        }
    //            points.insert(points.end(),c.getPoints().begin(), c.getPoints().end());

    loggerDataServer.log(printDuration(t0_retrievePoints, "retrievePoints_Thread"));

}

void
DataServer::retrievePoints_WithThreads(
        const std::vector<Client> &clients,
        short numOfThreads
) {
    auto t0_retrievePoints = CLOCK::now();  //  for logging, profiling, DBG// logging
    if (clients.empty()) return;

    unsigned long splitSize = clients.size() / numOfThreads;
    std::vector<Client> splits[numOfThreads];
    std::vector<std::thread> threadVec;
    for (int i = 0; i < numOfThreads; ++i) {
        splits[i] = (i == numOfThreads - 1) ?
                    std::vector<Client>(clients.begin() + i * splitSize, clients.end()) :
                    std::vector<Client>(clients.begin() + i * splitSize,
                                        clients.begin() + (i + 1) * splitSize);
        threadVec.emplace_back(&DataServer::retrievePoints_Thread, this, splits[i]);
    }
    for (auto &t:threadVec) t.join();

    loggerDataServer.log(
            printDuration(t0_retrievePoints, "retrievePoints_WithThreads"));

}


const std::vector<std::vector<Point> > &
DataServer::pickRandomPoints(
        const std::vector<Point> &points,
        int m
) {
    auto t0_rndPoints = CLOCK::now();     //  for logging, profiling, DBG

    if (points.empty() || m > points.size()) return randomPointsList;

    //    const Point &tinyRandomPoint = keysServer.tinyRandomPoint();
    //    std::vector<std::vector<Point> > randomPoints(DIM);
    //    randomPointsList.reserve(DIM);
    //    randomPointsList = std::vector<std::vector<Point>>(DIM);

    for (int dim = 0; dim < DIM; ++dim) {

        // for every dimension random m^dim points (m points for every 'slice') and one tiny-point
        //        randomPoints[dim].reserve(1 + std::pow(m, dim));
        randomPointsList[dim].reserve(1 + std::pow(m, dim));
        //        to avoid cases of null (zero) points being assigned to group and blowing up in size,
        //          we add a rep for null points to whom those points will be assigned
        //        randomPoints.push_back(keysServer.tinyRandomPoint());
        // todo consider not saving it in the randomPoints since it's already saved in the dataServ
        //        randomPoints[dim].emplace_back(tinyRandomPoint);
        randomPointsList[dim].emplace_back(tinyRandomPoint);

        // choose random indices
        std::vector<int> indices(points.size());
        for (int i = 0; i < indices.size(); ++i) indices[i] = i;
        auto rd = std::random_device{};
        //        auto rng = std::default_random_engine{rd()}; //for extra randomness. recommended.
        auto rng = std::default_random_engine{};
        std::shuffle(std::begin(indices), std::end(indices), rng);
        //todo shuffle only m instead of all ? check which is more efficient

        for (int i = 0; i < pow(m, dim + 1); ++i) {
            //            randomPoints[dim].emplace_back(points[indices[i]]);
            randomPointsList[dim].emplace_back(retrievedPoints[indices[i]]);
        }
    }

    loggerDataServer.log(printDuration(t0_rndPoints, "pickRandomPoints"));
    //    return randomPoints;
    return randomPointsList;

}

CmpDict
DataServer::createCmpDict(
        const std::vector<Point> &allPoints,
        const std::vector<std::vector<Point> > &randomPoints
) {
    auto t0_cmpDict = CLOCK::now();     //  for logging, profiling, DBG

    std::vector<
            std::unordered_map<
                    const Point,
                    std::unordered_map<
                            const Point,
                            helib::Ctxt> > >
            cmpDict(DIM);

    std::vector<CBit> res;

    for (short dim = 0; dim < DIM; ++dim) {
        cmpDict[dim].reserve(randomPoints[dim].size());
        for (const Point &rep : randomPoints[dim]) {
            for (const Point &point : allPoints) {
                //                if (!cmpDict[rep].empty() && !cmpDict[rep][point].isEmpty()))
                //                if (rep == point) continue; //this is checked inside isBigger function
                // todo in the future, for efficiency - check if exist
                res = rep.isBiggerThan(point, dim);
                /*
                //  nand = !(mu*nu) = 1-(mu*nu)     meaning both numbers were equal
                Ctxt nand(res[0]);              // nand = mu
                (nand *= res[1]).negate();      // nand = - (mu * nu)
                nand.addConstant(1l);      // nand = 1 - (mu * nu) = !(mu*nu)

                //  beq = mu OR nand = mu + nand - mu * nand
                Ctxt beq(res[0]);               // beq = mu
                (beq *= nand).negate();           // beq = - (mu * nand)
                (beq += res[0]) += nand;            // beq = mu + nand - (mu * nand) = mu OR nand

                cmpDict[dim][rep].emplace(point, nand);   //  rep > point or rep=point
                */
                cmpDict[dim][rep].emplace(point, res[0]);   //  rep > point
                cmpDict[dim][point].emplace(rep, res[1]);   //  rep < point
            }
            res = rep.isBiggerThan(tinyRandomPoint, dim);
            cmpDict[dim][rep].emplace(tinyRandomPoint, res[0]);   //  rep > point2
            cmpDict[dim][tinyRandomPoint].emplace(rep, res[1]);   //  rep < point2
        }
    }

    loggerDataServer.log(printDuration(t0_cmpDict, "createCmpDict"));
    return cmpDict;
}

void
DataServer::createCmpDict_Dim_Thread(short dim) {

    auto t0_cmpDict_thread = CLOCK::now();     //  for logging, profiling, DBG
    std::vector<CBit> res;

    cmpDict[dim].reserve(randomPointsList[dim].size());

    for (const Point &rep : randomPointsList[dim]) {
        for (const Point &point : retrievedPoints) {
            //                if (!cmpDict[rep].empty() && !cmpDict[rep][point].isEmpty()))
            //                if (rep == point) continue; //this is checked inside isBigger function
            // todo in the future, for efficiency - check if exist
            res = rep.isBiggerThan(point, dim);
            /*
            //  nand = !(mu*nu) = 1-(mu*nu)     meaning both numbers were equal
            Ctxt nand(res[0]);              // nand = mu
            (nand *= res[1]).negate();      // nand = - (mu * nu)
            nand.addConstant(1l);      // nand = 1 - (mu * nu) = !(mu*nu)

            //  beq = mu OR nand = mu + nand - mu * nand
            Ctxt beq(res[0]);               // beq = mu
            (beq *= nand).negate();           // beq = - (mu * nand)
            (beq += res[0]) += nand;            // beq = mu + nand - (mu * nand) = mu OR nand

            cmpDict[dim][rep].emplace(point, nand);   //  rep > point or rep=point
            */
            cmpDictLock.lock();
            cmpDict[dim][rep].emplace(point, res[0]);   //  rep > point
            cmpDict[dim][point].emplace(rep, res[1]);   //  rep < point
            cmpDictLock.unlock();
        }
        res = rep.isBiggerThan(tinyRandomPoint, dim);
        cmpDictLock.lock();
        cmpDict[dim][rep].emplace(tinyRandomPoint, res[0]);   //  rep > point2
        cmpDict[dim][tinyRandomPoint].emplace(rep, res[1]);   //  rep < point2
        cmpDictLock.unlock();
    }
    loggerDataServer.log(printDuration(t0_cmpDict_thread, "createCmpDict_Dim_Thread"));
}

void
DataServer::createCmpDict_WithThreads(
        short numOfThreads
) {
    auto t0_cmpDict_withThreads = CLOCK::now();     //  for logging, profiling, DBG

    numOfThreads = DIM; //fixme consider threading by a different param
    std::vector<std::thread> threadVec;

    //    cmpDict.reserve(DIM);

    for (short dim = 0; dim < DIM; ++dim) {

        threadVec.emplace_back(&DataServer::createCmpDict_Dim_Thread,
                               this,
                               dim);
    }
    for (auto &t:threadVec) t.join();

    loggerDataServer.log(
            printDuration(t0_cmpDict_withThreads, "createCmpDict_WithThreads"));
}

std::map<int, //DIM
        std::vector< //current slices for approp dimension
                Slice
        >
>
DataServer::splitIntoEpsNet(
        const std::vector<Point> &points,
        const std::vector<std::vector<Point> > &randomPoints,
        const CmpDict &cmpDict,
        const KeysServer &keysServer
) {
    auto t0_split = CLOCK::now();     //  for logging, profiling, DBG

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
            auto t0_itr_slice = CLOCK::now();     //  for logging, profiling, DBG
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

                /*
                                    cout << endl << endl;
                                    printNameVal(dim) << "Base Slice (prev) reps: ";
                                    for (auto const &baseRep: baseSlice.reps)
                                        printPoint(baseRep, keysServer);
                                    cout << endl << " ========== R: ";
                                    printPoint(R, keysServer);
                                    cout << " ========== " << endl;
                */

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
                    // if current poins is the current Rep no need to add it now - it will be added to the group later
                    //                        if (R==p) continue; //  R==p means R.id==p.id
                    //                        CBit isPointInPrevSlice(pointTuple.isIn);
                    CBit isPointInPrevSlice(pointTuple.second);

                    /*
                     if (p == R) continue;
                    // todo what about cases of p (non random point) that is equal to representative?
                    // make sure with adi and dan current solution makes sense
                    */
                    /*
                                            cout << endl;
                                            //                        if (7 == p.id) cout << "*******";// << endl;
                                            cout << " ~~~~~~ current point: { id=" << p.id
                                                 << " [" << p.pCoordinatesDBG[0] << "," << p.pCoordinatesDBG[1] << "] ";
                                            printPoint(p, keysServer);
                                            cout << "}";
                    */

                    CBit isInGroup = isRepInPrevSlice;
                    isInGroup *= isPointInPrevSlice; //fixme new   <--------------

                    // p < R
                    CBit pIsBelowCurrentRep(cmpDict[dim].at(R).at(p));
                    PpIsBelowCurrentRep = keysServer.decryptCtxt(pIsBelowCurrentRep);

                    CBit pIsAboveAllSmallerReps(pIsBelowCurrentRep); //todo other init

                    for (const Point &r: randomPoints[dim]) {
                        /*                            cout << "       --- other r: ";
                                                    printPoint(r, keysServer);*/
                        if ((r == R) || (p == r)) continue;
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
                    /*   if (PisInGroup) {
                           *//*
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
                                cout << "$$$$" << endl;*//*
                            cout << "\t\t$$$$$$";
                        }*/

                    newSlice.addPoint(pointIsInSlice, isInGroup);
                }

                slices[dim].emplace_back(newSlice);

                loggerDataServer.log(printDuration(t0_itr_rep, "Split Random Rep iteration"));
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
            cout << endl;
            loggerDataServer.log(printDuration(t0_itr_slice, "Split Slice iteration"));
        }
        cout << endl;
        loggerDataServer.log(printDuration(
                t0_itr_dim,
                "Split iteration for #" + std::to_string(dim) + " dimension"));
    }
    loggerDataServer.log(printDuration(t0_split, "splitIntoEpsNet"));

    /*       todo  checkout this function, from @file Ctxt.h, could be used for all iterarion of 2nd rep at once?
            // set out=prod_{i=0}^{n-1} v[j], takes depth log n and n-1 products
            // out could point to v[0], but having it pointing to any other v[i]
            // will make the result unpredictable.
            void totalProduct(Ctxt& out, const std::vector<Ctxt>& v);
            */
    return slices;
}

std::mutex slicesLock;

void
DataServer::splitIntoEpsNet_R_Thread(
        const Slice &baseSlice,
        const Point &R,
        int dim
) {
    auto t0_itr_rep = CLOCK::now();     //  for logging, profiling, DBG

    /*
            cout << endl << endl;
            printNameVal(dim) << "Base Slice (prev) reps: ";
            for (auto const &baseRep: baseSlice.reps)
                printPoint(baseRep, keysServer);
            cout << endl << " ========== R: ";
            printPoint(R, keysServer);
            cout << " ========== " << endl;
*/
    /**     for DBG  (todo remove)    **/
    long PisRepInPrevSlice;// = keysServer.decryptCtxt(isRepInPrevSlice);
    long PisInGroup;// = keysServer.decryptCtxt(isInGroup);
    long PpIsBelowCurrentRep;// = keysServer.decryptCtxt(pIsBelowCurrentRep);
    long PpIsAboveAllSmallerReps;// = keysServer.decryptCtxt(pIsAboveAllSmallerReps);
    long PotherRepIsAboveCurrentRep;// = keysServer.decryptCtxt(otherRepIsAboveCurrentRep);
    long PpIsAboveOtherSmallerRep;// = keysServer.decryptCtxt(pIsAboveOtherSmallerRep);
    long PpIsBelowCurrentRepAndAboveOtherRep;// = keysServer.decryptCtxt(pIsBelowCurrentRepAndAboveOtherRep);

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

    //                    for (const Point &p:baseSlice.retrievedPoints) {
    for (const PointTuple &pointTuple : baseSlice.pointTuples) {

        //                        const Point &p = pointTuple.point;
        const Point &p = pointTuple.first;
        // if current poins is the current Rep no need to add it now - it will be added to the group later
        //                        if (R==p) continue; //  R==p means R.id==p.id
        //                        CBit isPointInPrevSlice(pointTuple.isIn);
        CBit isPointInPrevSlice(pointTuple.second);

        /*
         if (p == R) continue;
        // todo what about cases of p (non random point) that is equal to representative?
        // make sure with adi and dan current solution makes sense
        */
        /*
                                cout << endl;
                                //                        if (7 == p.id) cout << "*******";// << endl;
                                cout << " ~~~~~~ current point: { id=" << p.id
                                     << " [" << p.pCoordinatesDBG[0] << "," << p.pCoordinatesDBG[1] << "] ";
                                printPoint(p, keysServer);
                                cout << "}";
        */

        CBit isInGroup = isRepInPrevSlice;
        isInGroup *= isPointInPrevSlice; //fixme new   <--------------

        // p < R
        CBit pIsBelowCurrentRep(cmpDict[dim].at(R).at(p));
        PpIsBelowCurrentRep = keysServer.decryptCtxt(pIsBelowCurrentRep);

        CBit pIsAboveAllSmallerReps(pIsBelowCurrentRep); //todo other init

        for (const Point &r: randomPointsList[dim]) {
            /*                            cout << "       --- other r: ";
                                        printPoint(r, keysServer);*/
            if ((r == R) || (p == r)) continue;
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
        /*   if (PisInGroup) {
               *//*
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
                    cout << "$$$$" << endl;*//*
                cout << "\t\t$$$$$$";
            }*/

        newSlice.addPoint(pointIsInSlice, isInGroup);
    }
    slicesLock.lock();
    slices[dim].emplace_back(newSlice);
    slicesLock.unlock();

    loggerDataServer.log(printDuration(t0_itr_rep, "Split Random Rep Thread"));
}

std::map<int, //DIM
        std::vector< //current slices for approp dimension
                Slice
        >
>
DataServer::splitIntoEpsNet_WithThreads() {
    auto t0_split = CLOCK::now();     //  for logging, profiling, DBG

    //    std::map<int, std::vector<Slice> > slices;
    if (retrievedPoints.empty()) return slices;        // sanity check

    // initialize base level of data
    Slice startingSlice;
    for (auto const &point:retrievedPoints)
        startingSlice.addPoint(point, cmpDict[0].at(point).at(tinyRandomPoint));
    slices[-1].push_back(startingSlice);

    for (int dim = 0; dim < DIM; ++dim) {
        auto t0_itr_dim = CLOCK::now();     //  for logging, profiling, DBG

        for (const Slice &baseSlice : slices[dim - 1]) {
            auto t0_itr_slice = CLOCK::now();     //  for logging, profiling, DBG
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

            std::vector<std::thread> threadVec;

            for (const Point &R : randomPointsList[dim]) {
                threadVec.emplace_back(&DataServer::splitIntoEpsNet_R_Thread,
                                       this,
                                       baseSlice,
                                       R,
                                       dim
                                       );

            }

            for (auto &t:threadVec) t.join();

            /*
            // todo handle tail retrievedPoints - retrievedPoints bigger than all the random retrievedPoints at current slice
            // init separate slice for tail retrievedPoints
            Slice tailSlice;
            for (const Point &p:baseSlice.retrievedPoints) {
                // init separate counter for tail retrievedPoints
                CBit pIsAboveAllReps = cmpDict[dim].at(p).at(tinyRandomPoint);
                for (const Point &R: randomPoints[dim])
                    pIsAboveAllReps *= cmpDict[dim].at(p).at(R);
                tailSlice.addPoint(p * pIsAboveAllReps, pIsAboveAllReps);
            }
            slices[dim].emplace_back(tailSlice);
             */
            cout << endl;
            loggerDataServer.log(printDuration(t0_itr_slice, "Split Slice iteration"));
        }


        cout << endl;
        loggerDataServer.log(printDuration(
                t0_itr_dim,
                "Split iteration for #" + std::to_string(dim) + " dimension"));
    }
    loggerDataServer.log(printDuration(t0_split, "splitIntoEpsNet"));

    /*       todo  checkout this function, from @file Ctxt.h, could be used for all iterarion of 2nd rep at once?
            // set out=prod_{i=0}^{n-1} v[j], takes depth log n and n-1 products
            // out could point to v[0], but having it pointing to any other v[i]
            // will make the result unpredictable.
            void totalProduct(Ctxt& out, const std::vector<Ctxt>& v);
            */
    return slices;
}

std::vector<std::tuple<Point, Slice> >
DataServer::calculateSlicesMeans(const std::vector<Slice> &slices,
                                 const KeysServer &keysServer) {
    auto t0_means = CLOCK::now();     //  for logging, profiling, DBG

    std::vector<std::tuple<Point, Slice> > slicesMeans;//(slices.size());

    for (const Slice &slice: slices) {

        std::vector<Point> points;
        points.reserve(slice.reps.size() + slice.points.size()); // preallocate memory
        points.insert(points.end(), slice.reps.begin(), slice.reps.end());
        points.insert(points.end(), slice.points.begin(), slice.points.end());
        Point sum(Point::addManyPoints(points, keysServer));

        const Point mean(keysServer.getQuotientPoint(sum, slice.counter, DIM));

        cout << "slice reps: ";
        printPoints(slice.reps, keysServer);
        cout << endl;
        cout << "slice points: ";
        printNonEmptyPoints(slice.points, keysServer);
        cout << endl;
        cout << "the currnt mean: ";
        printPoint(mean, keysServer);
        slicesMeans.emplace_back(mean, slice);
    }

    loggerDataServer.log(printDuration(t0_means, "calculateSlicesMeans"));

    return slicesMeans;
}

std::vector<Point> DataServer::collectMeans(const std::vector<std::tuple<Point, Slice>> &slices,
                                            const KeysServer &keysServer) {
    std::vector<Point> means;
    means.reserve(slices.size());
    for (auto const &slice:slices) means.push_back(std::get<0>(slice));
    return means;
}

std::vector<std::tuple<Point, Point, EncryptedNum> >
DataServer::collectMinimalDistancesAndClosestPoints(const std::vector<Point> &points,
                                                    const std::vector<Point> &means,
                                                    const KeysServer &keysServer) {
    auto t0_collectMinDist = CLOCK::now();

    std::vector<std::tuple<Point, Point, EncryptedNum> > minDistanceTuples;
    minDistanceTuples.reserve(points.size());
    for (const Point &point: points) {
        const std::pair<Point, EncryptedNum> &
                minDistFromMeans = point.findMinDistFromMeans(means, keysServer);
        minDistanceTuples.emplace_back(
                point,
                minDistFromMeans.first,
                minDistFromMeans.second);
    }

    loggerDataServer.log(
            printDuration(t0_collectMinDist, "collectMinimalDistancesAndClosestPoints"));

    return minDistanceTuples;
}

EncryptedNum DataServer::calculateThreshold(
        const std::vector<std::tuple<Point, Point, EncryptedNum>> &minDistanceTuples,
        const KeysServer &keysServer, int iterationNumber) {
    //  collect minimal distances
    EncryptedNum sum;
    helib::CtPtrs_vectorCt sum_wrapper(sum);
    std::vector<EncryptedNum> summandsVec(minDistanceTuples.size());
    for (auto const &tuple:minDistanceTuples) summandsVec.push_back(std::get<2>(tuple));
    helib::CtPtrMat_vectorCt summands_wrapper(summandsVec);

    //  sum all distance
    addManyNumbers(
            sum_wrapper,
            summands_wrapper
    );

    //  find average distance
    long num = long(NUMBER_OF_POINTS / pow(2, iterationNumber));
    printNameVal(pow(2, iterationNumber));
    printNameVal(num);
    EncryptedNum
            threshold =
            keysServer.getQuotient(
                    sum,
                    num);
    return threshold;
}

std::tuple<
        std::unordered_map<long, std::vector<std::pair<Point, CBit> > >,
        std::vector<std::pair<Point, CBit> >,
        std::vector<std::pair<Point, CBit> >
> DataServer::choosePointsByDistance(
        const std::vector<std::tuple<Point, Point, EncryptedNum>> &minDistanceTuples,
        std::vector<Point> means, EncryptedNum &threshold) {

    std::unordered_map<
            long, //mean index
            std::vector<std::pair<Point, CBit> > > groups;
    std::vector<std::pair<Point, CBit> > closest;
    std::vector<std::pair<Point, CBit> > farthest;

    for (auto const &tuple:minDistanceTuples) {
        Point point = std::get<0>(tuple);
        Point meanClosest = std::get<1>(tuple);
        EncryptedNum distance = std::get<2>(tuple);
        const helib::PubKey &public_key = point.public_key;

        //  check if distance within threshold margin
        helib::Ctxt mu(public_key), ni(public_key);
        helib::CtPtrs_vectorCt threshold_wrpr(threshold), distance_wrpr(distance);
        helib::compareTwoNumbers(mu, ni, threshold_wrpr, distance_wrpr); // fixme
        //  pick all points with distance bigger than avg
        farthest.emplace_back(point * mu, mu);
        //  pick all points with distance smaller than avg
        closest.emplace_back(point * ni, ni);

        for (int i = 0; i < means.size(); ++i) {
            //  check if the closest mean to the point is the current one
            helib::Ctxt muCid(public_key), niCid(public_key);
            helib::CtPtrs_vectorCt closestCid(meanClosest.cid), meanCid(means[i].cid);
            helib::compareTwoNumbers(muCid, niCid, closestCid, meanCid);

            //  check if point is within margin and her closest mean equals to the current one
            helib::Ctxt notMuCid(muCid);
            notMuCid.negate();
            notMuCid.addConstant(1l);
            //  result in isCloseToCurrentMean = !(muCid)
            helib::Ctxt isCloseToCurrentMean(niCid);
            isCloseToCurrentMean.negate();
            isCloseToCurrentMean.addConstant(1l);
            isCloseToCurrentMean *= notMuCid;
            //  result in isCloseToCurrentMean = !(muCid) && !(niCid)
            isCloseToCurrentMean *= ni;
            //  result in isCloseToCurrentMean = !(muCid) && !(niCid) && ni

            //  pick all points with distance smaller than avg, arrange by closest mean point
            groups[i].emplace_back(point * isCloseToCurrentMean, isCloseToCurrentMean);
        }
    }
    return {groups, closest, farthest};
}
