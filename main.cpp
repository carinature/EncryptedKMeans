//
// run a full protocol
//
#include <iostream>
#include <chrono>

#include "utils/aux.h"
#include "src/DataServer.h"

#include <bitset>
#include <src/coreset/run1meancore.h>

static Logger logger(log_debug, "loggerMain");
//static Logger loggerMain(log_debug, "loggerMain");

int main() {
    auto t0_main = CLOCK::now();

    // CT for coreset.csv result for multiple iterations
    auto t = time(nullptr);
    auto tm = *localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y_%m_%d_%H_%M_%S");
    string timestamp = oss.str();
    //    oss.flush();
    oss.clear();
    string prefix = IO_DIR + timestamp + "_";

    //    Logger logger;
    logger.log("Starting Protocol", log_trace);
    ////  Keys Server
    KeysServer keysServer;
    logger.log(printDuration(t0_main, "KeysServer Initialization"));
    DataServer dataServer(keysServer);

    ////    (generate data)
    const std::vector<ClientDevice> clients = generateDataClients(keysServer);

    ////    Retrieve Data from Clients
    const std::vector<Point> points = dataServer.retrievePoints_WithThreads(clients);
    cout << " ---   Points  ---" << endl;
    cout << printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    ////    Pick Random Representatives
    const Point &tinyRandomPoint = keysServer.tinyRandomPoint();

    const std::vector<std::vector<Point> >
            randomPoints = dataServer.pickRandomPoints(points);//, (1 / EPSILON)-1);

    cout << " ---   Random Points  ---" << endl;
    for (const auto &vec :randomPoints) cout << printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    ////  Create points-comparing dict
    ///     - for every 2 points (p1,p2) answers p1[dim]>p2[dim]
    const std::vector<
            std::unordered_map<
                    const Point,
                    std::unordered_map<
                            const Point,
                            helib::Ctxt> > >
            cmpDict = dataServer.createCmpDict_WithThreads(
            points,
            randomPoints
    );

    cout << " ---   The Dictionary  ---" << endl;
    for (short dim = 0; dim < DIM; ++dim) {
        cout << "    ======   ";
        printNameVal(dim);
        for (auto const&[point, cmpResults] : cmpDict[dim]) {
            cout << " --- ";
            cout << printPoint(point, keysServer);
            cout << " --- \n";
            for (auto const&[point2, val]: cmpResults) {
                cout << printPoint(point2, keysServer);
                long pVal = keysServer.decryptCtxt(val);
                printNameVal(pVal);
            }
            printNameVal(cmpResults.size()) << " --- --- ---" << endl;
        }
        printNameVal(cmpDict[dim].size()) << " === === ===" << endl;
    }
    cout << " --- --- --- --- ---" << endl;


    ////    Eps-Net - Split Data into Slices
    std::map<int, std::vector<Slice> >
            epsNet =
            dataServer.splitIntoEpsNet_WithThreads(
                    points,
                    randomPoints,
                    cmpDict
            );

    cout << " ---   Epsilon-Net  ---" << endl;
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "   ---   For dim " << dim << "  --- " << endl;
        for (Slice &cell: epsNet[dim]) {
            cell.printSlice(keysServer);
        }
        cout << "   ---     --- " << endl;
        cout << endl;
    }
    cout << " --- --- --- --- ---" << endl;

    ////    Calculate Eps-Net means
    std::vector<std::tuple<Point, Slice> >
            meanCellTuples = dataServer.calculateSlicesMeans_WithThreads(epsNet[DIM - 1]);

    cout << " ---   Means  ---" << endl;
    for (auto const &tup:meanCellTuples) {
        cout << "The Mean is: ";
        cout << printPoint(std::get<0>(tup), keysServer);
        cout << endl;
        std::get<1>(tup).printSlice(keysServer);
    }
    cout << " --- --- --- --- ---" << endl;

    ////    Find Minimal Distances from Means and the Closest Mean
    std::vector means = dataServer.collectMeans(meanCellTuples);
    std::vector<std::tuple<Point, Point, EncryptedNum> >
            minDistanceTuples =
            dataServer.collectMinimalDistancesAndClosestPoints_WithThreads(
                    points,
                    means);

    cout << " ---   Minimal Distances and Closest Means  ---" << endl;
    for (const auto &tuple:minDistanceTuples) {
        cout << " Point: ";
        cout << printPoint(std::get<0>(tuple), keysServer);
        cout << " Closest mean: ";
        cout << printPoint(std::get<1>(tuple), keysServer);
        long minDistance = keysServer.decryptNum(std::get<2>(tuple));
        printNameVal(minDistance);
    }
    cout << " --- --- --- --- ---" << endl;

    ////    Calculate Threshold
    EncryptedNum
            threshold =
            dataServer.calculateThreshold(
                    minDistanceTuples,
                    0);

    printNameVal(keysServer.decryptNum(threshold));

    std::tuple<
            std::unordered_map<long, std::vector<std::pair<Point, CBit> > >,
            //            std::vector<std::pair<Point, CBit> >,
            std::vector<std::pair<Point, CBit> >
    > groups = dataServer.choosePointsByDistance_WithThreads(
            minDistanceTuples,
            means,
            threshold
    );


    cout << " === === === === === ===" << endl;
    cout << " === FINAL PRINTOUT ===" << endl;
    cout << " === === === === === ===" << endl;
    cout << endl;
    cout << " ===   All Points  ===" << endl;
    cout << printPoints(points, keysServer);
    cout << " === === === === ===" << endl;
    cout << endl;
    cout << " ===   Random Points  ===" << endl;
    for (auto const &vec :randomPoints) cout << printPoints(vec, keysServer);
    cout << " === === === === ===" << endl;
    cout << endl;
    cout << " ===   Slices  ===" << endl;
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "   ---   For dim " << dim << "  --- " << endl;
        for (const Slice &slice: epsNet[dim]) slice.printSlice(keysServer);
        cout << endl;
    }
    cout << " === === === === ===" << endl;
    cout << endl;
    cout << " ===   Means  ===" << endl;
    for (auto const &tup:meanCellTuples) {
        cout << "The Mean is: ";
        cout << printPoint(std::get<0>(tup), keysServer);
        cout << endl;
        std::get<1>(tup).printSlice(keysServer);
    }
    cout << " === === === === ===" << endl;
    cout << endl;
    cout << " ===   Minimal Distances and Closest Means  ===" << endl;
    for (const auto &tuple:minDistanceTuples) {
        cout << " Point: ";
        cout << printPoint(std::get<0>(tuple), keysServer);
        cout << " Closest mean: ";
        cout << printPoint(std::get<1>(tuple), keysServer);
        long minDistance = keysServer.decryptNum(std::get<2>(tuple));
        printNameVal(minDistance);
    }
    cout << " === === === === ===" << endl;
    cout << endl;
    cout << " --- --- --- --- ---" << endl;

    printNameVal(keysServer.decryptNum(threshold));

    cout << " === Groups by Means === " << endl;
    std::unordered_map<long, std::vector<std::pair<Point, CBit> > >
            groups_by_means = std::get<0>(groups);
    for (auto const &[meanI, points] : groups_by_means) {
        cout << printPoint(means[meanI], keysServer);
        printNameVal(meanI) << "Close Points: ";
        std::vector<Point> pointsGroup;
        pointsGroup.reserve(points.size());
        for (auto const &pair: points) {
            cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
            cout << printPoint(pair.first, keysServer);
            pointsGroup.emplace_back(pair.first);
        }
        cout << endl;
        printNonEmptyPoints(pointsGroup, keysServer);
        cout << endl;

        //        decAndWriteToFile(pointsGroup, IO_DIR + to_string(meanI) + "_" + CHOSEN_FILE, keysServer);

    }
    cout << " === === === === === " << endl << endl;
    //
    //    cout << " === Closest Points === " << endl;
    //    for (auto const &pair : std::get<1>(groups)) {
    //        cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
    //        cout << printPoint(pair.first, keysServer);
    //    }
    //    cout << endl << " === === === === === " << endl << endl;

    cout << " === Farthest Points === " << endl;
    std::vector<std::pair<Point, CBit> > farthest = std::get<1>(groups);
    std::vector<Point> leftover;
    leftover.reserve(points.size());
    for (auto const &[point, isInFlag] : farthest) {
        cout << "-" << keysServer.decryptCtxt(isInFlag) << "-";
        cout << printPoint(point, keysServer);
        leftover.emplace_back(point);
    }
    leftover.shrink_to_fit();
    cout << endl << " === === === === === " << endl << endl;


    /**********   integration to coreset alg    ***********/
    //  for each mean-group in groups
    for (auto const &[meanI, points] : groups_by_means) {
        std::vector<std::vector<double>> pointsGroup;
        pointsGroup.reserve(points.size());
        std::vector<Point> forClean;
        // add all non-zero points
        //        for (auto const &pair: points) {
        //            const Point &currPoint(pair.first);
        for (auto const &[currPoint, isInFlag]: points) {
            std::vector<long> tempPoint = decryptPoint(currPoint, keysServer);
            std::vector<double> doublePoint(DIM);
            long checkNull = 0;
            for (const long &coordinate:decryptPoint(currPoint, keysServer)) {
                checkNull += coordinate;
                doublePoint.emplace_back(coordinate / CONVERSION_FACTOR);
            }
            //            if (checkNull) pointsGroup.emplace_back(doublePoint);
            pointsGroup.emplace_back(doublePoint);
            forClean.emplace_back(currPoint);
        }
        pointsGroup.shrink_to_fit();
        cout << "For Group of Points:\t";
        printNonEmptyPoints(forClean, keysServer);
        cout << endl;
        //      call coreset alg with mean-group
        cout << "Running 1-Mean Coreset Algorithm:" << endl;
        auto t0_itr_rep = CLOCK::now();     //  for logging, profiling, DBG
        runCoreset(pointsGroup, pointsGroup.size(), DIM, EPSILON);  // <|-------------
        logger.log(printDuration(t0_itr_rep, "runCoreset in iteration"));

    }

    string currPrefix = prefix + "iter_"; //fixme add iteration number
    //    string filename = currPrefix + "_coreset.csv";

    cout << "WRITE TO FILES" << endl;
    decAndWriteToFile(points, currPrefix + POINTS_FILE, keysServer);
    decAndWriteToFile(points, currPrefix + POINTS_COPY_FILE, keysServer);
    decAndWriteToFile(randomPoints[DIM - 1], currPrefix + RANDS_FILE, keysServer);
    decAndWriteToFile(means, currPrefix + MEANS_FILE, keysServer);
    decAndWriteToFile(leftover, currPrefix + LEFTOVER_FILE, keysServer);
    for (auto const &[meanI, points] : groups_by_means) {
        std::vector<Point> pointsGroup;
        pointsGroup.reserve(points.size());
        for (auto const &pair: points) pointsGroup.emplace_back(pair.first);
        decAndWriteToFile(pointsGroup,
                          currPrefix + to_string(meanI).append("_").append(CHOSEN_FILE),
                          keysServer);
    }

    //-----------------------------------------------------------

    logger.log(printDuration(t0_main, "Main"));
    logger.print_log(log_trace);//, false);

    return 0;
}

int main_without_threads() {
    auto t0_main = CLOCK::now();
    logger.log("Starting Protocol", log_trace);
    ////  Keys Server
    KeysServer keysServer;
    logger.log(printDuration(t0_main, "KeysServer Initialization"));
    DataServer dataServer(keysServer);

    ////    (generate data)
    const std::vector<ClientDevice> clients = generateDataClients(keysServer);

    ////    Retrieve Data from Clients
    const std::vector<Point> points = dataServer.retrievePoints(clients);
    cout << " ---   Points  ---" << endl;
    cout << printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    ////    Pick Random Representatives
    const Point &tinyRandomPoint = keysServer.tinyRandomPoint();

    const std::vector<std::vector<Point> >
            randomPoints = dataServer.pickRandomPoints(points);//, (1 / EPSILON)-1);

    cout << " ---   Random Points  ---" << endl;
    for (const auto &vec :randomPoints) cout << printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    ////  Create points-comparing dict
    ///     - for every 2 points (p1,p2) answers p1[dim]>p2[dim]
    const std::vector<
            std::unordered_map<
                    const Point,
                    std::unordered_map<
                            const Point,
                            helib::Ctxt> > >
            cmpDict = dataServer.createCmpDict(points, randomPoints);

    cout << " ---   The Dictionary  ---" << endl;
    for (short dim = 0; dim < DIM; ++dim) {
        cout << "    ======   ";
        printNameVal(dim);
        for (auto const&[point, map] : cmpDict[dim]) {
            cout << printPoint(point, keysServer);
            for (auto const&[point2, val]: map) {
                cout << printPoint(point2, keysServer);
                long pVal = keysServer.decryptCtxt(val);
                printNameVal(pVal);
            }
            printNameVal(map.size()) << " --- --- ---" << endl;
        }
        printNameVal(cmpDict[dim].size()) << " === === ===" << endl;
    }
    cout << " --- --- --- --- ---" << endl;


    ////    Eps-Net - Split Data into Slices
    std::map<int, std::vector<Slice> >
            epsNet = dataServer.splitIntoEpsNet(points, randomPoints, cmpDict);

    cout << " ---   Epsilon-Net  ---" << endl;
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "   ---   For dim " << dim << "  --- " << endl;
        for (Slice &cell: epsNet[dim]) {
            cell.printSlice(keysServer);
        }
        cout << "   ---     --- " << endl;
        cout << endl;
    }
    cout << " --- --- --- --- ---" << endl;

    ////    Calculate Eps-Net means
    std::vector<std::tuple<Point, Slice> >
            meanCellTuples = dataServer.calculateSlicesMeans(epsNet[DIM - 1]);

    cout << " ---   Means  ---" << endl;
    for (auto const &tup:meanCellTuples) {
        cout << "The Mean is: ";
        cout << printPoint(std::get<0>(tup), keysServer);
        cout << endl;
        std::get<1>(tup).printSlice(keysServer);
    }
    cout << " --- --- --- --- ---" << endl;

    ////    Find Minimal Distances from Means and the Closest Mean
    std::vector means = dataServer.collectMeans(meanCellTuples);
    std::vector<std::tuple<Point, Point, EncryptedNum> >
            minDistanceTuples =
            dataServer.collectMinimalDistancesAndClosestPoints(
                    points,
                    means);

    cout << " ---   Minimal Distances and Closest Means  ---" << endl;
    for (const auto &tuple:minDistanceTuples) {
        cout << " Point: ";
        cout << printPoint(std::get<0>(tuple), keysServer);
        cout << " Closest mean: ";
        cout << printPoint(std::get<1>(tuple), keysServer);
        long minDistance = keysServer.decryptNum(std::get<2>(tuple));
        printNameVal(minDistance);
    }
    cout << " --- --- --- --- ---" << endl;

    ////    Calculate Threshold
    EncryptedNum
            threshold =
            dataServer.calculateThreshold(
                    minDistanceTuples,
                    0);

    printNameVal(keysServer.decryptNum(threshold));

    std::tuple<
            std::unordered_map<long, std::vector<std::pair<Point, CBit> > >,
            std::vector<std::pair<Point, CBit> >,
            std::vector<std::pair<Point, CBit> >
    > groups = dataServer.choosePointsByDistance(
            minDistanceTuples,
            means,
            threshold
    );


    cout << " === === === === === ===" << endl;
    cout << " === FINAL PRINTOUT ===" << endl;
    cout << " === === === === === ===" << endl;
    cout << endl;
    cout << " ===   All Points  ===" << endl;
    cout << printPoints(points, keysServer);
    cout << " === === === === ===" << endl;
    cout << endl;
    cout << " ===   Random Points  ===" << endl;
    for (auto const &vec :randomPoints) cout << printPoints(vec, keysServer);
    cout << " === === === === ===" << endl;
    cout << endl;
    cout << " ===   Slices  ===" << endl;
    for (int dim = 0; dim < DIM; ++dim) {
        cout << "   ---   For dim " << dim << "  --- " << endl;
        for (const Slice &slice: epsNet[dim]) slice.printSlice(keysServer);
        cout << endl;
    }
    cout << " === === === === ===" << endl;
    cout << endl;
    cout << " ===   Means  ===" << endl;
    for (auto const &tup:meanCellTuples) {
        cout << "The Mean is: ";
        cout << printPoint(std::get<0>(tup), keysServer);
        cout << endl;
        std::get<1>(tup).printSlice(keysServer);
    }
    cout << " === === === === ===" << endl;
    cout << endl;
    cout << " ===   Minimal Distances and Closest Means  ===" << endl;
    for (const auto &tuple:minDistanceTuples) {
        cout << " Point: ";
        cout << printPoint(std::get<0>(tuple), keysServer);
        cout << " Closest mean: ";
        cout << printPoint(std::get<1>(tuple), keysServer);
        long minDistance = keysServer.decryptNum(std::get<2>(tuple));
        printNameVal(minDistance);
    }
    cout << " === === === === ===" << endl;
    cout << endl;
    cout << " --- --- --- --- ---" << endl;

    printNameVal(keysServer.decryptNum(threshold));

    cout << " === Groups by Means === " << endl;
    for (auto const &[meanI, points] : std::get<0>(groups)) {
        printNameVal(meanI) << "\tClose Points: ";
        for (auto const &pair: points) {
            cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
            cout << printPoint(pair.first, keysServer);
        }
        cout << endl;
    }
    cout << " === === === === === " << endl << endl;

    cout << " === Closest Points === " << endl;
    for (auto const &pair : std::get<1>(groups)) {
        cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
        cout << printPoint(pair.first, keysServer);
    }
    cout << endl << " === === === === === " << endl << endl;

    cout << " === Farthest Points === " << endl;
    for (auto const &pair : std::get<2>(groups)) {
        cout << "-" << keysServer.decryptCtxt(pair.second) << "-";
        cout << printPoint(pair.first, keysServer);
    }
    cout << endl << " === === === === === " << endl << endl;

#ifdef alt
    // Choose two random n-bit integers
    long pa = NTL::RandomBits_long(BIT_SIZE);
    long pb = NTL::RandomBits_long(BIT_SIZE + 1);

    long pMax = std::max(pa, pb);
    long pMin = std::min(pa, pb);
    bool pMu = pa > pb;
    bool pNi = pa < pb;

    //    auto pv = helib::PtrVector<helib::Ctxt>(pa);

    // Encrypt the individual bits
    for (long i = 0; i <= BIT_SIZE; i++) {
        //todo in future remove this `if`,
        // currently serves as an example for init/enc of 2 different-sized numbers
        if (i < BIT_SIZE)
            encryptionKey.Encrypt(enca[i], NTL::ZZX((pa > > i) & 1));
        encryptionKey.Encrypt(encb[i], NTL::ZZX((pb > > i) & 1));
        //        if (helib_bootstrap) { // put them at a lower level
        //            if (i < BIT_SIZE)
        //                enca[i].bringToSet(dataServer.getCtxtPrimes(5));
        //            encb[i].bringToSet(dataServer.getCtxtPrimes(5));
        //        }
    }

//
//    printNameVal(n);
//    std::vector<long> v0(n);
//    for (long i = 0; i < n; i++)
//        v0[i] = i;
//    helib::PtxtArray ptxtArray(dataServer.getContext(), v0);
//    helib::Ctxt c0(publicKey);
//    ptxtArray.encrypt(c0);
//
//    compareTwoNumbers(mu, ni, helib::CtPtrs_VecCt(c0), helib::CtPtrs_VecCt(c0), false, nullptr);
#elif alt2
    // Each bit of the binary number is encoded into a single ciphertext. Thus
    // for a 16 bit binary number, we will represent this as an array of 16
    // unique ciphertexts.
    // i.e. b0 = [0] [0] [0] ... [0] [0] [0]        ciphertext for bit 0
    //      b1 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 1
    //      b2 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 2
    // These 3 ciphertexts represent the 3-bit binary number 110b = 6

    // Note: several numbers can be encoded across the slots of each ciphertext
    // which would result in several parallel slot-wise operations.
    // For simplicity we place the same data into each slot of each ciphertext,
    // printing out only the back of each vector.
    // NB: fifteenOrLess4Four max is 15 bits. Later in the code we pop the MSB.
    long bitSize = 16;
    long outSize = 2 * bitSize;
    long a_data = NTL::RandomBits_long(bitSize);
    long b_data = NTL::RandomBits_long(bitSize);
    long c_data = NTL::RandomBits_long(bitSize);

    // Use a scratch ciphertext to populate vectors.
    const helib::PubKey &public_key = encryptionKey;
//    const helib::PubKey &public_key = dataServer.getSecKey();
    helib::Ctxt scratch(public_key);
    std::vector<helib::Ctxt> encrypted_a(bitSize, scratch);
    std::vector<helib::Ctxt> encrypted_b(bitSize, scratch);
    std::vector<helib::Ctxt> encrypted_c(bitSize, scratch);
    // Encrypt the data in binary representation.
    for (long i = 0; i < bitSize; ++i) {
        std::vector<long> a_vec(ea.size());
        std::vector<long> b_vec(ea.size());
        std::vector<long> c_vec(ea.size());
        // Extract the i'th bit of a,b,c.
        for (auto &slot : a_vec)
            slot = (a_data >> i) & 1;
        for (auto &slot : b_vec)
            slot = (b_data >> i) & 1;
        for (auto &slot : c_vec)
            slot = (c_data >> i) & 1;
        ea.encrypt(encrypted_a[i], public_key, a_vec);
        ea.encrypt(encrypted_b[i], public_key, b_vec);
        ea.encrypt(encrypted_c[i], public_key, c_vec);

        std::vector<long> slotsMin, slotsMax, slotsMu, slotsNi;
        // comparison only
        compareTwoNumbers(mu, ni,
                          helib::CtPtrs_vectorCt(encrypted_a),
                          helib::CtPtrs_vectorCt(encrypted_a),
                          false);
        ea.decrypt(mu, encryptionKey, slotsMu);
        ea.decrypt(ni, encryptionKey, slotsNi);
//        for (auto s : slotsMu) printNameVal(s);
//        for (auto s : slotsNi) printNameVal(s);

    }
    // Although in general binary numbers are represented here as
    // std::vector<helib::Ctxt> the binaryArith APIs for HElib use the PtrVector
    // wrappers instead, e.g. helib::CtPtrs_vectorCt. These are nothing more than
    // thin wrapper classes to standardise access to different vector types, such
    // as NTL::Vec and std::vector. They do not take ownership of the underlying
    // object but merely provide access to it.
    //
    // helib::CtPtrMat_vectorCt is a wrapper for
    // std::vector<std::vector<helib::Ctxt> >, used for representing a list of
    // encrypted binary numbers.
#endif

    //-----------------------------------------------------------

    logger.log(printDuration(t0_main, "Main"));
    logger.print_log(log_trace);//, false);

    return 0;
}



