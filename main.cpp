//
// run a full protocol
//
#include <iostream>
#include <chrono>
#include <fstream>

using std::cout;
using std::endl;

#include "properties.h"
#include "utils/Logger.h"
#include "utils/aux.h"
//#include "src/Client.h"
#include "src/DataServer.h"

//helib
#include <helib/binaryArith.h>
#include <helib/binaryCompare.h>
#include <bitset>

static Logger loggerMain(log_debug, "loggerMain");

int main() {
    auto t0_main = CLOCK::now();
    Logger logger;
    logger.log("Starting Protocol", log_trace);
    //  Keys Server
    KeysServer keysServer;
    DataServer dataServer(keysServer);
    logger.log(printDuration(t0_main, "KeysServer Initialization"));

    const std::vector<Client> clients = generateDataClients(keysServer);

    const std::vector<Point> points = DataServer::retrievePoints(clients);
    cout << " --- Points  ---" << endl;
    printPoints(points, keysServer);
    cout << " --- --- --- --- ---" << endl;

    const std::vector<std::vector<Point>> randomPoints = dataServer.pickRandomPoints(points, 0,
                                                                                     keysServer.tinyRandomPoint());
    cout << " --- Random Points  ---" << endl;
    for (auto vec :randomPoints)         printPoints(vec, keysServer);
    cout << " --- --- --- --- ---" << endl;

    //  create points-comparing dict - for every 2 points (p1,p2) answers p1[dim]>p2[dim]
    const std::vector<
            std::unordered_map<
                    const Point,
                    std::unordered_map<
                            const Point,
                            helib::Ctxt> > >
            cmpDict = DataServer::createCmpDict(points, randomPoints, keysServer.tinyRandomPoint());

//    const std::vector<
//            std::tuple<
//                    Point,
//                    std::vector<Point>,
//                    std::vector<Ctxt>
//            >
//    > groups = DataServer::split(points, randomPoints, cmpDict);
    std::map<int, //DIM
            std::vector< //current slices for approp dimension
                    Cell
            >
    >
    groups = DataServer::splitIntoEpsNet(points, randomPoints, cmpDict, keysServer);
    for (Cell & cell: groups[0]) {
        cell.printCell(keysServer);
    }

/*
    std::map<
            const Point,
            std::vector<
                    std::tuple<
                            const Point,
                            std::vector<Point>,
                            std::vector<Ctxt>
                    >
            >
    > allCells;
    //    allCells.reserve(pow((1 / epsilon), 2));*/
    /*    printNameVal(groups.size());
        cout << " --- Random Points  ---" << endl;
        for (auto g: groups) printPoint(std::get<0>(g), dataServer);
        cout << " --- Groups  ---" << endl;*//*
    for (std::tuple groupTuple: groups) {
        Point randomPoint = std::get<0>(groupTuple);
        std::vector<Point> group = std::get<1>(groupTuple);
        std::vector<Ctxt> groupSize = std::get<2>(groupTuple);
        *//*        cout << "on the left of random point: ";
                printPoint(randomPoint, dataServer);
                cout << " these points will be included: \n";
                printPoints(group, dataServer);
                cout << "group size: " << dataServer.decryptSize(groupSize) << endl;*//*
        std::vector<
                std::tuple<
                        Point,
                        std::vector<Point>,
                        std::vector<Ctxt>
                >
        > cells = DataServer::split(group, 1, keysServer);
        *//*        printNameVal(cells.size());
                cout << " --- Random Points For Cells ---" << endl;
                for (auto g: cells)  printPoint(std::get<0>(g), dataServer);
                cout << " --- Cells  ---" << endl;*//*
        for (std::tuple cellTuple: cells) {
            Point randomPointCell = std::get<0>(cellTuple);
            std::vector<Point> cell = std::get<1>(cellTuple);
            std::vector<Ctxt> cellSize = std::get<2>(cellTuple);
            cout << "on the LEFT of random point: ";
            printPoint(randomPoint, keysServer);
            cout << "and UNDER the (cell) random point: ";
            printPoint(randomPointCell, keysServer);
            cout << "these points will be included: \n";
            //            printPoints(cell, dataServer);
            printNonEmptyPoints(cell, keysServer);
            cout << "cell size: " << keysServer.decryptSize(cellSize) << endl << endl;
        }
        allCells.emplace(randomPoint, cells);
    }
*/
    //    for (auto &tup:allCells) {
    //        for
    //    }

    cout << " --- --- --- --- ---" << endl;

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
            encryptionKey.Encrypt(enca[i], NTL::ZZX((pa >> i) & 1));
        encryptionKey.Encrypt(encb[i], NTL::ZZX((pb >> i) & 1));
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
//    cout << "odin" << endl;
    // Encrypt the data in binary representation.
    for (long i = 0; i < bitSize; ++i) {
//        cout << "dva" << endl;
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
//        cout << "tri" << endl;

        std::vector<long> slotsMin, slotsMax, slotsMu, slotsNi;
        // comparison only
        compareTwoNumbers(mu, ni,
                          helib::CtPtrs_vectorCt(encrypted_a),
                          helib::CtPtrs_vectorCt(encrypted_a),
                          false);
//        cout << "chetiri" << endl;
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
    // std::vector<std::vector<helib::Ctxt>>, used for representing a list of
    // encrypted binary numbers.
#endif

    //-----------------------------------------------------------

    printDuration(t0_main, "Main");
    logger.print_log(log_trace);//, false);
    //    clientLogger.print_log(log_error, false);

    return 0;
}



