

#include "aux.h"
#include "src/KeysServer.h"
#include "src/Point.h"
#include "src/ClientDevice.h"

#include <sstream>      // std::stringstream
#include <random>


std::chrono::time_point<std::chrono::system_clock> NowTime() {
    return CLOCK::now();
}

std::string
printDuration(
        const std::chrono::time_point<std::chrono::system_clock> &t1,
        const std::string &funcName
) {
    auto t2 = CLOCK::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::string str =
            "\'" + funcName + "\' Finished in " + std::to_string(duration) + " milliseconds.\n";
    cout << str << endl;
    fcout << str << endl;
    return str;
}

std::vector<long>
decryptPoint(
        const Point &p,
        const KeysServer &keysServer
) {
    std::vector<long> pPoint(DIM);
    for (short dim = 0; dim < DIM; ++dim)
        pPoint[dim] = keysServer.decryptNum(p[dim]);
    return pPoint;
}

std::string
printPoint(
        const Point &p,
        const KeysServer &keysServer
) {
    //    cout << "( ";
    //    for (short dim = 0; dim < DIM; ++dim)
    //        cout << keysServer.decryptNum(p[dim]) << " ";
    ////    cout << ") ";
    //    cout << "id=" << p.id << " cid=" << keysServer.decryptNum(p.cid) << " ) ";
    std::stringstream tempstr;
    tempstr << "(";
    for (short dim = 0; dim < DIM - 1; ++dim)
        tempstr << keysServer.decryptNum(p[dim]) << ",";
    tempstr << keysServer.decryptNum(p[DIM - 1]) << "), ";
    return tempstr.str();
    //    cout << "(";
    //    for (short dim = 0; dim < DIM - 1; ++dim)
    //        cout << keysServer.decryptNum(p[dim]) << ",";
    //    cout << keysServer.decryptNum(p[DIM - 1]) << "), ";
    //    cout << "id=" << p.id << " cid=" << keysServer.decryptNum(p.cid) << " ) ";
}

std::string
printPoints(
        const std::vector<Point> &points,
        const KeysServer &keysServer
) {
    std::stringstream tempstr;
    tempstr << "   [ total of " << points.size() << " points ]   ";
    for (const Point &p:points)
        tempstr << printPoint(p, keysServer);
    return tempstr.str();
}

void printNonEmptyPoints(
        const std::vector<Point> &points,
        const KeysServer &keysServer
) {
    long arr[DIM], cnt = 0;
    for (const Point &p:points) {
        long sum = 0;
        for (short dim = 0; dim < DIM; ++dim) {
            arr[dim] = keysServer.decryptNum(p[dim]);
            sum += arr[dim];
        }
        if (sum) {
            ++cnt;
            cout << printPoint(p, keysServer);
        }
    }
    cout << " \t\t[ total of " << cnt << " points are not empty, out of " << points.size()
         << " ]    ";
}

void
decAndWriteToFile(
        const std::vector<Point> &points,
        const std::string &filename,
        const KeysServer &keysServer
) {
    std::vector<DecryptedPoint> decPoints;
    decPoints.reserve(points.size());
    for (const Point &p : points) decPoints.push_back(decryptPoint(p, keysServer));
    std::ofstream outputFileStream(filename);
    std::stringstream ss;
    long sum;
    for (const DecryptedPoint &p : decPoints) {
        sum = 0;
        ss.str(std::string());
        for (long coor : p) {
            sum += coor;
            double realCoor = double(coor) / CONVERSION_FACTOR;
            ss << realCoor << " ";
        }
        if (0 < sum) {
            outputFileStream << ss.rdbuf() << endl;
        }
    }
    //    outputFileStream.flush();
    outputFileStream.close();
}

/** @brief Generate random data.
    @returns a vector of clients with random points.
       client [0] stays empty
       client [1] has 1 point - {point0}
       client [2] has 2 points - {point0, point1}
       ...
       client [n] has n points -  {point0, point1, ... , pointN}
*/
std::vector<ClientDevice> generateDataClients(const KeysServer &keysServer) {
    //    std::uniform_real_distribution<double> dist(0, NUMBERS_RANGE);
    std::uniform_int_distribution<long> dist(1, NUMBERS_RANGE);
    long tempArr[DIM];
    std::vector<ClientDevice> clients(NUMBER_OF_CLIENTS, ClientDevice(keysServer));
    for (ClientDevice &client:clients) {
        for (int n = 0; n < NUMBER_OF_POINTS / NUMBER_OF_CLIENTS; ++n) {
            for (int dim = 0; dim < DIM; ++dim) tempArr[dim] = randomLongInRange(mt);
            client.encryptPoint(tempArr);
        }
    }
    /*
    //  another option
    for (int i = 1; i < clients.size(); ++i)
        for (int j = 0; j < i; ++j) {
            for (short dim = 0; dim < DIM; ++dim)
                arr[dim] = random() % NUMBERS_RANGE;
            clients[i].encryptPoint(arr);
         }
     */
    return clients;
}

void
Slice::printSlice(
        const KeysServer &keysServer
) const {
    cout << "For the Reps: ";
    cout << printPoints(reps, keysServer);
    cout << endl << " These " << keysServer.decryptSize(counter)
         << " Points will be included: ";//<< endl;
    printNonEmptyPoints(points, keysServer);
    //    cout << printPoints(points, keysServer);
    cout << "   ---     --- " << endl;
}

Slice
&Slice::addPoint(
        const Point &point,
        const Ctxt &isIncluded) {
    points.push_back(point);
    counter.push_back(isIncluded);
    pointTuples.emplace_back(point, isIncluded);
    return *this;
}

std::vector<Ctxt>
prefix(
        const std::vector<Ctxt> &v,
        long k) {
    std::vector<Ctxt> pref;
    pref.reserve(k);
    for (int i = 1; i <= k; ++i)
//        for (int i = k+1; i <= v.size(); ++i)
        pref.push_back(v[i-1]);
    return pref;
}

std::vector<Ctxt>
suffix(
        const std::vector<Ctxt> &v,
        long k) {
    std::vector<Ctxt> pref;
    pref.reserve(k);
    for (int i = k+1; i <= v.size(); ++i)
//        for (int i = 1; i <= k; ++i)
            pref.push_back(v[i-1]);
    return pref;
}

Ctxt
isEqual(
        const std::vector<Ctxt> &a,
        const std::vector<Ctxt> &b,
        long w
) {
    const helib::PubKey &public_key = a.front().getPubKey();
    Ctxt res(public_key);
    Ctxt temp(public_key);
    res.addConstant(1L);    //  res = 1

    for (int i = 1; i <= w; ++i) {
        temp = a[i-1];    //  a[i]
        temp += b[i-1];   //  a[i]+b[i]
        temp.addConstant(1L);   //  a[i]+b[i]+1
        res *= temp;
    }

    return res;
}

Ctxt
isGrt(
        const std::vector<Ctxt> &a,
        const std::vector<Ctxt> &b,
        long w,
        const helib::SecKey &secret_key
) {

    const helib::PubKey &public_key = a.front().getPubKey();

    //    cout << "AAA: ";
    helib::Ptxt<helib::BGV> plaintext_result(public_key);
    //    for (int bit = a.size() - 1; bit >= 0; --bit) {
    //        secret_key.Decrypt(plaintext_result, a[bit]);
    //        cout << plaintext_result[0] << endl;
    //    }
    //    cout<<endl;
    //    cout << "BBB: ";
    //    for (int bit = b.size() - 1; bit >= 0; --bit) {
    //        secret_key.Decrypt(plaintext_result, b[bit]);
    //        cout << plaintext_result[0] << endl;
    //    }
    //    cout<<endl;

    Ctxt res(public_key);
    Ctxt tempRes(public_key);
    Ctxt tempProduct(public_key);

    res.addConstant(1L);    //  res = 1
    res += b[w - 1];            //  res = 1+b[w]
    res *= a[w - 1];            //  res = ( 1+b[w] ) * a[w]
/*
    cout << " [ res = ( 1+b[w - 1] ) * a[w - 1] ] ";
    secret_key.Decrypt(plaintext_result, res);
    cout << plaintext_result[0] << endl;
*/

    for (int i = 1; i <= w - 1; ++i) {
//        cout << "-----------  " << i << "  ---------------\n";
        Ctxt temp(public_key);
        temp.addConstant(1L);   //  temp = 1
/*
        cout << " [ temp = 1 ] ";
        secret_key.Decrypt(plaintext_result, temp);
        cout << plaintext_result[0] << endl;
*/

        temp += b[i-1];           //  temp = 1+b[i]
/*
        cout << " [ temp = 1+b[i] ] ";
        secret_key.Decrypt(plaintext_result, temp);
        cout << plaintext_result[0] << endl;
*/

        temp *= a[i-1];           //  temp = ( 1+b[i] ) * a[i]
/*
        cout << " [ temp = ( 1+b[i] ) * a[i] ] ";
        secret_key.Decrypt(plaintext_result, temp);
        cout << plaintext_result[0] << endl;
*/

// fixme Adi said "suffix" in the article is works in a backwards way -
//  that the msb is in the lowest index and lsb in the highest
//  need to make sure
        temp *= isEqual(suffix(a, i), suffix(b, i), w - 1 - i);
/*
        cout << " [ temp *= isEqual(sfxA,sfxB,w-i) ] ";
        secret_key.Decrypt(plaintext_result, temp);
        cout << plaintext_result[0] << endl;
*/

        //fixme consider changing it to a sort of accumulation
        // since it's enough for an even number of temp-bits to be non zero to nullify
                res += temp;

//        tempProduct = temp;      // temp
///*
//        cout << " [ tempProduct = temp ] ";
//        secret_key.Decrypt(plaintext_result, tempProduct);
//        cout << plaintext_result[0] << endl;
//*/
//
//        tempProduct *= tempRes;    // tempRes * temp
///*
//        cout << " [ tempProduct = temp*temp ] ";
//        secret_key.Decrypt(plaintext_result, tempProduct);
//        cout << plaintext_result[0] << endl;
//*/
//
//        tempProduct.negate();       // - tempRes * temp
///*
//        cout << " [ tempProduct =  - tempRes * temp ] ";
//        secret_key.Decrypt(plaintext_result, tempProduct);
//        cout << plaintext_result[0] << endl;
//*/
//
//        tempProduct += temp;    // temp - tempRes * temp
///*
//        cout << " [ tempProduct = temp - tempRes * temp ] ";
//        secret_key.Decrypt(plaintext_result, tempProduct);
//        cout << plaintext_result[0] << endl;
//*/
//
//        tempProduct += tempRes;    // tempRes + temp - tempRes * temp
///*
//        cout << " [ tempProduct = tempRes + temp - tempRes * temp ] ";
//        secret_key.Decrypt(plaintext_result, tempProduct);
//        cout << plaintext_result[0] << endl;
//*/
//
////        tempRes += tempProduct; // tempRes = tempRes + temp - tempRes * temp
//        tempRes = tempProduct; // tempRes = tempRes + temp - tempRes * temp
///*
//        cout << " [ tempRes = tempRes + temp - tempRes * temp ] ";
//        secret_key.Decrypt(plaintext_result, tempRes);
//        cout << plaintext_result[0] << endl;
//*/
//
    }
///*
//    cout << " [ res before ] ";
//    secret_key.Decrypt(plaintext_result, res);
//    cout << plaintext_result[0] << endl;
//*/
//    tempProduct = tempRes;
//    tempProduct *= res;
//    tempProduct.negate();
//    tempProduct += tempRes;
////    tempProduct += res;
////    res = tempProduct;
//    res += tempProduct;
//
//
////    res += tempRes;
///*
//    cout << " [ res after ] ";
//    secret_key.Decrypt(plaintext_result, res);
//    cout << plaintext_result[0] << endl;
//*/

    return res;
}
