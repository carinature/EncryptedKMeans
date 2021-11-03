

#include "aux.h"
//#include "properties.h"
#include "src/KeysServer.h"
#include "src/Point.h"
#include "src/Client.h"

#include <sstream>      // std::stringstream
#include <random>

using std::cout;
using std::endl;

std::chrono::time_point<std::chrono::system_clock> NowTime() {
    return CLOCK::now();
}

std::string printDuration(const std::chrono::time_point<std::chrono::system_clock> &t1,
                          const std::string &funcName) {
    auto t2 = CLOCK::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::string str =
            "\'" + funcName + "\' Finished in " + std::to_string(duration) + " seconds.\n";
    cout << str << endl;
    fcout << str << endl; //todo remove?
    return str;
}

std::vector<long> decryptPoint(const Point &p, const KeysServer &keysServer) {
    std::vector<long> pPoint(DIM);
    for (short dim = 0; dim < DIM; ++dim)
        pPoint[dim] = keysServer.decryptNum(p[dim]);
    return pPoint;
}

//    void printPoints(std::vector<Point> & points){ todo not static or move to aux
void printPoint(const Point &p, const KeysServer &keysServer) {
    cout << "( ";
    for (short dim = 0; dim < DIM ; ++dim)
        cout << keysServer.decryptNum(p[dim]) << ",";
    cout << "id="<<p.id<< ",cid="<<keysServer.decryptNum(p.cid)<<" ) ";
}

void printPoints(const std::vector<Point> &points, const KeysServer &keysServer) {
    cout << "   [ total of " << points.size() << " points ]   ";
    for (const Point &p:points)         printPoint(p, keysServer);

}

void printNonEmptyPoints(const std::vector<Point> &points, const KeysServer &keysServer) {
    long arr[DIM], cnt = 0;
    for (const Point &p:points) {
        long sum = 0;
        for (short dim = 0; dim < DIM - 1; ++dim) {
            arr[dim] = keysServer.decryptNum(p[dim]);
            sum += arr[dim];
        }
        if (sum) {
            ++cnt;
            printPoint(p, keysServer);
        }
    }
    cout << " \t\t[ total of " << cnt << " points are not empty, out of " << points.size()
         << " ]    ";
}

/** @brief Generate random data.
    @returns a vector of clients with random points.
       client [0] stays empty
       client [1] has 1 point - {point0}
       client [2] has 2 points - {point0, point1}
       ...
       client [n] has n points -  {point0, point1, ... , pointN}
*/
std::vector<Client> generateDataClients(const KeysServer &keysServer) {
    std::random_device rd;
    //    std::mt19937 mt(rd());
    std::mt19937 mt;
    //    std::uniform_real_distribution<double> dist(0, NUMBERS_RANGE);
    std::uniform_int_distribution<long> dist(1, NUMBERS_RANGE);
    long tempArr[DIM];
    std::vector<Client> clients(NUMBER_OF_CLIENTS, Client(keysServer));
    for (Client &client:clients) {
        for (int dim = 0; dim < DIM; ++dim) tempArr[dim] = dist(mt);
        client.encryptPoint(tempArr);
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

void Slice::printSlice(const KeysServer &keysServer) const {
    cout << "For the Reps: ";
    printPoints(reps, keysServer);
    cout << endl << " These " << keysServer.decryptSize(counter)
         << " Points will be included: ";//<< endl;
    printNonEmptyPoints(points, keysServer);
    //    printPoints(points, keysServer);
    cout << "   ---     --- " << endl;
}

Slice &Slice::addPoint(const Point &point, const Ctxt &isIncluded) {
    points.push_back(point);
    counter.push_back(isIncluded);
    //        std::tuple<Point, CBit> pointTuple;
    pointTuples.emplace_back(point, isIncluded);
    return *this;
}
