

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

//    void printPoints(std::vector<Point> & points){ todo not static or move to aux
void printPoint(const Point &p, const KeysServer &keysServer) {
    cout << "( ";
    for (short dim = 0; dim < DIM - 1; ++dim)
        cout << keysServer.decryptNum(p[dim]) << ",";
    cout << keysServer.decryptNum(p[short(DIM - 1)]) << " ) ";
}

void printPoints(const std::vector<Point> &points, const KeysServer &keysServer) {
    cout << "   [ total of " << points.size() << " points ]   ";
    for (const Point &p:points) {
        cout << "( ";
        for (short dim = 0; dim < DIM - 1; ++dim)
            cout << keysServer.decryptNum(p[dim]) << ",";
        cout << keysServer.decryptNum(p[DIM - 1]) << " ) \t";
    }
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
            cout << "( ";
            for (short dim = 0; dim < DIM - 1; ++dim) {
                arr[dim] = keysServer.decryptNum(p[dim]);
                cout << arr[dim] << ",";
            }
            cout << keysServer.decryptNum(p[DIM - 1]) << " ) \t";
        }
    }
    cout << " \t\t[ total of " << cnt << " points are not empty, out of " << points.size() << " ]    ";
}

/*  @brief Generate random data.
 *  Returns a vector of clients with random points.
*      client [0] stays empty
*      client [1] has 1 point - {point0}
*      client [2] has 2 points - {point0, point1}
*      ...
*      client [n] has n points -  {point0, point1, ... , pointN}
*/  // todo move to aux ?
std::vector<Client> generateDataClients(const KeysServer &keysServer) {
    /*        //        Note that rand() is considered harmful, and is discouraged in C++14
            int uniquePointsNum = 3 + random() % NUMBER_OF_POINTS,
                    clientsNum = uniquePointsNum;
            printNameVal(NUMBER_OF_POINTS);
            printNameVal(uniquePointsNum);
            *//*        //  init coordinate arrays
                //        long arrs[uniquePointsNum][DIM];
                //        for (auto &arr: arrs)
                //            for (short dim = 0; dim < DIM; ++dim)
                //                arr[dim] = rand() % NUMBERS_RANGE;*//*
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
            }*/

    std::random_device rd;
//    std::mt19937 mt(rd());
    std::mt19937 mt;
    std::uniform_real_distribution<double> dist(0, NUMBERS_RANGE);
    long tempArr[DIM];
    std::vector<Client> clients(NUMBER_OF_CLIENTS, Client(keysServer));
    for (Client &client:clients) {
        for (int dim = 0; dim < DIM; ++dim) tempArr[dim] = dist(mt);
        client.encryptPoint(tempArr);
    }
    return clients;
}

void Slice::printSlice(const KeysServer &keysServer) const {
    cout << "For the Reps: ";
    printPoints(reps, keysServer);
    cout <<endl<< " These " << keysServer.decryptSize(included) << " Points will be included: " ;//<< endl;
    printNonEmptyPoints(includedPoints, keysServer);
    //    printPoints(includedPoints, keysServer);
    cout << "   ---     --- " << endl;
}
