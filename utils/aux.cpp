

#include "aux.h"
//#include "properties.h"
#include "src/KeysServer.h"
#include "src/Point.h"

#include <sstream>      // std::stringstream

using std::cout;
using std::endl;

std::chrono::time_point<std::chrono::system_clock> NowTime() {
    return std::chrono::high_resolution_clock::now();
}

std::string printDuration(const std::chrono::time_point<std::chrono::system_clock> &t1, const std::string &funcName) {
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::string str = "\'" + funcName + "\' Finished in " + std::to_string(duration) + " seconds.\n";
    cout << str << endl;
    fcout << str << endl; //todo remove?
    return str;
}

//    void printPoints(std::vector<Point> & points){ todo not static or move to aux
void printPoint(const Point &p, KeysServer &keysServer) {
    cout << "( ";
    for (short dim = 0; dim < DIM - 1; ++dim)
        cout << keysServer.decryptNum(p[dim]) << ",";
    cout << keysServer.decryptNum(p[short(DIM) - 1]) << " ) " << endl;
}
void printPoints(const std::vector<Point> &points, KeysServer &keysServer) {
    for (const Point &p:points) {
        cout << "( ";
        for (short dim = 0; dim < DIM - 1; ++dim)
            cout << keysServer.decryptNum(p[dim]) << ",";
        cout << keysServer.decryptNum(p[DIM - 1]) << " ) \t" ;
    }
    cout << endl<< "   --- total of " << points.size() << " points ---    " << endl;
}



