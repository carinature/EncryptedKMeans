

#ifndef ENCKMEAN_DATASERVER_H
#define ENCKMEAN_DATASERVER_H


#include "Client.h"

//move to the cpp file
#include <algorithm> //for the random shuffle
#include <random>

static Logger loggerDataServer(log_debug, "loggerDataServer");

using CmpDict =
const std::vector<
        std::unordered_map<
                const Point,
                std::unordered_map<
                        const Point,
                        helib::Ctxt
                >
        >
>;

class DataServer {

    //    const helib::PubKey encryptionKey;
    //     const helib::EncryptedArray ea;
    //    const helib::Ctxt scratch;

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
            tinyRandomPoint(keysServer.tinyRandomPoint()) {
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
    std::vector<std::vector<Point> >
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
    CmpDict
    createCmpDict(
            const std::vector<Point> &allPoints,
            const std::vector<std::vector<Point> > &randomPoints
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
            const std::vector<std::vector<Point> > &randomPoints,
            const CmpDict &cmpDict,
            const KeysServer &keysServer // for dbg todo remove
    );

    // TODO candidate for multithreading
    /**
     * @brief calculate cell-means
     * @param slices a list of Cells (each cell is a list of encrypted points and encrypted size)
     * @param keysServer
     * @return a list of slices and their corresponding means
     * @returns std::vector<std::tuple<Point, Slice> >
     * */
    static
    std::vector<std::tuple<Point, Slice> >
    calculateSlicesMeans(
            const std::vector<Slice> &slices,
            const KeysServer &keysServer
    );

    /**
     * @brief collect mean point from epsNet*/
    static
    const std::vector<Point>
    collectMeans(
            const std::vector<std::tuple<Point, Slice> > &slices,
            const KeysServer &keysServer
    ) {
        std::vector<Point> means;
        means.reserve(slices.size());
        for (auto const &slice:slices) means.push_back(std::get<0>(slice));
        return means;
    }

    //  collect all minimal distances into one place
    //      for each point
    //          calculate minimal distance from point to all means:
    //              tuple<point, mean, distance> point.mininmal-distance-from(means)
    //              for each mean
    //                  calculate distance from point to mean:
    //                      tuple<point, distance> point.distance-from(mean)
    //  calculate avg distance
    //  collect for each mean the points closeset to it
    //  pick all points with distance bigger than avg
    //



    static
    //    void
    std::vector<std::tuple<Point, Point, EncryptedNum> >
    collectDistancesFromMean(
            const std::vector<std::tuple<Point, Slice> > &slices,
            const std::vector<Point> &points,
            const KeysServer &keysServer
    ) {
        //  Collect Means
        const std::vector<Point> &means = collectMeans(slices, keysServer);

        for (const Point &point: points) {
            for (const Point &mean: means) {
                // calculate distance from mean
            }
            // caclculate minimal distances from all distances from means
        }

        return std::vector<std::tuple<Point, Point, EncryptedNum> >();
    }


    static
    //    void
    std::vector<std::tuple<Point, Point, EncryptedNum> >
    collectDistancesFromMeans(
            const std::vector<std::tuple<Point, Slice> > &slices,
            const std::vector<Point> &points,
            const KeysServer &keysServer
    ) {
        //  Collect Means
        const std::vector<Point> &means = collectMeans(slices, keysServer);

        for (const Point &point: points) {
            for (const Point &mean: means) {
                // calculate distance from mean
            }
            // caclculate minimal distances from all distances from means
        }

        return std::vector<std::tuple<Point, Point, EncryptedNum> >();
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

