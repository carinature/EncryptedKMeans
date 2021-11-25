
#ifndef ENCKMEAN_DATASERVER_H
#define ENCKMEAN_DATASERVER_H

#include "Client.h"


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

protected:
    //    const helib::PubKey &public_key;// = encryptionKey;
    const KeysServer &keysServer;
    const Point tinyRandomPoint;
public:
    /**
     * Constructor for \class{Client},
     * @param keysServer binds to the \class{KeysServer} responsible for the distributing the appropriate key
     * @brief simulates a mini-protocol between the keys server and the Data Server.
     * ks and ds will create a scratch that will allow to encrypt the results:
     * result bit in case of compare, and result num in case of add/multiplication.
     * */
    explicit DataServer(const KeysServer &keysServer) :
            keysServer(keysServer),
            tinyRandomPoint(keysServer.tinyRandomPoint())
//            ,
//            retrievedPoints(NUMBER_OF_POINTS)
//            ,
//            randomPointsList(DIM) //todo change name to randomPoints
//            ,
//            cmpDict(DIM, {{}})
    {
        //        dataServerLogger.log("DataServer()");
        cout << "DataServer()" << endl;
        cmpDict.resize(DIM);
        randomPointsList.resize(DIM);
        retrievedPoints.reserve(NUMBER_OF_POINTS);

    }

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

    std::vector<Point> retrievedPoints;
    std::mutex retrievedPointsLock;

    void
    retrievePoints_Thread(const std::vector<Client> &clients);

    void
    retrievePoints_WithThreads(
            const std::vector<Client> &clients
            , short numOfThreads = NUMBER_OF_THREADS
            );


    std::vector<std::vector<Point> > randomPointsList;

    /**
     * @brief request data from clients and conentrate into one list
     * @param points - all the points from all the clients in the data set
     * @param m - number of random representatives for each slice
     * @returns a list of #DIM lists - each containing m^d randomly chosen points
     * @return std::vector<Point>
     * */
//    std::vector<std::vector<Point> >
    const std::vector<std::vector<Point> > &
    pickRandomPoints(
            const std::vector<Point> &points,
            int m = 1 / EPSILON   //  -1
            //            const Point & tinyRandomPoint
    ) ;

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
    CmpDict
    createCmpDict(
            const std::vector<Point> &allPoints,
            const std::vector<std::vector<Point> > &randomPoints
            //            const Point & tinyRandomPoint
    );

//    CmpDict cmpDict;
std::vector<
        std::unordered_map<
                const Point,
                std::unordered_map<
                        const Point,
                        helib::Ctxt> > > cmpDict;
std::mutex cmpDictLock;

    void
    createCmpDict_WithThreads(short numOfThreads);

    // TODO candidate for multithreading
    /**
     * @brief Split into (1/eps) groups - each group is between 2 representative points.
     * @param points - a list of unordered points
     * @param randomPoints - a list of unordered points
     * @param dim - the index of coor by which the comparison is made. in other words - in what dimension the split is made.
     * @returns a list of pairs/tuples of a representative point and a list of the points in its Group (slice/slice).
     * */
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
    std::vector<Point>
    collectMeans(
            const std::vector<std::tuple<Point, Slice> > &slices,
            const KeysServer &keysServer
    );


    //  collect all minimal distances into one place:
    //      for each point
    //          calculate minimal distance from point to all means:
    //              tuple<point, mean, distance> point.mininmal-distance-from(means)
    //              for each mean
    //                  calculate distance from point to mean:
    //                      tuple<point, distance> point.distance-from(mean)
    /**
     * @brief collect all pairs of points and their closest mean point 
     * (not necessarily from the same slice) 
     * by calculating minimal distance from point to one of the means
     * @param points - all original points
     * @param means - all the means from the epsNet
     * @return tuples of [point, closest mean, minimal distance]
     * @returns
     * */
    // TODO candidate for multithreading
    static
    std::vector<std::tuple<Point, Point, EncryptedNum> >
    collectMinimalDistancesAndClosestPoints(
            const std::vector<Point> &points,
            //            const std::reference_wrapper<std::vector<Point>> &means,
            const std::vector<Point> &means,
            const KeysServer &keysServer
    );


    //  calculate avg distance
    /**
     * @brief calculates the avarage which will be used as a threshold for picking "closest" points
     * @param minDistanceTuples
     * @return Encrypted avrage
     * @returns EncryptedNum
     * */
    static
    EncryptedNum
    calculateThreshold(
            const std::vector<std::tuple<Point, Point, EncryptedNum> > &minDistanceTuples,
            const KeysServer &keysServer,
            int iterationNumber = 0
    );


    //  collect for each mean the points closest to it
    //  for each Point also includes a bit signifying if the point is included returns
    // TODO candidate for multithreading
    static
    std::tuple<
            std::unordered_map<long, std::vector<std::pair<Point, CBit> > >,
            std::vector<std::pair<Point, CBit> >,
            std::vector<std::pair<Point, CBit> >
    >
    choosePointsByDistance(
            const std::vector<std::tuple<Point, Point, EncryptedNum> > &minDistanceTuples,
            std::vector<Point> means,
            EncryptedNum &threshold
    );

    void createCmpDict_Dim_Thread(short dim);

    std::map<int, std::vector<Slice>> splitIntoEpsNet_WithThreads(const KeysServer &keysServer);
};


#endif //ENCKMEAN_DATASERVER_H

