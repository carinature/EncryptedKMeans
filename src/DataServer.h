
#ifndef ENCKMEAN_DATASERVER_H
#define ENCKMEAN_DATASERVER_H

#include <map>

#include "ClientDevice.h"
#include "tests/TestDataServer.h"

/**
 * CmpDict is a list (dictionary) of "comparisons";
 * it contains a list for each dimension.
 * in each list, for every pair of {randomPoint, dataPoint}
 * there is the result of (helib's) compare(randomPoint[dim], dataPoint[dim]) in
 * */
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

/**
 * DataServer is the server picking the core-set points returned to the customer
 * */
class DataServer {

    friend TestDataServer;
protected:
    //    const helib::PubKey &public_key;// = encryptionKey;
    const KeysServer &keysServer;
    const Point tinyRandomPoint;

private:
    /**
     * @var CmpDict contains the result of current iteration's compares.
     *  The compare are between iteration's data points and random representatives.
     *  CmpDict will be calculated for each iteration.
     *  @note comment for `using CmpDict`
     * */
    //    CmpDict cmpDict;
    std::vector<
            std::unordered_map<
                    const Point,
                    std::unordered_map<
                            const Point,
                            helib::Ctxt> > > cmpDict;
    std::mutex cmpDictLock;

    /**
     * @var The list of data points for the current iteration
     * */
    std::vector<Point> retrievedPoints;
    std::mutex retrievedPointsLock;

    /**
     * @var The list of data points randomly picked as representatives for the current iteration
     * */
    std::vector<std::vector<Point> > randomPointsList;

public:
    /**
     * Constructor for \class{ClientDevice},
     * @param keysServer binds to the \class{KeysServer} responsible for the distributing the appropriate key
     * @brief simulates a mini-protocol between the keys server and the Data Server.
     * ks and ds will create a scratch that will allow to encrypt the results:
     * result bit in case of compare, and result num in case of add/multiplication.
     * */
    explicit DataServer(const KeysServer &keysServer);

    void clearForNextIteration();

    /**
     * @brief A simulated retrieval of data from clients.
     * @param clients - a list of clients (chosen by the CA, to share a similar public key).
     * @returns a list of all the points.
    * @return std::vector<Point>
     * * */
    std::vector<Point>
    retrievePoints(
            const std::vector<ClientDevice> &clients
    );

    void
    retrievePoints_Thread(
            const std::vector<ClientDevice> &clients
    );

    std::vector<Point>
    retrievePoints_WithThreads(
            const std::vector<ClientDevice> &clients,
            short numOfThreads = NUMBER_OF_THREADS
    );


    /**
     * @brief request data from clients and concentrate into one list
     * @param points - all the points from all the clients in the data set
     * @param m - number of random representatives for each slice
     * @returns a list of #DIM lists - each containing m^d randomly chosen points
     * @return std::vector<Point>
     * */
    const std::vector<std::vector<Point> > &
    pickRandomPoints(
            const std::vector<Point> &points,
            int m = 1 / EPSILON
    );

    /**
     * @brief create a comparison dict:
     *   for each 2 points a,b returns the answer to a[dim]>b[dim]
     *   | where a and b are from 2 different groups (one is sub-group of the other),
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
    );

    void
    createCmpDict_Dim_Thread(
            short dim
    );

    CmpDict &
    createCmpDict_WithThreads(
            const std::vector<Point> &allPoints,
            const std::vector<std::vector<Point> > &randomPoints,
            int numOfThreads = NUMBER_OF_THREADS
    );

    /**
     * @brief Split into (1/eps) groups - each group is between 2 representative points.
     * @param points - a list of unordered points
     * @param randomPoints - a list of unordered points
     * @param dim - the coordinate by which the comparison is made. in other words - in what dimension the split is made.
     * @returns a list of pairs/tuples of a representative point and a list of the points in its Group (slice/slice).
     * */
    std::map<int, //DIM
            std::vector< //current slices corresponding to dimensions
                    Slice
            >
    >
    splitIntoEpsNet(
            const std::vector<Point> &points,
            const std::vector<std::vector<Point> > &randomPoints,
            const CmpDict &cmpDict
    );

    void
    splitIntoEpsNet_R_Thread(
            const Slice &baseSlice,
            const Point &R,
            int dim);

    std::map<int, //DIM
            std::vector<Slice> // slices corresponding to dimensions
    > splitIntoEpsNet_WithThreads(
            const std::vector<Point> &points,
            const std::vector<std::vector<Point> > &randomPoints,
            const CmpDict &cmpDict
    );

    /**
     * @brief calculate cell-means
     * @param slices a list of Cells (each cell is a list of encrypted points and encrypted size)
     * @param keysServer
     * @return a list of slices and their corresponding means
     * @returns std::vector<std::tuple<Point, Slice> >
     * */
    std::vector<std::tuple<Point, Slice>>
    calculateSlicesMeans(
            const std::vector<Slice> &slices
    );

    void
    calculateSliceMean_Slice_Thread(
            const Slice &slice
    ) const;

    std::vector<std::tuple<Point, Slice> >
    calculateSlicesMeans_WithThreads(
            const std::vector<Slice> &slices
    );

    /**
     * @brief collect mean point from epsNet
     * */
    std::vector<Point>
    collectMeans(
            const std::vector<std::tuple<Point, Slice>> &slices
    );

    //  collect all minimal distances into one place:
    //      for each point
    //          calculate minimal distance from point to all means:
    //              tuple<point, mean, distance> point.minimal-distance-from(means)
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
    std::vector<std::tuple<Point, Point, EncryptedNum>>
    collectMinimalDistancesAndClosestPoints(
            const std::vector<Point> &points,
                                            const std::vector<Point> &means);

    std::vector<std::tuple<Point, Point, EncryptedNum>>
    collectMinimalDistancesAndClosestPoints_WithThreads(
            const std::vector<Point> &points,
            const std::vector<Point> &means
    );

    /**
     * @brief calculate the average which will be used as a threshold for picking "closest" points
     * @param minDistanceTuples a list of tuples {point, closest_mean_point, distance}
     * @return Encrypted average of all minimal distances
     * @returns EncryptedNum
     * */
    EncryptedNum
    calculateThreshold(
            const std::vector<std::tuple<Point, Point, EncryptedNum>> &minDistanceTuples,
            int iterationNumber);


    /**
     * @brief pick points that are closest to current iteration means;
     *  meaning their distance from the closest mean is smaller than the 'threshold'
     *  (also calculated for this iteration).
     *  @note for each Point also includes a bit signifying if the point is included returns
     * @param minDistanceTuples a list of tuples {point, closest_mean_point, distance}
     *  @note collect for each mean the points closest to it
     * @param means the list of current iteration means
     * @param current iteration threshold for picking groups of closest points
     * @return a tuple of three lists:
     *  1. list of groups of points who are closest to the means - to be passed to the 1-mean alg
     *  2. list of all closest point
     *  3. list of "leftover" (furthest) points - the data set for the next iteration
     * @returns EncryptedNum
     * */
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

    std::tuple<
            std::unordered_map<long, std::vector<std::pair<Point, CBit>>>,
            std::vector<std::pair<Point, CBit>>
    >
    choosePointsByDistance_WithThreads_slower(
            const std::vector<std::tuple<Point, Point, EncryptedNum>> &minDistanceTuples,
            std::vector<Point> &means,
            EncryptedNum &threshold);

    std::tuple<
            std::unordered_map<long, std::vector<std::pair<Point, CBit>>>,
            std::vector<std::pair<Point, CBit>>
    >
    choosePointsByDistance_WithThreads(
            const std::vector<std::tuple<Point, Point, EncryptedNum>> &minDistanceTuples,
            std::vector<Point> &means,
            EncryptedNum &threshold);

};


#endif //ENCKMEAN_DATASERVER_H

