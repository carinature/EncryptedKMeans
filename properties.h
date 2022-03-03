

#ifndef ENCKMEAN_PROPERTIES_H
#define ENCKMEAN_PROPERTIES_H

#include <fstream>

using CLOCK = std::chrono::high_resolution_clock;
//using CLOCK = std::chrono::time_point<std::chrono::system_clock>;

static std::ofstream
fcout("/home/karina/CLionProjects/EncryptedKMeans/fcout");  //for DBG

//using std::string;
using std::cout;
using std::endl;
using std::cerr;

#include <helib/helib.h>

using DecryptedPoint = std::vector<long>;
using CBit = helib::Ctxt;
using EncNumber = NTL::Vec<CBit>;
/*
 * Parse JSON config file
 * */
#include <nlohmann/json.hpp>

using json = nlohmann::json;
static std::fstream config_json_file("/home/karina/CLionProjects/EncryptedKMeans/config.json");
static json jsonConfig = json::parse(config_json_file);

/*
 * Global Constatns & Data Properties
 * */
static const short NUMBER_OF_THREADS = jsonConfig["data_properties"]["number_of_threads"];
static const short NUMBER_OF_POINTS = jsonConfig["data_properties"]["number_of_points"];
//static const short NUMBER_OF_CLIENTS = NUMBER_OF_POINTS / NUMBER_OF_THREADS;
static const short NUMBER_OF_CLIENTS = jsonConfig["data_properties"]["number_of_clients"];
static const short NUMBER_OF_ITERATIONS = jsonConfig["data_properties"]["number_of_iterations"];
static const short DIM = jsonConfig["data_properties"]["DIM"];
static const short BIT_SIZE = jsonConfig["data_properties"]["bitSize"];
static const short OUT_SIZE = 2 * BIT_SIZE;
static const short NUMBERS_RANGE = pow(2, BIT_SIZE) - 1; // pow(2, BIT_SIZE)-0
static const short DISTANCE_BIT_SIZE =
        std::log2(NUMBER_OF_POINTS * (NUMBERS_RANGE * NUMBERS_RANGE)); //fixme make sure
static const short CID_BIT_SIZE = NUMBER_OF_POINTS; //fixme make sure
static const short range_lim = jsonConfig["data_properties"]["range_lim"];
static const std::string range_lim_comment = jsonConfig["data_properties"]["range_lim_comment"];
static const short low_limit = jsonConfig["data_properties"]["low_limit"];
static const short N_Threads = jsonConfig["data_properties"]["N_Threads"];
static const short decimal_digits = jsonConfig["data_properties"]["decimal_digits"];
static const short CONVERSION_FACTOR = pow(10, decimal_digits);
static const std::string conversion_factor_comment = jsonConfig["data_properties"]["conversion_factor_comment"];
static const double EPSILON = jsonConfig["data_properties"]["epsilon"];


/*
 * Flags
 * */
[[maybe_unused]] static const bool DBG = jsonConfig["flags"]["DBG"];
//[[maybe_unused]] static const bool VERBOSE = jsonConfig["flags"]["VERBOSE"];
//// in productopn sould be `#define`d and not `static const..`
#define VERBOSE false
//#define DBG true //#define DBG false
[[maybe_unused]] static const bool helib_bootstrap = jsonConfig["helib_flags"]["helib_bootstrap"];

/*
 * Data-Files Names
 * */
static const std::string IO_DIR = jsonConfig["files"]["io_dir"]; // todo USE
static const std::string POINTS_FILE = jsonConfig["files"]["points_file"];
static const std::string POINTS_COPY_FILE = jsonConfig["files"]["points_copy_file"];
static const std::string RANDS_FILE = jsonConfig["files"]["rands_file"];
static const std::string MEANS_FILE = jsonConfig["files"]["means_file"];
static const std::string CHOSEN_FILE = jsonConfig["files"]["chosen_file"];
static const std::string LEFTOVER_FILE = jsonConfig["files"]["leftover_file"];
static const std::string rands_bad_file = jsonConfig["files"]["rands_bad_file"];
static const std::string point_csv_file = jsonConfig["files"]["point_csv_file"];

static const unsigned long PLAINTEXT_PRIME_MODULUS = 2;
//static const unsigned long PLAINTEXT_PRIME_MODULUS = 4999;

//#define HELIB_DEBUG


#define CLIENT_CLASS_EXISTS true
#define ENC_NUM_IS_VEC true

//namespace helib {
//
//    typedef PtrVector <Ctxt> CtPtrs;
//
//    // CtPtrs_VecCt(NTL::Vec<Ctxt>)
//    typedef PtrVector_VecT <Ctxt> CtPtrs_VecCt;
//    // CtPtrs_vectorCt(std::vector<Ctxt>)
//    typedef PtrVector_vectorT <Ctxt> CtPtrs_vectorCt;
//    // CtPtrs_VecPt(NTL::Vec<Ctxt*>)
//    typedef PtrVector_VecPt <Ctxt> CtPtrs_VecPt;
//    // CtPtrs_vectorPt(std::vector<Ctxt*>)
//    typedef PtrVector_vectorPt <Ctxt> CtPtrs_vectorPt;
//
//    // struct PtrVector_VecT;    // constructed as PtrVector_VecT(NTL::Vec<T>)
//    // struct PtrVector_VecPt;   // constructed as PtrVector_VecPt(NTL::Vec<T*>)
//    // struct PtrVector_vectorT; // constructed as PtrVector_vectorT(std::vector<T>)
//    // struct PtrVector_vectorPt;// constructed PtrVector_vectorPt(std::vector<T*>)
//
//    // A slice of CtPtrs
//    typedef PtrVector_slice <Ctxt> CtPtrs_slice;
//
//}



#endif //ENCKMEAN_PROPERTIES_H
