

#ifndef ENCKMEAN_PROPERTIES_H
#define ENCKMEAN_PROPERTIES_H


#include <helib/helib.h>

using DecryptedPoint = std::vector<long>;     // typedef std::vector<long> DecryptedPoint;  // same as `typedef/using`
using CBit = helib::Ctxt;  // typedef helib::Ctxt CBit;   // same as `typedef/using`
using EncNumber = NTL::Vec<CBit>;  // typedef NTL::Vec<helib::Ctxt> EncNumber;   // same as `typedef/using`
using CLOCK = std::chrono::high_resolution_clock;
static std::ofstream fcout("/home/karina/CLionProjects/EncryptedKMeans/fcout");  //for DBG

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
static const short NUMBER_OF_POINTS = jsonConfig["data_properties"]["number_of_points"];
static const short NUMBER_OF_CLIENTS = NUMBER_OF_POINTS;
//static const short NUMBER_OF_CLIENTS = jsonConfig["data_properties"]["number_of_clients"];
static const short DIM = jsonConfig["data_properties"]["DIM"];
static const short range_lim = jsonConfig["data_properties"]["range_lim"];
static const short low_limit = jsonConfig["data_properties"]["low_limit"];
static const short BIT_SIZE = jsonConfig["data_properties"]["bitSize"];
static const short OUT_SIZE = 2 * BIT_SIZE;
static const short NUMBERS_RANGE = pow(2, BIT_SIZE);
static const short N_Threads = jsonConfig["data_properties"]["N_Threads"];
static const short decimal_digits = jsonConfig["data_properties"]["decimal_digits"];
static const short conversion_factor = pow(10, decimal_digits);
static const double EPSILON = jsonConfig["data_properties"]["epsilon"];

[[maybe_unused]] static const std::string range_lim_comment = jsonConfig["data_properties"]["range_lim_comment"];
[[maybe_unused]] static const std::string conversion_factor_comment = jsonConfig["data_properties"]["conversion_factor_comment"];

/*
 * Flags
 * */
[[maybe_unused]] static const bool DBG = jsonConfig["flags"]["DBG"];
[[maybe_unused]] static const bool VERBOSE = jsonConfig["flags"]["VERBOSE"];
//// in productopn sould be `#define`d and not `statc const..`
//#define VERBOSE true
//#define DBG true //#define DBG false
[[maybe_unused]] static const bool helib_bootstrap = jsonConfig["helib_flags"]["helib_bootstrap"];

/*
 * Data-Files Names
 * */
static const std::string points_file = jsonConfig["files"]["points_file"];
static const std::string points_copy_file = jsonConfig["files"]["points_copy_file"];
static const std::string chosen_file = jsonConfig["files"]["chosen_file"];
static const std::string leftover_file = jsonConfig["files"]["leftover_file"];
static const std::string means_file = jsonConfig["files"]["means_file"];
static const std::string rands_file = jsonConfig["files"]["rands_file"];
static const std::string rands_bad_file = jsonConfig["files"]["rands_bad_file"];
static const std::string point_csv_file = jsonConfig["files"]["point_csv_file"];

//#define HELIB_DEBUG


#define CLIENT_CLASS_EXISTS true

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
