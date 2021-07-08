

#ifndef ENCKMEAN_PROPERTIES_H
#define ENCKMEAN_PROPERTIES_H


#include <helib/helib.h>

//  todo relocate to aux or main ...?
using DecryptedPoint = std::vector<long>;     // typedef std::vector<long> DecryptedPoint;  // same as `typedef/using`
using Bit = helib::Ctxt;  // typedef helib::Ctxt Bit;   // same as `typedef/using`
using EncNumber = NTL::Vec<Bit>;  // typedef NTL::Vec<helib::Ctxt> EncNumber;   // same as `typedef/using`
static std::ofstream fcout("fcout");  //for DBG

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
static const short number_of_points = jsonConfig["data_properties"]["number_of_points"];
static const short dim = jsonConfig["data_properties"]["dim"];
static const short range_lim = jsonConfig["data_properties"]["range_lim"];
static const short low_limit = jsonConfig["data_properties"]["low_limit"];
static const short bit_size = jsonConfig["data_properties"]["bit_size"];
static const short N_Threads = jsonConfig["data_properties"]["N_Threads"];
static const short decimal_digits = jsonConfig["data_properties"]["decimal_digits"];
static const short conversion_factor = pow(10, decimal_digits);
static const double epsilon = jsonConfig["data_properties"]["epsilon"];

[[maybe_unused]] static const std::string range_lim_comment = jsonConfig["data_properties"]["range_lim_comment"];
[[maybe_unused]] static const std::string conversion_factor_comment = jsonConfig["data_properties"]["conversion_factor_comment"];

/*
 * Flags
 * */
[[maybe_unused]] static const bool DBG = jsonConfig["flags"]["DBG"];
[[maybe_unused]] static const bool VERBOSE = jsonConfig["flags"]["VERBOSE"];


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


#define DBG true //#define DBG false

//string files{};





#endif //ENCKMEAN_PROPERTIES_H
