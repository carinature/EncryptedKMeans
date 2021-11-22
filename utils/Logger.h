

#ifndef ENCRYPTEDKMEANS_LOGGER_H
#define ENCRYPTEDKMEANS_LOGGER_H

/**
 * @file Logger.h
 * */

#include <iostream>
#include <string>
#include <sstream>
#include <chrono>

enum LogLevel {
    log_default_level,
    log_trace,
    log_debug,
    log_info,
    log_warning,
    log_error,
    log_fatal,
};

/**
 * @class Logger
 * @brief Logger class handles logging (instead of stdout).
 * @param minLogLevel log all message from this level and up, discard all others.
 * @note The data will be logged to all levels below and up to msgLevel,
 *  and discarded from the those above
 *  (e.g msgLevel=info will log the data to dbg as well as info but not to error).
 * */
class Logger {
private:
    const LogLevel level = log_trace;
    std::ostringstream logs[log_fatal + 1];

    static std::string levelToString(LogLevel level);

public:
    const std::string name;

    explicit Logger(
            LogLevel minLogLevel = log_trace,
            std::string name = "General Logger");

    virtual ~Logger(); //note virtual? why

    /**
     * @brief Log messages to the wanted level and all those below.
     * @param msg The data to be loged .
     * @param msgLevel The maximum level to which the data is logged.
     * @note The data will be logged to all levels below and up to msgLevel,
     * and discarded from the those above
     * (e.g msgLevel=info will log the data to dbg as well as info but not to error).
     * */
    void log(const std::string &msg, LogLevel msgLevel = log_trace);

    /**
     * @brief Log messages to the wanted level and all those below.
     * @param msg The data to be loged .
     * @param msgLevel The maximum level to which the data is logged.
     * @note The data will be logged to all levels below and up to msgLevel, and
     *  discarded from the those above (e.g msgLevel=info will print the data from error,
     *  as well as info if all=true, but not from dbg).
     * */
    void print_log(LogLevel msgLevel = log_trace, bool all = false);

    friend std::ostream &operator<<(std::ostream &os, const Logger &logger);

    Logger &operator<<(const std::string &msg); //fixme change to work with ostream?

};


#endif //ENCRYPTEDKMEANS_LOGGER_H
