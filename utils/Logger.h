//
// Created by karina on 23/06/2021.
//

#ifndef ENCRYPTEDKMEANS_LOGGER_H
#define ENCRYPTEDKMEANS_LOGGER_H

#include <iostream>
#include <string>
#include <sstream>

/**
 * Logger Class
 * */
enum LogLevel {
    log_trace,
    log_debug,
    log_info,
    log_warning,
    log_error,
    log_fatal,
};

class Logger {
private:
    LogLevel level = log_trace;
    std::ostringstream logs[log_fatal + 1];

    static std::string levelToString(LogLevel level);

public:
    explicit Logger(LogLevel minLevel = log_trace);

    virtual ~Logger();

    void log(LogLevel msgLevel, const std::string &msg);

    void print_log(LogLevel msgLevel = log_trace, bool all = true);

    friend std::ostream &operator<<(std::ostream &os, const Logger &logger);
    Logger& operator<<(const std::string & msg); //fixme change to work with ostream?

};


#endif //ENCRYPTEDKMEANS_LOGGER_H
