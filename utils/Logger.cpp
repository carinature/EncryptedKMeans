

#include "Logger.h"
#include "../properties.h"

/**
 * Logger Class
 * */
Logger::Logger(LogLevel level) {
    this->level = level;
    auto timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    for (int l = level; l <= log_fatal; ++l) this->logs[l] << ctime(&timenow);
}

Logger::~Logger() {
//#if VERBOSE //todo consider for production version
    if (VERBOSE)
        for (int l = this->level; l <= log_fatal; ++l) {
            std::cout << " \n---------- " << levelToString(LogLevel(l)) << " ---------- " << std::endl;
            if (log_error <= l) {
                std::cerr << this->logs[l].str();
                std::cout << this->logs[l].str();
            }
            else std::cout << this->logs[l].str();
        }
    else std::cout << std::endl << "~Logger() " << std::endl;
//#else //VERBOSE==flase
//#endif //VERBOSE

}

void Logger::log(LogLevel msgLevel, const std::string &msg) {
    for (int l = this->level; l <= msgLevel; ++l) {
        this->logs[l] << "- " << msg << std::endl;
        if (log_error <= msgLevel) std::cerr << msg << std::endl;
    }
}

std::ostream &operator<<(std::ostream &os, const Logger &logger) {
    os << logger.logs;
    return os;
}

Logger &Logger::operator<<(const std::string &msg) {
    this->log(log_info, msg);
    return *this;
}


void Logger::print_log(LogLevel msgLevel, bool all) {
    if (all) {
        for (int l = std::max(this->level, msgLevel); l <= log_fatal; ++l) {
            std::cout << "  --- " << levelToString(LogLevel(l)) << " --- " << std::endl
                      << this->logs[l].str() << std::endl;
            if (log_error <= msgLevel) std::cerr << this->logs[l].str();

        }
    } else {
        std::cout << "  --- " << levelToString(msgLevel) << " --- " << std::endl
                  << this->logs[msgLevel].str() << std::endl;
        if (log_error <= msgLevel) std::cerr << this->logs[msgLevel].str();
    }
}

std::string Logger::levelToString(LogLevel level) {
    switch (level) {
        case log_trace:
            return "log_trace";
        case log_debug:
            return "log_debug";
        case log_info:
            return "log_info";
        case log_warning:
            return "log_warning";
        case log_error:
            return "log_error";
        default: //case log_fatal:
            return "log_fatal";
    }
}


