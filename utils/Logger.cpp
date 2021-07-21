

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
            cout << " \n---------- " << levelToString(LogLevel(l)) << " ---------- " << endl;
            if (log_error <= l) {
                std::cerr << this->logs[l].str();
                cout << this->logs[l].str();
            }
            else cout << this->logs[l].str();
        }
    else cout << endl << "~Logger() " << endl;
//#else //VERBOSE==flase
//#endif //VERBOSE

}

void Logger::log(LogLevel msgLevel, const std::string &msg) {
    for (int l = this->level; l <= msgLevel; ++l) {
        this->logs[l] << "- " << msg << endl;
        if (log_error <= msgLevel) std::cerr << msg << endl;
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
            cout << "  --- " << levelToString(LogLevel(l)) << " --- " << endl
                      << this->logs[l].str() << endl;
            if (log_error <= msgLevel) std::cerr << this->logs[l].str();

        }
    } else {
        cout << "  --- " << levelToString(msgLevel) << " --- " << endl
                  << this->logs[msgLevel].str() << endl;
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


