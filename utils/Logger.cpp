
#include "Logger.h"

#include <utility>
#include "properties.h"

Logger::Logger(LogLevel minLogLevel, std::string name) :
        name(std::move(name)),
        level(minLogLevel) {
    auto timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    for (int l = minLogLevel; l <= log_fatal; ++l) this->logs[l] << ctime(&timeNow);

    auto t = time(nullptr);
    auto tm = *localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y_%m_%d_%H_%M_%S");
    timestamp = oss.str();
    logfilename = IO_DIR + timestamp + "_logfile";
    logfile.open(logfilename);
}

Logger::~Logger() {

#if VERBOSE //todo consider for production version
    if (VERBOSE)
        for (int l = this->level; l <= log_fatal; ++l) {
            cout << " \n---------- " << levelToString(LogLevel(l)) << " ---------- " << endl;
            if (log_error <= l) cerr << this->logs[l].str();
            cout << this->logs[l].str();
        }
#else //VERBOSE==flase
    cout << endl << "~Logger() for "<<name << endl; //fixme uncomment
#endif //VERBOSE


//fixme file not created and/or is empty
    std::ofstream logfiletemp(logfilename);
    logfiletemp.open(logfilename);

    logfiletemp << logs[log_trace].str();
    logfiletemp << logs[log_trace].rdbuf();

//    cout << "++++++++++++++++++++++==="
//         << logfiletemp.rdbuf()//fixme this part is never printed
//         << "===++++++++++++++++++++++";

    logfiletemp.close();
}

void Logger::log(const std::string &msg, LogLevel msgLevel) {
    for (int l = this->level; l <= msgLevel; ++l) {
        this->logs[l] << "- " << msg << endl;
        if (log_error <= msgLevel) std::cerr << msg << endl;
    }
}

void Logger::print_log(LogLevel msgLevel, bool all) {
    cout << " --- print_log ---\n";
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

std::ostream &operator<<(std::ostream &os, const Logger &logger) {
    os << logger.logs;
    return os;
}

Logger &Logger::operator<<(const std::string &msg) {
    this->log(msg, log_info);
    return *this;
}



