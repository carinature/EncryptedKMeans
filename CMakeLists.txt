cmake_minimum_required(VERSION 3.19)
project(EncryptedKMeans)

## Use -std=c++17 as default.
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_SOURCE_DIR "/home/karina/CLionProjects/EncryptedKMeans")
# Includes Directory
set(BASEPATH "${CMAKE_SOURCE_DIR}")
include_directories("${BASEPATH}")

# FetchContent can be used to automatically download the repository as a dependency at configure time.
include(FetchContent)

# JSON Parser
FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.7.3)
FetchContent_GetProperties(json)
if (NOT json_POPULATED)
    FetchContent_Populate(json)
    add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()

# HElib
find_package(helib REQUIRED)

add_executable(EncryptedKMeans
        main.cpp
        properties.h
        utils/aux.cpp
        utils/Logger.cpp
        src/DataServer.cpp
        src/KeysServer.cpp
        src/Client.cpp
#        src/Point.cpp src/Point.h
        )
target_link_libraries(EncryptedKMeans PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(EncryptedKMeans PUBLIC helib)

# Testing
add_executable(Tests
        properties.h
        utils/aux.cpp
        utils/Logger.cpp
        src/DataServer.cpp
        src/KeysServer.cpp
        src/Client.cpp
#        tests
        tests/TestClient.cpp
        tests/test.cpp
        tests/TestDataServer.cpp
#        src/Point.cpp src/Point.h
        )
target_link_libraries(Tests PUBLIC helib)
target_link_libraries(Tests PRIVATE nlohmann_json::nlohmann_json)

