cmake_minimum_required(VERSION 3.19)
project(EncryptedKMeans)

## Use -std=c++17 as default.
set(CMAKE_CXX_STANDARD 17)

find_package(helib REQUIRED)

add_executable(EncryptedKMeans
        main.cpp
        properties.h
        utils/aux.cpp
#        aux.h
        src/DataServer.cpp
#        DataServer.h
        src/KeysServer.cpp
#        KeysServer.h
        src/Client.cpp
#        Client.h
        utils/Logger.cpp utils/Logger.h)

target_link_libraries(EncKMean helib)

