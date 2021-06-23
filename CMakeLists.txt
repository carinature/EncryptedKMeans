cmake_minimum_required(VERSION 3.19)
project(EncKMean)

## Use -std=c++17 as default.
set(CMAKE_CXX_STANDARD 17)

find_package(helib REQUIRED)

add_executable(EncKMean DataServer.cpp aux.cpp aux.h DataServer.h KeysServer.cpp KeysServer.h Client.cpp Client.h properties.h main.cpp)

target_link_libraries(EncKMean helib)

