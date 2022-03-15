# CMake generated Testfile for 
# Source directory: /home/carina/CLionProjects/EncryptedKMeans/cmake-build-debug/_deps/json-src/test/cmake_import
# Build directory: /home/carina/CLionProjects/EncryptedKMeans/cmake-build-debug/_deps/json-build/test/cmake_import
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(cmake_import_configure "/snap/clion/180/bin/cmake/linux/bin/cmake" "-G" "Ninja" "-A" "" "-Dnlohmann_json_DIR=/home/carina/CLionProjects/EncryptedKMeans/cmake-build-debug/_deps/json-build" "/home/carina/CLionProjects/EncryptedKMeans/cmake-build-debug/_deps/json-src/test/cmake_import/project")
set_tests_properties(cmake_import_configure PROPERTIES  FIXTURES_SETUP "cmake_import" _BACKTRACE_TRIPLES "/home/carina/CLionProjects/EncryptedKMeans/cmake-build-debug/_deps/json-src/test/cmake_import/CMakeLists.txt;1;add_test;/home/carina/CLionProjects/EncryptedKMeans/cmake-build-debug/_deps/json-src/test/cmake_import/CMakeLists.txt;0;")
add_test(cmake_import_build "/snap/clion/180/bin/cmake/linux/bin/cmake" "--build" ".")
set_tests_properties(cmake_import_build PROPERTIES  FIXTURES_REQUIRED "cmake_import" _BACKTRACE_TRIPLES "/home/carina/CLionProjects/EncryptedKMeans/cmake-build-debug/_deps/json-src/test/cmake_import/CMakeLists.txt;8;add_test;/home/carina/CLionProjects/EncryptedKMeans/cmake-build-debug/_deps/json-src/test/cmake_import/CMakeLists.txt;0;")
