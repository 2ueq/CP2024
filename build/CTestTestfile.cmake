# CMake generated Testfile for 
# Source directory: /home/coups/CP2024
# Build directory: /home/coups/CP2024/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTests "UnitTests_catch2")
set_tests_properties(UnitTests PROPERTIES  _BACKTRACE_TRIPLES "/home/coups/CP2024/CMakeLists.txt;21;add_test;/home/coups/CP2024/CMakeLists.txt;0;")
subdirs("_deps/catch2-build")
