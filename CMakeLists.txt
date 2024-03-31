# test framework repository
include(FetchContent)
FetchContent_Declare(catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.1.0
)
FetchContent_MakeAvailable(Catch2)

# main (irrelevant now; temporary settings here)
cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 14)
project(Test_of_Catch2 CXX)
add_executable(a.out main.cpp)
target_compile_options(a.out PRIVATE -Wall -Wextra)

# unit tests
enable_testing()
add_executable(UnitTests_catch2 UnitTest.cpp)
target_compile_options(UnitTests_catch2 PRIVATE -Wall -Wextra)
target_link_libraries(UnitTests_catch2 Catch2::Catch2WithMain)
add_test(UnitTest UnitTests_catch2)
