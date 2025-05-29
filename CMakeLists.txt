cmake_minimum_required(VERSION 3.16)
project(KPC)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-O3 -g -march=native)

find_package(TBB REQUIRED)

add_executable(main main.cpp)
add_executable(R2 R2.cpp)
add_executable(R3 R3.cpp)
target_link_libraries(R3 TBB::tbb)
add_executable(DegenCol DegenCol.cpp)
add_executable(baseline baseline.cpp)

add_executable(base2 base2.cpp)
