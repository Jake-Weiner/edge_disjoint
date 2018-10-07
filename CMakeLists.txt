cmake_minimum_required (VERSION 3.0)

project (Lagrangian_Relax)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 -g")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

include_directories("/home/jake/PhD/Edge_Disjoint/c++/boost_1_67_0/Build/include")
link_directories("/home/jake/PhD/Edge_Disjoint/c++/boost_1_67_0/Build/lib")

add_executable(main main.cpp)


