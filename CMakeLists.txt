cmake_minimum_required (VERSION 3.0)

project (Lagrangian_Relax)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 -g")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

include_directories("/home/jake/PhD/Edge_Disjoint/c++/boost_1_67_0/Build/include")
include_directories("/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/LAPSO_Req")
link_directories("/home/jake/PhD/Edge_Disjoint/c++/boost_1_67_0/Build/lib")

add_executable(main_test main_test.cpp ED.cpp LaPSO.cpp CpuTimer.cpp anyoption.cpp)


