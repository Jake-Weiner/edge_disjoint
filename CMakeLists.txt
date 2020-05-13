cmake_minimum_required (VERSION 3.0)

project (Lagrangian_Relax)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 -g -O0 -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories("/home/jake/PhD/Edge_Disjoint/c++/boost_1_67_0/Build/include")
include_directories("/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/LAPSO_Req")
include_directories("/home/jake/PhD/Coding_Libraries/cplex/include")
include_directories("/home/jake/PhD/Coding_Libraries/concert/include")

link_directories("/home/jake/PhD/Coding_Libraries/cplex/lib/x86-64_linux/static_pic")
link_directories("/home/jake/PhD/Coding_Libraries/concert/lib/x86-64_linux/static_pic")
link_directories("/home/jake/PhD/Edge_Disjoint/c++/boost_1_67_0/Build/lib")

FILE(GLOB_RECURSE SOURCES RELATIVE MyCppSources src/*.cpp)

add_executable(main_test main_test.cpp ED.cpp LaPSO.cpp CpuTimer.cpp anyoption.cpp VolVolume.cpp prep_mip.cpp djikstra.cpp)

target_link_libraries(main_test concert ilocplex cplex m pthread dl)
