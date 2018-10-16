cmake_minimum_required (VERSION 3.0)

project (test)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 -g")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")


add_executable(main non_dominating.cpp)


