cmake_minimum_required(VERSION 3.6)
project(HTmag)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}") #  -fopenmp -march=native -O3
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3")

file(GLOB_RECURSE HEADERS include/*.hpp)
file(GLOB_RECURSE SOURCE source/*.cpp)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}/release)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_CURRENT_SOURCE_DIR}/debug)

add_executable(HTmag demo.cpp ${HEADERS} ${SOURCE})

include(../matplotlib-cpp/matplotlibcpp.cmake)

target_link_libraries(HTmag PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../gsl_wrapper/libgsl_wrapper-${CMAKE_BUILD_TYPE}.a)
target_link_libraries(HTmag PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../eos/libeos-${CMAKE_BUILD_TYPE}.a)

find_package(Boost COMPONENTS filesystem system iostreams REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(HTmag PRIVATE ${Boost_LIBRARIES})

find_path( CLN_INCLUDE_DIR "cln/cln.h" DOC "Directory of cln.h" )
find_library( CLN_LIB "cln" DOC "Path to libcln" )

target_link_libraries(HTmag PRIVATE gsl gslcblas ${CLN_LIB})