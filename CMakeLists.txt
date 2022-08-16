cmake_minimum_required(VERSION 2.8.3)
set(CMAKE_CXX_STANDARD 14)
project(libnbspline)

# set default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-deprecated-declarations")

set(BUILD_SHARED_LIBS ON)

file(GLOB INCLUDE_FILES "include/*.h")

find_package(Eigen3 REQUIRED)
find_package(PythonLibs REQUIRED)
include_directories(
    ${EIGEN3_INCLUDE_DIRS} 
    ${PYTHON_INCLUDE_DIRS}
    matplotlib-cpp
    include)

include(CTest)

add_executable(functionality src/test_functionality.cpp)
add_test(functionality functionality)
target_link_libraries(functionality ${PYTHON_LIBRARIES})