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

# Add nbspline library
add_library(nbspline src/nbspline.cpp)
target_include_directories (nbspline PUBLIC 
    ${PROJECT_SOURCE_DIR} 
    ${PROJECT_SOURCE_DIR}/include)

add_executable(functionality_1d src/test_functionality_1d.cpp)
add_test(functionality_1d functionality_1d)
target_link_libraries(functionality_1d 
    ${PYTHON_LIBRARIES}
    nbspline)

add_executable(functionality_3d src/test_functionality_3d.cpp)
target_link_libraries(functionality_3d 
    ${PYTHON_LIBRARIES}
    nbspline)