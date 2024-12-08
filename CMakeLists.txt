cmake_minimum_required(VERSION 3.16)
project(NonObtuseTriangulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CGAL REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(opt_triangulation src/main.cpp src/TriangulationOptimizer.cpp src/InputParser.cpp src/OutputGenerator.cpp)

target_link_libraries(opt_triangulation CGAL::CGAL)
