cmake_minimum_required(VERSION 2.8)
project(SubgraphMatching)

include_directories(
        ${PROJECT_SOURCE_DIR}/
)

set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -ffast-math -g -Wall -march=native -pthread")

add_subdirectory(graph)
add_subdirectory(matching)
add_subdirectory(utility)
