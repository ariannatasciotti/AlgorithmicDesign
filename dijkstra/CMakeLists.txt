cmake_minimum_required (VERSION 2.6)
project (dijkstra)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O4")

FILE (GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.c"
	               "${PROJECT_SOURCE_DIR}/include/*.h")

include_directories("${PROJECT_SOURCE_DIR}/include")



add_executable (dijkstra_test ${SOURCES})