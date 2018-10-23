cmake_minimum_required (VERSION 2.8.11)
project (rtinaweekend)

if(UNIX AND NOT APPLE)
    set(CMAKE_CXX_FLAGS "-std=c++0x")
    add_definitions(-DLINUX)
endif ()

add_subdirectory(src)

