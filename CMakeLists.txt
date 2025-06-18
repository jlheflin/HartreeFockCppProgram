cmake_minimum_required(VERSION 3.6)
project(MyProgram)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
set(CMAKE_BUILD_TYPE Debug)
# Include headers
include_directories(include)

#Add source files
set(SOURCES
    src/main.cpp
    src/classes.cpp
    src/functions.cpp
)

add_executable(program ${SOURCES})
target_include_directories(program PRIVATE include)
