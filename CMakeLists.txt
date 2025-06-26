cmake_minimum_required(VERSION 3.5...4.0)
project(MyProgram)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
set(CMAKE_BUILD_TYPE Debug)
option(BUILD_SHARED_LIBS "Enable building shared libraries" ON)

include(FetchContent)
include(GNUInstallDirs)

set(CMAKE_INSTALL_LIBDIR "lib")
set(CMAKE_INSTALL_INCLUDEDIR "include")
set(FETCHCONTENT_FULLY_DISCONNECTED OFF)

# include(cmake/fetch_parallelzone.cmake)
# include(cmake/fetch_tensorwrapper.cmake)

find_package(Eigen3 CONFIG REQUIRED)
# find_package(nlohmann_json CONFIG REQUIRED)
# find_package(cereal CONFIG REQUIRED)
set(SOURCES
    src/main.cpp
    src/classes.cpp
    src/functions.cpp
)

add_executable(program ${SOURCES})
target_include_directories(program PRIVATE include)
target_link_libraries(program PRIVATE
    Eigen3::Eigen
    # nlohmann_json::nlohmann_json
    # cereal::cereal parallelzone
)
install(TARGETS program
    RUNTIME DESTINATION bin
)
