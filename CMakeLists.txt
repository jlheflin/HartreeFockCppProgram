cmake_minimum_required(VERSION 3.5...4.0)
project(MyProgram)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
set(CMAKE_BUILD_TYPE Debug)

include(GNUInstallDirs)

set(CMAKE_INSTALL_LIBDIR "lib")
set(CMAKE_INSTALL_INCLUDEDIR "include")

find_package(Eigen3 CONFIG REQUIRED)
if(Eigen3_FOUND)
    message(STATUS "Eigen3 found!")
endif()
find_package(Libint2 CONFIG REQUIRED)
if(Libint2_FOUND)
    message(STATUS "Libint2 found!")
endif()

file(GLOB_RECURSE src_files CONFIGURE_DEPENDS src/*.cpp)
file(GLOB_RECURSE inc_files CONFIGURE_DEPENDS include/*.hpp)

add_executable(program ${src_files} ${inc_files})
target_include_directories(program PRIVATE include)
target_link_libraries(program PRIVATE
    Eigen3::Eigen
    Libint2::cxx
)

install(TARGETS program
    RUNTIME DESTINATION bin
)
