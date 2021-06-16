cmake_minimum_required(VERSION 3.14)

project(tatami
    VERSION 1.0.0
    DESCRIPTION "A C++ interface for assorted matrix representations"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(tatami INTERFACE)

target_include_directories(tatami INTERFACE include/)

add_subdirectory(tests)
