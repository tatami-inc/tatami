cmake_minimum_required(VERSION 3.24)

project(tatami
    VERSION 4.1.0
    DESCRIPTION "A C++ interface for assorted matrix representations"
    LANGUAGES CXX)

add_library(tatami INTERFACE)
add_library(tatami::tatami ALIAS tatami)

target_compile_features(tatami INTERFACE cxx_std_17)

include(GNUInstallDirs)
target_include_directories(tatami INTERFACE 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/tatami_tatami>"
)

# Dependencies
option(TATAMI_FETCH_EXTERN "Automatically fetch tatami's external dependencies." ON)
if(TATAMI_FETCH_EXTERN)
    add_subdirectory(extern)
else()
    find_package(ltla_subpar 0.5.0 CONFIG REQUIRED)
    find_package(ltla_sanisizer 0.2.0 CONFIG REQUIRED)
endif()

target_link_libraries(tatami INTERFACE ltla::subpar ltla::sanisizer)

# Building the test-related machinery, if we are compiling this library directly.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(TATAMI_TESTS "Build tatami's test suite." ON)
else()
    option(TATAMI_TESTS "Build tatami's test suite." OFF)
endif()
if(TATAMI_TESTS)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif() 
endif()

# Installing for find_package.
include(CMakePackageConfigHelpers)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tatami_tatami)

install(TARGETS tatami 
    EXPORT tatamiTargets)

install(EXPORT tatamiTargets
    FILE tatami_tatamiTargets.cmake
    NAMESPACE tatami::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatamiConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/tatami_tatamiConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatamiConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatamiConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami)
