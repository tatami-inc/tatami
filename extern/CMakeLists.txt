cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  byteme 
  GIT_REPOSITORY https://github.com/LTLA/byteme
  GIT_TAG 1deeedd74cc531bd1f0838f319be04c6c41dfe32
)

FetchContent_MakeAvailable(byteme)
