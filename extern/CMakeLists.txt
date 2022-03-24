cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  byteme 
  GIT_REPOSITORY https://github.com/clusterfork/buffin
  GIT_TAG 55730cfa3b01eb4d741e99b4069ea99ef2cf359f
)

FetchContent_MakeAvailable(byteme)
