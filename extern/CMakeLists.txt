cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  buffin
  GIT_REPOSITORY https://github.com/clusterfork/buffin
  GIT_TAG 0d46932ee74c4687bf0070a4dc81b91c65ac9be7
)

FetchContent_MakeAvailable(buffin)
