#ifndef TATAMI_TEST_TEMP_FILE_PATH_HPP
#define TATAMI_TEST_TEMP_FILE_PATH_HPP

#if __has_include(<filesystem>)

#include <filesystem>

namespace tatami_test {

inline auto temp_file_path(std::string base) {
    auto path = std::filesystem::temp_directory_path();
    path.append(base);
    return path;
}

}

#else

#include <experimental/filesystem>

namespace tatami_test {

inline auto temp_file_path(std::string base) {
    auto path = std::experimental::filesystem::temp_directory_path();
    path.append(base);
    return path;
}

}

#endif

#endif
