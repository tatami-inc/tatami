#ifndef TATAMI_TEST_TEMP_FILE_PATH_HPP
#define TATAMI_TEST_TEMP_FILE_PATH_HPP

#if __has_include(<filesystem>)

#include <filesystem>

namespace tatami_test {

inline auto temp_file_path(const std::string& base) {
    auto path = std::filesystem::temp_directory_path();
    path.append(base);
    return path;
}

inline void remove_file_path(const std::string& path) {
    if (std::filesystem::exists(path)) {
        std::filesystem::remove_all(path);
    }
}

}

#else

#include <experimental/filesystem>

namespace tatami_test {

inline auto temp_file_path(const std::string& base) {
    auto path = std::experimental::filesystem::temp_directory_path();
    path.append(base);
    return path;
}

inline void remove_file_path(const std::string& path) {
    if (std::experimental::filesystem::exists(path)) {
        std::experimental::filesystem::remove_all(path);
    }
}

}

#endif

#endif
