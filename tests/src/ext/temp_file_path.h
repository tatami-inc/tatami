#ifndef TEMP_FILE_PATH
#define TEMP_FILE_PATH

#if __has_include(<filesystem>)

#include <filesystem>

inline auto temp_file_path(std::string base) {
    auto path = std::filesystem::temp_directory_path();
    path.append(base);
    return path;
}

#else

#include <experimental/filesystem>

inline auto temp_file_path(std::string base) {
    auto path = std::experimental::filesystem::temp_directory_path();
    path.append(base);
    return path;
}

#endif

#endif
