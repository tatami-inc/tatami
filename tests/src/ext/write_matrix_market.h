#ifndef WRITE_MATRIX_MARKET_H
#define WRITE_MATRIX_MARKET_H

#include <fstream>
#include <string>

template<class U, class V, class W>
void write_matrix_market(std::string path, size_t nr, size_t nc, const U& vals, const V& rows, const W& cols) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("failed to open output file");
    }

    out << "%%MatrixMarket matrix coordinate integer general\n";
    out << nr << " " << nc << " " << vals.size();

    for (size_t i = 0; i < vals.size(); ++i) {
        out << "\n" << rows[i] + 1 << " " << cols[i] + 1 << " " << vals[i];
    }
    out << std::endl;
    out.close();
    return;
}

#if __has_include(<filesystem>)

#include <filesystem>

inline auto temp_file_path(std::string base) {
    auto path = std::filesystem::temp_directory_path();
    path += base;
    return path;
}

#else

#include <experimental/filesystem>

inline auto temp_file_path(std::string base) {
    auto path = std::filesystem::temp_directory_path();
    path += base;
    return path;
}

#endif

#endif
