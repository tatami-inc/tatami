#ifndef WRITE_MATRIX_MARKET_H
#define WRITE_MATRIX_MARKET_H

#include <fstream>
#include <string>

typedef std::vector<int> IntVec;

template<class Stream, class U, class V, class W>
void write_matrix_market(Stream& stream, size_t nr, size_t nc, const U& vals, const V& rows, const W& cols) {
    stream << "%%MatrixMarket matrix coordinate integer general\n";
    stream << nr << " " << nc << " " << vals.size();

    for (size_t i = 0; i < vals.size(); ++i) {
        stream << "\n" << rows[i] + 1 << " " << cols[i] + 1 << " " << vals[i];
    }
    stream << "\n";

    return;
}

template<class U, class V, class W>
void write_matrix_market(std::string path, size_t nr, size_t nc, const U& vals, const V& rows, const W& cols) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("failed to open output file");
    }
    write_matrix_market(out, nr, nc, vals, rows, cols);
    out.close();
    return;
}

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
