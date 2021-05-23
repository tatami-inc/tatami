#include <fstream>
#include <deque>
#include <iostream>
#include "vector/standard_vector.hpp"
#include "data_frame/data_frame.hpp"
#include "matrix/dense_matrix.hpp"
#include "summarized_experiment/summarized_experiment.hpp"

using int32_deque = bioc::standard_vector<double, std::deque<int32_t> >;
using double_deque = bioc::standard_vector<double, std::deque<double> >;
using double_deque_dense = bioc::dense_matrix<double, std::deque<double> >;

template<typename T>
inline T lshift (char x, int shift) {
    return static_cast<T>(static_cast<unsigned char>(x)) << shift;
}

// Binary file reader for int32's.
std::deque<int32_t> read_binary_int32(const std::string& fname) {
    std::deque<int32_t> output;
    std::ifstream input(fname, std::ios::binary);

    std::vector<char> buffer(sizeof(int32_t));
    while (input.peek() != EOF) {
        input.read(buffer.data(), sizeof(int32_t));
        
        // Format is guaranteed to be big endian.
        uint32_t placeholder = lshift<uint32_t>(buffer[3], 0)
                             | lshift<uint32_t>(buffer[2], 8)
                             | lshift<uint32_t>(buffer[1], 16)
                             | lshift<uint32_t>(buffer[0], 24);

        output.push_back(*reinterpret_cast<int32_t*>(&placeholder));
    }

    return output;
}

// Binary reader for doubles.
std::deque<double> read_binary_double(const std::string& fname) {
    std::deque<double> output;
    std::ifstream input(fname, std::ios::binary);

    std::vector<char> buffer(sizeof(double));
    while (input.peek() != EOF) {
        input.read(buffer.data(), sizeof(double));

        // Format is guaranteed to be big endian.
        uint64_t placeholder = lshift<uint64_t>(buffer[7], 0)
                             | lshift<uint64_t>(buffer[6], 8)
                             | lshift<uint64_t>(buffer[5], 16)
                             | lshift<uint64_t>(buffer[4], 24)
                             | lshift<uint64_t>(buffer[3], 64)
                             | lshift<uint64_t>(buffer[2], 40)
                             | lshift<uint64_t>(buffer[1], 48)
                             | lshift<uint64_t>(buffer[0], 56);

        output.push_back(*reinterpret_cast<double*>(&placeholder));
    }

    return output;
}

//int main() {
//    auto stuff0 = read_binary_int32("colData/stuff");
//    std::shared_ptr<bioc::numeric_vector> stuff(new int32_deque(std::move(stuff0)));
//    std::shared_ptr<bioc::data_frame> coldata(new bioc::data_frame(stuff->length()));
//    coldata->set_column("stuff", stuff);
//
//    auto vals0 = read_binary_double("assays/counts/values");
//    std::shared_ptr<bioc::numeric_matrix> vals(new double_deque_dense(10, 26, std::move(vals0)));
//    std::cout << stuff->length() << "\t" << stuff->get(1) << std::endl;
//
//    std::shared_ptr<bioc::summarized_experiment> se (new bioc::summarized_experiment(vals->nrow(), vals->ncol()));
//    se->set_column_data(coldata);
//
//    std::vector<double> buffer(vals -> nrow());
//    std::cout << vals->nrow() << "\t" << vals->ncol() << "\t" << vals->get_column(2, buffer.data())[0] << std::endl;
//
//    return 0;
//}

class SummarizedExperiment {
public:
   SummarizedExperiment(int nr, int nc, std::string assay) : se(new bioc::summarized_experiment(nr, nc)) {
        auto vals0 = read_binary_double(assay);
        std::shared_ptr<bioc::numeric_matrix> vals(new double_deque_dense(nr, nc, std::move(vals0)));
        se->set_assay("WHEE", vals);
        return;
   }

   int nrow() const {
       return se->nrow();
   }

   int ncol() const {
       return se->ncol();
   }
private:
    std::shared_ptr<bioc::summarized_experiment> se;
};

#include <emscripten/bind.h>

EMSCRIPTEN_BINDINGS(my_class_example) {
    emscripten::class_<SummarizedExperiment>("SummarizedExperiment")
        .constructor<int, int, std::string>()
        .function("nrow", &SummarizedExperiment::nrow)
        .function("ncol", &SummarizedExperiment::ncol)
        ;
}
