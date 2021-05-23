#ifndef DATA_FRAME_H
#define DATA_FRAME_H

#include <memory>
#include <unordered_map>
#include "../vector/vector.hpp"

namespace bioc {

class data_frame {
public:
    data_frame(size_t nr) : nrows(nr) {} 

    // Standard getters and setters.
    size_t ncol() const {
        return columns.size();
    }

    size_t nrow() const {
        return nrows;
    }

    std::shared_ptr<vector> get_column(const std::string& name) {
        return columns[name];
    }

    void set_column(const std::string& name, const std::shared_ptr<vector>& col) {
        columns[name] = col;
        return;
    }

    void remove_column(const std::string& name) {
        columns.erase(name);
        return;
    }

    bool has_column(const std::string& name) {
        return columns.find(name) != columns.end();
    }

    // Specialist methods for row names.
    bool get_row_names (std::vector<std::string>& rn) const {
        rn.clear();
        if (is_named) {
            rn.reserve(nrows);
            for (const auto& x : row_names) {
                rn[x.second] = x.first;
            }
        }
        return is_named;
    }

    bool has_row_name (const std::string& name) const {
        return is_named && row_names.find(name) != row_names.end();
    }

    void set_row_names(const std::vector<std::string>& rn) {
        if (rn.size()!=nrows) {
            throw std::runtime_error("row names must be of length equal to the number of rows");
        }
        row_names.clear();

        size_t counter = 0;
        is_named = true;

        for (const auto& x : rn) {
            if (row_names.count(x) > 0) {
                throw std::runtime_error("row names must be unique");
            }
            row_names[x] = counter;
            ++counter;
        }
        return;
    }

    void remove_row_names() {
        row_names.clear();
        is_named = false;
        return;
    }

//    // Subsetting.
//    std::shared_ptr<data_frame> subset_columns(const std::vector<std::string>& cols, bool clone=false) {
//        std::shared_ptr<data_frame> output(this->clone());
//        for (const auto& x : cols) {
//            output->remove_column(x);
//        }
//        return ouput;
//    }
//
//    std::shared_ptr<data_frame> subset_rows(const std::vector<int>& rows, bool clone=false) {
//        std::shared_ptr<data_frame> output(this->clone());
//        for (const auto& x : cols) {
//            output->remove_column(x);
//        }
//        return ouput;
//    }
//
//    std::shared_ptr<data_frame> subset_rows(const std::vector<std::string>& rn, bool clone=false) {
//
//    }
//
//    std::shared_ptr<data_frame> subset(const std::vector<std::string>& rn, bool clone=false) {
//
//    }
//
//    std::shared_ptr<data_frame> subset_rows(const std::vector<std::string>& rn, bool clone=false) {
//
//    }
//
//    // Cloning.

protected:
    std::unordered_map<std::string, std::shared_ptr<vector> > columns;
    std::unordered_map<std::string, size_t> row_names; 
    size_t nrows;
    bool is_named=false;
};

}
#endif
