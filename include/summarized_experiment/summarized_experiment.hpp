#ifndef SUMMARIZED_EXPERIMENT_H
#define SUMMARIZED_EXPERIMENT_H

#include <memory>
#include "../data_frame/data_frame.hpp"
#include "../matrix/matrix.hpp"
#include <unordered_map>

namespace bioc {

class summarized_experiment {
public:
    summarized_experiment(size_t nr, size_t nc) : row_data(new data_frame(nr)), column_data(new data_frame(nc)) {} 

    // Standard getters and setters.
    size_t ncol() const {
        return row_data->nrow();
    }

    size_t nrow() const {
        return column_data->nrow();
    }

    // Getters and setters for the assay data.
    std::shared_ptr<matrix> get_assay(const std::string& name) {
        return assays[name];
    }

    void set_assay(const std::string& name, std::shared_ptr<matrix> value) {
        assays[name] = value;
        return;
    }

    void remove_assay(const std::string& name) {
        assays.erase(name);
        return;
    }

    void get_assay_names(std::vector<std::string>& names) const {
        names.clear();
        names.reserve(assays.size());
        for (const auto& x : assays) {
            names.push_back(x.first);
        }
        return;
    }

    bool has_assay(const std::string& name) const {
        return assays.find(name) != assays.end();
    }
    
    // Getters and setters for the row/column data.
    std::shared_ptr<data_frame> get_column_data() const {
        return column_data;
    }

    void set_column_data(std::shared_ptr<data_frame> value) {
        column_data = value;
        return;
    }

    std::shared_ptr<data_frame> get_row_data() const {
        return row_data;
    }

    void set_row_data(std::shared_ptr<data_frame> value) {
        row_data = value;
        return;
    }

    // Specialist methods for the names.
    bool get_column_names(std::vector<std::string>& names) const { 
        return column_data->get_row_names(names);
    }

    bool get_row_names(std::vector<std::string>& names) const { 
        return row_data->get_row_names(names);
    }

    void set_column_names(const std::vector<std::string>& names) { 
        column_data->set_row_names(names);
        return;
    }

    void set_row_names(const std::vector<std::string>& names) { 
        row_data->set_row_names(names);
        return;
    }

    void remove_column_names() { 
        column_data->remove_row_names();
        return;
    }

    void remove_row_names() {
        row_data->remove_row_names();
        return;
    }

    bool has_column_name(const std::string& name) const{ 
        return column_data->has_row_name(name);
    }

    bool has_row_name(const std::string& name) const{ 
        return row_data->has_row_name(name);
    }

private:
    std::unordered_map<std::string, std::shared_ptr<matrix> > assays;
    std::shared_ptr<data_frame> row_data;
    std::shared_ptr<data_frame> column_data;
};

}

#endif
