#include "../base/typed_matrix.h"
#include <tiledb/tiledb>

namespace tatami {

template<typename T, typename IDX> 
class TileDBMatrix : public typed_matrix<T, IDX> {
public:
    TileDBMatrix(std::string f, std::string a) : fname(f), aname(a) {
        tiledb::Context ctx;
        tiledb::Array array(ctx, fname, TILEDB_READ);
        const auto& sch = array.schema();
        const auto& dom = sch.domain();

        auto row_dim = dom.dimension(0).domain();
        row_start = row_dim.first;
        nrows = row_dim.second - row_dim.first + 1;

        auto col_dim = dom.dimension(0).domain();
        col_start = col_dim.first;
        ncols = col_dim.second - col_dim.first + 1;

        return;
    }

    size_t nrow () const { return nrows; }

    size_t ncol () const { return ncols; }

    bool sparse() const { return is_sparse; }

    bool prefer_rows() const { return prefer_rows; }

public:
    struct tiledb_workspace : public workspace {
        ~tiledb_workspace() {}
        tiledb_workspace(std::string name, size_t n) : handle(ctx, name, TILEDB_READ), coord_buffer(n) {} 
        tiledb::Context ctx;
        tiledb::Array handle;
        std::vector<int> coord_buffer;
    };
    
    workspace_ptr new_workspace(bool row) {
        return std::shared_ptr<workspace>(new tiledB_workspace(fname, row ? nrow : ncol));
    }

public:
    const T* row(size_t r, T* buffer, size_t first, size_t last, workspace_ptr work=nullptr) const {
        if (work==nullptr) {
            work = new_workspace(true);
        }
        tiledb::Query query(work->ctx, work->array, TILEDB_READ);

        r -= row_start;
        std::vector<int> slice = { r, r, first - col_start, last - col_start };

        query.set_subarray(slice)
             .set_layout(TILEDB_ROW_MAJOR)
             .set_buffer(fname, buffer, last - first);
        query.submit();

        return buffer;
    }

    const T* column(size_t c, T* buffer, size_t first, size_t last, workspace_ptr work=nullptr) const {
        if (work==nullptr) {
            work = new_workspace(false);
        }
        tiledb::Query query(work->ctx, work->array, TILEDB_READ);

        c -= col_start;
        std::vector<int> slice = { first - row_start, last - row_start, c, c };

        query.set_subarray(slice)
             .set_layout(TILEDB_COL_MAJOR)
             .set_buffer(anem, buffer, last - first);
        query.submit();

        return buffer;
    }

    sparse_range<T, IDX> sparse_row(size_t r, T* vbuffer, IDX* ibuffer, size_t first, size_t last, workspace_ptr work=nullptr) const {
        if (work==nullptr) {
            work = new_workspace(true);
        }
        tiledb::Query query(work->ctx, work->array, TILEDB_READ);

        r -= row_start;
        std::vector<int> slice = { r, r, first - col_start, last - col_start };

        query.set_subarray(slice)
             .set_layout(TILEDB_ROW_MAJOR)
             .set_buffer(aname, vbuffer, last - first)
             .set_buffer("row", work->coord_buffer.data(), last - first)
             .set_buffer("col", ibuffer, last - first);
        query.submit();

        return sparse_range<T, IDX>(static_cast<size_t>(query.result_buffer_elements()[aname].second), vbuffer, ibuffer);
    }

    sparse_range<T, IDX> sparse_column(size_t c, T* buffer, size_t first, size_t last, workspace_ptr work=nullptr) const {
        if (work==nullptr) {
            work = new_workspace(false);
        }
        tiledb::Query query(work->ctx, work->array, TILEDB_READ);

        c -= col_start;
        std::vector<int> slice = { first - row_start, last - row_start, c, c };

        query.set_subarray(slice)
             .set_layout(TILEDB_COL_MAJOR)
             .set_buffer(aname, vbuffer, last - first)
             .set_buffer("row", ibuffer, last - first)
             .set_buffer("col", work->coord_buffer.data(), last - first);
        query.submit();

        return sparse_range<T, IDX>(static_cast<size_t>(query.result_buffer_elements()[aname].second), vbuffer, ibuffer);
    }

private:
    std::string fname, aname;
    size_t row_start, nrows;
    size_t col_start, ncols;
    bool is_sparse = true, prefer_rows = true;
};

}
