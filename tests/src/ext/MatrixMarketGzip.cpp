#include <gtest/gtest.h>

#include "tatami/ext/MatrixMarket.hpp"

/* Copied verbatim from https://www.zlib.net/zpipe.c */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "zlib.h"

#define CHUNK 16384

int def(FILE *source, FILE *dest, int level) {
    int ret, flush;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    /* allocate deflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 | 16, 8, Z_DEFAULT_STRATEGY);
    if (ret != Z_OK)
        return ret;

    /* compress until end of file */
    do {
        strm.avail_in = fread(in, 1, CHUNK, source);
        if (ferror(source)) {
            (void)deflateEnd(&strm);
            return Z_ERRNO;
        }
        flush = feof(source) ? Z_FINISH : Z_NO_FLUSH;
        strm.next_in = in;

        /* run deflate() on input until output buffer not full, finish
           compression if all of source has been read in */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = deflate(&strm, flush);    /* no bad return value */
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            have = CHUNK - strm.avail_out;
            if (fwrite(out, 1, have, dest) != have || ferror(dest)) {
                (void)deflateEnd(&strm);
                return Z_ERRNO;
            }
        } while (strm.avail_out == 0);
        assert(strm.avail_in == 0);     /* all input will be used */

        /* done when last data in file processed */
    } while (flush != Z_FINISH);
    assert(ret == Z_STREAM_END);        /* stream will be complete */

    /* clean up and return */
    (void)deflateEnd(&strm);
    return Z_OK;
}

#include <limits>
#include <random>
#include "write_matrix_market.h"

TEST(MatrixMarketGzip, GzipTest) {
    std::vector<int> vals, rows, cols;
    size_t NR = 2000, NC = 1000;
    std::mt19937_64 rng(1234567890);
    std::uniform_real_distribution<> unif(0.0, 1.0);

    for (size_t i = 0; i < NC; ++i) {
        for (size_t j = 0; j < NR; ++j) {
            if (unif(rng) < 0.05) {
                rows.push_back(j);
                cols.push_back(i);
                vals.push_back(unif(rng) * 100);
            }
        }
    }

    // Writing to file first.
    std::string path = "testing-MM-gzip.mtx";
    write_matrix_market(path, NR, NC, vals, rows, cols);

    FILE * ihandle = std::fopen(path.c_str(), "rb");
    std::string gzname = path + ".gz";
    FILE * ohandle = std::fopen(gzname.c_str(), "wb");
    def(ihandle, ohandle, Z_DEFAULT_COMPRESSION);
    std::fclose(ihandle);
    std::fclose(ohandle);

    auto loaded = tatami::MatrixMarket::load_layered_sparse_matrix(path.c_str());
    auto loaded_gz = tatami::MatrixMarket::load_layered_sparse_matrix_gzip(gzname.c_str());

    auto indptrs = tatami::compress_sparse_triplets<false>(NR, NC, vals, rows, cols);
    auto ref = std::shared_ptr<tatami::NumericMatrix>(new tatami::CompressedSparseColumnMatrix<double, int,
                                                                                               decltype(vals),
                                                                                               decltype(rows),
                                                                                               decltype(indptrs)
                                                                                              >(NR, NC, std::move(vals), std::move(rows), std::move(indptrs)));

    EXPECT_EQ(ref->nrow(), loaded.matrix->nrow());
    EXPECT_EQ(ref->nrow(), loaded_gz.matrix->nrow());
    EXPECT_EQ(ref->ncol(), loaded.matrix->ncol());
    EXPECT_EQ(ref->ncol(), loaded_gz.matrix->ncol());
    EXPECT_TRUE(loaded.matrix->sparse());
    EXPECT_TRUE(loaded_gz.matrix->sparse());
    EXPECT_FALSE(loaded.matrix->prefer_rows());
    EXPECT_FALSE(loaded_gz.matrix->prefer_rows());

    for (size_t i = 0; i < NC; ++i) {
        int adjusted = loaded.permutation[i];
        EXPECT_EQ(adjusted, loaded_gz.permutation[i]);

        auto stuff = loaded.matrix->column(adjusted);
        EXPECT_EQ(stuff, loaded_gz.matrix->column(adjusted));
        EXPECT_EQ(stuff, ref->column(i));
    }
}
