#ifndef HDF5_CUSTOM_LOCK_H
#define HDF5_CUSTOM_LOCK_H

extern std::mutex hdf5_lock;

template<class Function>
static void hdf5_serialize(Function f) {
    std::lock_guard<std::mutex> thing(hdf5_lock);
    f();
}

#define TATAMI_HDF5_PARALLEL_LOCK hdf5_serialize

#endif
