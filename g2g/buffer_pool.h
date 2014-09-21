#ifndef __THREAD_BUFFER_POOL_H__
#define __THREAD_BUFFER_POOL_H__

#include <mkl.h>
#include <stdexcept>
#include <iostream>
#include <cstring>

static const int ALIGNMENT = 64;

template<class T>
class ThreadBufferPool {
public:
    ThreadBufferPool(int buffers, int buffer_elements) : pieces(buffers) {
        aligned_buffer_size = buffer_elements * sizeof(T);
        aligned_buffer_size = (aligned_buffer_size + ALIGNMENT - (aligned_buffer_size % ALIGNMENT));
        pool = (T *) mkl_malloc(buffers * aligned_buffer_size, ALIGNMENT);
        total_size = pieces * aligned_buffer_size;
        requested = 0;
    }
    int buffer_size() const{
        return aligned_buffer_size;
    }
    ~ThreadBufferPool() {
        mkl_free(pool);
    }
    void reset() {
        requested = 0;
    }
    T * get_pool() {
        if(requested >= pieces){
            throw new std::runtime_error("Thread pool has received more requests than can be satisfied");
        }
        T * ret = (T *)((intptr_t)pool + requested * aligned_buffer_size);
        requested++;
        return ret;
    }
private:
    T * pool;
    int pieces, aligned_buffer_size, total_size, requested;
};

#endif
