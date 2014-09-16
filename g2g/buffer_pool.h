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
        aligned_buffer_elements = ((buffer_elements + ALIGNMENT - 1) / ALIGNMENT) * ALIGNMENT;
        aligned_buffer_size = aligned_buffer_elements * sizeof(T);
        pool = (T *) mkl_malloc(buffers * aligned_buffer_size, ALIGNMENT);
        total_size = pieces * aligned_buffer_size;
        requested = 0;
    }
    ThreadBufferPool(const ThreadBufferPool & other) {
        copy(other);
    }
    ThreadBufferPool & operator=(const ThreadBufferPool & other){
        if(this != &other) {
            mkl_free(pool); 
            pool = NULL;
            copy(other);
        }
        return *this;
    }
    ~ThreadBufferPool() {
        mkl_free(pool);
        pool = NULL;
    }
    void reset() {
        requested = 0;
    }
    T * get_pool() {
        if(requested >= pieces){
            throw new std::runtime_error("Thread pool has received more requests than can be satisfied");
        }
        T * ret = pool + requested * aligned_buffer_elements;
        requested++;
        return ret;
    }
private:
    void copy(const ThreadBufferPool & other) {
        total_size = other.total_size; pieces = other.pieces;
        aligned_buffer_elements = other.aligned_buffer_elements;
        aligned_buffer_size = other.aligned_buffer_size;
        requested = other.requested;

        pool = (T *) mkl_malloc(total_size, ALIGNMENT);
        memcpy(pool, other.pool, total_size);

    }

    T * pool;
    int pieces, aligned_buffer_elements, aligned_buffer_size, total_size, requested;
};

#endif
