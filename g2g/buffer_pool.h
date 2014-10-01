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
    ThreadBufferPool(int buffers, int buffer_elements) : _pieces(buffers) {
        _buffer_size = buffer_elements * sizeof(T);
        _pool = (T *) mkl_malloc(buffers * _buffer_size, ALIGNMENT);
        _total_size = _pieces * _buffer_size;
        _skip = _pos = 0;
    }
    int buffer_size() const{
        return _buffer_size;
    }
    ~ThreadBufferPool() {
        mkl_free(_pool);
    }
    void reset() {
        _pos = 0;
    }
    T * pool_start() const {
        return this->_pool;
    }
    void set_skip(int skip) {
        _skip = skip;
    }
    inline T * get_pool() {
        if(_pos >= _total_size){
            throw new std::runtime_error("Thread pool has received more requests than can be satisfied");
        }
        if(_skip == 0) {
            throw new std::runtime_error("Skip was not set");
        }
        T * ret = (T *)((unsigned long)_pool + _pos * sizeof(T));
        _pos += _skip;
        return ret;
    }
private:
    T * _pool;
    int _pieces, _buffer_size, _total_size, _skip, _pos;
};

#endif
