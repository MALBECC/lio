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
  ThreadBufferPool(int elements) {
    int bytes = elements * sizeof(T);
    bytes = bytes + ALIGNMENT - (bytes % ALIGNMENT);
    _pool = (T *) mkl_malloc(bytes, ALIGNMENT);
    _total_size = bytes;
    _pos = 0;
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
  inline T * get_pool(int elements) {
    if(_pos >= _total_size){
      throw new std::runtime_error("Thread pool has received more requests than can be satisfied");
    }
    T * ret = _pool + _pos;
    _pos += elements;
    return ret;
  }
private:
  T * _pool;
  int _buffer_size, _total_size, _pos;
};

#endif
