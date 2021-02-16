#ifndef _DATATYPE_H
#define _DATATYPE_H
#include <cstddef>
#include <vector>
#include <omp.h>
#include <cstdint>

typedef unsigned idx_t;

#ifndef INLINE
# if __GNUC__ && !__GNUC_STDC_INLINE__
#  define INLINE extern inline
# else
#  define INLINE inline
# endif
#endif

template<typename T>
size_t get_bytes(const std::vector<T> & vec)
{return vec.capacity()*sizeof(T) + sizeof(std::vector<T>);}

#endif //_DATATYPE_H
