//
// Created by arrouan on 11/08/15.
//

#ifndef AEVOL_MAKE_UNIQUE_H
#define AEVOL_MAKE_UNIQUE_H

#if __cplusplus == 201103L
#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
  return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}
#endif

#endif //AEVOL_MAKE_UNIQUE_H
