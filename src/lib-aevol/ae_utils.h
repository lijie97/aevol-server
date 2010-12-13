#ifndef __AE_UTILS_H__
#define __AE_UTILS_H__

#include <inttypes.h>

class utils
{
  public :
    static inline int32_t mod( int32_t a, int32_t b );
    static inline int32_t min( int32_t a, int32_t b );
    static inline int32_t max( int32_t a, int32_t b );
    static inline void exchange( int32_t &a, int32_t &b );
    static inline int16_t hamming( const char* str1, const char* str2 );

  protected :
    utils( const utils &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
  
};


inline int32_t utils::mod( int32_t a, int32_t b )
{
  while ( a < 0 ) a += b;
  
  return (a % b);
}

inline int32_t utils::min( int32_t a, int32_t b )
{
  return ((a < b)? a : b);
}

inline int32_t utils::max( int32_t a, int32_t b )
{
  return ((a > b)? a : b);
}

inline void utils::exchange( int32_t &a, int32_t &b )
{
  int32_t tmp = a;
  a = b;
  b = tmp;
}

#endif // __AE_UTILS_H__
