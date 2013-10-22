// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
// 
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// 
//*****************************************************************************


#ifndef __AE_UTILS_H__
#define __AE_UTILS_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <assert.h>



// =================================================================
//                            Project Files
// =================================================================




// =================================================================
//                          Class declarations
// =================================================================






class ae_utils
{
  public :
    static inline int32_t mod( int32_t a, int32_t b );
    static inline int32_t min( int32_t a, int32_t b );
    static inline int32_t max( int32_t a, int32_t b );
    static inline void    exchange( int32_t &a, int32_t &b );
    static inline int16_t hamming( const char* str1, const char* str2 );

  protected :
    ae_utils( const ae_utils &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
  
};


inline int32_t ae_utils::mod( int32_t a, int32_t b )
{
  assert( b > 0 );
  
  while ( a < 0 )  a += b;
  while ( a >= b ) a -= b;
  
  return a;
}

inline int32_t ae_utils::min( int32_t a, int32_t b )
{
  return ((a < b)? a : b);
}

inline int32_t ae_utils::max( int32_t a, int32_t b )
{
  return ((a > b)? a : b);
}

inline void ae_utils::exchange( int32_t &a, int32_t &b )
{
  int32_t tmp = a;
  a = b;
  b = tmp;
}

#endif // __AE_UTILS_H__
