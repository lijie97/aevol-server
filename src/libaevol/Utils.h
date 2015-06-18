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


#ifndef AEVOL_UTILS_H__
#define AEVOL_UTILS_H__


// =================================================================
//                              Libraries
// =================================================================
#include <cinttypes>
#include <cassert>
#include <cstdlib>

#include <string>
#include <iostream>



using std::string;
using std::cout;
using std::endl;


// =================================================================
//                            Project Files
// =================================================================




// =================================================================
//                          Class declarations
// =================================================================

namespace aevol {

class Utils
{
 public :
  static inline int32_t mod(int32_t a, int32_t b);
  static inline int64_t mod(int64_t a, int64_t b);
  static inline int32_t min(int32_t a, int32_t b);
  static inline int32_t max(int32_t a, int32_t b);
  static inline void    exchange(int32_t &a, int32_t &b);
  static inline int16_t hamming(const char* str1, const char* str2);

  static inline void ExitWithMsg(const string& msg,
                                 const string& file, int line);

 protected :
  Utils( const Utils &model )
  {
    printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  };
};


inline int32_t Utils::mod(int32_t a, int32_t b)
{
  assert(b > 0);
  
  while (a < 0)  a += b;
  while (a >= b) a -= b;
  
  return a;
}

inline int64_t Utils::mod(int64_t a, int64_t b)
{
  assert(b > 0);
  
  while (a < 0)  a += b;
  while (a >= b) a -= b;
  
  return a;
}

inline int32_t Utils::min(int32_t a, int32_t b)
{
  return ((a < b)? a : b);
}

inline int32_t Utils::max(int32_t a, int32_t b)
{
  return ((a > b)? a : b);
}

inline void Utils::exchange(int32_t &a, int32_t &b)
{
  int32_t tmp = a;
  a = b;
  b = tmp;
}

/**
 * Print an error message and exit (with error status)
 *
 * \param msg message to be printed
 * \param file should be __FILE__
 * \param line should be __LINE__
 */
inline void Utils::ExitWithMsg(const string& msg,
                               const string& file, int line) {
  cout << file << ":" << line << ": error: " << msg << endl;
  exit(EXIT_FAILURE);
}

} // namespace aevol
#endif // AEVOL_UTILS_H__
