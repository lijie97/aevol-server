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


#ifndef __AE_TIME_H__
#define __AE_TIME_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================






class Time
{
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  Time(void) = delete; //< Default ctor
  Time(const Time&) = delete; //< Copy ctor
  Time(Time&&) = delete; //< Move ctor

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~Time(void) = delete;

  // =================================================================
  //                        Accessors: getters
  // =================================================================
  static inline int64_t get_time(void) { return time; };

  // =================================================================
  //                        Accessors: setters
  // =================================================================
  static inline void set_time(int64_t t) { time = t; };

  // =================================================================
  //                              Operators
  // =================================================================

  // =================================================================
  //                            Public Methods
  // =================================================================
  static inline void plusplus(void) { time++; };

  // =================================================================
  //                           Public Attributes
  // =================================================================





 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  static int64_t time;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================

// =====================================================================
//                           Setters' definitions
// =====================================================================

// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

int64_t get_time(void);

} // namespace aevol

#endif // __AE_TIME_H__
