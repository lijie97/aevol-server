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


#ifndef __AE_HABITAT_H__
#define __AE_HABITAT_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <zlib.h>


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================






class Habitat
{
 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Habitat(void); //< Default ctor
  Habitat(const Habitat&) = delete; //< Copy ctor
  Habitat(Habitat&&) = delete; //< Move ctor
  Habitat(gzFile backup_file);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Habitat(void) = default; //< Destructor

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  double compound_amount(void) const {return compound_amount_;};

  // ==========================================================================
  //                                 Setters
  // ==========================================================================
  void set_compound_amount(double compound_amount) {
    compound_amount_ = compound_amount;
  };

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void save(gzFile backup_file) const;





 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================

  // Amount of secreted compound currently present in the grid cell 
  double compound_amount_;
};


// ============================================================================
//                           Getters' definitions
// ============================================================================

// ============================================================================
//                           Setters' definitions
// ============================================================================

// ============================================================================
//                          Operators' definitions
// ============================================================================

// ============================================================================
//                       Inline functions' definition
// ============================================================================

} // namespace aevol

#endif // __AE_HABITAT_H__
