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


#ifndef AEVOL_HABITAT_H__
#define AEVOL_HABITAT_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <memory>

#include <zlib.h>

#include "PhenotypicTargetHandler.h"


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
  Habitat(const Habitat&, bool share_phenotypic_target);
  Habitat(gzFile backup_file,
          std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler_);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Habitat(void) = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void ApplyVariation();
  void save(gzFile backup_file,
            bool skip_phenotypic_target = false) const;
  void load(gzFile backup_file,
            std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler);

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  double compound_amount(void) const {return compound_amount_;};
  const PhenotypicTarget& phenotypic_target() const {
    return phenotypic_target_handler_->phenotypic_target();
  }
  const PhenotypicTargetHandler& phenotypic_target_handler() const {
    return *phenotypic_target_handler_;
  }
  PhenotypicTargetHandler& phenotypic_target_handler_nonconst() const {
    return *phenotypic_target_handler_;
  }

  // ==========================================================================
  //                                 Setters
  // ==========================================================================
  void set_compound_amount(double compound_amount) {
    compound_amount_ = compound_amount;
  };

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  // Amount of secreted compound currently present in the grid cell 
  double compound_amount_;

  /** Handler for the phenotypic target and its "evolution" over time */
  std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler_;
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

#endif // AEVOL_HABITAT_H__
