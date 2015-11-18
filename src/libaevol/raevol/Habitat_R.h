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


#ifndef AEVOL_HABITAT_R_H__
#define AEVOL_HABITAT_R_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <memory>

#include <zlib.h>

#include "PhenotypicTargetHandler_R.h"
#include "Habitat.h"


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================






class Habitat_R : public virtual Habitat
{
 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Habitat_R(void); //< Default ctor
  Habitat_R(const Habitat_R&) = delete; //< Copy ctor
  Habitat_R(Habitat_R&&) = delete; //< Move ctor
  Habitat_R(const Habitat_R&, bool share_phenotypic_target);
  Habitat_R(gzFile backup_file,
            PhenotypicTargetHandler_R* phenotypic_target_handler_);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Habitat_R(void) = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  virtual void ApplyVariation();
  virtual void load(gzFile backup_file,
                    PhenotypicTargetHandler_R* phenotypic_target_handler);
  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  virtual const PhenotypicTarget& phenotypic_target() const {
    Utils::ExitWithDevMsg("You should not call a phenotypic target without age id in RAevol", __FILE__, __LINE__);
  }

  const PhenotypicTarget_R& phenotypic_target(  int8_t age ) const {
    return *phenotypic_targets_.at(age);
  }

  int8_t number_of_phenotypic_targets() const {
    return phenotypic_targets_.size();
  }

  virtual const PhenotypicTargetHandler_R& phenotypic_target_handler() const {
    return *(dynamic_cast<PhenotypicTargetHandler_R*> (phenotypic_target_handler_));
  }
  virtual PhenotypicTargetHandler_R& phenotypic_target_handler_nonconst() const {
    return *(dynamic_cast<PhenotypicTargetHandler_R*> (phenotypic_target_handler_));;
  }

  // ==========================================================================
  //                                 Setters
  // ==========================================================================
  // This function keep only the last element of the vector
  void resetPhenotypicTargets();
  void addEnv( int8_t env_id );

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  std::vector<PhenotypicTarget_R*> phenotypic_targets_;
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

#endif // AEVOL_HABITAT_R_H__
