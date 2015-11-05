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


#ifndef AEVOL_PHENOTYPIC_TARGET_HANDLER_R_H__
#define AEVOL_PHENOTYPIC_TARGET_HANDLER_R_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <memory>
#include <list>

#include "PhenotypicTarget.h"
#include "Gaussian.h"
#include "ae_enums.h"
#include "JumpingMT.h"
#include "AbstractFuzzy.h"

using std::list;


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================




/**
 * Manages a phenotypic target and its "evolution" over time
 *
 * Handles a phenotypic target, the variation and/or noise that may be applied
 * to it as well as the set of possible phenotypic targets and the rules that
 * define how and when we switch from one to another
 */
class PhenotypicTargetHandler_R : public virtual PhenotypicTargetHandler
{
 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  PhenotypicTargetHandler_R(void); //< Default ctor
  PhenotypicTargetHandler_R(const PhenotypicTargetHandler_R&); //< Copy ctor
  PhenotypicTargetHandler_R(PhenotypicTargetHandler_R&&) = delete; //< Move ctor
  PhenotypicTargetHandler_R(gzFile backup_file);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~PhenotypicTargetHandler_R(void); //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void BuildPhenotypicTarget();
  void ApplyVariation( const Habitat_R& habitat );

  void save(gzFile backup_file) const;
  void load(gzFile backup_file);

  // ==========================================================================
  //                                 Getters
  // ==========================================================================


  // ==========================================================================
  //                                 Setters
  // ==========================================================================


 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================


  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  std::vector<<PhenotypicTarget> > phenotypic_target_models_;
};

// ============================================================================
//                       Inline functions' definition
// ============================================================================

} // namespace aevol

#endif // AEVOL_PHENOTYPIC_TARGET_HANDLER_R_H__
