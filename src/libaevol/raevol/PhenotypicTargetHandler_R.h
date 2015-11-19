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
//#include <list>
#include <vector>

#include "PhenotypicTargetHandler.h"
#include "PhenotypicTarget_R.h"
#include "../Utils.h"
//#include "Habitat_R.h"

//using std::list;


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================
class Habitat_R;



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
  void ApplyVariation( Habitat_R& habitat );
  void BuildPhenotypicTargets();
  void BuildPhenotypicTarget( int8_t id);
  virtual void BuildPhenotypicTarget() {
      Utils::ExitWithDevMsg("You should not call a phenotypic target without age id in RAevol", __FILE__, __LINE__);
  }
  void print_geometric_areas();
  void save(gzFile backup_file) const;
  void load(gzFile backup_file);

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  PhenotypicTarget_R* model_pointer(int8_t env_id) const{
    return phenotypic_target_models_.at(env_id);
  }

  // ==========================================================================
  //                                 Setters
  // ==========================================================================
  void set_gaussians(const std::vector<std::list<Gaussian>>& gaussians_list) {
    env_gaussians_list_ = gaussians_list;
  }

  void set_switch_probability(double p) {
    env_switch_probability_ = p;
  }

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================


  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  std::vector<PhenotypicTarget_R*> phenotypic_target_models_;
  std::vector<std::list<Gaussian>> env_gaussians_list_;
  double env_switch_probability_;
};

// ============================================================================
//                       Inline functions' definition
// ============================================================================

} // namespace aevol

#endif // AEVOL_PHENOTYPIC_TARGET_HANDLER_R_H__
