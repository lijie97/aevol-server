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




// ============================================================================
//                                   Includes
// ============================================================================
#include "Habitat_R.h"
#include "PhenotypicTargetHandler_R.h"

#include <iostream>


using std::cout;
using std::endl;


namespace aevol {


//##############################################################################
//                                                                             #
//                                Class Habitat                                #
//                                                                             #
//##############################################################################

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
Habitat_R::Habitat_R(void) : phenotypic_targets_(0){
  compound_amount_ = 0.0;
  phenotypic_target_handler_ = new PhenotypicTargetHandler_R();
}

Habitat_R::Habitat_R(const Habitat_R& rhs, bool share_phenotypic_target) : phenotypic_targets_(0){
  assert(share_phenotypic_target);
  compound_amount_ = rhs.compound_amount_;
  phenotypic_target_handler_ = rhs.phenotypic_target_handler_;
}

Habitat_R::Habitat_R(gzFile backup_file,
                 PhenotypicTargetHandler_R*
                    phenotypic_target_handler) : phenotypic_targets_(0) {
  load(backup_file, phenotypic_target_handler);
}

// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                           Setters' definitions
// ============================================================================
  void Habitat_R::resetPhenotypicTargets() {
    PhenotypicTarget_R* last_env = phenotypic_targets_.back();
    int8_t size = phenotypic_targets_.size();
    phenotypic_targets_.clear();
    phenotypic_targets_.resize(1);
    phenotypic_targets_.at(0) = last_env;
    phenotypic_targets_.reserve(size);
  }

  void Habitat_R::addEnv( int8_t env_id ) {
    PhenotypicTarget_R* env_to_add = phenotypic_target_handler().model_pointer( env_id );
    phenotypic_targets_.push_back(env_to_add);
  }



// ============================================================================
//                                   Methods
// ============================================================================
void Habitat_R::ApplyVariation() {
  dynamic_cast<PhenotypicTargetHandler_R*>(phenotypic_target_handler_)->ApplyVariation( *this );
}

void Habitat_R::load(gzFile backup_file, 
                     PhenotypicTargetHandler_R* phenotypic_target_handler) {
  gzread(backup_file, &compound_amount_, sizeof(compound_amount_));
  if (phenotypic_target_handler == NULL)
    phenotypic_target_handler_ = new PhenotypicTargetHandler_R(backup_file);
  else
    phenotypic_target_handler_ = phenotypic_target_handler;
}
// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
