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

Habitat_R::Habitat_R(const Habitat_R& rhs, bool share_phenotypic_target) {
  assert(share_phenotypic_target);
  phenotypic_targets_ = rhs.phenotypic_targets_;
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
    //printf("Taille de l'habitat après reset : %d\n", phenotypic_targets_.size());
  }

  void Habitat_R::initializePhenotypicTargets(int nb_indiv_age) {
    phenotypic_targets_.clear();
    //printf("Taille de l'habitat après le clear dans initialize... : %d\n", phenotypic_targets_.size());
    phenotypic_targets_.reserve(nb_indiv_age);
    PhenotypicTarget_R* env_to_add;
    for (int i = 0; i < nb_indiv_age; ++i) {
      PhenotypicTarget_R* env_to_add = phenotypic_target_handler().model_pointer( 0 );
      phenotypic_targets_.push_back(env_to_add);
    }
    //printf("Taille de l'habitat avant applyvariation : %d\n", phenotypic_targets_.size());
    ApplyVariation();  
  }

  void Habitat_R::addEnv( int8_t env_id ) {
    PhenotypicTarget_R* env_to_add = phenotypic_target_handler().model_pointer( env_id );
    phenotypic_targets_.push_back(env_to_add);
  }

  void Habitat_R::changeEnv( int8_t ind, int8_t env_id ) {
    PhenotypicTarget_R* new_env = phenotypic_target_handler().model_pointer( env_id );
    phenotypic_targets_.at(ind) = new_env;
  }



// ============================================================================
//                                   Methods
// ============================================================================
void Habitat_R::ApplyVariation() {
  printf("Appel au apply_variation de habitat_R\n");
  dynamic_cast<PhenotypicTargetHandler_R*>(phenotypic_target_handler_)->ApplyVariation( *this );
}
void Habitat_R::save(gzFile backup_file,
                     bool skip_phenotypic_target /*= false*/) const {
  printf("Appel a la sauvegarde de Habitat_R\n");
  gzwrite(backup_file, &compound_amount_, sizeof(compound_amount_));
  if (not skip_phenotypic_target) {
    phenotypic_target_handler_->save(backup_file);
  }
  // We have to save phenotypic_targets_ but we cannot since its a vector of pointer
  // Thus we save the ids
  int8_t nb_env = phenotypic_targets_.size();
  gzwrite(backup_file, &nb_env, sizeof(nb_env));
  int8_t id = 0;
  for (PhenotypicTarget_R* env: phenotypic_targets_) {
    id = env->get_id();
    gzwrite(backup_file, &id, sizeof(id));
  }
}

void Habitat_R::load(gzFile backup_file, 
                     PhenotypicTargetHandler_R* phenotypic_target_handler) {
  printf("Appel au chargement de Habitat_R\n");
  gzread(backup_file, &compound_amount_, sizeof(compound_amount_));
  if (phenotypic_target_handler == NULL)
    phenotypic_target_handler_ = new PhenotypicTargetHandler_R(backup_file);
  else
    phenotypic_target_handler_ = phenotypic_target_handler;

  int8_t nb_env = 0;
  gzread(backup_file, &nb_env, sizeof(nb_env));
  int8_t id = 0;
  PhenotypicTarget_R* env_to_add = NULL;
  for (int8_t i = 0 ; i < nb_env ; i++) {
    gzread(backup_file, &id, sizeof(id));
    env_to_add = this->phenotypic_target_handler().model_pointer( id );
    phenotypic_targets_.push_back(env_to_add);
  }
}
// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
