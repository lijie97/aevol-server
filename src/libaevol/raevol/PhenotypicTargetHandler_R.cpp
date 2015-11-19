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
#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#include "PhenotypicTargetHandler_R.h"
#include "Habitat_R.h"
#include "ExpSetup.h"
#include "HybridFuzzy.h"
#include "Utils.h"

#include <iostream>


using std::cout;
using std::endl;


namespace aevol {


//##############################################################################
//                                                                             #
//                        Class PhenotypicTargetHandler                        #
//                                                                             #
//##############################################################################

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
PhenotypicTargetHandler_R::PhenotypicTargetHandler_R() : PhenotypicTargetHandler(), phenotypic_target_models_(0) {
  env_switch_probability_ = 0.1;
}

PhenotypicTargetHandler_R::PhenotypicTargetHandler_R(
    const PhenotypicTargetHandler_R& rhs) : PhenotypicTargetHandler(rhs) {
  phenotypic_target_models_ = rhs.phenotypic_target_models_;
  env_gaussians_list_ = rhs.env_gaussians_list_;
  env_switch_probability_ = rhs.env_switch_probability_;
}

PhenotypicTargetHandler_R::PhenotypicTargetHandler_R(gzFile backup_file) {
  load(backup_file);
}

// ============================================================================
//                                 Destructor
// ============================================================================
PhenotypicTargetHandler_R::~PhenotypicTargetHandler_R() {
  delete cur_noise_;
}

// ============================================================================
//                                   Methods
// ============================================================================

void PhenotypicTargetHandler_R::ApplyVariation( Habitat_R& habitat ) {
  switch (var_method_) {
    case NO_VAR :
      return;
    case AUTOREGRESSIVE_MEAN_VAR :
      Utils::ExitWithDevMsg("Not implemented yet", __FILE__, __LINE__);
      //ApplyAutoregressiveMeanVariation();
      break;
    case AUTOREGRESSIVE_HEIGHT_VAR :
      Utils::ExitWithDevMsg("Not implemented yet", __FILE__, __LINE__);
      //ApplyAutoregressiveHeightVariation();
      break;
    case SWITCH_IN_A_LIST : {
      // Yoram : reprise du code que j'avais rajouté dans Raevol 3
      // Pour l'instant les signaux ne sont pas gérés, mais j'aimerais ne pas avoir à les gérer manuellement

      // A security in order to preserve the program from an infinite loop : while( id_new_env == id_old_env )
      int8_t nb_env_in_list = phenotypic_target_models_.size();
      int8_t last_age = habitat.number_of_phenotypic_targets();
      //printf("last_age = %d\n", last_age);
      if ( nb_env_in_list <= 1 )
      {
        break;
      }

      //reset the vector of phenotypic targets keeping only the last environment
      habitat.resetPhenotypicTargets();

      // Shortcuts used :
      PhenotypicTarget_R actual_env = habitat.phenotypic_target(1);
      int8_t id_old_env = actual_env.get_id();
      int8_t id_new_env = 0;

      //Special case for the first env that may change also :
      if ( var_prng_->random() < env_switch_probability_)
      {         
        //we have to change to a new env that have an id different from the old one
        while( id_new_env == id_old_env )
        {
          id_new_env = var_prng_->random(nb_env_in_list);
        }
        //The environment has changed
        id_old_env = id_new_env;
        habitat.changeEnv(0,id_new_env);
      }

      // At each age we have to add the environment of this age to habitat
      for (int8_t i = 1; i < last_age ; i++)
      {
        id_new_env = id_old_env;

        // if we have to change of environment :
        if ( var_prng_->random() < env_switch_probability_)
        {         
          //we have to change to a new env that have an id different from the old one
          while( id_new_env == id_old_env )
          {
            id_new_env = var_prng_->random(nb_env_in_list);
          }
          //The environment has changed
          id_old_env = id_new_env;
        }
        habitat.addEnv(id_new_env);
      }

      break; }
    default :
      Utils::ExitWithDevMsg("Unknown variation method", __FILE__, __LINE__);
      break;
  }

  // Phenotypic target has changed, recompute its area
  // Yoram : il faudra vérifier ce que ça fait et si on en a besoin
  //phenotypic_target_->ComputeArea();
}

void PhenotypicTargetHandler_R::BuildPhenotypicTargets() {
  // First of all we have to know how many models do we have :
  int8_t nb_models = env_gaussians_list_.size();
  //debug
  //printf("PhenotypicTargetHandler_R::BuildPhenotypicTargets : we have %d env\n", nb_models);
  for (int8_t i = 0; i < nb_models ; i++) {
    BuildPhenotypicTarget(i);
  }
}

void PhenotypicTargetHandler_R::InitPhenotypicTargets() {
  // First of all we have to know how many models do we have :
  int8_t nb_models = env_gaussians_list_.size();
  //debug
  //printf("PhenotypicTargetHandler_R::BuildPhenotypicTargets : we have %d env\n", nb_models);
  for (int8_t i = 0; i < nb_models ; i++) {
    phenotypic_target_models_.push_back(new PhenotypicTarget_R( i ));
  }
}



void PhenotypicTargetHandler_R::BuildPhenotypicTarget( int8_t id) {
  // NB : Extreme points (at abscissa X_MIN and X_MAX) will be generated, we need to erase the list first
  PhenotypicTarget_R* phenotypic_target = phenotypic_target_models_.at(id);
  phenotypic_target->fuzzy()->reset();

  // Generate sample points from gaussians
  if (not env_gaussians_list_.empty()) {
    for (int16_t i = 0; i <= sampling_; i++) {
      Point new_point = Point(
          X_MIN + (double) i * (X_MAX - X_MIN) / (double) sampling_, 0.0);
      for (const Gaussian& g: env_gaussians_list_.at(id))
        new_point.y += g.compute_y(new_point.x);
      phenotypic_target->fuzzy()->add_point(new_point.x, new_point.y);
    }

    if (FuzzyFactory::fuzzyFactory->get_fuzzy_flavor() == 1) {
      HybridFuzzy* fuz = (HybridFuzzy*) phenotypic_target->fuzzy();

      for (int i = 1; i < fuz->get_pheno_size(); i++) {
        if (fuz->get_points()[i] == 0.0) {
          int minL = i - 1;
          int maxL = i + 1;
          int dist = 1;

          while (fuz->get_points()[maxL] == 0.0) {
            maxL++;
            dist++;
          }
          double inc = 0.0;
          if (fuz->get_points()[maxL] > fuz->get_points()[minL]) {
            inc = (fuz->get_points()[maxL] - fuz->get_points()[minL]) / dist;
          } else {
            inc = (fuz->get_points()[minL] - fuz->get_points()[maxL]) / dist;
            minL = maxL;
          }

          for (int j = i; j < maxL; j++) {
            fuz->get_points()[j] = fuz->get_points()[minL] + inc;
            inc += inc;
          }

        }
      }
    }
  }
  // Add lower and upper bounds
  phenotypic_target->fuzzy()->clip(AbstractFuzzy::min, Y_MIN);
  phenotypic_target->fuzzy()->clip(AbstractFuzzy::max, Y_MAX);

  // Simplify (get rid of useless points)
  phenotypic_target->fuzzy()->simplify();

  // Compute areas (total and by feature)
  phenotypic_target->ComputeArea();
}

void PhenotypicTargetHandler_R::print_geometric_areas() {
  double area = 0.0;
  for (int8_t i = 0; i < phenotypic_target_models_.size() ; i++) {
    area = phenotypic_target_models_.at(i)->fuzzy()->get_geometric_area();
    printf("Entire geometric area of the phenotypic target %d: %f\n", i,
           area);
  }
}

void PhenotypicTargetHandler_R::save(gzFile backup_file) const {
  printf("Appel a la sauvegarde de PhenotypicTargetHandler_R\n");
  PhenotypicTargetHandler::save(backup_file);
  
  // Sauvegarde en plus
  gzwrite(backup_file, &env_switch_probability_, sizeof(env_switch_probability_));

  // Save gaussians :
  int8_t nb_gaussian_list = env_gaussians_list_.size();
  int8_t nb_gaussians = 0;
  gzwrite(backup_file, &nb_gaussian_list, sizeof(nb_gaussian_list));
  for (const std::list<Gaussian>& gaussian_list: env_gaussians_list_) {
    nb_gaussians = gaussian_list.size();
    gzwrite(backup_file, &nb_gaussians, sizeof(nb_gaussians));
    for (const Gaussian & g: gaussian_list) {
      g.save(backup_file);
    }
  }

  // Save segmentation :
  int8_t nb_models = phenotypic_target_models_.size();
  gzwrite(backup_file, &nb_models, sizeof(nb_models));
  for (PhenotypicTarget_R* model: phenotypic_target_models_) {
    model->save( backup_file);
  }
}

void PhenotypicTargetHandler_R::load(gzFile backup_file) {
  printf("Appel au chargement de PhenotypicTargetHandler_R\n");
  PhenotypicTargetHandler::load(backup_file);
  // Chargement en plus
  gzread(backup_file, &env_switch_probability_, sizeof(env_switch_probability_));

  //Load gaussians :
  int8_t nb_gaussian_list = 0;
  int8_t nb_gaussians = 0;
  std::list<Gaussian> current_gaussians;
  gzread(backup_file, &nb_gaussian_list, sizeof(nb_gaussian_list));
  for( int8_t i = 0; i<nb_gaussian_list; i++) {
    env_gaussians_list_.push_back( std::list<Gaussian>());
    current_gaussians = env_gaussians_list_.back();
    gzread(backup_file, &nb_gaussians, sizeof(nb_gaussians));
    for (int8_t i = 0 ; i < nb_gaussians ; i++) {
      current_gaussians.push_back(Gaussian(backup_file));
    }
  }

  // Now that gaussians are loader we can build our PhenotypicTargets
  InitPhenotypicTargets();

  //load segmentation :
  int8_t nb_models = 0;
  gzread(backup_file, &nb_models, sizeof(nb_models));
  for (int8_t i = 0 ; i < nb_models ; i++) {
    phenotypic_target_models_.at(i)->load( backup_file );
  }

  BuildPhenotypicTargets();
}

// ============================================================================
//                            Non inline accessors
// ============================================================================



} // namespace aevol
