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
PhenotypicTargetHandler_R::PhenotypicTargetHandler_R() : 
PhenotypicTargetHandler(), phenotypic_target_models_(0), phenotypic_targets_(0) {
  env_switch_probability_ = 0.1;
}

PhenotypicTargetHandler_R::PhenotypicTargetHandler_R(
    const PhenotypicTargetHandler_R& rhs) : PhenotypicTargetHandler(rhs) {
  phenotypic_target_models_ = rhs.phenotypic_target_models_;
  phenotypic_targets_       = rhs.phenotypic_targets_;
  env_gaussians_list_       = rhs.env_gaussians_list_;
  env_signals_list_         = rhs.env_signals_list_;
  signals_models_           = rhs.signals_models_;
  env_switch_probability_   = rhs.env_switch_probability_;
  _nb_indiv_age             = rhs._nb_indiv_age;
}

PhenotypicTargetHandler_R::PhenotypicTargetHandler_R(gzFile backup_file) {
  load(backup_file);
  //printf("Number of env %ld\n",phenotypic_targets_.size());
}

// ============================================================================
//                                 Destructor
// ============================================================================
PhenotypicTargetHandler_R::~PhenotypicTargetHandler_R() {
  for(Protein_R* prot : signals_models_) {
    delete prot;
  }
  signals_models_.clear();
}

// ============================================================================
//                                   Methods
// ============================================================================
void PhenotypicTargetHandler_R::ApplyVariation() {
  hasChanged_ = false;
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
      //printf("ApplyVariation SWITCH_IN_A_LIST\n");

      // A security in order to preserve the program from an infinite loop : while( id_new_env == id_old_env )
      int8_t nb_env_in_list = phenotypic_target_models_.size();
      int8_t last_age = phenotypic_targets_.size();
      //printf("last_age = %d\n", last_age);
      if ( nb_env_in_list <= 1 ) {
        break;
      } else if (deterministic && nb_env_in_list == 2) {

        if (!init_2_env) {
          init_2_env = true;

          phenotypic_targets_.clear();
          phenotypic_targets_.resize(_nb_indiv_age);
          //printf("Reserver %d\n",_nb_indiv_age);


          int8_t half = (int) _nb_indiv_age / 2;

          //printf("Half is %d %d \n",half, _nb_indiv_age);

          for (int i = 0; i < half; i++) {
            //printf("Add env %d\n",i);
            addEnv(i,0);
          }

          for (int i = 0; i < half; i++) {
            //printf("Add env %d\n",half+i);
            addEnv(half+i,1);
          }
          phenotypic_targets_.resize(_nb_indiv_age);
        }
        //printf("env is ok now or not\n");
        break;
      }

      int* list_of_old_target_id = new int[phenotypic_targets_.size()];
      int nb_old_env = phenotypic_targets_.size();
      int i = 0;

//      printf("Nb env targers %d\n",nb_old_env);

      for (int i = 0; i < phenotypic_targets_.size(); i++) {

//        printf("Target at %d is %d\n",i,phenotypic_targets_[i]->get_id());
        list_of_old_target_id[i] = phenotypic_targets_[i]->get_id();

      }

      //reset the vector of phenotypic targets keeping only the last environment
      ResetPhenotypicTargets();

  //    printf("Start creating new env\n");

      // Shortcuts used :
      int8_t id_old_env = phenotypic_targets_.at(0)->get_id();
      int8_t id_new_env = 0;

    //  printf("Computing first env\n");
      //Special case for the first env that may change also :
      if ( var_prng_->random() < env_switch_probability_) {         
        //we have to change to a new env that have an id different from the old one
        while( id_new_env == id_old_env ) {
          id_new_env = var_prng_->random(nb_env_in_list);
        }
        //The environment has changed
        id_old_env = id_new_env;
        changeEnv(0,id_new_env);
      }

      // At each age we have to add the environment of this age
      for (int8_t i = 1; i < last_age ; i++) {

        id_new_env = id_old_env;

        // if we have to change of environment :
        double env_chang = var_prng_->random();

        if ( env_chang < env_switch_probability_) {
          //we have to change to a new env that have an id different from the old one
          while( id_new_env == id_old_env ) {
            id_new_env = var_prng_->random(nb_env_in_list);
          }
          //The environment has changed
          id_old_env = id_new_env;
        }
        addEnv(i,id_new_env);
      //  printf("Computing first env %d : %d\n",i,id_new_env);
      }

      //printf("Nb indiv age %d -- %ld\n",_nb_indiv_age,phenotypic_targets_.size());

      phenotypic_targets_.resize(_nb_indiv_age);

      //printf("AX Nb indiv age %d -- %ld\n",_nb_indiv_age,phenotypic_targets_.size());

      //printf("Computing env changed or not\n");

      if (nb_old_env == phenotypic_targets_.size()) {
        i=0;
        for (PhenotypicTarget_R* target : phenotypic_targets_) {
          if (list_of_old_target_id[i] != target->get_id()) {
            hasChanged_ = true;
            break;
          }
          i++;
        }
      } else {
        hasChanged_ = true;
      }


//      printf("BX Nb indiv age %d -- %ld\n",_nb_indiv_age,phenotypic_targets_.size());
/*
      i = 0;
      for (PhenotypicTarget_R* target : phenotypic_targets_) {

        printf("Env Model %d at %d\n", i, target->get_id());
        i++;
      }
*/

      break; }
    default :
      Utils::ExitWithDevMsg("Unknown variation method", __FILE__, __LINE__);
      break;
  }
}

void PhenotypicTargetHandler_R::InitPhenotypicTargetsAndModels( int8_t nb_indiv_age ) {
  InitPhenotypicTargetsModels();
  BuildPhenotypicTargetsModels();
  InitPhenotypicTargets(nb_indiv_age);
}

void PhenotypicTargetHandler_R::print_geometric_areas() {
  double area = 0.0;
  for (int8_t i = 0; i < phenotypic_target_models_.size() ; i++) {
    area = phenotypic_target_models_.at(i)->fuzzy()->get_geometric_area();
    printf("Entire geometric area of the phenotypic target %d: %f\n", i, area);
  }
}

void PhenotypicTargetHandler_R::save(gzFile backup_file) const {
  //printf("Appel a la sauvegarde de PhenotypicTargetHandler_R\n");
  PhenotypicTargetHandler::save(backup_file);
  
  // Sauvegarde en plus
  gzwrite(backup_file, &env_switch_probability_, sizeof(env_switch_probability_));
  //printf("save nb age %d \n",_nb_indiv_age);

  gzwrite(backup_file, &_nb_indiv_age, sizeof(_nb_indiv_age));

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

  // Save signals_models :
  int8_t nb_signals_models = signals_models_.size();
  gzwrite(backup_file, &nb_signals_models, sizeof(nb_signals_models)); 
  for (Protein_R* p: signals_models_) {
    p->save(backup_file);
  }

  // We do not need to save signals_models_list as its just a copy

  // Save env_signals_list :
  int8_t nb_signal_list = env_signals_list_.size();
  int8_t nb_signals = 0;
  gzwrite(backup_file, &nb_signal_list, sizeof(nb_signal_list));
  for (const std::list<int8_t>& signals_list: env_signals_list_) {
    nb_signals = signals_list.size();
    gzwrite(backup_file, &nb_signals, sizeof(nb_signals));
    for (const int8_t & id: signals_list) {
      gzwrite(backup_file, &id, sizeof(id));
    }
  }

  // Save segmentation :
  int8_t nb_models = phenotypic_target_models_.size();
  gzwrite(backup_file, &nb_models, sizeof(nb_models));
  for (PhenotypicTarget_R* model: phenotypic_target_models_) {
    model->save( backup_file);
  }

  // We have to save phenotypic_targets_ but we cannot since its a vector of pointer
  // Thus we save the ids
  int8_t nb_env = phenotypic_targets_.size();
  gzwrite(backup_file, &nb_env, sizeof(nb_env));
  //printf("Nb env %d\n",nb_env);

  int8_t id = 0;
  for (PhenotypicTarget_R* env: phenotypic_targets_) {
    id = env->get_id();
    gzwrite(backup_file, &id, sizeof(id));
  }
}

void PhenotypicTargetHandler_R::load(gzFile backup_file) {
  //printf("Appel au chargement de PhenotypicTargetHandler_R\n");
  PhenotypicTargetHandler::load(backup_file);
  // Chargement en plus
  gzread(backup_file, &env_switch_probability_, sizeof(env_switch_probability_));

  gzread(backup_file, &_nb_indiv_age, sizeof(_nb_indiv_age));

  //Load gaussians :
  int8_t nb_gaussian_list = 0;
  int8_t nb_gaussians = 0;
  gzread(backup_file, &nb_gaussian_list, sizeof(nb_gaussian_list));
  //printf("Loading %d gaussians list\n", nb_gaussian_list);
  for( int8_t i = 0; i<nb_gaussian_list; i++) {
    env_gaussians_list_.push_back( std::list<Gaussian>());
    gzread(backup_file, &nb_gaussians, sizeof(nb_gaussians));
    //printf("There are %d gaussian in gaussians list %d\n", nb_gaussians, i);
    for (int8_t j = 0 ; j < nb_gaussians ; j++) {
      env_gaussians_list_.back().push_back(Gaussian(backup_file));
      //printf("Nb gaussians in current_gaussians : %d\n", env_gaussians_list_.back().size());
      /*printf("Gaussian %d. Height = %f, Mean = %f, width = %f\n",j, 
       env_gaussians_list_.back().back().get_height(),
       env_gaussians_list_.back().back().get_mean(),
       env_gaussians_list_.back().back().get_width()
        );*/
    }
  }

  // Load signals_models
  int8_t nb_signals_models = 0;
  gzread(backup_file, &nb_signals_models, sizeof(nb_signals_models)); 
  for (int8_t i = 0; i< nb_signals_models; i++) {
    signals_models_.push_back(new Protein_R(backup_file));
  }

  // Rebuild signals_models_list
  for(Protein_R* prot : signals_models_) {
    signals_models_list_.push_back(prot);
  }

  // Load env_signals_list
  int8_t nb_signal_list = 0;
  int8_t nb_signals = 0;
  int8_t id = 0;
  gzread(backup_file, &nb_signal_list, sizeof(nb_signal_list));
  for( int8_t i = 0; i<nb_signal_list; i++) {
    env_signals_list_.push_back( std::list<int8_t>());
    gzread(backup_file, &nb_signals, sizeof(nb_signals));
    for (int8_t j = 0 ; j < nb_signals ; j++) {
      gzread(backup_file, &id, sizeof(id));
      env_signals_list_.back().push_back(id);
    }
  }

  // Now that gaussians are loaded we can build our PhenotypicTargetsModels
  InitPhenotypicTargetsModels();

  //load segmentation :
  int8_t nb_models = 0;
  gzread(backup_file, &nb_models, sizeof(nb_models));
  for (int8_t i = 0 ; i < nb_models ; i++) {
    phenotypic_target_models_.at(i)->load( backup_file );
  }

  BuildPhenotypicTargetsModels();
  //Debug
  //print_geometric_areas();

  // We load the phenotypic targets after having built the models
  int8_t nb_env = 0;
  gzread(backup_file, &nb_env, sizeof(nb_env));
  id = 0;
  PhenotypicTarget_R* env_to_add = NULL;

  phenotypic_targets_.resize(nb_env);

  for (int8_t i = 0 ; i < nb_env ; i++) {
    gzread(backup_file, &id, sizeof(id));
    //printf("Restore %d : %d\n",i,id);
    addEnv(i,id);
  }

  //printf("Number of env XX %ld\n",phenotypic_targets_.size());
}

// ============================================================================
//                              Protected Methods
// ============================================================================
void PhenotypicTargetHandler_R::InitPhenotypicTargetsModels() {
  // First of all we have to know how many models do we have :
  int8_t nb_models = env_gaussians_list_.size();
  //debug
  //printf("PhenotypicTargetHandler_R::InitPhenotypicTargets : we have %d env\n", nb_models);

  phenotypic_target_models_.resize(nb_models);

  for (int8_t i = 0; i < nb_models ; i++) {
    phenotypic_target_models_[i] = new PhenotypicTarget_R( i );
  }

  phenotypic_target_models_.resize(nb_models);
}

void PhenotypicTargetHandler_R::BuildPhenotypicTargetsModels() {
  // First of all we have to know how many models do we have :
  int8_t nb_models = env_gaussians_list_.size();
  //debug
  //printf("PhenotypicTargetHandler_R::BuildPhenotypicTargets : we have %d env\n", nb_models);
  for (int8_t i = 0; i < nb_models ; i++) {
    BuildPhenotypicTargetModel(i);
  }
}

void PhenotypicTargetHandler_R::BuildPhenotypicTargetModel( int8_t id) {
  //printf("Appel a BuildPhenotypicTargetModel avec id = %d\n", id);
  // NB : Extreme points (at abscissa X_MIN and X_MAX) will be generated, we need to erase the list first

  PhenotypicTarget_R* phenotypic_target = phenotypic_target_models_.at(id);
  phenotypic_target->fuzzy()->reset();

  //printf("On a %d gaussiennes\n", env_gaussians_list_.at(id).size());

  // Generate sample points from gaussians
  if (not env_gaussians_list_.at(id).empty()) {
    for (int16_t i = 0; i <= sampling_; i++) {
      Point new_point = Point(
          X_MIN + (double) i * (X_MAX - X_MIN) / (double) sampling_, 0.0);
      int gi = 0;
      for (const Gaussian& g: env_gaussians_list_.at(id)) {
        gi++;
        new_point.y += g.compute_y(new_point.x);
      }
      phenotypic_target->fuzzy()->add_point(new_point.x, new_point.y);
    }

    if (FuzzyFactory::fuzzyFactory->get_fuzzy_flavor() == 1) {
      HybridFuzzy* fuz = (HybridFuzzy*) phenotypic_target->fuzzy();

      for (int i = 1; i < fuz->get_pheno_size(); i++) {
        if (fuz->points()[i] == 0.0) {
          int minL = i - 1;
          int maxL = i + 1;
          int dist = 1;

          while (fuz->points()[maxL] == 0.0) {
            maxL++;
            dist++;
          }
          double inc = 0.0;
          if (fuz->points()[maxL] > fuz->points()[minL]) {
            inc = (fuz->points()[maxL] - fuz->points()[minL]) / dist;
          } else {
            inc = (fuz->points()[minL] - fuz->points()[maxL]) / dist;
            minL = maxL;
          }

          for (int j = i; j < maxL; j++) {
            fuz->points()[j] = fuz->points()[minL] + inc;
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
  double area = phenotypic_target->fuzzy()->get_geometric_area();
    //printf("Entire geometric area of the phenotypic target %d: %f\n", id, area);

  //Add its signals to the env
  // build a temporary real signals list (not a id list)
  std::list<Protein_R*> signals_list;
  for(int8_t & signal_id : env_signals_list_.at(id))
  {
    signals_list.push_back(signals_models_.at(signal_id));
  }
  phenotypic_target->set_signals(signals_list);
}

void PhenotypicTargetHandler_R::ResetPhenotypicTargets() {
  PhenotypicTarget_R* last_env = phenotypic_targets_.back();
  int8_t size = phenotypic_targets_.size();
  phenotypic_targets_.clear();
  phenotypic_targets_.resize(1);
  phenotypic_targets_.at(0) = last_env;
  phenotypic_targets_.resize(size);
  //printf("Taille de l'habitat après reset : %d\n", phenotypic_targets_.size());
}

void PhenotypicTargetHandler_R::InitPhenotypicTargets(int8_t nb_indiv_age) {
  phenotypic_targets_.clear();
  //printf("Taille de l'habitat après le clear dans initialize... : %d\n", phenotypic_targets_.size());
  //phenotypic_targets_.resize(nb_indiv_age);
  PhenotypicTarget_R* env_to_add;
  phenotypic_targets_.resize(nb_indiv_age);

  for (int i = 0; i < nb_indiv_age; ++i) {
    addEnv(i,0);
  }

  _nb_indiv_age = nb_indiv_age;
  phenotypic_targets_.resize(nb_indiv_age);
  //printf("NB AGE AT INIT %d %d\n",_nb_indiv_age,nb_indiv_age);
  //printf("Taille de l'habitat avant applyvariation : %d\n", phenotypic_targets_.size());
  ApplyVariation();  
}

void PhenotypicTargetHandler_R::addEnv( int time, int8_t env_id ) {
  assert(env_id >= 0 && env_id <= phenotypic_target_models_.size());
  //printf("Add env %d : %d\n",time,phenotypic_target_models_[env_id]->get_id());

  phenotypic_targets_[time] = phenotypic_target_models_.at(env_id);
}

void PhenotypicTargetHandler_R::changeEnv( int8_t ind, int8_t env_id ) {
  assert(env_id >= 0 && env_id <= phenotypic_target_models_.size());
  phenotypic_targets_.at(ind) = phenotypic_target_models_.at(env_id);
}

void PhenotypicTargetHandler_R::set_single_env(int8_t id) {
  phenotypic_targets_.clear();
  phenotypic_targets_.reserve(_nb_indiv_age);


  for (int i = 0; i < _nb_indiv_age; i++) {
    addEnv(i,id);
  }
}

void PhenotypicTargetHandler_R::set_two_env(int8_t id_1, int8_t id_2) {
  phenotypic_targets_.clear();
  phenotypic_targets_.reserve(_nb_indiv_age);

  int8_t half = (int) _nb_indiv_age / 2;

  for (int i = 0; i < half; i++) {
    addEnv(i,id_1);
  }

  for (int i = 0; i < half; i++) {
    addEnv(half+i,id_2);
  }
}
// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
