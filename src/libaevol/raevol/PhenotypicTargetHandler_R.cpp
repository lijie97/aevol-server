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
PhenotypicTargetHandler_R::PhenotypicTargetHandler_R() {
  // The phenotypic target
#if __cplusplus == 201103L
  phenotypic_target_ = make_unique<PhenotypicTarget>();
#else
  phenotypic_target_ = std::make_unique<PhenotypicTarget>();
#endif

  // Sampling
  sampling_ = 0;

  // Variation
  var_prng_   = NULL;
  var_method_ = NO_VAR;
  var_sigma_  = 0.0;
  var_tau_    = 0;

  // Noise
  cur_noise_          = NULL;
  noise_method_       = NO_NOISE;
  noise_prng_         = NULL;
  noise_prob_         = 0.0;
  noise_alpha_        = 0.0;
  noise_sigma_        = 0.0;
  noise_sampling_log_ = 8;
}

PhenotypicTargetHandler_R::PhenotypicTargetHandler_R(
    const PhenotypicTargetHandler_R& rhs) {
  // ------------------------------------------------ Current Phenotypic Target
#if __cplusplus == 201103L
  phenotypic_target_ = make_unique<PhenotypicTarget>(*(rhs.phenotypic_target_));
#else
  phenotypic_target_ = std::make_unique<PhenotypicTarget>(*(rhs.phenotypic_target_));
#endif

  // ---------------------------------------------------------------- Gaussians
  initial_gaussians_ = rhs.initial_gaussians_;
  current_gaussians_ = rhs.current_gaussians_;

  // ----------------------------------------------------------------- Sampling
  sampling_ = rhs.sampling_;

  // ---------------------------------------------------------------- Variation
  var_method_ = rhs.var_method_;
  var_prng_ = rhs.var_prng_;
  var_sigma_ = rhs.var_sigma_;
  var_tau_ = rhs.var_tau_;

  // -------------------------------------------------------------------- Noise
  cur_noise_ = rhs.cur_noise_;
  noise_prng_ = rhs.noise_prng_;
  noise_method_ = rhs.noise_method_;
  noise_alpha_ = rhs.noise_alpha_;
  noise_sigma_ = rhs.noise_sigma_;
  noise_prob_ = rhs.noise_prob_;
  noise_sampling_log_ = rhs.noise_sampling_log_;
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
      if ( nb_env_in_list <= 1 )
      {
        break;
      }

      //reset the vector of phenotypic targets keeping only the last environment
      habitat.resetPhenotypicTargets();

      // Shortcuts used :
      PhenotypicTarget_R actual_env = habitat.phenotypic_target(0);
      int8_t id_old_env = actual_env.get_id();
      int8_t id_new_env = 0;



      // At each age we have to add the environment of this age to habitat
      for (int8_t i = 0; i < last_age ; i++)
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

void PhenotypicTargetHandler_R::save(gzFile backup_file) const {
  // --------------------------------------------------------------------------
  //  Write phenotypic target segmentation
  phenotypic_target_->SaveSegmentation(backup_file);

  // --------------------------------------------------------------------------
  //  Write current gaussians (initial gaussians will be stored later if
  // necessary)
  int8_t nb_gaussians = current_gaussians_.size();
  gzwrite(backup_file, &nb_gaussians, sizeof(nb_gaussians));
  for (const Gaussian & g: current_gaussians_)
    g.save(backup_file);

  // --------------------------------------------------------------------------
  //  Write sampling
  gzwrite(backup_file, &sampling_, sizeof(sampling_));

  // --------------------------------------------------------------------------
  //  Write variation data
  int8_t tmp_var_method = var_method_;
  gzwrite(backup_file, &tmp_var_method,  sizeof(tmp_var_method));

  if (var_method_ != NO_VAR) {
    var_prng_->save(backup_file);
    gzwrite(backup_file, &var_sigma_, sizeof(var_sigma_));
    gzwrite(backup_file, &var_tau_,   sizeof(var_tau_));
  }

  // --------------------------------------------------------------------------
  //  Write noise data
  int8_t tmp_noise_method = noise_method_;
  gzwrite(backup_file, &tmp_noise_method, sizeof(tmp_noise_method));

  if (noise_method_ != NO_NOISE) {
    int8_t tmp_save_cur_noise = (cur_noise_ != NULL);
    gzwrite(backup_file, &tmp_save_cur_noise,  sizeof(tmp_save_cur_noise));
    if (tmp_save_cur_noise) cur_noise_->save(backup_file);

    noise_prng_->save(backup_file);
    gzwrite(backup_file, &noise_alpha_,  sizeof(noise_alpha_));
    gzwrite(backup_file, &noise_sigma_,  sizeof(noise_sigma_));
    gzwrite(backup_file, &noise_prob_,   sizeof(noise_prob_));
    gzwrite(backup_file, &noise_sampling_log_, sizeof(noise_sampling_log_));
  }

  // ---------------------------------------------------------------
  //  If needed, keep a copy of the initial state of the gaussians
  // ---------------------------------------------------------------
  if (var_method_ != NO_VAR || noise_method_ != NO_NOISE) {
    size_t nb_gaussians = initial_gaussians_.size();
    gzwrite(backup_file, &nb_gaussians, sizeof(nb_gaussians));
    for (const Gaussian & g: initial_gaussians_)
      g.save(backup_file);
  }
}

void PhenotypicTargetHandler_R::load(gzFile backup_file) {
  // --------------------------------------------------------------------------
  //  Retrieve phenotypic target segmentation
#if __cplusplus == 201103L
  phenotypic_target_ = make_unique<PhenotypicTarget>();
#else
  phenotypic_target_ = std::make_unique<PhenotypicTarget>();
#endif

  phenotypic_target_->LoadSegmentation(backup_file);

  // --------------------------------------------------------------------------
  //  Retrieve current gaussians
  int8_t nb_gaussians;
  gzread(backup_file, &nb_gaussians, sizeof(nb_gaussians));
  for (int8_t i = 0 ; i < nb_gaussians ; i++)
    current_gaussians_.push_back(Gaussian(backup_file));

  // --------------------------------------------------------------------------
  //  Retrieve sampling
  gzread(backup_file, &sampling_, sizeof(sampling_));

  // --------------------------------------------------------------------------
  //  Retrieve variation data
  int8_t tmp_var_method;
  gzread(backup_file, &tmp_var_method, sizeof(tmp_var_method));
  var_method_ = (PhenotypicTargetVariationMethod) tmp_var_method;

  if (var_method_ != NO_VAR) {
    var_prng_ = std::make_shared<JumpingMT>(backup_file);
    gzread(backup_file, &var_sigma_, sizeof(var_sigma_));
    gzread(backup_file, &var_tau_, sizeof(var_tau_));
  }

  // --------------------------------------------------------------------------
  //  Retrieve noise data
  int8_t tmp_noise_method;
  gzread(backup_file, &tmp_noise_method, sizeof(tmp_noise_method));
  noise_method_ = (PhenotypicTargetNoiseMethod) tmp_noise_method;

  if (noise_method_ != NO_NOISE) {
    int8_t tmp_cur_noise_saved;
    gzread(backup_file, &tmp_cur_noise_saved,  sizeof(tmp_cur_noise_saved));
    if (tmp_cur_noise_saved)
      cur_noise_ = FuzzyFactory::fuzzyFactory->create_fuzzy(backup_file);

    noise_prng_ = std::make_shared<JumpingMT>(backup_file);
    gzread(backup_file, &noise_alpha_, sizeof(noise_alpha_));
    gzread(backup_file, &noise_sigma_, sizeof(noise_sigma_));
    gzread(backup_file, &noise_prob_,  sizeof(noise_prob_));
    gzread(backup_file, &noise_sampling_log_, sizeof(noise_sampling_log_));
  }

  // --------------------------------------------------------------------------
  //  If needed, retrieve a copy of the initial state of the gaussians
  if (var_method_ != NO_VAR || noise_method_ != NO_NOISE) {
    size_t nb_gaussians;
    gzread(backup_file, &nb_gaussians, sizeof(nb_gaussians));
    for (size_t i = 0 ; i < nb_gaussians ; i++)
      initial_gaussians_.emplace_back(backup_file);
  }

  // --------------------------------------------------------------------------
  //  Build the phenotypic target
  PhenotypicTargetHandler::BuildPhenotypicTarget();
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
