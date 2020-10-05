//
// Created by arrouan on 18/08/2020.
//

#include "SIMD_PhenotypicTargetHandler_R.h"
#include "ExpSetup.h"
#include "Protein_7.h"

namespace aevol {

SIMD_PhenotypicTargetHandler_R::SIMD_PhenotypicTargetHandler_R(PhenotypicTargetHandler_R *handler, ExpSetup* exp_s) {

  var_method_ = handler->var_method();

  if (var_method_ == NO_VAR)
    var_method_ = SWITCH_IN_A_LIST;

  env_gaussians_list_.resize(handler->env_gaussians_list_.size());
  nb_env_ = env_gaussians_list_.size();

  int i = 0;
  for (auto gaussian_env : handler->env_gaussians_list_) {
    for (auto gauss : gaussian_env) {
      env_gaussians_list_[i].emplace_back(gauss);
    }
    handler->env_signals_list_[i];
    i++;
  }

  env_signals_list_.resize(handler->env_signals_list_.size());
  i=0;
  for (auto env_signal : handler->env_signals_list_) {
    for (auto e_signal : env_signal) {
      env_signals_list_[i].push_back(e_signal);
    }
    i++;

  }

  i = 0;
  signals_models_.resize(handler->signals_models_.size());
  for (auto signal_protein : handler->signals_models_) {
    Protein_7* prot = new Protein_7(signal_protein);
    signals_models_[i] = prot;
    i++;
  }

  env_switch_probability_ = handler->env_switch_probability_;
  nb_indiv_age_ = handler->_nb_indiv_age;

  var_prng_ = handler->var_prng_;

  sampling_ = handler->sampling();

  targets_fuzzy_by_id_ = new Vector_Fuzzy*[nb_env_];

  for (int env_id = 0; env_id < nb_env_; env_id++) {
    targets_fuzzy_by_id_[env_id] = new Vector_Fuzzy();

    if (not env_gaussians_list_.at(env_id).empty()) {
      for (int16_t i = 0; i <= sampling_; i++) {
        Point new_point = Point(
            X_MIN + (double) i * (X_MAX - X_MIN) / (double) sampling_, 0.0);
        int gi = 0;
        for (const Gaussian& g: env_gaussians_list_.at(env_id)) {
          gi++;
          new_point.y += g.compute_y(new_point.x);
        }
        targets_fuzzy_by_id_[env_id]->add_point(new_point.x, new_point.y);
      }
    }
    // Add lower and upper bounds
    targets_fuzzy_by_id_[env_id]->clip(AbstractFuzzy::min, Y_MIN);
    targets_fuzzy_by_id_[env_id]->clip(AbstractFuzzy::max, Y_MAX);

    // Simplify (get rid of useless points)
    targets_fuzzy_by_id_[env_id]->simplify();
  }

  if (var_method_ == SWITCH_IN_A_LIST) {
    targets_fuzzy_ = new Vector_Fuzzy*[nb_indiv_age_];
    list_env_id_ = new int16_t[nb_indiv_age_];

    if (nb_env_ <= 1) {
      for (int age = 0; age < nb_indiv_age_; age++) {
        targets_fuzzy_[age] = targets_fuzzy_by_id_[0];
        list_env_id_[age] = 0;
      }
    } else {
      for (int age = 0; age < nb_indiv_age_; age++) {
        targets_fuzzy_[age] = targets_fuzzy_by_id_[handler->phenotypic_targets_[age]->get_id()];
        list_env_id_[age] = handler->phenotypic_targets_[age]->get_id();
        printf("Init Env at age %d is %d\n",age,handler->phenotypic_targets_[age]->get_id());
      }
    }

    nb_eval_ = exp_s->get_list_eval_step()->size();
  }


}

void SIMD_PhenotypicTargetHandler_R::ApplyVariation() {

printf("ApplyVar\n");
  switch (var_method_) {
  case NO_VAR :
    return;
  case SWITCH_IN_A_LIST : {
    if ( nb_env_ <= 1 ) {
      break;
    }

    int16_t* list_of_old_target_id = list_env_id_;
    list_env_id_ = new int16_t[nb_indiv_age_];

    // Shortcuts used
    int16_t id_new_env = list_of_old_target_id[nb_indiv_age_-1];
    int16_t id_old_env = list_of_old_target_id[nb_indiv_age_-1];

    hasChanged_ = false;

    for (int16_t i = 0; i < nb_indiv_age_ ; i++) {
        // if we have to change of environment :
        double env_chang = var_prng_->random();

        if (env_chang < env_switch_probability_) {
          //we have to change to a new env that have an id different from the old one
          while (id_new_env == id_old_env) {
            id_new_env = var_prng_->random(nb_env_);
          }
          //The environment has changed
          id_old_env = id_new_env;
        }

        list_env_id_[i] = id_new_env;
        targets_fuzzy_[i] = targets_fuzzy_by_id_[id_new_env];
        printf("ENV at Age %d is %d\n",i,id_new_env);
        if (list_env_id_[i] != list_of_old_target_id[i])
          hasChanged_ = true;
    }
    break; }
  case ONE_AFTER_ANOTHER:
    break;
  default :
    Utils::ExitWithDevMsg("Unknown variation method", __FILE__, __LINE__);
    break;
  }
}

}
