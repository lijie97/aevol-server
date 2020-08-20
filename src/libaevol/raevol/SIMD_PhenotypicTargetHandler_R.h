//
// Created by arrouan on 18/08/2020.
//

#ifndef AEVOL_SIMD_PHENOTYPICTARGETHANDLER_R_H
#define AEVOL_SIMD_PHENOTYPICTARGETHANDLER_R_H

#include "population/SIMD_Individual.h"
#include "PhenotypicTargetHandler_R.h"
#include "phenotype/Gaussian.h"
#include "raevol/Protein_R.h"
#include "Vector_Fuzzy.h"

#include <list>
#include <vector>

namespace aevol {
class pProtein;

class SIMD_PhenotypicTargetHandler_R {
 public:
  SIMD_PhenotypicTargetHandler_R(PhenotypicTargetHandler_R* handler, ExpSetup* exp_s);

  void ApplyVariation();

  std::vector<pProtein*> signals_models_;
  std::vector<std::list<int16_t>> env_signals_list_;

  Vector_Fuzzy** targets_fuzzy_;
  PhenotypicTargetVariationMethod var_method_;

  int16_t* list_env_id_;

  int16_t nb_indiv_age_;
  int16_t nb_eval_;
  int16_t nb_env_;
 protected:
  std::vector<std::list<Gaussian>> env_gaussians_list_;


  double env_switch_probability_;
  int16_t sampling_;


  Vector_Fuzzy** targets_fuzzy_by_id_;

  std::shared_ptr<JumpingMT> var_prng_;

  bool hasChanged_ = false;
};

}

#endif //AEVOL_SIMD_PHENOTYPICTARGETHANDLER_R_H
