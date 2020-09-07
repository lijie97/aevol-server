//
// Created by arrouan on 27/07/17.
//

#ifndef AEVOL_INDIVIDUAL_7_H
#define AEVOL_INDIVIDUAL_7_H

#include <vector>
#include <map>
#include <set>

#include "DnaFactory.h"
#include "Observable.h"
#include "PhenotypicTargetHandler.h"
#include "Stats.h"
#include "Vector_Fuzzy.h"
#include "ae_enums.h"

namespace aevol {

class ExpManager;
class Dna_7;
class Abstract_Metadata;
class Map_Metadata;

class Individual_7 : public Observable {
 public:
  Individual_7(ExpManager* exp_m, double w_max, DnaFactory* dna_factory);

  Individual_7(ExpManager* exp_m, Individual_7* clone, DnaFactory* dna_factory);

  ~Individual_7();

#ifdef PHENOTYPE_VECTOR
  double phenotype[PHENOTYPE_VECTOR_SIZE];
  double delta[PHENOTYPE_VECTOR_SIZE];
#else
  Vector_Fuzzy* phenotype;
  Vector_Fuzzy* delta;
#endif
  double fitness;
  double metaerror;

  Dna_7* dna_;

  int32_t indiv_id;
  int32_t parent_id;

  int32_t usage_count_ = 1;

  ExpManager* exp_m_;
  DnaFactory* dna_factory_;

  int global_id = -1;

  double w_max_;

  /** Variables for Tree mgmt **/
  int32_t nb_genes_activ     = 0;
  int32_t nb_genes_inhib     = 0;
  int32_t nb_func_genes      = 0;
  int32_t nb_non_func_genes  = 0;
  int32_t nb_coding_RNAs     = 0;
  int32_t nb_non_coding_RNAs = 0;
  /** END of Variables for Tree Mgmt **/

  void reset_stats() {
    nb_genes_activ     = 0;
    nb_genes_inhib     = 0;
    nb_func_genes      = 0;
    nb_non_func_genes  = 0;
    nb_coding_RNAs     = 0;
    nb_non_coding_RNAs = 0;
  }

  void rebuild_index();

  Abstract_Metadata* metadata_;
};



}
#endif //RAEVOL_CUDA_SIMD_INDIVIDUAL_H
