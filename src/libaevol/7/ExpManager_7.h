//
// Created by arrouan on 07/09/2020.
//

#ifndef AEVOL_EXPMANAGER_7_H
#define AEVOL_EXPMANAGER_7_H

#include "Stats.h"
#include "DnaFactory.h"
#include "Observable.h"
#include "PhenotypicTargetHandler.h"
#include "Vector_Fuzzy.h"

namespace aevol {

constexpr const char* PROM_SEQ_LEAD = "0101011001110010010110";
constexpr const char* PROM_SEQ_LAG  = "1010100110001101101001";

constexpr const char* SHINE_DAL_SEQ_LEAD = "011011000";
constexpr const char* SHINE_DAL_SEQ_LAG  = "100100111";

constexpr const char* PROTEIN_END_LEAD = "001";
constexpr const char* PROTEIN_END_LAG  = "110";

const int32_t PHENOTYPE_VECTOR_SIZE            = 120000;
constexpr const double D_PHENOTYPE_VECTOR_SIZE = 120000.0;


class ExpManager;
class Stats_7;

class ExpManager_7 : public Observable{
 public:
  ExpManager_7(ExpManager* exp_m);

  ~ExpManager_7();

  void run_a_step(double w_max, double selection_pressure,bool optim_prom = false);

  void do_mutation(int indiv_id);

  void start_stop_RNA(int indiv_id);
  void opt_prom_compute_RNA(int indiv_id);
  void compute_RNA(int indiv_id);
  void start_protein(int indiv_id);
  void compute_protein(int indiv_id);
  void translate_protein(int indiv_id, double w_max);
  void compute_phenotype(int indiv_id);
  void compute_fitness(int indiv_id, double selection_pressure);


  void check_result();
  void check_dna();
  void check_struct();
  void check_individual(int indiv_id, int x, int y);
  bool standalone() const { return standalone_; }


  void build_phenotypic_target(PhenotypicTargetHandler* phenotypic_target_handler);


  void set_stats(Stats* stats) { stats_ = stats; }

  Individual_7** internal_simd_struct;
  Individual_7** prev_internal_simd_struct;
  Individual_7* best_indiv;

  int32_t* next_generation_reproducer_;

  int32_t nb_indivs_;
  int32_t nb_clones_;

  static bool standalone_simd;//= true;
  int rna_grain_size = 32;
  int protein_grain_size = 32;

  DnaFactory* dna_factory_;

 private:
  ExpManager* exp_m_;
  int* dna_size;
#ifdef PHENOTYPE_VECTOR
  double* target;
#else
  Vector_Fuzzy* target;
#endif
  bool standalone_;
  bool first_gener_ = true;

  Stats_7* stats_best = nullptr;
  Stats_7* stats_mean = nullptr;



  long apply_mutation[1024];


  Stats* stats_;

  void selection(int indiv_id);

  void check_selection(int indiv_id);



};
}

#endif //AEVOL_EXPMANAGER_7_H
