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
// ****************************************************************************


#ifndef AEVOL_EXPMANAGER_7_H
#define AEVOL_EXPMANAGER_7_H

#include "Stats.h"
#include "DnaFactory.h"
#include "Observable.h"
#include "PhenotypicTargetHandler.h"
#include "raevol/SIMD_PhenotypicTargetHandler_R.h"
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

  void setup_individuals(double w_max, double selection_pressure);
  void run_a_step(double w_max, double selection_pressure);

  void do_mutation(int indiv_id);

  void start_stop_RNA(int indiv_id);
  void opt_prom_compute_RNA(int indiv_id);
  void compute_RNA(int indiv_id);
  void start_protein(int indiv_id);
  void compute_protein(int indiv_id);
  void translate_protein(int indiv_id, double w_max);
  void compute_phenotype(int indiv_id);
  void compute_fitness(int indiv_id, double selection_pressure, int env_id = -1);

  void write_stat();

  void check_result();
  void check_dna();
  void check_struct();
  void check_individual(int indiv_id, int x, int y);
  static bool standalone() { return standalone_simd; }


  void build_phenotypic_target(PhenotypicTargetHandler* phenotypic_target_handler);
#ifdef __REGUL
  void compute_network(int indiv_id, double selection_pressure);
  void update_network(int indiv_id, double selection_pressure);
  void evaluate_network(int indiv_id, double selection_pressure, int env_id);
  void finalize_network(int indiv_id, double selection_pressure);
  void solve_network(int indiv_id, double selection_pressure);
  void update_phenotype( int indiv_id );
#endif

  void set_stats(Stats* stats) { stats_ = stats; }

  Individual_7** current_individuals;
  Individual_7** previous_individuals;
  Individual_7* best_indiv;

  int32_t nb_indivs_;
  int32_t nb_clones_;

  static bool standalone_simd;//= true;
  int rna_grain_size = 32;
  int protein_grain_size = 32;

  DnaFactory* dna_factory_;

  long cumulate_size = 0;
  long cumulate_diff = 0;

  double* fitness_sum_tab_;
#ifdef __REGUL
  Vector_Fuzzy** targets;
SIMD_PhenotypicTargetHandler_R* phenotypic_target_handler_;
#else
#ifdef PHENOTYPE_VECTOR
  double* target;
#else
  Vector_Fuzzy* target;
#endif
#endif

 private:
  ExpManager* exp_m_;
  int* dna_size;

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
