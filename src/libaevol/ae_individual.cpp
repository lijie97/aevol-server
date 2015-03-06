// Aevol - An in silico experimental evolution platform
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

#include <cstdio>
#include <cmath>
#include <sys/stat.h>
#include <algorithm>
#include <cassert>
#include <list>
#include <algorithm>

#include "ae_codon.h"
#include "ae_exp_setup.h"
#include "ae_exp_manager.h"
#include "ae_grid_cell.h"
#include "ae_genetic_unit.h"
#include "ae_population.h"
#include "ae_vis_a_vis.h"
#include "ae_utils.h"
#include "fuzzy.h"

#ifdef __NO_X
  #ifdef __REGUL
    #include "ae_individual_R.h"
  #else
    #include "ae_individual.h"
  #endif
#elif defined __X11
  #ifdef __REGUL
    #include "ae_individual_R_X11.h"
  #else
    #include "ae_individual_X11.h"
  #endif
#endif
namespace aevol {

/// TODO
const char* ae_individual::get_strain_name(void) const {
  return _strain_name;
}

/// TODO
int32_t ae_individual::get_age(void) const {
  return _age;
}

/// TODO
int32_t ae_individual::get_id() const {
  return _id;
}

/// TODO
double* ae_individual::get_dist_to_target_by_segment() const {
  return _dist_to_target_by_segment;
}

/*!
  Get the individual's rank in the population (1 for the worst indiv, POP_SIZE for the best)

  Warning: be sure you call sort_individuals() before using get_rank_in_population
*/
int32_t ae_individual::get_rank() const {
  return _rank;
}

/// TODO
ae_exp_manager* ae_individual::get_exp_m() const {
  return _exp_m;
}

/// TODO
ae_jumping_mt* ae_individual::get_mut_prng() const {
  return _mut_prng;
}

/// TODO
ae_jumping_mt* ae_individual::get_stoch_prng() const {
  return _stoch_prng;
}

/*!
  Returns the number of genetic units
*/
int16_t ae_individual::get_nb_genetic_units() const {
  return _genetic_unit_list.size();
}

/// Get the number of plasmids. That is, the number of genetic units
/// minus one DNA-based genetic unit.
int32_t ae_individual::get_nb_plasmids() const {
  return _genetic_unit_list.size() - 1;
}

/// TODO
int32_t ae_individual::get_amount_of_dna() const {
  int32_t amount = 0;
  for (const auto& gen_unit: _genetic_unit_list)
    amount += gen_unit.get_dna()->get_length();
  return amount;
}

/// Return the list of genetic units.
const std::list<ae_genetic_unit>& ae_individual::get_genetic_unit_list_std() const {
  return _genetic_unit_list;
}

std::list<ae_genetic_unit>& ae_individual::get_genetic_unit_list_std_nonconst() {
  return _genetic_unit_list;
}

/// Remove all the elements from the GU list except the firt and the
/// last ones. If the GU list has less that 2 elements, do nothing.
void ae_individual::drop_nested_genetic_units() {
  if (_genetic_unit_list.size() <= 2)
    return;

  _genetic_unit_list.erase(std::next(_genetic_unit_list.begin()),
                           std::prev(_genetic_unit_list.end()));
  assert(_genetic_unit_list.size() <= 2);
}

/// Returns genetic unit number `num_unit` (0 for main chromosome)
const ae_genetic_unit& ae_individual::get_genetic_unit(int16_t num_unit) const {
  assert(num_unit < static_cast<int32_t>(_genetic_unit_list.size()));
  auto it = _genetic_unit_list.begin();
  std::advance(it, num_unit);
  return *it;
}

/// Returns genetic unit number `num_unit` (0 for main chromosome) as
/// a non-constant reference. To be used when the purpose if to alter
/// the individual.
ae_genetic_unit& ae_individual::get_genetic_unit_nonconst(int16_t num_unit) {
  assert(num_unit < static_cast<int32_t>(_genetic_unit_list.size()));
  auto it = _genetic_unit_list.begin();
  std::advance(it, num_unit);
  return *it;
}

/// TODO
double ae_individual::get_dist_to_target_by_feature(ae_env_axis_feature feature) const {
  assert(_distance_to_target_computed);

  return _dist_to_target_by_feature[feature];
}

/// TODO
double ae_individual::get_fitness() const {
  assert(_fitness_computed);

  return _fitness;
}

/// TODO
double ae_individual::get_fitness_by_feature(ae_env_axis_feature feature) const {
  assert(_fitness_computed);

  return _fitness_by_feature[feature];
}

/// TODO
ae_grid_cell* ae_individual::get_grid_cell() const {
  return _grid_cell;
}

/// TODO
bool ae_individual::get_placed_in_population () const {
  return _placed_in_population;
}

/*!
  Returns the sequence of genetic unit number <num_unit> (0 for main chromosome)
*/
const char* ae_individual::get_genetic_unit_sequence(int16_t num_unit) const {
  return get_genetic_unit(num_unit).get_sequence();
}

/*!
  Returns the sequence length of genetic unit number <num_unit> (0 for main chromosome)
*/
int32_t ae_individual::get_genetic_unit_seq_length(int16_t num_unit) const {
  return get_genetic_unit(num_unit).get_seq_length();
}

/// TODO
Fuzzy* ae_individual::get_phenotype_activ() const {
  return _phenotype_activ;
}

/// TODO
Fuzzy* ae_individual::get_phenotype_inhib() const {
  return _phenotype_inhib;
}

/// TODO
ae_phenotype* ae_individual::get_phenotype() const {
  return _phenotype;
}

/// TODO
ae_replication_report* ae_individual::get_replic_report() const {
  return _replic_report;
}

/// TODO
const std::list<ae_protein*>& ae_individual::get_protein_list() const {
  return _protein_list;
}

/// TODO
const std::list<ae_rna*>& ae_individual::get_rna_list() const {
  return _rna_list;
}

/// TODO
int32_t ae_individual::get_total_genome_size() {
  if (! _statistical_data_computed) compute_statistical_data();
  return _total_genome_size;
}

/// TODO
int32_t ae_individual::get_nb_coding_RNAs() {
  if (! _statistical_data_computed) compute_statistical_data();
  return _nb_coding_RNAs;
}

/// TODO
int32_t ae_individual::get_nb_non_coding_RNAs() {
  if (! _statistical_data_computed) compute_statistical_data();
  return _nb_non_coding_RNAs;
}

/// TODO
double ae_individual::get_overall_size_coding_RNAs() {
  if (! _statistical_data_computed) compute_statistical_data();
  return _overall_size_coding_RNAs;
}

/// TODO
double ae_individual::get_av_size_coding_RNAs() {
  if (! _statistical_data_computed) compute_statistical_data();
  if (_nb_coding_RNAs != 0)
  {
    return _overall_size_coding_RNAs / _nb_coding_RNAs;
  }
  else return 0.0;
}

/// TODO
double ae_individual::get_overall_size_non_coding_RNAs() {
  if (! _statistical_data_computed) compute_statistical_data();
  return _overall_size_non_coding_RNAs;
}

/// TODO
double ae_individual::get_av_size_non_coding_RNAs() {
  if (! _statistical_data_computed) compute_statistical_data();
  if (_nb_non_coding_RNAs != 0)
  {
    return _overall_size_non_coding_RNAs / _nb_non_coding_RNAs;
  }
  else return 0.0;
}

/// TODO
int32_t ae_individual::get_nb_genes_activ() {
  if (! _statistical_data_computed) compute_statistical_data();
  return _nb_genes_activ;
}

/// TODO
int32_t ae_individual::get_nb_genes_inhib() {
  if (! _statistical_data_computed) compute_statistical_data();
  return _nb_genes_inhib;
}

/// TODO
int32_t ae_individual::get_nb_functional_genes() {
  if (! _statistical_data_computed) compute_statistical_data();
  return _nb_functional_genes;
}

/// TODO
int32_t ae_individual::get_nb_non_functional_genes() {
  if (! _statistical_data_computed) compute_statistical_data();
  return _nb_non_functional_genes;
}

/// TODO
double ae_individual::get_overall_size_functional_genes() {
  if (! _statistical_data_computed) compute_statistical_data();
  return _overall_size_functional_genes;
}

/// TODO
double ae_individual::get_av_size_functional_genes() {
  if (! _statistical_data_computed) compute_statistical_data();
  if (_nb_functional_genes != 0)
  {
    return _overall_size_functional_genes / _nb_functional_genes;
  }
  else return 0.0;
}

/// TODO
double ae_individual::get_overall_size_non_functional_genes() {
  if (! _statistical_data_computed) compute_statistical_data();
  return _overall_size_non_functional_genes;
}

/// TODO
double ae_individual::get_av_size_non_functional_genes() {
  if (! _statistical_data_computed) compute_statistical_data();
  if (_nb_non_functional_genes != 0)
  {
    return _overall_size_non_functional_genes / _nb_non_functional_genes;
  }
  else return 0.0;
}

/// TODO
int32_t ae_individual::get_nb_bases_in_0_CDS() {
  if (! _statistical_data_computed) compute_statistical_data();
  if (! _non_coding_computed) compute_non_coding();
  return _nb_bases_in_0_CDS;
}

/// TODO
int32_t ae_individual::get_nb_bases_in_0_functional_CDS() {
  if (! _statistical_data_computed) compute_statistical_data();
  if (! _non_coding_computed) compute_non_coding();
  return _nb_bases_in_0_functional_CDS;
}

/// TODO
int32_t ae_individual::get_nb_bases_in_0_non_functional_CDS() {
  if (! _statistical_data_computed) compute_statistical_data();
  if (! _non_coding_computed) compute_non_coding();
  return _nb_bases_in_0_non_functional_CDS;
}

/// TODO
int32_t ae_individual::get_nb_bases_in_0_RNA() {
  if (! _statistical_data_computed) compute_statistical_data();
  if (! _non_coding_computed) compute_non_coding();
  return _nb_bases_in_0_RNA;
}

/// TODO
int32_t ae_individual::get_nb_bases_in_0_coding_RNA() {
  if (! _statistical_data_computed) compute_statistical_data();
  if (! _non_coding_computed) compute_non_coding();
  return _nb_bases_in_0_coding_RNA;
}

/// TODO
int32_t ae_individual::get_nb_bases_in_0_non_coding_RNA() {
  if (! _statistical_data_computed) compute_statistical_data();
  if (! _non_coding_computed) compute_non_coding();
  return _nb_bases_in_0_non_coding_RNA;
}

/// TODO
int32_t ae_individual::get_nb_bases_in_neutral_regions() {
  if (! _statistical_data_computed) compute_statistical_data();
  if (! _non_coding_computed) compute_non_coding();
  return _nb_bases_in_neutral_regions;
}

/// TODO
int32_t ae_individual::get_nb_neutral_regions() {
  if (! _statistical_data_computed) compute_statistical_data();
  if (! _non_coding_computed) compute_non_coding();
  return _nb_neutral_regions;
}

/// TODO
double ae_individual::get_modularity() {
  printf("\n  WARNING : modularity measure not yet implemented.\n");
  //~ if (_modularity < 0) compute_modularity();
  //~ return _modularity;
  return 0;
}

/// TODO
double ae_individual::get_w_max() const {
  return _w_max;
}

// ------------------------------------------------------------- Mutation rates
/// TODO
double ae_individual::get_point_mutation_rate() const {
  return _mut_params->get_point_mutation_rate();
}

/// TODO
double ae_individual::get_small_insertion_rate() const {
  return _mut_params->get_small_insertion_rate();
}

/// TODO
double ae_individual::get_small_deletion_rate() const {
  return _mut_params->get_small_deletion_rate();
}

/// TODO
int16_t ae_individual::get_max_indel_size() const {
  return _mut_params->get_max_indel_size();
}

// ---------------------------------- Rearrangement rates (without alignements)
/// TODO
double ae_individual::get_duplication_rate() const {
  return _mut_params->get_duplication_rate();
}

/// TODO
double ae_individual::get_deletion_rate() const {
  return _mut_params->get_deletion_rate();
}

/// TODO
double ae_individual::get_translocation_rate() const {
  return _mut_params->get_translocation_rate();
}

/// TODO
double ae_individual::get_inversion_rate() const {
  return _mut_params->get_inversion_rate();
}

// ------------------------------------- Rearrangement rates (with alignements)
/// TODO
double ae_individual::get_neighbourhood_rate() const {
  return _mut_params->get_neighbourhood_rate();
}

/// TODO
double ae_individual::get_duplication_proportion() const {
  return _mut_params->get_duplication_proportion();
}

/// TODO
double ae_individual::get_deletion_proportion() const {
  return _mut_params->get_deletion_proportion();
}

/// TODO
double ae_individual::get_translocation_proportion() const {
  return _mut_params->get_translocation_proportion();
}

/// TODO
double ae_individual::get_inversion_proportion() const {
  return _mut_params->get_inversion_proportion();
}

// ---------------------------------------------------------------- Transfer
bool ae_individual::get_with_4pts_trans() const {
  return _mut_params->get_with_4pts_trans();
}

bool ae_individual::get_with_HT() const {
  return _mut_params->get_with_HT();
}

bool ae_individual::get_repl_HT_with_close_points() const {
  return _mut_params->get_repl_HT_with_close_points();
}

double ae_individual::get_HT_ins_rate() const {
  return _mut_params->get_HT_ins_rate();
}

double ae_individual::get_HT_repl_rate() const {
  return _mut_params->get_HT_repl_rate();
}

double ae_individual::get_repl_HT_detach_rate() const {
  return _mut_params->get_repl_HT_detach_rate();
}


// ---------------------------------------------------------------- Alignements
bool ae_individual::get_with_alignments(void) const {
 return _mut_params->get_with_alignments();
}

ae_align_fun_shape ae_individual::get_align_fun_shape() const {
 return _mut_params->get_align_fun_shape();
}

double ae_individual::get_align_sigm_lambda() const {
 return _mut_params->get_align_sigm_lambda();
}

int16_t ae_individual::get_align_sigm_mean() const {
 return _mut_params->get_align_sigm_mean();
}

int16_t ae_individual::get_align_lin_min() const {
 return _mut_params->get_align_lin_min();
}

int16_t ae_individual::get_align_lin_max() const {
 return _mut_params->get_align_lin_max();
}

int16_t ae_individual::get_align_max_shift() const {
 return _mut_params->get_align_max_shift();
}

int16_t ae_individual::get_align_w_zone_h_len() const {
 return _mut_params->get_align_w_zone_h_len();
}

int16_t ae_individual::get_align_match_bonus() const {
 return _mut_params->get_align_match_bonus();
}

int16_t ae_individual::get_align_mismatch_cost() const {
 return _mut_params->get_align_mismatch_cost();
}

/// TODO
bool ae_individual::get_with_stochasticity() const {
  return _with_stochasticity;
}

void ae_individual::set_allow_plasmids(bool allow_plasmids) {
  _allow_plasmids=allow_plasmids;
}

// Genome size constraints
/// TODO
int32_t ae_individual::get_min_genome_length() const {
  return _min_genome_length;
}

/// TODO
int32_t ae_individual::get_max_genome_length() const {
  return _max_genome_length;
}

// Plasmids settings
/// TODO
bool ae_individual::get_allow_plasmids() const {
  return _allow_plasmids;
}

/*!
  \brief Return the _int_probes

  \return _int_probes
*/
int32_t* ae_individual::get_int_probes () const {
  return _int_probes;
}

/*!
  \brief Return the _double_probes

  \return _double_probes
*/
double* ae_individual::get_double_probes (void) const {
  return _double_probes;
}


// =====================================================================
//                           Setters' definitions
// =====================================================================
void ae_individual::set_strain_name(char* name) {
  assert(name && strlen(name) < INT8_MAX); // Conservative, could be <=
  int8_t name_len = strlen(name);
  delete [] _strain_name;
  _strain_name = new char[name_len+1];
  memcpy(_strain_name, name, name_len+1);
}

/// TODO
void ae_individual::set_id(int32_t id) {
  _id = id;

  if (_replic_report != NULL)
  {
    _replic_report->set_id(id);
  }
}

/// TODO
void ae_individual::set_rank(int32_t rank) {
  _rank = rank;

  if (_replic_report != NULL)
  {
    _replic_report->set_rank(rank);
  }
}

/// TODO
void ae_individual::set_placed_in_population(bool placed_in_population) {
  _placed_in_population = placed_in_population;
}

/*!
  Set the individual's replication report
  To be used by post-treatment only
*/
void ae_individual::set_replication_report(ae_replication_report * rep) {
  _replic_report = rep;
}

/// TODO
void ae_individual::set_w_max(double w_max) {
  _w_max = w_max;
}



// Genome size constraints
/// TODO
void ae_individual::set_min_genome_length(int32_t min_genome_length) {
  _min_genome_length = min_genome_length;
}

/// TODO
void ae_individual::set_max_genome_length(int32_t max_genome_length) {
  _max_genome_length = max_genome_length;
}


void ae_individual::set_point_mutation_rate(double point_mutation_rate) {
  _mut_params->set_point_mutation_rate(point_mutation_rate);
}

void ae_individual::set_small_insertion_rate(double small_insertion_rate) {
  _mut_params->set_small_insertion_rate(small_insertion_rate);
}

void ae_individual::set_small_deletion_rate(double small_deletion_rate) {
  _mut_params->set_small_deletion_rate(small_deletion_rate);
}

void ae_individual::set_max_indel_size(int16_t max_indel_size) {
  _mut_params->set_max_indel_size(max_indel_size);
}

void ae_individual::set_duplication_rate(double duplication_rate) {
  _mut_params->set_duplication_rate(duplication_rate);
}

void ae_individual::set_deletion_rate(double deletion_rate) {
  _mut_params->set_deletion_rate(deletion_rate);
}

void ae_individual::set_translocation_rate(double translocation_rate) {
  _mut_params->set_translocation_rate(translocation_rate);
}

void ae_individual::set_inversion_rate(double inversion_rate) {
  _mut_params->set_inversion_rate(inversion_rate);
}

void ae_individual::set_neighbourhood_rate(double neighbourhood_rate) {
  _mut_params->set_neighbourhood_rate(neighbourhood_rate);
}

void ae_individual::set_duplication_proportion(double duplication_proportion) {
  _mut_params->set_duplication_proportion(duplication_proportion);
}

void ae_individual::set_deletion_proportion(double deletion_proportion) {
  _mut_params->set_deletion_proportion(deletion_proportion);
}

void ae_individual::set_translocation_proportion(double translocation_proportion) {
  _mut_params->set_translocation_proportion(translocation_proportion);
}

void ae_individual::set_inversion_proportion(double inversion_proportion) {
  _mut_params->set_inversion_proportion(inversion_proportion);
}

void ae_individual::set_with_4pts_trans(bool with_4pts_trans) {
        _mut_params->set_with_4pts_trans(with_4pts_trans);
}

void ae_individual::set_with_alignments(bool with_alignments) {
        _mut_params->set_with_alignments(with_alignments);
}

void ae_individual::set_with_HT(bool with_HT) {
        _mut_params->set_with_HT(with_HT);
}

void ae_individual::set_repl_HT_with_close_points(bool repl_HT_with_close_points) {
        _mut_params->set_repl_HT_with_close_points(repl_HT_with_close_points);
}

void ae_individual::set_HT_ins_rate(double HT_ins_rate) {
        _mut_params->set_HT_ins_rate(HT_ins_rate);
}

void ae_individual::set_HT_repl_rate(double HT_repl_rate) {
        _mut_params->set_HT_repl_rate(HT_repl_rate);
}

void ae_individual::set_repl_HT_detach_rate(double repl_HT_detach_rate) {
        _mut_params->set_repl_HT_detach_rate(repl_HT_detach_rate);
}


void ae_individual::set_with_stochasticity(bool with_stoch) {
  _with_stochasticity = with_stoch;
}

void ae_individual::set_stoch_prng(ae_jumping_mt* prng) {
  _stoch_prng = prng;
}

void ae_individual::set_mut_prng(ae_jumping_mt* prng) {
  _mut_prng = prng;
}

/*!
  \brief Change the _int_probes

  \param int_probes 5 int32_t* that constitute a probe
*/
void ae_individual::set_int_probes (int32_t* int_probes) {
  _int_probes = int_probes;
}

/*!
  \brief Change the _double_probes

  \param double_probes 5 double* that constitute a probe
*/
void ae_individual::set_double_probes (double* double_probes) {
  _double_probes = double_probes;
}

// =====================================================================
//                       functions' definition
// =====================================================================

void ae_individual::reset_dist_to_target_by_segment(double* dist_to_target_by_segment) {
  if (_dist_to_target_by_segment != NULL) delete [] _dist_to_target_by_segment;
  _dist_to_target_by_segment = dist_to_target_by_segment;
}

void ae_individual::renew_dist_to_target_by_feature() {
  if (_dist_to_target_by_feature != NULL) delete [] _dist_to_target_by_feature;
  _dist_to_target_by_feature = new double [NB_FEATURES];
}


void ae_individual::renew_fitness_by_feature() {
  if (_fitness_by_feature != NULL) delete [] _fitness_by_feature;
  _fitness_by_feature = new double [NB_FEATURES];
}

void ae_individual::do_transcription_translation_folding() {
  if (_transcribed == true && _translated == true && _folded == true) return;

  do_transcription();
  do_translation();
  do_folding();

  make_protein_list();
}

#ifdef DEBUG
  void ae_individual::assert_promoters()
  {
    // Perform assertion for each genetic unit
    for (auto& gen_unit: _genetic_unit_list)
      gen_unit.assert_promoters();
  }

  void ae_individual::assert_promoters_order()
  {
    // Perform assertion for each genetic unit
    for (auto& gen_unit: _genetic_unit_list)
      gen_unit.assert_promoters_order();
  }
#endif

//                             Constructors

ae_individual::ae_individual(ae_exp_manager* exp_m,
                              ae_jumping_mt* mut_prng,
                              ae_jumping_mt* stoch_prng,
                              ae_params_mut* param_mut,
                              double w_max,
                              int32_t min_genome_length,
                              int32_t max_genome_length,
                              bool allow_plasmids,
                              int32_t id,
                              char* strain_name,
                              int32_t age) {
  // Experiment manager
  _exp_m = exp_m;

  // PRNGs
  _mut_prng   = mut_prng;
  _stoch_prng = stoch_prng;

  // Replication Report
  _replic_report = NULL;

  // ID and rank of the indiv ; name and "age" of the strain
  set_id(id);
  _rank = -1; // TODO: UNRANKED
  _age = age;
  _strain_name = new char[strlen(strain_name)+1];
  strcpy(_strain_name, strain_name);

  _phenotype_activ  = NULL;
  _phenotype_inhib  = NULL;
  _phenotype        = NULL;

  _dist_to_target_by_segment  = NULL;
  _dist_to_target_by_feature  = new double [NB_FEATURES];
  for (int i=0; i<NB_FEATURES; i++)
  {
    _dist_to_target_by_feature[i]=0;
  }
  _fitness_by_feature         = new double [NB_FEATURES];
  for (int i=0; i<NB_FEATURES; i++)
  {
    _fitness_by_feature[i]=0;
  }

  _fitness = 0.0;
  // When using structured population, this is the cell the individual is in
  _grid_cell = NULL;

  // The chromosome and plasmids (if allowed)

  // Generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for (int8_t i = 0 ; i < 5 ; i++)
  {
    _int_probes[i]    = 0;
    _double_probes[i] = 0.0;
  }

  // Mutation rates etc...
  _mut_params = new ae_params_mut(*param_mut);

  // Artificial chemistry
  _w_max = w_max;

  // Genome size constraints
  _min_genome_length = min_genome_length;
  _max_genome_length = max_genome_length;

  // Plasmids settings
  _allow_plasmids         = allow_plasmids;


  // --------------------------------------------------
  // "State" of the individual
  // --------------------------------------------------
  _evaluated                    = false;
  _transcribed                  = false;
  _translated                   = false;
  _folded                       = false;
  _phenotype_computed           = false;
  _distance_to_target_computed  = false;
  _fitness_computed             = false;
  _statistical_data_computed    = false;
  _non_coding_computed          = false;
  _modularity_computed          = false;
  _placed_in_population         = false;



  // ----------------------------------------
  // Statistical data
  // ----------------------------------------
  // Genome, RNAs and genes size and stuff
  _total_genome_size  = 0;

  _nb_coding_RNAs     = 0;
  _nb_non_coding_RNAs = 0;

  _overall_size_coding_RNAs     = 0.0;
  _overall_size_non_coding_RNAs = 0.0;

  _nb_genes_activ           = 0;
  _nb_genes_inhib           = 0;
  _nb_functional_genes      = 0;
  _nb_non_functional_genes  = 0;

  _overall_size_functional_genes      = 0.0;
  _overall_size_non_functional_genes  = 0.0;

  // Coding / non-coding
  _nb_bases_in_0_CDS                = 0;
  _nb_bases_in_0_functional_CDS     = 0;
  _nb_bases_in_0_non_functional_CDS = 0;
  _nb_bases_in_0_RNA                = 0;
  _nb_bases_in_0_coding_RNA         = 0;
  _nb_bases_in_0_non_coding_RNA     = 0;
  _nb_bases_in_neutral_regions      = 0;
  _nb_neutral_regions               = 0;

  _modularity = 0.0;
}

/*!
  This constructor retreives an individual from a backup file.

  Since this generation has already been processed, no unnecessary calculation (e.g. fitness) will be done.
  No transcription, translation or other process of that kind is performed.
*/
ae_individual::ae_individual(ae_exp_manager* exp_m, gzFile backup_file) {
  _exp_m = exp_m;

  // Retrieve the name and "age" of the strain
  int8_t strain_string_len;
  gzread(backup_file, &strain_string_len, sizeof(strain_string_len));
  _strain_name = new char[strain_string_len+1];
  gzread(backup_file, _strain_name, strain_string_len+1);
  gzread(backup_file, &_age, sizeof(_age));

  // Retrieve the PRNGs
  #ifdef DISTRIBUTED_PRNG
    _mut_prng   = new ae_jumping_mt(backup_file);
    _stoch_prng = new ae_jumping_mt(backup_file);
  #else
    if (exp_m == NULL)
    {
      // Detached mode
      _mut_prng   = NULL;
      _stoch_prng = NULL;
    }
    else
    {
      _mut_prng   = exp_m->get_pop()->get_mut_prng();
      _stoch_prng = exp_m->get_pop()->get_stoch_prng();
    }
  #endif

  // Retreive id and rank
  gzread(backup_file, &_id,    sizeof(_id));
  gzread(backup_file, &_rank,  sizeof(_rank));

  // Retrieve spatial coordinates
  int16_t x;
  int16_t y;
  gzread(backup_file, &x, sizeof(x));
  gzread(backup_file, &y, sizeof(y));
  _placed_in_population = false;
  if (_exp_m != NULL)
  {
    set_grid_cell(_exp_m->get_grid_cell(x, y));
  }

  // Retreive generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  gzread(backup_file, _int_probes,     5 * sizeof(*_int_probes));
  gzread(backup_file, _double_probes,  5 * sizeof(*_double_probes));

  // Retrieve mutational parameters
  _mut_params = new ae_params_mut(backup_file);

  // ------------------------------------------------- Phenotypic stochasticity
  gzread(backup_file, &_with_stochasticity, sizeof(_with_stochasticity));

  // Retrieve artificial chemistry parameters
  gzread(backup_file, &_w_max, sizeof(_w_max));

  // Retrieve genome size constraints
  gzread(backup_file, &_min_genome_length, sizeof(_min_genome_length));
  gzread(backup_file, &_max_genome_length, sizeof(_max_genome_length));

  // Retrieve plasmids settings
  int8_t tmp_allow_plasmids;
  gzread(backup_file, &tmp_allow_plasmids, sizeof(tmp_allow_plasmids));
  _allow_plasmids = tmp_allow_plasmids ? 1 : 0;

  // Retreive genetic units
  int16_t nb_gen_units;
  gzread(backup_file, &nb_gen_units,  sizeof(nb_gen_units));

  for (int16_t i = 0 ; i < nb_gen_units ; i++)
    _genetic_unit_list.emplace_back(this, backup_file);

  // --------------------------------------------------------------------------------------------
  // No more data to retreive, the following are only structural initializations (no data is set)
  // --------------------------------------------------------------------------------------------

  // Create empty fuzzy sets for activation and inhibition
  _phenotype_activ  = NULL;
  _phenotype_inhib  = NULL;
  _phenotype        = NULL;

  _dist_to_target_by_segment  = NULL;
  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];

  for (int8_t i = 0 ; i < NB_FEATURES ; i++)
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }


  // Replication report, protein list and rna list initialization
  _replic_report = NULL;  // We are reloading from a backup, the replication report for the loaded generation
                          // is already in the tree file corresponding to this generation and needs not be re-written.
                          // NB : If the replication report is needed in future development, it will have to be
                          // loaded from the tree file.

  // Initialize the computational state of the individual
  _evaluated                    = false;
  _transcribed                  = false;
  _translated                   = false;
  _folded                       = false;
  _phenotype_computed           = false;
  _distance_to_target_computed  = false;
  _fitness_computed             = false;
  _statistical_data_computed    = false;
  _non_coding_computed          = false;
  _modularity_computed          = false;


  // Initialize statistical data
  _total_genome_size                  = 0;
  _nb_coding_RNAs                     = 0;
  _nb_non_coding_RNAs                 = 0;
  _overall_size_coding_RNAs           = 0;
  _overall_size_non_coding_RNAs       = 0;
  _nb_genes_activ                     = 0;
  _nb_genes_inhib                     = 0;
  _nb_functional_genes                = 0;
  _nb_non_functional_genes            = 0;
  _overall_size_functional_genes      = 0;
  _overall_size_non_functional_genes  = 0;

  _nb_bases_in_0_CDS                  = -1;
  _nb_bases_in_0_functional_CDS       = -1;
  _nb_bases_in_0_non_functional_CDS   = -1;
  _nb_bases_in_0_RNA                  = -1;
  _nb_bases_in_0_coding_RNA           = -1;
  _nb_bases_in_0_non_coding_RNA       = -1;
  _nb_bases_in_neutral_regions        = -1;
  _nb_neutral_regions                 = -1;

  _modularity = -1;

  //evaluate();
}

// Copy constructor
ae_individual::ae_individual(const ae_individual &model, bool replication_report_copy /* = FALSE */) {
  _exp_m = model._exp_m;

  // PRNGs
  #ifdef DISTRIBUTED_PRNG
    _mut_prng   = new ae_jumping_mt(*(model._mut_prng));
    _stoch_prng = new ae_jumping_mt(*(model._stoch_prng));
  #else
    _mut_prng   = model._mut_prng;
    _stoch_prng = model._stoch_prng;
  #endif

  int strain_string_len = strlen(model._strain_name);
  _strain_name = new char[strain_string_len+1];
  memcpy(_strain_name, model._strain_name, strain_string_len+1);
  _age  = model._age;

  _id   = model._id;
  _rank = model._rank;

  _evaluated                    = false;//model._evaluated;
  _transcribed                  = false;//model._transcribed;
  _translated                   = false;//model._translated;
  _folded                       = false;//model._folded;
  _phenotype_computed           = model._phenotype_computed;

  _with_stochasticity = model._with_stochasticity;

  // Artificial chemistry parameters
  _w_max = model._w_max;

  // The distance to target and what results from it depend on the environment
  // and must hence be recomputed with the (possibly different) environment.
  _distance_to_target_computed  = false;
  _fitness_computed             = false;
  _statistical_data_computed    = model._statistical_data_computed;

  _non_coding_computed          = model._non_coding_computed;
  _modularity_computed          = model._modularity_computed;

  _placed_in_population         = model._placed_in_population;
  _grid_cell = model._grid_cell;

  // Copy model's genetic units
  // Should actually use ae_genetic_unit copy ctor which is disabled.
  for (const auto& gu: model._genetic_unit_list)
    _genetic_unit_list.emplace_back(this, gu);

  // Copy phenotype
  if (_phenotype_computed)
  {
    _phenotype_activ  = new Fuzzy(*(model._phenotype_activ));
    _phenotype_inhib  = new Fuzzy(*(model._phenotype_inhib));
    _phenotype        = new ae_phenotype(this, *(model._phenotype));
  }
  else
  {
    _phenotype_activ  = NULL;
    _phenotype_inhib  = NULL;
    _phenotype        = NULL;
  }


  // Copy fitness-related stuff
  _dist_to_target_by_segment  = NULL;
  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];

  for (int8_t i = 0 ; i < NB_FEATURES ; i++)
  {
    _dist_to_target_by_feature[i] = model._dist_to_target_by_feature[i];
    _fitness_by_feature[i] = model._fitness_by_feature[i];
  }

  _fitness = model._fitness;


  // Copy statistical data
  _total_genome_size                  = model._total_genome_size;
  _nb_coding_RNAs                     = model._nb_coding_RNAs;
  _nb_non_coding_RNAs                 = model._nb_non_coding_RNAs;
  _overall_size_coding_RNAs           = model._overall_size_coding_RNAs;
  _overall_size_non_coding_RNAs       = model._overall_size_non_coding_RNAs;
  _nb_genes_activ                     = model._nb_genes_activ;
  _nb_genes_inhib                     = model._nb_genes_inhib;
  _nb_functional_genes                = model._nb_functional_genes;
  _nb_non_functional_genes            = model._nb_non_functional_genes;
  _overall_size_functional_genes      = model._overall_size_functional_genes;
  _overall_size_non_functional_genes  = model._overall_size_non_functional_genes;

  _nb_bases_in_0_CDS                  = model._nb_bases_in_0_CDS;
  _nb_bases_in_0_functional_CDS       = model._nb_bases_in_0_functional_CDS;
  _nb_bases_in_0_non_functional_CDS   = model._nb_bases_in_0_non_functional_CDS;
  _nb_bases_in_0_RNA                  = model._nb_bases_in_0_RNA;
  _nb_bases_in_0_coding_RNA           = model._nb_bases_in_0_coding_RNA;
  _nb_bases_in_0_non_coding_RNA       = model._nb_bases_in_0_non_coding_RNA;
  _nb_bases_in_neutral_regions        = model._nb_bases_in_neutral_regions;
  _nb_neutral_regions                 = model._nb_neutral_regions;

  _modularity = model._modularity;


  // Create a new replication report to store mutational events
  if (replication_report_copy && _exp_m->get_output_m()->get_record_tree() && (_exp_m->get_output_m()->get_tree_mode() == NORMAL) && (model._replic_report != NULL))
  {
    _replic_report = new ae_replication_report(*model._replic_report);
    _replic_report->set_indiv(this);
    // TODO: remove this after checking it is the old way
    //_exp_m->get_output_m()->get_tree()->set_replic_report(_id, _replic_report);
  }
  else
  {
    _replic_report = NULL;
  }

  // Generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for (int8_t i = 0 ; i < 5 ; i++)
  {
    _int_probes[i]    = model._int_probes[i];
    _double_probes[i] = model._double_probes[i];
  }

  // Mutation rates etc...
  _mut_params = new ae_params_mut(*(model._mut_params));


  // Genome size constraints
  _min_genome_length = model._min_genome_length;
  _max_genome_length = model._max_genome_length;

  // Plasmids settings
  _allow_plasmids         = model._allow_plasmids;

  evaluate();
}

/*!
  This constructor creates a new individual with the same genome as it's parent.
  The location of promoters will be copied but no further process will be performed.

  The phenotype and the fitness are not set, neither is the statistical data.
*/
ae_individual::ae_individual(ae_individual* const parent, int32_t id, ae_jumping_mt* mut_prng, ae_jumping_mt* stoch_prng) {
  _exp_m = parent->_exp_m;

  // PRNGs
  _mut_prng   = mut_prng;
  _stoch_prng = stoch_prng;

  int strain_string_len = strlen(parent->_strain_name);
  _strain_name = new char[strain_string_len+1];
  memcpy(_strain_name, parent->_strain_name, strain_string_len+1);
  _age  = parent->_age + 1;

  _id   = id;
  _rank = -1;

  _evaluated                    = false;
  _transcribed                  = false;
  _translated                   = false;
  _folded                       = false;
  _phenotype_computed           = false;
  _distance_to_target_computed  = false;
  _fitness_computed             = false;
  _statistical_data_computed    = false;
  _non_coding_computed          = false;
  _modularity_computed          = false;

  _placed_in_population = false;
  _grid_cell = NULL;

  _with_stochasticity = parent->_with_stochasticity;

  // Artificial chemistry parameters
  _w_max = parent->_w_max;

  // Create new genetic units with their DNA copied from here
  // NOTE : The RNA lists (one per genetic unit) will also be copied so that we don't
  // need to look for promoters on the whole genome
  for (auto& gu: parent->_genetic_unit_list)
    _genetic_unit_list.emplace_back(this, &gu);

  _phenotype_activ  = NULL;
  _phenotype_inhib  = NULL;
  _phenotype        = NULL;

  // Initialize all the fitness-related stuff
  _dist_to_target_by_segment  = NULL;
  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];

  for (int8_t i = 0 ; i < NB_FEATURES ; i++)
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }

  // Create a new replication report to store mutational events
  if (_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL)
  {
    _replic_report = new ae_replication_report(this, parent);

    // TODO: remove this after checking it is the old way
    //_exp_m->get_output_m()->get_tree()->set_replic_report(_id, _replic_report);
  }
  else
  {
    _replic_report = NULL;
  }

  // Generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for (int8_t i = 0 ; i < 5 ; i++)
  {
    _int_probes[i]    = parent->_int_probes[i];
    _double_probes[i] = parent->_double_probes[i];
  }

  // Mutation rates etc...
  _mut_params = new ae_params_mut(*(parent->_mut_params));

  // Genome size constraints
  _min_genome_length = parent->_min_genome_length;
  _max_genome_length = parent->_max_genome_length;

  // Plasmids settings
  _allow_plasmids         = parent->_allow_plasmids;

  // Initialize statistical data
  _total_genome_size                  = 0;
  _nb_coding_RNAs                     = 0;
  _nb_non_coding_RNAs                 = 0;
  _overall_size_coding_RNAs           = 0;
  _overall_size_non_coding_RNAs       = 0;
  _nb_genes_activ                     = 0;
  _nb_genes_inhib                     = 0;
  _nb_functional_genes                = 0;
  _nb_non_functional_genes            = 0;
  _overall_size_functional_genes      = 0;
  _overall_size_non_functional_genes  = 0;

  _nb_bases_in_0_CDS                  = -1;
  _nb_bases_in_0_functional_CDS       = -1;
  _nb_bases_in_0_non_functional_CDS   = -1;
  _nb_bases_in_0_RNA                  = -1;
  _nb_bases_in_0_coding_RNA           = -1;
  _nb_bases_in_0_non_coding_RNA       = -1;
  _nb_bases_in_neutral_regions        = -1;
  _nb_neutral_regions                 = -1;

  _modularity = -1;

  //evaluate();
}


// =================================================================
//                             Destructors
// =================================================================
ae_individual::~ae_individual() {
  #ifdef DISTRIBUTED_PRNG
    delete _mut_prng;
    delete _stoch_prng;
  #endif

  delete [] _strain_name;

  // The _replic_report pointer is destroyed, but not the report itself,
  // it will be deleted later, when the tree is written on disk and emptied.

  // Proteins and RNAs are recycled, don't delete them.

  // When the unit is destoyed, its dna is destroyed too, thus the pointer
  // to the ae_dna_replication_report is destroyed. But the
  // dna_replic_report object itself is not deleted, its address is
  // still contained in the global replic_report object in the tree.

  delete _phenotype_activ;
  delete _phenotype_inhib;
  delete _phenotype;

  if (_dist_to_target_by_segment != NULL) delete [] _dist_to_target_by_segment;
  delete [] _dist_to_target_by_feature;

  delete [] _fitness_by_feature;

  // Generic probes
  delete [] _int_probes;
  delete [] _double_probes;

  delete _mut_params;

  /*if(_replic_report!= NULL)
  {
    delete _replic_report;
  }*/
}

// =================================================================
//                        Non-inline Accessors
// =================================================================
void ae_individual::set_exp_m(ae_exp_manager* exp_m) {
  _exp_m = exp_m;

  // Update pointer to exp_manager in each GU
  for (auto& gen_unit: _genetic_unit_list)
    gen_unit.set_exp_m(_exp_m);
}

/// TODO
void ae_individual::set_grid_cell(ae_grid_cell* grid_cell) {
  _placed_in_population = true;
  _grid_cell = grid_cell;
  if (_grid_cell->get_individual() != this)
  {
    _grid_cell->set_individual(this);
  }
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_individual::compute_phenotype() {
  if (_phenotype_computed) return; // Phenotype has already been computed, nothing to do.
  _phenotype_computed = true;

  // Make sure the transcription, translation and folding stages have taken place
  do_transcription_translation_folding();


  // We will use two fuzzy sets :
  //   * _phenotype_activ for the proteins realising a set of functions
  //   * _phenotype_inhib for the proteins inhibiting a set of functions
  // The phenotype will then be given by the sum of these 2 fuzzy sets
  _phenotype_activ = new Fuzzy();
  _phenotype_inhib = new Fuzzy();

  for (const auto& gen_unit: _genetic_unit_list) {
    _phenotype_activ->add(*gen_unit.get_activ_contribution());
    _phenotype_inhib->add(*gen_unit.get_inhib_contribution());
  }

  _phenotype_activ->clip(Fuzzy::max,   Y_MAX);
  _phenotype_inhib->clip(Fuzzy::min, - Y_MAX);

  _phenotype = new ae_phenotype();
  _phenotype->add(*_phenotype_activ);
  _phenotype->add(*_phenotype_inhib);
  _phenotype->clip(Fuzzy::min, Y_MIN);
  _phenotype->simplify();
}

void ae_individual::compute_distance_to_target(Environment* envir) {
// Compute the areas between the phenotype and the environment for each environmental segment
// If the environment is not segmented, the total area is computed
  if (_distance_to_target_computed)
    return; // _distance_to_target has already been computed, nothing to do.

  _distance_to_target_computed = true;

  if (not _phenotype_computed)
    compute_phenotype();

  // Compute the difference between the (whole) phenotype and the environment
  Fuzzy* delta = new Fuzzy(*_phenotype);
  delta->sub(*envir);

  ae_env_segment** segments = envir->get_segments();
  _dist_to_target_by_segment = new double [envir->get_nb_segments()];
  for (size_t i = 0 ; i < envir->get_nb_segments() ; i++)
  {
    _dist_to_target_by_segment[i] = 0;
  }

  // TODO : We should take into account that we compute the areas in order (from the leftmost segment, rightwards)
  //   => We shouldn't parse the whole list of points on the left of the segment we are considering (we have
  //      already been through them!)

  for (size_t i = 0 ; i < envir->get_nb_segments() ; i++)
  {
    _dist_to_target_by_segment[i] = delta->get_geometric_area(segments[i]->start, segments[i]->stop);
    _dist_to_target_by_feature[segments[i]->feature] += _dist_to_target_by_segment[i];
  }

  delete delta;
}

/*!
  Computation of a "proper" fitness value (one that increases when the individual is fitter)

  Computation of a "proper" fitness value (one that increases when the individual is fitter)
  The behaviour of this function depends on many parameters and most notably on whether it is
  a "composite" fitness or not, and on the selection scheme.
*/
void ae_individual::compute_fitness(Environment* envir) {
  if (_fitness_computed) return; // Fitness has already been computed, nothing to do.
  _fitness_computed = true;

  #ifdef NORMALIZED_FITNESS
    for (int8_t i = 0 ; i < NB_FEATURES ; i++)
    {
      if (envir->get_area_by_feature(i)==0.)
      {
        _fitness_by_feature[i] = 0.;
      }
      else
      {
        _fitness_by_feature[i] =  (envir->get_area_by_feature(i) - _dist_to_target_by_feature[i]) / envir->get_area_by_feature(i);
        if ((_fitness_by_feature[i] < 0.) && (i != METABOLISM)) // non-metabolic fitness can NOT be lower than zero (we do not want individual to secrete a negative quantity of public good)
        {
          _fitness_by_feature[i] = 0.;
        }
      }
    }

    if ((! _placed_in_population) || (! _exp_m->get_with_secretion()))
    {
      _fitness = _fitness_by_feature[METABOLISM];
    }
    else
    {
      _fitness =  _fitness_by_feature[METABOLISM] * (1 + _exp_m->get_secretion_contrib_to_fitness() * (_grid_cell->get_compound_amount() - _exp_m->get_secretion_cost() * _fitness_by_feature[SECRETION]));
    }

    if (_exp_m->get_selection_scheme() == FITNESS_PROPORTIONATE) // Then the exponential selection is integrated inside the fitness value
    {
      _fitness = exp(-_exp_m->get_selection_pressure() * (1 - _fitness));
    }
  #else
    for (int8_t i = 0 ; i < NB_FEATURES ; i++)
    {
      if (i == SECRETION)
      {
        _fitness_by_feature[SECRETION] =  exp(- _exp_m->get_selection_pressure() * _dist_to_target_by_feature[SECRETION])
                                        - exp(- _exp_m->get_selection_pressure() * envir->get_area_by_feature(SECRETION));

        if (_fitness_by_feature[i] < 0)
        {
          _fitness_by_feature[i] = 0;
        }
      }
      else
      {
        _fitness_by_feature[i] = exp(- _exp_m->get_selection_pressure() * _dist_to_target_by_feature[i]);
      }
    }

    // Calculate combined, total fitness here!
    // Multiply the contribution of metabolism and the amount of compound in the environment
    if ((! _placed_in_population) || (! _exp_m->get_with_secretion()))
    {
      _fitness =  _fitness_by_feature[METABOLISM] ;
    }
    else
    {
      _fitness =  _fitness_by_feature[METABOLISM]
                  *  (1 + _exp_m->get_secretion_contrib_to_fitness() * _grid_cell->get_compound_amount()
                         - _exp_m->get_secretion_cost() * _fitness_by_feature[SECRETION]);
    }
  #endif
}


void ae_individual::clear_everything_except_dna_and_promoters() {
  _evaluated                    = false;
  _transcribed                  = false;
  _translated                   = false;
  _folded                       = false;
  _phenotype_computed           = false;
  _distance_to_target_computed  = false;
  _fitness_computed             = false;
  _statistical_data_computed    = false;
  _non_coding_computed          = false;
  _modularity_computed          = false;

  for (auto& gen_unit: _genetic_unit_list)
    gen_unit.reset_expression();

  if (_phenotype_activ != NULL)
  {
    delete _phenotype_activ;
    _phenotype_activ = NULL;
  }

  if (_phenotype_inhib != NULL)
  {
    delete _phenotype_inhib;
    _phenotype_inhib = NULL;
  }

  if (_phenotype != NULL)
  {
    delete _phenotype;
    _phenotype = NULL;
  }

  // Initialize all the fitness-related stuff
  if (_dist_to_target_by_segment != NULL)
    {
      delete [] _dist_to_target_by_segment;
      _dist_to_target_by_segment  = NULL;
    }

  for (int8_t i = 0 ; i < NB_FEATURES ; i++)
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }

  // For each RNA / individual / genetic_unit delete proteins it knows
  // Deleting the protein itself is made only once

  for (const auto& gen_unit: _genetic_unit_list) {
    auto rna_list = gen_unit.get_rna_list();
    for (const auto& rna: rna_list[LEADING])
      rna->clear_transcribed_proteins();
    for (const auto& rna: rna_list[LAGGING])
      rna->clear_transcribed_proteins();
    // TODO vld: convert to STL the next 2 lines
    // (gen_unit->get_protein_list()[LEADING])->erase(true);
    // (gen_unit->get_protein_list()[LAGGING])->erase(true);
  }



  // Initialize statistical data
  _total_genome_size                  = 0;
  _nb_coding_RNAs                     = 0;
  _nb_non_coding_RNAs                 = 0;
  _overall_size_coding_RNAs           = 0;
  _overall_size_non_coding_RNAs       = 0;
  _nb_genes_activ                     = 0;
  _nb_genes_inhib                     = 0;
  _nb_functional_genes                = 0;
  _nb_non_functional_genes            = 0;
  _overall_size_functional_genes      = 0;
  _overall_size_non_functional_genes  = 0;

  _nb_bases_in_0_CDS                  = -1;
  _nb_bases_in_0_functional_CDS       = -1;
  _nb_bases_in_0_non_functional_CDS   = -1;
  _nb_bases_in_0_RNA                  = -1;
  _nb_bases_in_0_coding_RNA           = -1;
  _nb_bases_in_0_non_coding_RNA       = -1;
  _nb_bases_in_neutral_regions        = -1;
  _nb_neutral_regions                 = -1;

  _modularity = -1;
}

void ae_individual::reevaluate(Environment* envir /*= NULL*/) {
  // useful for post-treatment programs that replay mutations
  // on a single individual playing the role of the successive
  // ancestors

  clear_everything_except_dna_and_promoters();

  if (envir == NULL) envir = _exp_m->get_env();
  evaluate(envir);
}


void ae_individual::add_GU(char * &sequence, int32_t length) {
  clear_everything_except_dna_and_promoters();
  _genetic_unit_list.emplace_back(this, sequence, length);
}

/// Overloaded version to prevent the use of ae_genetic_unit disabled
/// copy ctor. Forwards arguments to ae_genetic_unit's ctor.
void ae_individual::add_GU(ae_individual* indiv,
                           int32_t chromosome_length,
                           ae_jumping_mt* prng) {
  clear_everything_except_dna_and_promoters();
  _genetic_unit_list.emplace_back(indiv, chromosome_length, prng);
}

void ae_individual::remove_GU(int16_t num_unit)
{
  clear_everything_except_dna_and_promoters();
  auto it = _genetic_unit_list.begin();
  std::advance(it, num_unit);
  _genetic_unit_list.erase(it);
}



void ae_individual::do_transcription() {
  if (_transcribed)
    return; // Transcription has already been performed, nothing to do.
  _transcribed = true;

  for (auto& gen_unit: _genetic_unit_list) {
    gen_unit.do_transcription();
    {
      const auto& rna_list = gen_unit.get_rna_list();
      {
        const auto& lead = rna_list[LEADING];
        _rna_list.insert(_rna_list.end(), lead.begin(), lead.end());
      }
      {
        const auto& lagg = rna_list[LAGGING];
        _rna_list.insert(_rna_list.end(), lagg.begin(), lagg.end());
      }
    }
  }
}

void ae_individual::do_translation() {
  if (_translated)
    return; // ARNs have already been translated, nothing to do.
  _translated = true;

  if (not _transcribed)
    do_transcription();

  for (auto& gen_unit: _genetic_unit_list) {
    gen_unit.do_translation();
    // Create proxies
    const auto& lead = gen_unit.get_protein_list(LEADING);
    const auto& lagg = gen_unit.get_protein_list(LAGGING);

    _protein_list.insert(_protein_list.end(), lead.begin(), lead.end());
    _protein_list.insert(_protein_list.end(), lagg.begin(), lagg.end());
  }
}

void ae_individual::do_folding() {
  if (_folded)
    return; // Proteins have already been folded, nothing to do.
  _folded = true;

  if (not _translated)
    do_translation();

  for (auto& gen_unit: _genetic_unit_list)
    gen_unit.compute_phenotypic_contribution();
}

void ae_individual::evaluate(Environment* envir /*= NULL*/) {
  if (_evaluated == true) return; // Individual has already been evaluated, nothing to do.
  _evaluated = true;

  if (envir == NULL) envir = _exp_m->get_env();

  // ----------------------------------------------------------------------
  // 1) Transcription - Translation - Folding
  // ----------------------------------------------------------------------
  do_transcription_translation_folding();

  // ----------------------------------------------------------------------
  // 2) Compute phenotype and compare it to the environment => fitness
  // ----------------------------------------------------------------------
  compute_phenotype();
  compute_distance_to_target(envir);
  compute_fitness(envir);

  if (_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
  {
    for (auto& gen_unit: _genetic_unit_list) {
      gen_unit.compute_distance_to_target(envir);
      gen_unit.compute_fitness(envir);
    }
  }
}



void ae_individual::inject_GU(ae_individual* donor) {
  // Add the GU at the end of the list
  _genetic_unit_list.emplace_back(this, donor->_genetic_unit_list.back());
}

void ae_individual::inject_2GUs(ae_individual* partner) {
  // We swap GUs from the end of the list.

  // TODO vld: As far as I understood the old code (47b27578), the
  // elements were not swapped but appended to the end of the other GU
  // lists. Error?  Original author (Dule, commit 47b27578), asked for
  // clarification by e-mail on 2015-02-23.

  const auto& gu_list_back_it = std::prev(_genetic_unit_list.end()); // initial last cell from _genetic_unit_list
  _genetic_unit_list.splice(_genetic_unit_list.end(), partner->_genetic_unit_list, std::prev(partner->_genetic_unit_list.end()));
  partner->_genetic_unit_list.splice(partner->_genetic_unit_list.end(), _genetic_unit_list, gu_list_back_it);
}

void ae_individual::compute_statistical_data() {
  if (_statistical_data_computed) return; // Statistical data has already been computed, nothing to do.
  _statistical_data_computed = true;

  if (not _phenotype_computed)
  {
    compute_phenotype();
  }

  for (const auto& gen_unit: _genetic_unit_list) {
    _total_genome_size                 += gen_unit.get_dna()->get_length();
    _nb_coding_RNAs                    += gen_unit.get_nb_coding_RNAs();
    _nb_non_coding_RNAs                += gen_unit.get_nb_non_coding_RNAs();
    _overall_size_coding_RNAs          += gen_unit.get_overall_size_coding_RNAs();
    _overall_size_non_coding_RNAs      += gen_unit.get_overall_size_non_coding_RNAs();
    _nb_genes_activ                    += gen_unit.get_nb_genes_activ();
    _nb_genes_inhib                    += gen_unit.get_nb_genes_inhib();
    _nb_functional_genes               += gen_unit.get_nb_functional_genes();
    _nb_non_functional_genes           += gen_unit.get_nb_non_functional_genes();
    _overall_size_functional_genes     += gen_unit.get_overall_size_functional_genes();
    _overall_size_non_functional_genes += gen_unit.get_overall_size_non_functional_genes();
  }

  if (_replic_report != NULL)
  {
    // Finalize statistical data in the replication report
    _replic_report->signal_end_of_replication();
  }
}

void ae_individual::compute_non_coding() {
  if (_non_coding_computed) return; // _non_coding has already been computed, nothing to do.
  _non_coding_computed = true;

  _nb_bases_in_0_CDS                  = 0;
  _nb_bases_in_0_functional_CDS       = 0;
  _nb_bases_in_0_non_functional_CDS   = 0;
  _nb_bases_in_0_RNA                  = 0;
  _nb_bases_in_0_coding_RNA           = 0;
  _nb_bases_in_0_non_coding_RNA       = 0;
  _nb_bases_in_neutral_regions        = 0;
  _nb_neutral_regions                 = 0;

  for (auto& gen_unit: _genetic_unit_list) {
    _nb_bases_in_0_CDS                  += gen_unit.get_nb_bases_in_0_CDS();
    _nb_bases_in_0_functional_CDS       += gen_unit.get_nb_bases_in_0_functional_CDS();
    _nb_bases_in_0_non_functional_CDS   += gen_unit.get_nb_bases_in_0_non_functional_CDS();
    _nb_bases_in_0_RNA                  += gen_unit.get_nb_bases_in_0_RNA();
    _nb_bases_in_0_coding_RNA           += gen_unit.get_nb_bases_in_0_coding_RNA();
    _nb_bases_in_0_non_coding_RNA       += gen_unit.get_nb_bases_in_0_non_coding_RNA();
    _nb_bases_in_neutral_regions        += gen_unit.get_nb_bases_in_neutral_regions();
    _nb_neutral_regions                 += gen_unit.get_nb_neutral_regions();
  }
}

void ae_individual::save(gzFile backup_file) const {
  // Write the name and "age" of the strain
  int8_t strain_string_len = strlen(_strain_name);
  gzwrite(backup_file, &strain_string_len, sizeof(strain_string_len));
  gzwrite(backup_file, _strain_name, strain_string_len+1);
  gzwrite(backup_file, &_age, sizeof(_age));

  #ifdef DISTRIBUTED_PRNG
    // Write the PRNG's state
    _mut_prng->save(backup_file);
    _stoch_prng->save(backup_file);
  #endif

  // Write id and rank
  gzwrite(backup_file, &_id,   sizeof(_id));
  gzwrite(backup_file, &_rank, sizeof(_rank));

  // Write the position of the individual
  int16_t x = _grid_cell->get_x();
  int16_t y = _grid_cell->get_y();
  gzwrite(backup_file, &x, sizeof(x));
  gzwrite(backup_file, &y, sizeof(y));

  // Write generic probes
  gzwrite(backup_file, _int_probes,    5 * sizeof(*_int_probes));
  gzwrite(backup_file, _double_probes, 5 * sizeof(*_double_probes));

  // Write mutational parameters
  _mut_params->save(backup_file);

  // ------------------------------------------------- Phenotypic stochasticity
  gzwrite(backup_file, &_with_stochasticity, sizeof(_with_stochasticity));

  // Write artificial chemistry parameters
  gzwrite(backup_file, &_w_max, sizeof(_w_max));

  // Write genome size constraints
  gzwrite(backup_file, &_min_genome_length, sizeof(_min_genome_length));
  gzwrite(backup_file, &_max_genome_length, sizeof(_max_genome_length));

  // Write plasmids settings
  int8_t tmp_allow_plasmids = _allow_plasmids;
  gzwrite(backup_file, &tmp_allow_plasmids, sizeof(tmp_allow_plasmids));

  // Write genetic units
  int16_t nb_gen_units = _genetic_unit_list.size();
  gzwrite(backup_file, &nb_gen_units,  sizeof(nb_gen_units));

  for (const auto& gen_unit: _genetic_unit_list)
    gen_unit.save(backup_file);
}

int32_t ae_individual::get_nb_terminators() {
  int32_t nb_term = 0;
  for (auto& gen_unit: _genetic_unit_list)
    nb_term += gen_unit.get_nb_terminators();
  return nb_term;
}


/*!
  \brief Compute reproduction statistics and statistics about the offsprings of the current individual

  * Make nb_children replications of the current individual.
  * For each replication, determine if the offsprings is neutral, beneficial or deleterious by comparison of fitness with the current individual (the parent)
  * If statistics about offsprings are required (offsprings_statistics != NULL), fitness mean, fitness variance, size mean, size variance, functional gene number mean, functional gene number variance fo the nb_children offsprings are computed
  * If information about each children are required (replication_file != NULL), fitness, genome_size, nb of functional genes, number of coding bases, number of transcribed but not translated bases, number of non transcribed bases of each offsprings are written in replication_file

  \param nb_children              number of replications made to have the statistics
  \param reproduction_statistics  statistics about the replications (proportion of neutral offsprings, proportion of beneficial offsprings, proportion of deleterious offsprings)
  \param offsprings_statistics    statistics about the nb_children offsprings (fitness mean, fitness variance, size mean, size variance, functional gene number mean,
                                    functional gene number variance) compute if not null
  \param replication_file         file with information about each children of the current individual (fitness, genome_size, nb of functional genes, number of coding bases,
                                    number of transcribed but not translated bases, number of non transcribed bases) if not null
*/
void ae_individual::compute_experimental_f_nu(int32_t nb_children, double* reproduction_statistics, double* offsprings_statistics, FILE* replication_file) {
  double initial_fitness = get_fitness();

  if(reproduction_statistics != NULL)
  {
    reproduction_statistics[0] = 0; // proportion of neutral offsprings
    reproduction_statistics[1] = 0; // proportion of beneficial offsprings
    reproduction_statistics[2] = 0; // proportion of deleterious offsprings
  }
  else
  {
    printf("%s:%d: error: reproduction_statistics was not initialized\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  if(offsprings_statistics != NULL)
  {
    offsprings_statistics[0] = 0; // offspring fitness mean
    offsprings_statistics[1] = 0; // offspring fitness variance
    offsprings_statistics[2] = 0; // offspring size mean
    offsprings_statistics[3] = 0; // offspring size variance
    offsprings_statistics[4] = 0; // offspring functional gene number mean
    offsprings_statistics[5] = 0; // offspring functional gene number variance
  }

  // ------------------------------------------
  //      Simulate fitness degradation
  // ------------------------------------------

  double fitness_child = 0.0;
  double metabolic_error_child = 0.0;

  // replicate this individual to create 'nb_children' children
  ae_individual * child = NULL;

  int32_t genome_size = 0;
  int32_t nb_functional_genes = 0;
  int32_t nb_bases_in_0_functional_CDS = 0;
  int32_t nb_bases_in_0_coding_RNA = 0;

  for (int i = 0; i < nb_children; i++)
  {
    child = _exp_m->get_exp_s()->get_sel()->do_replication(this, _id);
    fitness_child = child->get_fitness();
    metabolic_error_child = child->get_dist_to_target_by_feature(METABOLISM);

    if (fabs(initial_fitness - fitness_child) < 1e-10*std::max(initial_fitness, fitness_child))
    {
      reproduction_statistics[0] += 1;
    }
    else if (fitness_child > initial_fitness)
    {
      reproduction_statistics[1] += 1;
    }
    else
    {
      reproduction_statistics[2] += 1;
    }

    genome_size = child->get_total_genome_size();
    nb_functional_genes = child->get_nb_functional_genes();
    nb_bases_in_0_functional_CDS = child->get_nb_bases_in_0_functional_CDS();
    nb_bases_in_0_coding_RNA = child->get_nb_bases_in_0_coding_RNA();

    if(offsprings_statistics != NULL)
    {
      offsprings_statistics[0] += fitness_child;
      offsprings_statistics[1] += pow(fitness_child,2);
      offsprings_statistics[2] += (double)genome_size;
      offsprings_statistics[3] += pow((double)genome_size,2);
      offsprings_statistics[4] += (double)nb_functional_genes;
      offsprings_statistics[5] += pow((double)nb_functional_genes,2);
    }

    if(replication_file != NULL)
    {
      fprintf(replication_file, "%le %le %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 "\n", fitness_child, metabolic_error_child, genome_size, nb_functional_genes, genome_size-nb_bases_in_0_functional_CDS, nb_bases_in_0_functional_CDS-nb_bases_in_0_coding_RNA, nb_bases_in_0_coding_RNA);
    }

    delete child;
  }

  //compute Fv
  reproduction_statistics[0] /= (double) nb_children;
  reproduction_statistics[1] /= (double) nb_children;
  reproduction_statistics[2] /= (double) nb_children;

  if(offsprings_statistics != NULL)
  {
    offsprings_statistics[0] /= (double) nb_children;
    offsprings_statistics[1] /= (double) nb_children;
    offsprings_statistics[2] /= (double) nb_children;
    offsprings_statistics[3] /= (double) nb_children;
    offsprings_statistics[4] /= (double) nb_children;
    offsprings_statistics[5] /= (double) nb_children;

    offsprings_statistics[1] -= pow(offsprings_statistics[0],2);
    offsprings_statistics[3] -= pow(offsprings_statistics[2],2);
    offsprings_statistics[5] -= pow(offsprings_statistics[4],2);
  }
}


/// Compute reproduction theoretical proportion of neutral offsprings.
///
/// Compute the theoretical proportion of neutral offsprings given the
/// Carole's formula, based on the mutations and rearrangement rates
/// and not on multiple replications.
///
/// \return theoretical proportion of neutral offsprings
double ae_individual::compute_theoritical_f_nu() {
  // We first have to collect information about genome structure.
  // Abbreviations are chosen according to Carole's formula.
  // Please notice that compared to the formula we have the beginning
  // and ends of neutral regions instead of 'functional regions'
  ae_genetic_unit& chromosome = _genetic_unit_list.front();
  int32_t L       = chromosome.get_dna()->get_length();
  int32_t N_G     = chromosome.get_nb_neutral_regions(); // which is not exactly Carole's original definition
  int32_t* b_i    = chromosome.get_beginning_neutral_regions();
  int32_t* e_i    = chromosome.get_end_neutral_regions();
  int32_t lambda  = chromosome.get_nb_bases_in_neutral_regions();
  int32_t l       = L - lambda; // nb bases in 'functional regions'

  int32_t* lambda_i = NULL;  // nb bases in ith neutral region
  if (N_G > 0) // all the chromosome may be functional
  {
    lambda_i = new int32_t [ N_G ];

    for (int32_t i = 0; i < N_G - 1; i++)
    {
      lambda_i[i] = e_i[i] - b_i[i] + 1;
    }
    if (b_i[N_G-1] > e_i[N_G-1]) // last neutral region is overlapping on the beginning of chromosome
    {
      lambda_i[N_G-1] = (e_i[N_G-1] + L) - b_i[N_G-1] + 1;
    }
    else // no overlap
    {
      lambda_i[N_G-1] = e_i[N_G-1] - b_i[N_G-1] + 1;
    }
  }

  // we now compute the probabilities of neutral reproduction for
  // each type of mutation and rearrangement and update Fv
  double Fv = 1;

  // mutation + insertion + deletion
  double nu_local_mutation = 1 - ((double) l)/L;
  Fv  = pow(1 - get_point_mutation_rate()  * (1 - nu_local_mutation), L);
  Fv *= pow(1 - get_small_insertion_rate() * (1 - nu_local_mutation), L);
  Fv *= pow(1 - get_small_deletion_rate()  * (1 - nu_local_mutation), L);

  // inversion ~ two local mutations
  double nu_inversion = nu_local_mutation * nu_local_mutation;
  Fv *= pow(1 - get_inversion_rate()       * (1 - nu_inversion)     , L);

  // translocation ~ inversion + insertion (mathematically)
  Fv *= pow(1 - get_translocation_rate()   * (1 - nu_inversion * nu_local_mutation), L);

  // long deletion
  double nu_deletion = 0; // if N_G == 0, a deletion is always not neutral
  for (int32_t i = 0; i < N_G; i++) { nu_deletion += lambda_i[i] * (lambda_i[i] + 1); }
  nu_deletion /= ((double) 2*L*L);
  Fv *= pow(1 - get_deletion_rate()        * (1 - nu_deletion)      , L);

  // duplication ~ big deletion + insertion
  Fv *= pow(1 - get_duplication_rate()     * (1 - nu_deletion * nu_local_mutation), L);

  if (lambda_i != NULL) delete [] lambda_i;

  return Fv;
}

/// Remove the bases that are not in coding RNA.
///
/// Remove the bases that are not in coding RNA and test at each loss
/// that fitness is not changed.
void ae_individual::remove_non_coding_bases() {
  for (auto& gen_unit: _genetic_unit_list)
    gen_unit.remove_non_coding_bases();

  _non_coding_computed = false;
  assert(get_nb_bases_in_0_coding_RNA() == 0);
}

/// Double the bases that are not in coding RNA.
///
/// Double the bases that are not in coding RNA by addition of random
/// bases and test at each addition that fitness is not changed.
void ae_individual::double_non_coding_bases(void) {
  _statistical_data_computed = false;
  _non_coding_computed = false;
  _total_genome_size = 0;
  int32_t initial_non_coding_base_nb = get_nb_bases_in_0_coding_RNA();

  for (auto& gen_unit: _genetic_unit_list)
    gen_unit.double_non_coding_bases();

  _statistical_data_computed = false;
  _non_coding_computed = false;
  _total_genome_size = 0;

  assert(get_nb_bases_in_0_coding_RNA()==2*initial_non_coding_base_nb);
}

// =================================================================
//                           Protected Methods
// =================================================================

// TODO vld: refactor make_protein_list and make_rna_list
void ae_individual::make_protein_list() {
  // Clear list
  _protein_list.clear();

  // Make a copy of each genetic unit's protein list
  for (const auto& gen_unit: _genetic_unit_list)
  {
    // Creates proxies
    const auto& lead = gen_unit.get_protein_list(LEADING);
    const auto& lagg = gen_unit.get_protein_list(LAGGING);

    _protein_list.insert(_protein_list.end(), lead.begin(), lead.end());
    _protein_list.insert(_protein_list.end(), lagg.begin(), lagg.end());
  }
}

void ae_individual::make_rna_list() {
  // Clear list
  _rna_list.clear();

  // Make a copy of each genetic unit's rna list
  for (const auto& gen_unit: _genetic_unit_list)
  {
    // Create proxies
    const auto& rna_list = gen_unit.get_rna_list();
    const auto& lead = rna_list[LEADING];
    const auto& lagg = rna_list[LAGGING];

    _rna_list.insert(_rna_list.end(), lead.begin(), lead.end());
    _rna_list.insert(_rna_list.end(), lagg.begin(), lagg.end());
  }
}

} // namespace aevol
