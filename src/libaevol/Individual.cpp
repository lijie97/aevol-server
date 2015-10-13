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

#include "Codon.h"
#include "ExpSetup.h"
#include "ExpManager.h"
#include "GridCell.h"
#include "GeneticUnit.h"
#include "VisAVis.h"
#include "Utils.h"

#ifdef __NO_X
  #ifdef __REGUL
    #include "ae_individual_R.h"
  #else
    #include "Individual.h"
  #endif
#elif defined __X11
  #ifdef __REGUL
    #include "ae_individual_R_X11.h"
  #else
    #include "Individual_X11.h"
  #endif
#endif
namespace aevol {

// ===========================================================================
//                             Constructors
// ===========================================================================
/**
 * // TODO <david.parsons@inria.fr>
 */
Individual::Individual(ExpManager * exp_m,
                       std::shared_ptr<JumpingMT> mut_prng,
                       std::shared_ptr<JumpingMT> stoch_prng,
                       std::shared_ptr<MutationParams> param_mut,
                       double w_max,
                       int32_t min_genome_length,
                       int32_t max_genome_length,
                       bool allow_plasmids,
                       int32_t id,
                       const char* strain_name,
                       int32_t age)
{
  // Experiment manager
  _exp_m = exp_m;

  // PRNGs
  _mut_prng   = mut_prng;
  _stoch_prng = stoch_prng;

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
  _fitness_by_feature = new double [NB_FEATURES];
  for (int i=0; i<NB_FEATURES; i++)
  {
    _dist_to_target_by_feature[i]=0;
    _fitness_by_feature[i]=0;
  }

  _fitness = 0.0;

  // TODO <david.parsons@inria.fr> Deal with cell coordinates on the grid
  // When using structured population, this is the cell the individual is in
  // x = y = -1;

  // The chromosome and plasmids (if allowed)
  // TODO <david.parsons@inria.fr> ???

  // Generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for (int8_t i = 0 ; i < 5 ; i++)
  {
    _int_probes[i]    = 0;
    _double_probes[i] = 0.0;
  }

  // Mutation rates etc...
  _mut_params = param_mut;

  // Artificial chemistry
  _w_max = w_max;

  // Genome size constraints
  _min_genome_length = min_genome_length;
  _max_genome_length = max_genome_length;

  // Plasmids settings
  _allow_plasmids = allow_plasmids;


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
  _placed_in_population         = false;



  // ----------------------------------------
  // Statistical data
  // ----------------------------------------
  _modularity = 0.0;
}

/**
 * This constructor retrieves an individual from a backup file.
 *
 * Since this generation has already been processed, no unnecessary calculation
 * (e.g. fitness) will be done.
 * No transcription, translation or other process of that kind is performed.
 */
Individual::Individual(ExpManager* exp_m, gzFile backup_file)
{
  _exp_m = exp_m;

  // Retrieve the name and "age" of the strain
  int8_t strain_string_len;
  gzread(backup_file, &strain_string_len, sizeof(strain_string_len));
  _strain_name = new char[strain_string_len+1];
  gzread(backup_file, _strain_name, strain_string_len+1);
  gzread(backup_file, &_age, sizeof(_age));

  // Retrieve the PRNGs
  if (exp_m == NULL)
  {
    // Detached mode
    _mut_prng   = NULL;
    _stoch_prng = NULL;
  }
  else
  {
    // TODO: => prngs as parameters
    _mut_prng   = exp_m->world()->get_mut_prng();
    _stoch_prng = exp_m->world()->get_stoch_prng();
    assert(_mut_prng);
    assert(_stoch_prng);
  }

  // Retrieve id and rank
  gzread(backup_file, &_id,    sizeof(_id));
  gzread(backup_file, &_rank,  sizeof(_rank));

  // Retrieve spatial coordinates
  // gzread(backup_file, &x, sizeof(x));
  // gzread(backup_file, &y, sizeof(y));
  _placed_in_population = false;

  // Retrieve generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  gzread(backup_file, _int_probes,     5 * sizeof(*_int_probes));
  gzread(backup_file, _double_probes,  5 * sizeof(*_double_probes));

  // Retrieve mutational parameters
  _mut_params = std::make_shared<MutationParams>(backup_file);

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

  // Retrieve genetic units
  int16_t nb_gen_units;
  gzread(backup_file, &nb_gen_units,  sizeof(nb_gen_units));

  for (int16_t i = 0 ; i < nb_gen_units ; i++)
    _genetic_unit_list.emplace_back(this, backup_file);

  // --------------------------------------------------------------------------
  // No more data to retrieve, the following are only structural
  // initializations (no data is set)
  // --------------------------------------------------------------------------

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

  // Initialize the computational state of the individual
  _evaluated                    = false;
  _transcribed                  = false;
  _translated                   = false;
  _folded                       = false;
  _phenotype_computed           = false;
  _distance_to_target_computed  = false;
  _fitness_computed             = false;

  _modularity = -1;
}

/**
 * Copy constructor
 */
Individual::Individual(const Individual& other) {
  _exp_m = other._exp_m;

  // PRNGs
  _mut_prng   = other._mut_prng;
  _stoch_prng = other._stoch_prng;

  int strain_string_len = strlen(other._strain_name);
  _strain_name = new char[strain_string_len+1];
  memcpy(_strain_name, other._strain_name, strain_string_len+1);
  _age  = other._age;

  _id   = other._id;
  _rank = other._rank;

  _evaluated                    = other._evaluated;
  _transcribed                  = other._transcribed;
  _translated                   = other._translated;
  _folded                       = other._folded;
  _phenotype_computed           = other._phenotype_computed;

  _with_stochasticity = other._with_stochasticity;

  // Artificial chemistry parameters
  _w_max = other._w_max;

  _distance_to_target_computed  = other._distance_to_target_computed;
  _fitness_computed             = other._fitness_computed;
  _placed_in_population         = other._placed_in_population;

  // Copy genetic units from other
  // Should actually use GeneticUnit copy ctor which is disabled.
  for (const auto& gu: other._genetic_unit_list)
    _genetic_unit_list.emplace_back(this, gu);

  // Copy phenotype
  if (_phenotype_computed) {
    _phenotype_activ  = new Fuzzy(*(other._phenotype_activ));
    _phenotype_inhib  = new Fuzzy(*(other._phenotype_inhib));
    _phenotype        = new Phenotype(*(other._phenotype));
  }
  else {
    _phenotype_activ  = NULL;
    _phenotype_inhib  = NULL;
    _phenotype        = NULL;
  }


  // Copy fitness-related stuff
  _dist_to_target_by_segment  = NULL;
  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];

  for (int8_t i = 0 ; i < NB_FEATURES ; i++) {
    _dist_to_target_by_feature[i] = other._dist_to_target_by_feature[i];
    _fitness_by_feature[i] = other._fitness_by_feature[i];
  }

  _fitness = other._fitness;

  // Copy statistical data
  metrics_ = other.metrics_ ?
           new Metrics(*other.metrics_) :
           nullptr;
  nc_metrics_ = other.nc_metrics_ ?
              new NonCodingMetrics(*other.nc_metrics_) :
              nullptr;

  _modularity = other._modularity;

  // Generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for (int8_t i = 0 ; i < 5 ; i++) {
    _int_probes[i]    = other._int_probes[i];
    _double_probes[i] = other._double_probes[i];
  }

  // Mutation rates etc...
  _mut_params = std::make_shared<MutationParams>(*(other._mut_params));

  // Genome size constraints
  _min_genome_length = other._min_genome_length;
  _max_genome_length = other._max_genome_length;

  // Plasmids settings
  _allow_plasmids = other._allow_plasmids;
}

/**
 * Reproduction constructor
 *
 * This constructor creates a new individual with the same genome as it's
 * parent. The location of promoters will be copied but no further process will
 * be performed.
 *
 * The phenotype and the fitness are not set, neither is the statistical data.
*/
Individual::Individual(const Individual * parent, int32_t id,
                             std::shared_ptr<JumpingMT> mut_prng,
                             std::shared_ptr<JumpingMT> stoch_prng)
{
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

  _placed_in_population = false;
  // x = y = -1;

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

  // Generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for (int8_t i = 0 ; i < 5 ; i++)
  {
    _int_probes[i]    = parent->_int_probes[i];
    _double_probes[i] = parent->_double_probes[i];
  }

  // Mutation rates etc...
  _mut_params = std::make_shared<MutationParams>(*(parent->_mut_params));

  // Genome size constraints
  _min_genome_length = parent->_min_genome_length;
  _max_genome_length = parent->_max_genome_length;

  // Plasmids settings
  _allow_plasmids = parent->_allow_plasmids;

  // Initialize statistical data
  _modularity = -1;
}

Individual* Individual::CreateIndividual(ExpManager * exp_m,
                                         gzFile backup_file) {
  Individual* indiv = NULL;
#ifdef __NO_X
    #ifndef __REGUL
      indiv = new Individual(exp_m, backup_file);
    #else
      indiv = new ae_individual_R(exp_m, backup_file);
    #endif
  #elif defined __X11
    #ifndef __REGUL
      indiv = new Individual_X11(exp_m, backup_file);
    #else
      indiv = new ae_individual_R_X11(exp_m, backup_file);
    #endif
  #endif

  return indiv;
}

/*!
  \brief Create a clone

  \param dolly original individual to be cloned
  \param id ID of the clone
  \return clone of dolly (not evaluated)
*/
Individual* Individual::CreateClone(const Individual* dolly, int32_t id) {
  Individual* indiv = new Individual(*dolly);
  indiv->set_id(id);
  return indiv;
}


// =================================================================
//                             Destructor
// =================================================================
Individual::~Individual()
{
  delete [] _strain_name;

  // Proteins and RNAs are recycled, don't delete them.

  delete _phenotype_activ;
  delete _phenotype_inhib;
  delete _phenotype;

  if (_dist_to_target_by_segment != NULL) delete [] _dist_to_target_by_segment;
  delete [] _dist_to_target_by_feature;

  delete [] _fitness_by_feature;

  // Generic probes
  delete [] _int_probes;
  delete [] _double_probes;

  delete metrics_;
  delete nc_metrics_;
}

// =================================================================
//                        Non-inline Accessors
// =================================================================
void Individual::set_exp_m(ExpManager * exp_m) {
  _exp_m = exp_m;

  // Update pointer to exp_manager in each GU
  for (auto& gen_unit: _genetic_unit_list)
    gen_unit.set_exp_m(_exp_m);
}

/// TODO
void Individual::set_grid_cell(GridCell * grid_cell)
{
  _grid_cell = grid_cell;
  _placed_in_population = true;
  // x = grid_cell->get_x();
  // y = grid_cell->get_y();
  if (grid_cell->get_individual() != this)
    grid_cell->set_individual(this);
}

/// TODO
const char*Individual::get_strain_name() const {
  return _strain_name;
}

/// TODO
int32_t Individual::get_age() const {
  return _age;
}

/// TODO
int32_t Individual::get_id() const {
  return _id;
}

/// TODO
double*Individual::get_dist_to_target_by_segment() const {
  return _dist_to_target_by_segment;
}

/*!
  Get the individual's rank in the population (1 for the worst indiv, POP_SIZE for the best)

  Warning: be sure you call sort_individuals() before using get_rank_in_population
*/
int32_t Individual::get_rank() const {
  return _rank;
}

/// TODO
ExpManager *Individual::get_exp_m() const {
  return _exp_m;
}

/// TODO
std::shared_ptr<JumpingMT> Individual::get_mut_prng() const {
  return _mut_prng;
}

/// TODO
std::shared_ptr<JumpingMT> Individual::get_stoch_prng() const {
  return _stoch_prng;
}

/*!
  Returns the number of genetic units
*/
int16_t Individual::get_nb_genetic_units() const {
  return _genetic_unit_list.size();
}

/// Get the number of plasmids. That is, the number of genetic units
/// minus one DNA-based genetic unit.
int32_t Individual::get_nb_plasmids() const {
  return _genetic_unit_list.size() - 1;
}

/// TODO
int32_t Individual::get_amount_of_dna() const {
  int32_t amount = 0;
  for (const auto& gen_unit: _genetic_unit_list)
    amount += gen_unit.get_dna()->length();
  return amount;
}

/// Return the list of genetic units.
const std::list<GeneticUnit>&Individual::get_genetic_unit_list() const {
  return _genetic_unit_list;
}

std::list<GeneticUnit>&Individual::get_genetic_unit_list_nonconst() {
  return _genetic_unit_list;
}

/// Remove all the elements from the GU list except the firt and the
/// last ones. If the GU list has less that 2 elements, do nothing.
void Individual::drop_nested_genetic_units() {
  if (_genetic_unit_list.size() <= 2)
    return;

  _genetic_unit_list.erase(std::next(_genetic_unit_list.begin()),
                           std::prev(_genetic_unit_list.end()));
  assert(_genetic_unit_list.size() <= 2);
}

/// Returns genetic unit number `num_unit` (0 for main chromosome)
const GeneticUnit& Individual::get_genetic_unit(int16_t num_unit) const {
  assert(num_unit < static_cast<int32_t>(_genetic_unit_list.size()));
  auto it = _genetic_unit_list.cbegin();
  std::advance(it, num_unit);
  return *it;
}

/// Returns genetic unit number `num_unit` (0 for main chromosome) as
/// a non-constant reference. To be used when the purpose if to alter
/// the individual.
GeneticUnit& Individual::get_genetic_unit_nonconst(int16_t num_unit) {
  assert(num_unit < static_cast<int32_t>(_genetic_unit_list.size()));
  auto it = _genetic_unit_list.begin();
  std::advance(it, num_unit);
  return *it;
}

/// TODO
double Individual::get_dist_to_target_by_feature(PhenotypicFeature feature) const {
  assert(_distance_to_target_computed);

  return _dist_to_target_by_feature[feature];
}

/// TODO
double Individual::get_fitness() const {
  assert(_fitness_computed);

  return _fitness;
}

/// TODO
double Individual::get_fitness_by_feature(PhenotypicFeature feature) const {
  assert(_fitness_computed);

  return _fitness_by_feature[feature];
}

/// TODO
GridCell *Individual::get_grid_cell() const {
  return _grid_cell;
}

/// TODO
const Habitat& Individual::habitat() const {
  return _grid_cell->habitat();
}

/// TODO
bool Individual::get_placed_in_population () const {
  return _placed_in_population;
}

/*!
  Returns the sequence of genetic unit number <num_unit> (0 for main chromosome)
*/
const char*Individual::get_genetic_unit_sequence(int16_t num_unit) const {
  return get_genetic_unit(num_unit).get_sequence();
}

/*!
  Returns the sequence length of genetic unit number <num_unit> (0 for main chromosome)
*/
int32_t Individual::get_genetic_unit_seq_length(int16_t num_unit) const {
  return get_genetic_unit(num_unit).get_seq_length();
}

/// TODO
Fuzzy* Individual::get_phenotype_activ() const {
  return _phenotype_activ;
}

/// TODO
Fuzzy* Individual::get_phenotype_inhib() const {
  return _phenotype_inhib;
}

/// TODO
Phenotype* Individual::get_phenotype() const {
  return _phenotype;
}

const PhenotypicTarget& Individual::phenotypic_target() const {
  return _grid_cell->phenotypic_target();
}

/// TODO
const std::list<Protein*>& Individual::get_protein_list() const {
  return _protein_list;
}

/// TODO
const std::list<const Rna*>& Individual::get_rna_list() const {
  return _rna_list;
}

/// TODO
int32_t Individual::get_total_genome_size() const {
  assert(metrics_ != nullptr);
  return metrics_->total_genome_size();
}

/// TODO
int16_t Individual::get_nb_coding_RNAs() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_coding_RNAs();
}

/// TODO
int16_t Individual::get_nb_non_coding_RNAs() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_non_coding_RNAs();
}

/// TODO
int32_t Individual::get_overall_size_coding_RNAs() const {
  assert(metrics_ != nullptr);
  return metrics_->overall_size_coding_RNAs();
}

/// TODO
double Individual::get_av_size_coding_RNAs() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_coding_RNAs() == 0 ?
         0.0 :
         metrics_->overall_size_coding_RNAs() /
         metrics_->nb_coding_RNAs();
}

/// TODO
int32_t Individual::get_overall_size_non_coding_RNAs() const {
  assert(metrics_ != nullptr);
  return metrics_->overall_size_non_coding_RNAs();
}

/// TODO
double Individual::get_av_size_non_coding_RNAs() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_non_coding_RNAs() == 0 ?
         0.0 :
         metrics_->overall_size_non_coding_RNAs() /
         metrics_->nb_non_coding_RNAs();
}

/// TODO
int16_t Individual::get_nb_genes_activ() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_genes_activ();
}

/// TODO
int16_t Individual::get_nb_genes_inhib() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_genes_inhib();
}

/// TODO
int16_t Individual::get_nb_functional_genes() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_functional_genes();
}

/// TODO
int16_t Individual::get_nb_non_functional_genes() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_non_functional_genes();
}

/// TODO
int32_t Individual::get_overall_size_functional_genes() const {
  assert(metrics_ != nullptr);
  return metrics_->overall_size_functional_genes();
}

/// TODO
double Individual::get_av_size_functional_genes() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_functional_genes() == 0 ?
         0.0 :
         metrics_->overall_size_functional_genes() /
         metrics_->nb_functional_genes();
}

/// TODO
int32_t Individual::get_overall_size_non_functional_genes() const {
  assert(metrics_ != nullptr);
  return metrics_->overall_size_non_functional_genes();
}

/// TODO
double Individual::get_av_size_non_functional_genes() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_non_functional_genes() == 0 ?
         0.0 :
         metrics_->overall_size_non_functional_genes() /
         metrics_->nb_non_functional_genes();
}

/// TODO
int32_t Individual::get_nb_bases_in_0_CDS() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_0_CDS();
}

/// TODO
int32_t Individual::get_nb_bases_in_0_functional_CDS() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_0_functional_CDS();
}

/// TODO
int32_t Individual::get_nb_bases_in_0_non_functional_CDS() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_0_non_functional_CDS();
}

/// TODO
int32_t Individual::get_nb_bases_in_0_RNA() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_0_RNA();
}

/// TODO
int32_t Individual::get_nb_bases_in_0_coding_RNA() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_0_coding_RNA();
}

/// TODO
int32_t Individual::get_nb_bases_in_0_non_coding_RNA() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_0_non_coding_RNA();
}

/// TODO
int32_t Individual::get_nb_bases_in_neutral_regions() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_neutral_regions();
}

/// TODO
int32_t Individual::get_nb_neutral_regions() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_neutral_regions();
}

/// TODO
double Individual::get_modularity() {
  printf("\n  WARNING : modularity measure not yet implemented.\n");
  //~ if (modularity_ < 0) compute_modularity();
  //~ return modularity_;
  return 0;
}

/// TODO
double Individual::get_w_max() const {
  return _w_max;
}

// ------------------------------------------------------------- Mutation rates
/// TODO
double Individual::get_point_mutation_rate() const {
  return _mut_params->get_point_mutation_rate();
}

/// TODO
double Individual::get_small_insertion_rate() const {
  return _mut_params->get_small_insertion_rate();
}

/// TODO
double Individual::get_small_deletion_rate() const {
  return _mut_params->get_small_deletion_rate();
}

/// TODO
int16_t Individual::get_max_indel_size() const {
  return _mut_params->get_max_indel_size();
}

// ---------------------------------- Rearrangement rates (without alignements)
/// TODO
double Individual::get_duplication_rate() const {
  return _mut_params->get_duplication_rate();
}

/// TODO
double Individual::get_deletion_rate() const {
  return _mut_params->get_deletion_rate();
}

/// TODO
double Individual::get_translocation_rate() const {
  return _mut_params->get_translocation_rate();
}

/// TODO
double Individual::get_inversion_rate() const {
  return _mut_params->get_inversion_rate();
}

// ------------------------------------- Rearrangement rates (with alignements)
/// TODO
double Individual::get_neighbourhood_rate() const {
  return _mut_params->get_neighbourhood_rate();
}

/// TODO
double Individual::get_duplication_proportion() const {
  return _mut_params->get_duplication_proportion();
}

/// TODO
double Individual::get_deletion_proportion() const {
  return _mut_params->get_deletion_proportion();
}

/// TODO
double Individual::get_translocation_proportion() const {
  return _mut_params->get_translocation_proportion();
}

/// TODO
double Individual::get_inversion_proportion() const {
  return _mut_params->get_inversion_proportion();
}

// ---------------------------------------------------------------- Transfer
bool Individual::get_with_4pts_trans() const {
  return _mut_params->get_with_4pts_trans();
}

bool Individual::get_with_HT() const {
  return _mut_params->get_with_HT();
}

bool Individual::get_repl_HT_with_close_points() const {
  return _mut_params->get_repl_HT_with_close_points();
}

double Individual::get_HT_ins_rate() const {
  return _mut_params->get_HT_ins_rate();
}

double Individual::get_HT_repl_rate() const {
  return _mut_params->get_HT_repl_rate();
}

double Individual::get_repl_HT_detach_rate() const {
  return _mut_params->get_repl_HT_detach_rate();
}


// ---------------------------------------------------------------- Alignements
bool Individual::get_with_alignments() const {
  return _mut_params->get_with_alignments();
}

AlignmentFunctionShape Individual::get_align_fun_shape() const {
  return _mut_params->get_align_fun_shape();
}

double Individual::get_align_sigm_lambda() const {
  return _mut_params->get_align_sigm_lambda();
}

int16_t Individual::get_align_sigm_mean() const {
  return _mut_params->get_align_sigm_mean();
}

int16_t Individual::get_align_lin_min() const {
  return _mut_params->get_align_lin_min();
}

int16_t Individual::get_align_lin_max() const {
  return _mut_params->get_align_lin_max();
}

int16_t Individual::get_align_max_shift() const {
  return _mut_params->get_align_max_shift();
}

int16_t Individual::get_align_w_zone_h_len() const {
  return _mut_params->get_align_w_zone_h_len();
}

int16_t Individual::get_align_match_bonus() const {
  return _mut_params->get_align_match_bonus();
}

int16_t Individual::get_align_mismatch_cost() const {
  return _mut_params->get_align_mismatch_cost();
}

/// TODO
bool Individual::get_with_stochasticity() const {
  return _with_stochasticity;
}

void Individual::set_allow_plasmids(bool allow_plasmids) {
  _allow_plasmids=allow_plasmids;
}

// Genome size constraints
/// TODO
int32_t Individual::get_min_genome_length() const {
  return _min_genome_length;
}

/// TODO
int32_t Individual::get_max_genome_length() const {
  return _max_genome_length;
}

// Plasmids settings
/// TODO
bool Individual::get_allow_plasmids() const {
  return _allow_plasmids;
}

/*!
  \brief Return the _int_probes

  \return _int_probes
*/
int32_t*Individual::get_int_probes () const {
  return _int_probes;
}

/*!
  \brief Return the _double_probes

  \return _double_probes
*/
double*Individual::get_double_probes (void) const {
  return _double_probes;
}


// =====================================================================
//                           Setters' definitions
// =====================================================================
void Individual::set_strain_name(char* name) {
  assert(name && strlen(name) < INT8_MAX); // Conservative, could be <=
  int8_t name_len = strlen(name);
  delete [] _strain_name;
  _strain_name = new char[name_len+1];
  memcpy(_strain_name, name, name_len+1);
}

/// TODO
void Individual::set_id(int32_t id) {
  _id = id;
}

/// TODO
void Individual::set_rank(int32_t rank) {
  _rank = rank;
}

/// TODO
void Individual::set_placed_in_population(bool placed_in_population) {
  _placed_in_population = placed_in_population;
}

/// TODO
void Individual::set_w_max(double w_max) {
  _w_max = w_max;
}



// Genome size constraints
/// TODO
void Individual::set_min_genome_length(int32_t min_genome_length) {
  _min_genome_length = min_genome_length;
}

/// TODO
void Individual::set_max_genome_length(int32_t max_genome_length) {
  _max_genome_length = max_genome_length;
}


void Individual::set_point_mutation_rate(double point_mutation_rate) {
  _mut_params->set_point_mutation_rate(point_mutation_rate);
}

void Individual::set_small_insertion_rate(double small_insertion_rate) {
  _mut_params->set_small_insertion_rate(small_insertion_rate);
}

void Individual::set_small_deletion_rate(double small_deletion_rate) {
  _mut_params->set_small_deletion_rate(small_deletion_rate);
}

void Individual::set_max_indel_size(int16_t max_indel_size) {
  _mut_params->set_max_indel_size(max_indel_size);
}

void Individual::set_duplication_rate(double duplication_rate) {
  _mut_params->set_duplication_rate(duplication_rate);
}

void Individual::set_deletion_rate(double deletion_rate) {
  _mut_params->set_deletion_rate(deletion_rate);
}

void Individual::set_translocation_rate(double translocation_rate) {
  _mut_params->set_translocation_rate(translocation_rate);
}

void Individual::set_inversion_rate(double inversion_rate) {
  _mut_params->set_inversion_rate(inversion_rate);
}

void Individual::set_neighbourhood_rate(double neighbourhood_rate) {
  _mut_params->set_neighbourhood_rate(neighbourhood_rate);
}

void Individual::set_duplication_proportion(double duplication_proportion) {
  _mut_params->set_duplication_proportion(duplication_proportion);
}

void Individual::set_deletion_proportion(double deletion_proportion) {
  _mut_params->set_deletion_proportion(deletion_proportion);
}

void Individual::set_translocation_proportion(double translocation_proportion) {
  _mut_params->set_translocation_proportion(translocation_proportion);
}

void Individual::set_inversion_proportion(double inversion_proportion) {
  _mut_params->set_inversion_proportion(inversion_proportion);
}

void Individual::set_with_4pts_trans(bool with_4pts_trans) {
  _mut_params->set_with_4pts_trans(with_4pts_trans);
}

void Individual::set_with_alignments(bool with_alignments) {
  _mut_params->set_with_alignments(with_alignments);
}

void Individual::set_with_HT(bool with_HT) {
  _mut_params->set_with_HT(with_HT);
}

void Individual::set_repl_HT_with_close_points(bool repl_HT_with_close_points) {
  _mut_params->set_repl_HT_with_close_points(repl_HT_with_close_points);
}

void Individual::set_HT_ins_rate(double HT_ins_rate) {
  _mut_params->set_HT_ins_rate(HT_ins_rate);
}

void Individual::set_HT_repl_rate(double HT_repl_rate) {
  _mut_params->set_HT_repl_rate(HT_repl_rate);
}

void Individual::set_repl_HT_detach_rate(double repl_HT_detach_rate) {
  _mut_params->set_repl_HT_detach_rate(repl_HT_detach_rate);
}


void Individual::set_with_stochasticity(bool with_stoch) {
  _with_stochasticity = with_stoch;
}

void Individual::set_stoch_prng(std::shared_ptr<JumpingMT> prng) {
  _stoch_prng = prng;
}

void Individual::set_mut_prng(std::shared_ptr<JumpingMT> prng) {
  _mut_prng = prng;
}

/*!
  \brief Change the _int_probes

  \param int_probes 5 int32_t* that constitute a probe
*/
void Individual::set_int_probes (int32_t* int_probes) {
  _int_probes = int_probes;
}

/*!
  \brief Change the _double_probes

  \param double_probes 5 double* that constitute a probe
*/
void Individual::set_double_probes (double* double_probes) {
  _double_probes = double_probes;
}

// =====================================================================
//                       functions' definition
// =====================================================================

void Individual::reset_dist_to_target_by_segment(double* dist_to_target_by_segment) {
  if (_dist_to_target_by_segment != NULL) delete [] _dist_to_target_by_segment;
  _dist_to_target_by_segment = dist_to_target_by_segment;
}

void Individual::renew_dist_to_target_by_feature() {
  if (_dist_to_target_by_feature != NULL) delete [] _dist_to_target_by_feature;
  _dist_to_target_by_feature = new double [NB_FEATURES];
}


void Individual::renew_fitness_by_feature() {
  if (_fitness_by_feature != NULL) delete [] _fitness_by_feature;
  _fitness_by_feature = new double [NB_FEATURES];
}

void Individual::do_transcription_translation_folding() {
  if (_transcribed == true && _translated == true && _folded == true) return;

  do_transcription();
  do_translation();
  do_folding();

  make_protein_list();
}

#ifdef DEBUG
  void Individual::assert_promoters()
  {
    // Perform assertion for each genetic unit
    // for (auto& gen_unit: _genetic_unit_list)
    //   gen_unit.assert_promoters();
  }

  void Individual::assert_promoters_order()
  {
    // Perform assertion for each genetic unit
    // for (auto& gen_unit: _genetic_unit_list)
    //   gen_unit.assert_promoters_order();
  }
#endif

// =================================================================
//                            Public Methods
// =================================================================
void Individual::compute_phenotype() {
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

  _phenotype = new Phenotype();
  _phenotype->add(*_phenotype_activ);
  _phenotype->add(*_phenotype_inhib);
  _phenotype->clip(Fuzzy::min, Y_MIN);
  _phenotype->simplify();
}

void Individual::compute_distance_to_target(const PhenotypicTarget& target) {
// Compute the areas between the phenotype and the target for each segment
// If the target is not segmented, the total area is computed
  if (_distance_to_target_computed)
    return; // _distance_to_target has already been computed, nothing to do.

  _distance_to_target_computed = true;

  if (not _phenotype_computed)
    compute_phenotype();

  // Compute the difference between the (whole) phenotype and the target
  Fuzzy* delta = new Fuzzy(*_phenotype);
  delta->sub(target);

  PhenotypicSegment ** segments = target.segments();
  _dist_to_target_by_segment = new double [target.nb_segments()];
  for (size_t i = 0 ; i < static_cast<size_t>(target.nb_segments()) ; i++) {
    _dist_to_target_by_segment[i] = 0;
  }

  // TODO : We should take into account that we compute the areas in order (from the leftmost segment, rightwards)
  //   => We shouldn't parse the whole list of points on the left of the segment we are considering (we have
  //      already been through them!)

  for (size_t i = 0 ; i < static_cast<size_t>(target.nb_segments()) ; i++) {
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
void Individual::compute_fitness(const PhenotypicTarget& target) {
  if (_fitness_computed) return; // Fitness has already been computed, nothing to do.
  _fitness_computed = true;

#ifdef NORMALIZED_FITNESS
    for (int8_t i = 0 ; i < NB_FEATURES ; i++)
    {
      if (envir->get_area_by_feature(i)==0.)
      {
        fitness_by_feature_[i] = 0.;
      }
      else
      {
        fitness_by_feature_[i] =  (envir->get_area_by_feature(i) - dist_to_target_by_feature_[i]) / envir->get_area_by_feature(i);
        if ((fitness_by_feature_[i] < 0.) && (i != METABOLISM)) // non-metabolic fitness can NOT be lower than zero (we do not want individual to secrete a negative quantity of public good)
        {
          fitness_by_feature_[i] = 0.;
        }
      }
    }

    if ((! _placed_in_population) || (! exp_m_->get_with_secretion()))
    {
      fitness_ = fitness_by_feature_[METABOLISM];
    }
    else
    {
      fitness_ =  fitness_by_feature_[METABOLISM] * (1 + exp_m_->get_secretion_contrib_to_fitness() * (_grid_cell->compound_amount() - exp_m_->get_secretion_cost() * fitness_by_feature_[SECRETION]));
    }

    if (exp_m_->get_selection_scheme() == FITNESS_PROPORTIONATE) // Then the exponential selection is integrated inside the fitness value
    {
      fitness_ = exp(-exp_m_->get_selection_pressure() * (1 - fitness_));
    }
  #else
  for (int8_t i = 0 ; i < NB_FEATURES ; i++)
  {
    if (i == SECRETION)
    {
      _fitness_by_feature[SECRETION] =  exp(- _exp_m->get_selection_pressure() * _dist_to_target_by_feature[SECRETION])
                                        - exp(- _exp_m->get_selection_pressure() * target.area_by_feature(SECRETION));

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
  // Multiply the contribution of metabolism and the amount of compound in the
  // habitat
  if ((! _placed_in_population) || (! _exp_m->get_with_secretion()))
  {
    _fitness =  _fitness_by_feature[METABOLISM] ;
  }
  else
  {
    _fitness =  _fitness_by_feature[METABOLISM]
                *  (1 + _exp_m->get_secretion_contrib_to_fitness() * get_grid_cell()->compound_amount()
                    - _exp_m->get_secretion_cost() * _fitness_by_feature[SECRETION]);
  }
#endif
}


void Individual::clear_everything_except_dna_and_promoters() {
  _evaluated                    = false;
  _transcribed                  = false;
  _translated                   = false;
  _folded                       = false;
  _phenotype_computed           = false;
  _distance_to_target_computed  = false;
  _fitness_computed             = false;

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

  for (auto& gen_unit: _genetic_unit_list)
    gen_unit.clear_transcribed_proteins();

  // Clear RNA and proteins
  _rna_list.clear();
  _protein_list.clear();

  // Reset statistical data
  delete metrics_;
  metrics_ = nullptr;
  delete nc_metrics_;
  nc_metrics_ = nullptr;

  _modularity = -1;
}

void Individual::reevaluate() {
  // useful for post-treatment programs that replay mutations
  // on a single individual playing the role of the successive
  // ancestors

  clear_everything_except_dna_and_promoters();
  Evaluate();
}


void Individual::add_GU(char * &sequence, int32_t length) {
  clear_everything_except_dna_and_promoters();
  _genetic_unit_list.emplace_back(this, sequence, length);
}

/// Overloaded version to prevent the use of GeneticUnit disabled
/// copy ctor. Forwards arguments to GeneticUnit's ctor.
void Individual::add_GU(Individual * indiv,
                           int32_t chromosome_length,
                           std::shared_ptr<JumpingMT> prng) {
  clear_everything_except_dna_and_promoters();
  _genetic_unit_list.emplace_back(indiv, chromosome_length, prng);
}

void Individual::remove_GU(int16_t num_unit)
{
  clear_everything_except_dna_and_promoters();
  auto it = _genetic_unit_list.begin();
  std::advance(it, num_unit);
  _genetic_unit_list.erase(it);
}



void Individual::do_transcription() {
  if (_transcribed)
    return; // Transcription has already been performed, nothing to do.
  _transcribed = true;

  for (auto& gen_unit: _genetic_unit_list) {
    gen_unit.do_transcription();
    const auto& rna_list = gen_unit.get_rna_list();
    for (auto& strand: {LEADING, LAGGING}) {
      for (auto& rna: rna_list[strand])
        _rna_list.push_back(&rna);
    }
  }
}

void Individual::do_translation() {
  if (_translated)
    return; // ARNs have already been translated, nothing to do.
  _translated = true;

  if (not _transcribed)
    do_transcription();

  for (auto& gen_unit: _genetic_unit_list) {
    gen_unit.do_translation();
    // append all proteins from `gen_unit` to `protein_list_`
    for (auto& strand_id: {LEADING, LAGGING}) {
      auto& strand = gen_unit.get_protein_list(strand_id);
      for (auto& p: strand)
        _protein_list.push_back(&p);
    }
  }
}

void Individual::do_folding() {
  if (_folded)
    return; // Proteins have already been folded, nothing to do.
  _folded = true;

  if (not _translated)
    do_translation();

  for (auto& gen_unit: _genetic_unit_list)
    gen_unit.compute_phenotypic_contribution();
}

void Individual::Evaluate() {
  EvaluateInContext(_grid_cell->habitat());
}

void Individual::EvaluateInContext(const Habitat& habitat) {
  if (_evaluated == true) return; // Individual has already been evaluated, nothing to do.
  _evaluated = true;

  // ----------------------------------------------------------------------
  // Transcription - Translation - Folding
  // ----------------------------------------------------------------------
  do_transcription_translation_folding();

  // ----------------------------------------------------------------------
  // Compute phenotype and compare it to the target => fitness
  // ----------------------------------------------------------------------
  compute_phenotype();
  compute_distance_to_target(habitat.phenotypic_target());
  compute_fitness(habitat.phenotypic_target());

  if (_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
    for (auto& gen_unit: _genetic_unit_list) {
      gen_unit.compute_distance_to_target(habitat.phenotypic_target());
      gen_unit.compute_fitness(habitat.phenotypic_target());
    }
}



void Individual::inject_GU(Individual * donor) {
  // Add the GU at the end of the list
  _genetic_unit_list.emplace_back(this, donor->_genetic_unit_list.back());
}

void Individual::inject_2GUs(Individual * partner) {
  // We swap GUs from the end of the list.

  // TODO vld: As far as I understood the old code (47b27578), the
  // elements were not swapped but appended to the end of the other GU
  // lists. Error?  Original author (Dule, commit 47b27578), asked for
  // clarification by e-mail on 2015-02-23.

  const auto& gu_list_back_it = std::prev(_genetic_unit_list.end()); // initial last cell from _genetic_unit_list
  _genetic_unit_list.splice(_genetic_unit_list.end(), partner->_genetic_unit_list, std::prev(partner->_genetic_unit_list.end()));
  partner->_genetic_unit_list.splice(partner->_genetic_unit_list.end(), _genetic_unit_list, gu_list_back_it);
}

void Individual::compute_statistical_data() {
  if (metrics_ != nullptr) return; // Statistical data has already been computed,
                                 // nothing to do.

  metrics_ = new Metrics();

  if (not _phenotype_computed)
    compute_phenotype();

  for (const auto& gen_unit : _genetic_unit_list) {
    metrics_->Accumulate(gen_unit);
  }
}

void Individual::compute_non_coding() {
  if (nc_metrics_ != nullptr) return; // NC stats have already been computed,
                                    // nothing to do.
  nc_metrics_ = new NonCodingMetrics();

  for (auto& gen_unit: _genetic_unit_list) {
    gen_unit.compute_non_coding();
    nc_metrics_->Accumulate(gen_unit);
  }
}

void Individual::save(gzFile backup_file) const {
  // Write the name and "age" of the strain
  int8_t strain_string_len = strlen(_strain_name);
  gzwrite(backup_file, &strain_string_len, sizeof(strain_string_len));
  gzwrite(backup_file, _strain_name, strain_string_len+1);
  gzwrite(backup_file, &_age, sizeof(_age));

  // Write id and rank
  gzwrite(backup_file, &_id,   sizeof(_id));
  gzwrite(backup_file, &_rank, sizeof(_rank));

  // Write the position of the individual
  // gzwrite(backup_file, &x, sizeof(x));
  // gzwrite(backup_file, &y, sizeof(y));

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

int32_t Individual::get_nb_terminators() {
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
void Individual::compute_experimental_f_nu(int32_t nb_children, double* reproduction_statistics, double* offsprings_statistics, FILE* replication_file) {
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
  Individual * child = NULL;

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
double Individual::compute_theoritical_f_nu() {
  // We first have to collect information about genome structure.
  // Abbreviations are chosen according to Carole's formula.
  // Please notice that compared to the formula we have the beginning
  // and ends of neutral regions instead of 'functional regions'
  GeneticUnit& chromosome = _genetic_unit_list.front();
  int32_t L       = chromosome.get_dna()->length();
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
void Individual::remove_non_coding_bases() {
  for (auto& gen_unit: _genetic_unit_list)
    gen_unit.remove_non_coding_bases();

  // Delete the obsolete stats
  delete metrics_;
  metrics_ = NULL;
  delete nc_metrics_;
  nc_metrics_ = NULL;

#ifdef DEBUG
    compute_statistical_data();
    compute_non_coding();
    assert(get_nb_bases_in_0_coding_RNA() == 0);
  #endif
}

/// Double the bases that are not in coding RNA.
///
/// Double the bases that are not in coding RNA by addition of random
/// bases and test at each addition that fitness is not changed.
void Individual::double_non_coding_bases(void) {
  metrics_->total_genome_size_ = 0;
  int32_t initial_non_coding_base_nb = get_nb_bases_in_0_coding_RNA();

  for (auto& gen_unit: _genetic_unit_list)
    gen_unit.double_non_coding_bases();

  // Delete the obsolete stats
  delete metrics_;
  metrics_ = NULL;
  delete nc_metrics_;
  nc_metrics_ = NULL;

#ifdef DEBUG
    compute_statistical_data();
    compute_non_coding();
    assert(get_nb_bases_in_0_coding_RNA() == 2 * initial_non_coding_base_nb);
  #endif
}

// =================================================================
//                           Protected Methods
// =================================================================

// TODO vld: refactor make_protein_list and make_rna_list
void Individual::make_protein_list() {
  // Clear list
  _protein_list.clear();

  // Make a copy of each genetic unit's protein list
  for (auto& gen_unit: _genetic_unit_list)
  {
    // append all proteins from `gen_unit` to `protein_list_`
    for (auto& strand_id: {LEADING, LAGGING}) {
      auto& strand = gen_unit.get_protein_list(strand_id);
      for (auto& p: strand)
        _protein_list.push_back(&p);
    }
  }
}

void Individual::make_rna_list() {
  // Clear list
  _rna_list.clear();

  // Make a copy of each genetic unit's rna list
  for (const auto& gen_unit: _genetic_unit_list)
  {
    // Create proxies
    const auto& rna_list = gen_unit.get_rna_list();
    const auto& lead = rna_list[LEADING];
    const auto& lagg = rna_list[LAGGING];

    // append pointers to rna material to local rna_list_
    for (auto& strand: {LEADING, LAGGING})
      for (auto& rna: rna_list[strand])
        _rna_list.push_back(&rna);
  }
}

} // namespace aevol
