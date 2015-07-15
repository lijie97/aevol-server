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




// =================================================================
//                              Includes
// =================================================================
#include "GeneticUnit.h"

#include <cassert>
#include <list>
#include <algorithm>

#include "ExpManager.h"
#include "ExpSetup.h"
#include "Codon.h"
#include "Mutation.h"
#include "ae_enums.h"

#ifdef __REGUL
  #include "ae_individual_R.h"
#else
#include "Individual.h"
#endif

#include "Fuzzy.h"
#include "PointMutation.h"
#include "SmallDeletion.h"
#include "SmallInsertion.h"
#include "Duplication.h"
#include "Deletion.h"
#include "Translocation.h"
#include "Inversion.h"
#include "InsertionHT.h"
#include "ReplacementHT.h"

namespace aevol {

// =================================================================
//                       Miscellaneous Functions
// =================================================================
int compare_prot_pos(const void* pos, const void* prot) // This function has to be a plain int
// to comply with the definition of bsearch()
{
  if (((Protein *)prot)->get_shine_dal_pos() == *(int32_t*)pos) return 0;
  else return 1;
}

//##############################################################################
//                                                                             #
//                              Class genetic_unit                             #
//                                                                             #
//##############################################################################
// =====================================================================
//                          Accessors' definitions
// =====================================================================
ExpManager* GeneticUnit::get_exp_m() const
{
  return _exp_m;
}

Individual* GeneticUnit::get_indiv() const
{
  return _indiv;
}

Dna* GeneticUnit::get_dna() const
{
  assert(_dna->get_length() != 0);
  return _dna;
}

const Promoters2Strands& GeneticUnit::get_rna_list() const {
  return _rna_list;
}

std::list<Protein>& GeneticUnit::get_protein_list(Strand strand) {
  return _protein_list[strand];
}

void GeneticUnit::clear_protein_list(Strand strand) {
  _protein_list[strand].clear();
}

Fuzzy* GeneticUnit::get_activ_contribution() const
{
  return _activ_contribution;
}

Fuzzy* GeneticUnit::get_inhib_contribution() const
{
  return _inhib_contribution;
}

Fuzzy* GeneticUnit::get_phenotypic_contribution() const
{
  assert(_phenotypic_contribution != NULL);
  return _phenotypic_contribution;
}

/*!
  Returns the DNA sequence
*/
const char* GeneticUnit::get_sequence( void ) const
{
  return _dna->get_data();
}

/*!
  Returns the DNA sequence length
*/
int32_t GeneticUnit::get_seq_length( void ) const
{
  return _dna->get_length();
}

int32_t GeneticUnit::get_nb_coding_RNAs( void ) const
{
  return _nb_coding_RNAs;
}

int32_t GeneticUnit::get_nb_non_coding_RNAs( void ) const
{
  return _nb_non_coding_RNAs;
}

double GeneticUnit::get_overall_size_coding_RNAs( void ) const
{
  return _overall_size_coding_RNAs;
}

double GeneticUnit::get_av_size_coding_RNAs( void ) const
{
  if ( _nb_coding_RNAs != 0 )
  {
    return _overall_size_coding_RNAs / _nb_coding_RNAs;
  }
  else return 0.0;
}

double GeneticUnit::get_overall_size_non_coding_RNAs( void ) const
{
  return _overall_size_non_coding_RNAs;
}

double GeneticUnit::get_av_size_non_coding_RNAs( void ) const
{
  if ( _nb_non_coding_RNAs != 0 )
  {
    return _overall_size_non_coding_RNAs / _nb_non_coding_RNAs;
  }
  else return 0.0;
}

int32_t GeneticUnit::get_nb_genes_activ( void ) const
{
  return _nb_genes_activ;
}

int32_t GeneticUnit::get_nb_genes_inhib( void ) const
{
  return _nb_genes_inhib;
}

int32_t GeneticUnit::get_nb_functional_genes( void ) const
{
  return _nb_fun_genes;
}

int32_t GeneticUnit::get_nb_non_functional_genes( void ) const
{
  return _nb_non_fun_genes;
}

double GeneticUnit::get_overall_size_functional_genes( void ) const
{
  return _overall_size_fun_genes;
}

double GeneticUnit::get_av_size_functional_genes( void ) const
{
  if ( _nb_fun_genes != 0 )
  {
    return _overall_size_fun_genes / _nb_fun_genes;
  }
  else return 0.0;
}

double GeneticUnit::get_overall_size_non_functional_genes( void ) const
{
  return _overall_size_non_fun_genes;
}

double GeneticUnit::get_av_size_non_functional_genes( void ) const
{
  if ( _nb_non_fun_genes != 0 )
  {
    return _overall_size_non_fun_genes / _nb_non_fun_genes;
  }
  else return 0.0;
}

int32_t GeneticUnit::get_nb_bases_in_0_CDS() const
{
  assert (_non_coding_computed );
  return _nb_bases_in_0_CDS;
}

int32_t GeneticUnit::get_nb_bases_in_0_functional_CDS() const
{
  assert (_non_coding_computed );
  return _nb_bases_in_0_functional_CDS;
}

int32_t GeneticUnit::get_nb_bases_in_0_non_functional_CDS() const
{
  assert (_non_coding_computed );
  return _nb_bases_in_0_non_functional_CDS;
}

int32_t GeneticUnit::get_nb_bases_in_0_RNA() const
{
  assert (_non_coding_computed );
  return _nb_bases_in_0_RNA;
}

int32_t GeneticUnit::get_nb_bases_in_0_coding_RNA() const
{
  assert (_non_coding_computed );
  return _nb_bases_in_0_coding_RNA;
}

int32_t GeneticUnit::get_nb_bases_in_0_non_coding_RNA() const
{
  assert (_non_coding_computed );
  return _nb_bases_in_0_non_coding_RNA;
}

int32_t GeneticUnit::get_nb_bases_non_essential() const
{
  assert (_non_coding_computed );
  return _nb_bases_non_essential;
}

int32_t GeneticUnit::get_nb_bases_non_essential_including_nf_genes() const
{
  assert (_non_coding_computed );
  return _nb_bases_non_essential_including_nf_genes;
}

int32_t GeneticUnit::get_nb_bases_in_neutral_regions() const
{
  assert (_non_coding_computed );
  return _nb_bases_in_neutral_regions;
}

int32_t GeneticUnit::get_nb_neutral_regions() const
{
  assert (_non_coding_computed );
  return _nb_neutral_regions;
}

int32_t* GeneticUnit::get_beginning_neutral_regions() const
{
  assert (_non_coding_computed );
  return _beginning_neutral_regions;
}

int32_t* GeneticUnit::get_end_neutral_regions() const
{
  assert (_non_coding_computed );
  return _end_neutral_regions;
}

double GeneticUnit::get_modularity() const
{
  return _modularity;
}

double GeneticUnit::get_dist_to_target_by_feature( PhenotypicFeature feature ) const
{
  assert(_distance_to_target_computed);

  return _dist_to_target_by_feature[feature];
}

double GeneticUnit::get_fitness() const
{
  assert(_fitness_computed);

  return _fitness;
}

double GeneticUnit::get_fitness_by_feature( PhenotypicFeature feature ) const
{
  assert( _fitness_computed );

  return _fitness_by_feature[feature];
}

int32_t GeneticUnit::get_min_gu_length( void ) const
{
  return _min_gu_length;
}

int32_t GeneticUnit::get_max_gu_length( void ) const
{
  return _max_gu_length;
}

void GeneticUnit::set_min_gu_length( int32_t min_gu_length )
{
  _min_gu_length = min_gu_length;
}

void GeneticUnit::set_max_gu_length( int32_t max_gu_length )
{
  _max_gu_length = max_gu_length;
}

void GeneticUnit::set_exp_m( ExpManager * exp_m )
{
  _exp_m = exp_m;
}

// =====================================================================
//                       functions' definition
// =====================================================================
void GeneticUnit::print_rnas() const {
  print_rnas(_rna_list);
}

/* static */ void GeneticUnit::print_rnas(const Promoters2Strands& rnas) {
  print_rnas(rnas[LEADING], LEADING);
  print_rnas(rnas[LAGGING], LAGGING);
}

/* static */ void GeneticUnit::print_rnas(const Promoters1Strand& rnas, Strand strand) {
  printf( "  %s ( %" PRId32 " )\n", strand == LEADING ? " LEADING " : "LAGGING", static_cast<int32_t>(rnas.size()));
  for (auto& rna: rnas) {
    assert( rna.get_strand() == strand );
    printf( "    Promoter on %s at %" PRId32 "\n", strand == LEADING ? " LEADING " : "LAGGING", rna.get_promoter_pos());
  }
}

bool GeneticUnit::is_start( Strand strand, int32_t index ) const
{
  return ( get_codon( strand, index ) == CODON_START );
}

bool GeneticUnit::is_stop( Strand strand, int32_t index ) const
{
  return ( get_codon( strand, index ) == CODON_STOP );
}

void GeneticUnit::remove_all_promoters( void )
{
  _rna_list[LEADING].clear();
  _rna_list[LAGGING].clear();
}

void GeneticUnit::move_all_promoters_after( int32_t pos, int32_t delta_pos )
{
  move_all_leading_promoters_after( pos, delta_pos );
  move_all_lagging_promoters_after( pos, delta_pos );
}

void GeneticUnit::extract_promoters_included_in(int32_t pos_1,
                                                int32_t pos_2,
                                                Promoters2Strands& extracted_promoters)
{
  assert(pos_1 >= 0);
  assert(pos_1 < pos_2);
  assert(pos_2 <= _dna->get_length());

  if ( pos_2 - pos_1 < PROM_SIZE )
    return;
  extract_leading_promoters_starting_between(pos_1, pos_2 - PROM_SIZE + 1, extracted_promoters[LEADING]);
  extract_lagging_promoters_starting_between(pos_1 + PROM_SIZE - 1, pos_2, extracted_promoters[LAGGING]);
}

void GeneticUnit::extract_promoters_starting_between(int32_t pos_1,
                                                     int32_t pos_2,
                                                     Promoters2Strands& extracted_promoters)
{
  extract_leading_promoters_starting_between(pos_1, pos_2, extracted_promoters[LEADING]);
  extract_lagging_promoters_starting_between(pos_1, pos_2, extracted_promoters[LAGGING]);
}


/*!
  \brief  Remove those promoters that would be broken if the chromosome was cut at pos.

  Remove promoters that include BOTH the base before AND after pos (marked X in the cartoon below).
  If the genome is smaller than the size of a promoter, all the promoters will be removed.

  \verbatim
     -------------------------------------------------------
    |   |   |   |   | X | X |   |   |   |   |   |   |   |   |
     -------------------------------------------------------
    ^                   ^
    0                  pos
  \endverbatim
*/
void GeneticUnit::remove_promoters_around( int32_t pos )
{
  if ( _dna->get_length() >= PROM_SIZE )
  {
    remove_leading_promoters_starting_between( Utils::mod(pos - PROM_SIZE + 1, _dna->get_length()), pos );
    remove_lagging_promoters_starting_between( pos, Utils::mod(pos + PROM_SIZE - 1, _dna->get_length()) );
  }
  else
  {
    remove_all_promoters();
  }
}


/*!
  \brief  Remove those promoters that would be broken if the sequence [pos_1 ; pos_2[ was deleted.

  Remove promoters that     * include BOTH the base before AND after pos_1 (marked X in the cartoon below).
                            * include BOTH the base before AND after pos_2 (marked Y in the cartoon below).
                            * are completely contained between pos_1 and pos_2.
  If the remaining sequence, i.e. [pos_2 ; pos_1[ is smaller than the size of a promoter, all the promoters will be removed.

  \verbatim
     -------------------------------------------------------
    |   |   |   |   | X | X |   |   |   | Y | Y |   |   |   |
     -------------------------------------------------------
    ^                   ^                   ^
    0                 pos_1               pos_2
  \endverbatim
*/
void GeneticUnit::remove_promoters_around( int32_t pos_1, int32_t pos_2 )
{
  if ( Utils::mod(pos_1 - pos_2, _dna->get_length()) >= PROM_SIZE )
  {
    remove_leading_promoters_starting_between( Utils::mod(pos_1 - PROM_SIZE + 1, _dna->get_length()), pos_2 );
    remove_lagging_promoters_starting_between( pos_1, Utils::mod(pos_2 + PROM_SIZE - 1, _dna->get_length()) );
  }
  else
  {
    remove_all_promoters();
  }
}


/// Look for promoters that are astride pos and add them to the list
/// of promoters (_rna_list).
///
/// Look for promoters that include BOTH the base before AND after pos (marked X in the cartoon below).
/// If the genome is smaller than the size of a promoter, no search is performed.
///
/// \verbatim
///    -------------------------------------------------------
///   |   |   |   |   | X | X |   |   |   |   |   |   |   |   |
///    -------------------------------------------------------
///   ^                   ^
///   0                  pos
/// \endverbatim
void GeneticUnit::look_for_new_promoters_around(int32_t pos)
{
  assert(pos >= 0 && pos <= _dna->get_length());

  if (_dna->get_length() >= PROM_SIZE)
  {
    look_for_new_leading_promoters_starting_between(
        Utils::mod(pos - PROM_SIZE + 1, _dna->get_length()),
        pos);
    look_for_new_lagging_promoters_starting_between(
        pos,
        Utils::mod(pos + PROM_SIZE - 1, _dna->get_length()));
  }
}


/// Look for promoters that contain at least 1 base lying in [pos_1 ;
/// pos_2[ and add them to the list of promoters (_rna_list).
///
/// Look for promoters that   * include BOTH the base before AND after pos_1 (marked X in the cartoon below).
///                           * include BOTH the base before AND after pos_2 (marked Y in the cartoon below).
///                           * are completely contained between pos_1 and pos_2.
/// If the genome is smaller than the size of a promoter, no search is performed.
///
/// \verbatim
///    -------------------------------------------------------
///   |   |   |   |   | X | X |   |   |   | Y | Y |   |   |   |
///    -------------------------------------------------------
///   ^                   ^                   ^
///   0                 pos_1               pos_2
/// \endverbatim
void GeneticUnit::look_for_new_promoters_around( int32_t pos_1, int32_t pos_2 )
{
  //~ if ( Utils::mod( pos_1 - pos_2, _dna->get_length()) == PROM_SIZE - 1 )
  //~ {
  //~ // We have to look at every possible position on the genome.
  //~ locate_promoters();
  //~ }
  /*else*/ if ( _dna->get_length() >= PROM_SIZE )
  {
    look_for_new_leading_promoters_starting_between( Utils::mod(pos_1 - PROM_SIZE + 1, _dna->get_length()), pos_2 );
    look_for_new_lagging_promoters_starting_between( pos_1, Utils::mod(pos_2 + PROM_SIZE - 1, _dna->get_length()) );
  }
}

void GeneticUnit::copy_promoters_starting_between(int32_t pos_1,
                                                  int32_t pos_2,
                                                  Promoters2Strands& new_promoter_lists)
{
  copy_leading_promoters_starting_between(pos_1, pos_2, new_promoter_lists[LEADING]);
  copy_lagging_promoters_starting_between(pos_1, pos_2, new_promoter_lists[LAGGING]);
}

void GeneticUnit::copy_promoters_included_in(int32_t pos_1,
                                             int32_t pos_2,
                                             Promoters2Strands& new_promoter_lists )
{
  if ( Utils::mod( pos_2 - pos_1 - 1, _dna->get_length() ) + 1 >= PROM_SIZE )
  {
    copy_leading_promoters_starting_between( pos_1, Utils::mod( pos_2 - PROM_SIZE + 1, _dna->get_length() ), new_promoter_lists[LEADING] );
    copy_lagging_promoters_starting_between( Utils::mod( pos_1 + PROM_SIZE - 1, _dna->get_length() ), pos_2, new_promoter_lists[LAGGING] );
  }
}

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
/*!
  \brief Create a new genetic unit for indiv with a random DNA sequence of length length

  Promoters will be looked for on the whole sequence but no further process
  will be performed.
*/
GeneticUnit::GeneticUnit(Individual * indiv,
                         int32_t length,
                         std::shared_ptr<JumpingMT> prng)
{
  _indiv = indiv;
  _exp_m = indiv->get_exp_m();

  _transcribed                        = false;
  _translated                         = false;
  _phenotypic_contributions_computed  = false;
  _non_coding_computed                = false;
  _distance_to_target_computed        = false;
  _fitness_computed                   = false;

  _min_gu_length = -1;
  _max_gu_length = -1;

  _dna = new Dna( this, length, prng );

  // Create empty fuzzy sets for the phenotypic contributions
  _activ_contribution = new Fuzzy();
  _inhib_contribution = new Fuzzy();
  _phenotypic_contribution = NULL;
  // NB : _phenotypic_contribution is only an indicative value,
  //      it is not used for the whole phenotype computation

  // _dist_to_target_per_segment depends on the segmentation of the environment
  // and will hence be newed at evaluation time
  _dist_to_target_per_segment = NULL;

  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }

  // Look for promoters
  locate_promoters();

  init_statistical_data();
}

/// Create a new genetic unit for `indiv` with sequence `seq` of size
/// `length` [and containing promoters `prom_list`]
///
///  Promoters will be looked for if prom_list is not provided (this
///  may take some time).
///
/// WARNING:
///   seq will be used directly which means the caller must not delete it
///   The same goes for prom_list if it is provided.
GeneticUnit::GeneticUnit(Individual * indiv,
                         char* seq,
                         int32_t length,
                         const Promoters2Strands& prom_list /* = {{},{}} */)
{
  _exp_m = indiv->get_exp_m();
  _indiv = indiv;

  _transcribed                        = false;
  _translated                         = false;
  _phenotypic_contributions_computed  = false;
  _non_coding_computed                = false;
  _distance_to_target_computed        = false;
  _fitness_computed                   = false;

  _min_gu_length = -1;
  _max_gu_length = -1;

  _dna = new Dna( this, seq, length );

  if (not prom_list[LEADING].empty() and not prom_list[LAGGING].empty()) { // if not default `prom_list`
    // Copy rna lists
    _rna_list = prom_list;
    Dna::set_GU(_rna_list, this);
  } else {
    for (auto& strand: {LEADING, LAGGING})
      assert(_rna_list[strand].empty());

    // Look for promoters
    locate_promoters();
  }

  // Create empty fuzzy sets for the phenotypic contributions
  _activ_contribution = new Fuzzy();
  _inhib_contribution = new Fuzzy();
  _phenotypic_contribution = NULL;
  // NB : _phenotypic_contribution is only an indicative value,
  //      it is not used for the whole phenotype computation

  // Initialize all the fitness-related stuff
  _dist_to_target_per_segment = NULL;

  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }

  init_statistical_data();
}





/*!
  \brief Copy constructor.

  Copies the DNA and recomputes all the rest.
  It is slower than copying as much as possible and regenerate only what is necessary but it works whatever the state of the model GU.
*/
GeneticUnit::GeneticUnit(Individual* indiv, const GeneticUnit& model)
{
  _exp_m = indiv->get_exp_m();
  _indiv = indiv;

  _transcribed                        = false;
  _translated                         = false;
  _phenotypic_contributions_computed  = false;
  _non_coding_computed                = false;
  _distance_to_target_computed        = model._distance_to_target_computed;
  _fitness_computed                   = model._fitness_computed;

  _min_gu_length = model._min_gu_length;
  _max_gu_length = model._max_gu_length;

  // Copy DNA
  _dna = new Dna( this, *(model._dna) );

  // Create empty fuzzy sets for the phenotypic contributions
  _activ_contribution = new Fuzzy();
  _inhib_contribution = new Fuzzy();
  _phenotypic_contribution = NULL;
  // NB : _phenotypic_contribution is only an indicative value, not used for the whole phenotype computation
  _dist_to_target_per_segment = NULL;
  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }

  // Compute everything
  init_statistical_data();
  locate_promoters();
  do_transcription();
  do_translation();
  compute_phenotypic_contribution();
}

/**
 * Reproduction constructor
 *
 * Create a new genetic unit copying the DNA sequence and the promoter list
 * from the provided `parent`
 */
GeneticUnit::GeneticUnit(Individual* indiv, const GeneticUnit* parent)
{
  _exp_m = indiv->get_exp_m();
  _indiv = indiv;

  _transcribed                        = false;
  _translated                         = false;
  _phenotypic_contributions_computed  = false;
  _non_coding_computed                = false;
  _distance_to_target_computed        = false;
  _fitness_computed                   = false;

  _min_gu_length = parent->_min_gu_length;
  _max_gu_length = parent->_max_gu_length;

  // Copy DNA
  _dna = new Dna( this, parent->_dna );

  // Copy promoter list (_rna_list)
  // Note that the length of the RNA will have to be recomputed (do_transcription)
  for (auto& strand: {LEADING, LAGGING}) {
    for (auto& rna: parent->_rna_list[strand]) {
#ifndef __REGUL
      _rna_list[strand].emplace_back(this, rna);
#else
      _rna_list[strand].emplace_back(this, (dynamic_cast<ae_rna_R>(rna)));
#endif
    }
  }

  // Create empty fuzzy sets for the phenotypic contributions
  _activ_contribution = new Fuzzy();
  _inhib_contribution = new Fuzzy();
  _phenotypic_contribution = NULL;
  // NB : _phenotypic_contribution is only an indicative value, not used for the whole phenotype computation

  // Initialize all the fitness-related stuff
  _dist_to_target_per_segment = NULL;

  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }

  init_statistical_data();
}

GeneticUnit::GeneticUnit( Individual * indiv, gzFile backup_file )
{
  _exp_m = indiv->get_exp_m();
  _indiv = indiv;

  _transcribed                        = false;
  _translated                         = false;
  _phenotypic_contributions_computed  = false;
  _non_coding_computed                = false;
  _distance_to_target_computed        = false;
  _fitness_computed                   = false;

  _dna = new Dna( this, backup_file );

  gzread( backup_file, &_min_gu_length, sizeof(_min_gu_length) );
  gzread( backup_file, &_max_gu_length, sizeof(_max_gu_length) );

  // Create empty fuzzy sets for the phenotypic contributions
  _activ_contribution = new Fuzzy();
  _inhib_contribution = new Fuzzy();
  _phenotypic_contribution = NULL;
  // NB : _phenotypic_contribution is only an indicative value, not used for the whole phenotype computation

  // Initialize all the fitness-related stuff
  _dist_to_target_per_segment = NULL;

  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }


  // Look for promoters
  locate_promoters();

  init_statistical_data();
}

/*!
  \brief Create a new genetic unit for indiv with a sequence saved in a text file

  Promoters will be looked for on the whole sequence but no further process
  will be performed.
*/
GeneticUnit::GeneticUnit( Individual * indiv, char* organism_file_name )
{
  _exp_m = indiv->get_exp_m();
  _indiv = indiv;

  _transcribed                        = false;
  _translated                         = false;
  _phenotypic_contributions_computed  = false;
  _non_coding_computed                = false;
  _distance_to_target_computed        = false;
  _fitness_computed                   = false;

  _dna = new Dna( this, organism_file_name );

  // Create empty fuzzy sets for the phenotypic contributions
  _activ_contribution = new Fuzzy();
  _inhib_contribution = new Fuzzy();
  _phenotypic_contribution = NULL;
  // NB : _phenotypic_contribution is only an indicative value,
  //      it is not used for the whole phenotype computation

  // Initialize all the fitness-related stuff
  _dist_to_target_per_segment = NULL;

  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }

  // Look for promoters
  locate_promoters();

  init_statistical_data();
}

// =================================================================
//                             Destructors
// =================================================================
GeneticUnit::~GeneticUnit( void )
{
  delete _dna;
  delete _activ_contribution;
  delete _inhib_contribution;
  if ( _phenotypic_contribution != NULL ) delete _phenotypic_contribution;

  if ( _dist_to_target_per_segment != NULL ) delete [] _dist_to_target_per_segment;

  assert( _dist_to_target_by_feature != NULL );
  delete [] _dist_to_target_by_feature;
  assert( _fitness_by_feature != NULL );
  delete [] _fitness_by_feature;

  delete [] _beginning_neutral_regions;
  delete [] _end_neutral_regions;
}

// =================================================================
//                            Public Methods
// =================================================================
/// Look for promoters in the genome and create a new Rna in the
/// corresponding strand's RNA list
void GeneticUnit::locate_promoters( void )
{

  // TODO vld 2015-04-14: make it return the generated rna-list rather
  // than alter the current one?

  int8_t dist; // Hamming distance of the sequence from the promoter consensus

  // Empty RNA list
  for (auto& strand: _rna_list) {
    strand.clear();
  }

  if ( _dna->get_length() < PROM_SIZE )
    return;
  for ( int32_t i = 0 ; i < _dna->get_length() ; i++ ) {
#ifndef __REGUL
    if (is_promoter(LEADING, i, dist)) // dist takes the hamming distance of the sequence from the consensus
      _rna_list[LEADING].emplace_back(this, LEADING, i, dist);
    if (is_promoter(LAGGING, _dna->get_length() - i - 1, dist))
      _rna_list[LAGGING].emplace_back(this, LAGGING, _dna->get_length() - i - 1, dist);
#else
    if (is_promoter(LEADING, i, dist))
      _rna_list[LEADING].emplace_back(this, LEADING, i, dist);
    if (is_promoter(LAGGING, _dna->get_length() - i - 1, dist))
      _rna_list[LAGGING].emplace_back(this, LAGGING, _dna->get_length() - i - 1, dist);
#endif
  }
}

void GeneticUnit::do_transcription( void )
{
  if ( _transcribed ) return;
  _transcribed = true;

  int32_t transcript_start  = -1;
  int32_t genome_length     = _dna->get_length();

  // If the genome is not long enough to bear a promoter and a terminator,
  // we set all its RNAs to a length of -1
  // (NB but a terminator can share code with the promoter, making it
  // possible for the genome to be no longer than the promoter)
  if (genome_length < PROM_SIZE) {
    for (auto& strand: _rna_list)
      for (auto& rna: strand)
        rna.set_transcript_length(-1);
    return;
  }

  for (auto& strand_id: {LEADING, LAGGING}) {
    auto& strand = _rna_list[strand_id];
    for (auto rna = strand.begin(); rna != strand.end(); ++rna) {
      transcript_start = rna->get_first_transcribed_pos();
      rna->set_transcript_length(-1);

      int32_t i;
      for (i = 0 ; i < genome_length ; ++i) {
        if (is_terminator(strand_id, transcript_start + (strand_id == LEADING ? i : -i))) {
          // Found terminator => set transcript's length
          rna->set_transcript_length(i + TERM_SIZE);

          // Deduce the length of all the RNAs that share the same terminator
          // These are the RNAs whose promoter is entirely (and strictly) included
          // between the promoter and the terminator of the RNA we have just treated.
          // They are hence the RNAs whose promoter starts at most i bases after the
          // current rna's promoter
          for (auto rna2 = std::next(rna); rna2 != strand.end(); ++rna2) {
            // We know that if rna2 is after rna then:
            // - rna_2.pos > rna.pos for LEADING strand
            // - rna_2.pos < rna.pos for LAGGING strand
            // because the list is sorted.

            auto delta_pos = abs(rna2->get_promoter_pos() - rna->get_promoter_pos());
            if (delta_pos <= i) {
              rna2->set_transcript_length(i - delta_pos + TERM_SIZE);
              // Step forward in RNA list
              rna = rna2;
            } else {
              // The promoter of rna_2 is after (or contains a part of) the terminator of rna,
              // we will need to search its own terminator
              break;
            }
          }
          // Terminator found for this RNA, nothing else to do (for this RNA)
          break;
        }
      }

      if (i == genome_length)
        // We have searched the whole genome and found no terminator for this promoter.
        // We consider that no RNA can actually be produced, hence we set the transcript
        // length to -1. This will prevent the search for coding sequences downstream of this promoter.
        // However, we do not destroy the Rna object, it must still be kept in memory and
        // transmitted to the offspring in case a mutation recreates a terminator.
        rna->set_transcript_length( -1 );
    }
  }
}

void GeneticUnit::do_translation()
{
  if (_translated)
    return;
  _translated = true;
  if (not _transcribed)
    do_transcription();

  int32_t transcript_start  = -1;
  int32_t transcript_length = -1;
  int32_t genome_length     = _dna->get_length();

  for (auto strand: {LEADING, LAGGING}) {
    for (auto& rna: _rna_list[strand]) {
      transcript_start  = rna.get_first_transcribed_pos();
      transcript_length = rna.get_transcript_length();

      // Try every position where a translation process could occur
      // Minimum number of bases needed is SHINE_DAL_SIZE + SHINE_START_SPACER + 3 * CODON_SIZE
      // (3 codons for START + STOP + at least one amino-acid)
      for (int32_t i = 0;
           transcript_length - i >= SHINE_DAL_SIZE + SHINE_START_SPACER + 3 * CODON_SIZE;
           ++i) {
        if (is_shine_dalgarno(strand, Utils::mod(transcript_start
                                                    + (strand == LEADING ? i : -i), genome_length))
            and is_start(strand,
                         Utils::mod(transcript_start
                                       + (strand == LEADING ? 1 : -1)
                                       * (i + SHINE_DAL_SIZE + SHINE_START_SPACER), genome_length))) {
          // We found a translation initiation, we can now build the
          // protein until we find a STOP codon or until we reach the
          // end of the transcript (in which case the protein is not
          // valid)

          // First of all, we will check whether this CDS has already
          // been translated (because it is present on another RNA).
          // In that case, we don't need to tranlate it again, we only
          // need to increase the protein's concentration according to
          // the promoter transcription level
          int32_t shine_dal_pos = Utils::mod(transcript_start +
                                                 (strand == LEADING ? i : -i),
                                             genome_length);
          auto& protein_strand = _protein_list[strand];
          auto protein = find_if(protein_strand.begin(),
                                 protein_strand.end(),
                                 [shine_dal_pos](Protein& p)
                                 { return p.get_shine_dal_pos() == shine_dal_pos; });

          if (protein != protein_strand.end()) {
            protein->add_RNA(&rna);
            rna.add_transcribed_protein(&*protein);
          }
          else {
            // Build codon list and make new protein when stop found
            int32_t j = i + SHINE_DAL_SIZE + SHINE_START_SPACER + CODON_SIZE; // next codon to examine

            std::list<Codon*> codon_list;

            while (transcript_length - j >= CODON_SIZE) {
              auto codon = new Codon(_dna,
                                        strand,
                                        Utils::mod(transcript_start + (strand == LEADING ? j : -j),
                                                      genome_length));

              if (codon->is_stop()) {
                if (not codon_list.empty()) { // at least one amino-acid
                  // The protein is valid, create the corresponding object
                  protein_strand.emplace_back(this, codon_list, strand,
                                              shine_dal_pos, &rna,
                                              get_indiv()->get_w_max());
                  auto& protein = protein_strand.back();
                  codon_list.clear(); // has been copied into `protein`
                  rna.add_transcribed_protein(&protein);

                  if (protein.get_is_functional()) {
                    _nb_fun_genes++;
                    _overall_size_fun_genes += protein.get_length() * CODON_SIZE;

                    if (protein.get_height() > 0) _nb_genes_activ++;
                    else                           _nb_genes_inhib++;
                  }
                  else {
                    _nb_non_fun_genes++;
                    _overall_size_non_fun_genes += (strand == LEADING) ?
                                                   protein.get_length() * CODON_SIZE :
                                                   ( protein.get_length() + 2 ) * CODON_SIZE;
                  }
                }
                delete codon;
                codon = nullptr;
                break;
              }
              else {
                codon_list.push_back(codon);
                codon = nullptr; // don't delete codon: recycled into the list
              }
              j += CODON_SIZE;
            }
            if (not codon_list.empty())
              for (auto& c: codon_list)
                delete c;
          }
        }
      }

      // Statistics
      if (not rna.get_transcribed_proteins().empty()) { // coding RNA
        _nb_coding_RNAs++;
        _overall_size_coding_RNAs += rna.get_transcript_length();
      }
      else { // non-coding RNA
        _nb_non_coding_RNAs++;
        _overall_size_non_coding_RNAs += rna.get_transcript_length();
      }
    }
  }
}

void GeneticUnit::compute_phenotypic_contribution( void )
{
  if ( _phenotypic_contributions_computed ) return;
  _phenotypic_contributions_computed = true;
  if ( ! _translated ) do_translation();

  for (const auto& strand: _protein_list) // two strands: LEADING & LAGGING
    for (const auto& prot: strand)
      if (prot.get_is_functional())
        ((prot.get_height() > 0) ? _activ_contribution: _inhib_contribution)
            ->add_triangle(prot.get_mean(),
                           prot.get_width(),
                           prot.get_height() * prot.get_concentration());
  // if ( prot->get_height() > 0 )
  //   _activ_contribution->add_triangle(prot->get_mean(),
  //                                     prot->get_width(),
  //                                     prot->get_height() * prot->get_concentration() );
  // else
  //   _inhib_contribution->add_triangle(prot->get_mean(),
  //                                     prot->get_width(),
  //                                     prot->get_height() * prot->get_concentration() );

  // It is not necessary to add a lower bound to _activ_contribution as there can be no negative y
  // The same goes for the upper bound for _inhib_contribution
  _activ_contribution->clip(Fuzzy::max,   Y_MAX );
  _inhib_contribution->clip(Fuzzy::min, - Y_MAX );
  _activ_contribution->simplify();
  _inhib_contribution->simplify();

  if ( _exp_m->get_output_m()->get_compute_phen_contrib_by_GU() )
  {
    _phenotypic_contribution = new Phenotype();
    _phenotypic_contribution->add( *_activ_contribution );
    _phenotypic_contribution->add( *_inhib_contribution );
    _phenotypic_contribution->simplify();
  }
}

/*!
  \brief Compute the areas between the phenotype and the environment for each environmental segment.

  If the environment is not segmented, the total area is computed
*/
void GeneticUnit::compute_distance_to_target(const PhenotypicTarget& target) {
  if ( _distance_to_target_computed ) return; // _distance_to_target has already been computed, nothing to do.
  _distance_to_target_computed = true;

  compute_phenotypic_contribution();

  // Compute the difference between the (whole) phenotype and the environment
  Fuzzy* delta = new Fuzzy( *_phenotypic_contribution );
  delta->sub( target );

  PhenotypicSegment ** segments = target.segments();

  // TODO <david.parsons@inria.fr> We should take into account that we compute the areas in order (from the leftmost segment, rightwards)
  //   => We shouldn't parse the whole list of points on the left of the segment we are considering (we have
  //      already been through them!)

  if ( _dist_to_target_per_segment == NULL )
  {
    _dist_to_target_per_segment = new double [target.nb_segments()]; // Can not be allocated in constructor because number of segments is then unknow
  }
  for (int8_t i = 0 ; i < target.nb_segments() ; i++)
  {
    _dist_to_target_per_segment[i] = delta->get_geometric_area( segments[i]->start, segments[i]->stop );
    _dist_to_target_by_feature[segments[i]->feature] += _dist_to_target_per_segment[i];
  }

  delete delta;
}

/*!
  \brief Compute a "proper" fitness value (one that increases when the individual is fitter).

  The behaviour of this function depends on many parameters and most notably on whether it is
  a "composite" fitness or not, and on the selection scheme.
*/
void GeneticUnit::compute_fitness(const PhenotypicTarget& target)
{
  if ( _fitness_computed ) return; // Fitness has already been computed, nothing to do.
  _fitness_computed = true;

#ifdef NORMALIZED_FITNESS

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    if (target.get_area_by_feature(i)==0.)
    {
      _fitness_by_feature[i] = 0.;
    }
    else
    {
      _fitness_by_feature[i] =  ( target.get_area_by_feature(i) - _dist_to_target_by_feature[i] ) / target.get_area_by_feature(i);
      if ( (_fitness_by_feature[i] < 0.) && (i != METABOLISM) ) // non-metabolic fitness can NOT be lower than zero (we do not want individual to secrete a negative quantity of public good)
      {
        _fitness_by_feature[i] = 0.;
      }
    }
  }

  if ((! _indiv->get_placed_in_population()) || (! _exp_m->get_with_secretion() ))
  {
    _fitness = _fitness_by_feature[METABOLISM];
  }
  else
  {
    _fitness =  _fitness_by_feature[METABOLISM] * ( 1 + _exp_m->get_secretion_contrib_to_fitness() * ( _indiv->get_grid_cell()->get_compound_amount() - _exp_m->get_secretion_cost() * _fitness_by_feature[SECRETION] ) );
  }

  if ( _exp_m->get_selection_scheme() == FITNESS_PROPORTIONATE ) // Then the exponential selection is integrated inside the fitness value
  {
    _fitness = exp( -_exp_m->get_selection_pressure() * (1 - _fitness) );
  }

#else

    for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
    {
      if ( i == SECRETION )
      {
        _fitness_by_feature[SECRETION] =  exp( -_exp_m->get_selection_pressure() * _dist_to_target_by_feature[SECRETION] )
                                          - exp( -_exp_m->get_selection_pressure() * target.area_by_feature(SECRETION) );

        if ( _fitness_by_feature[i] < 0 )
        {
          _fitness_by_feature[i] = 0;
        }
      }
      else
      {
        _fitness_by_feature[i] = exp( -_exp_m->get_selection_pressure() * _dist_to_target_by_feature[i] );
      }
    }

    // Calculate combined, total fitness here!
    // Multiply the contribution of metabolism and the amount of compound in the environment
    if ((! _indiv->get_placed_in_population()) || (! _exp_m->get_with_secretion() ))
    {
      _fitness = _fitness_by_feature[METABOLISM] ;
    }
    else
    {
      _fitness =  _fitness_by_feature[METABOLISM] *
                  ( 1 + _exp_m->get_secretion_contrib_to_fitness() * _indiv->get_grid_cell()->compound_amount()
                    - _exp_m->get_secretion_cost() * _fitness_by_feature[SECRETION] );
    }

#endif

  }


  void GeneticUnit::reset_expression( void )
  {
    // useful if the DNA sequence has changed (cf post-treatment programs
    // which replay mutations)

    _transcribed                        = false;
    _translated                         = false;
    _phenotypic_contributions_computed  = false;
    _non_coding_computed                = false;
    _distance_to_target_computed        = false;
    _fitness_computed                   = false;


    // I do not erase the RNA lists, because they were updated
    // during the mutations (cf Dna::undergo_this_mutation)
    // TODO : Reinitialize _transcribed proteins ?

    if ( _activ_contribution != NULL )
    {
      delete _activ_contribution;
      _activ_contribution = new Fuzzy();
    }

    if ( _inhib_contribution != NULL )
    {
      delete _inhib_contribution;
      _inhib_contribution = new Fuzzy();
    }

    if ( _phenotypic_contribution != NULL )
    {
      delete _phenotypic_contribution; // Not re-created now, will be conditionally allocated in compute_phenotypic_contribution
      _phenotypic_contribution = NULL;
    }

    init_statistical_data();
  }


  void GeneticUnit::print_coding_rnas() {
    for (int strand_id = LEADING; strand_id <= LAGGING; ++strand_id) {
      auto& strand = _rna_list[strand_id];
      printf("  %s \n", StrandName[strand_id]);
      for (const auto& rna: strand)
        if (rna.is_coding())
          printf("Promoter at %" PRId32 ", last transcribed position at %" PRId32 "\n", rna.get_promoter_pos(), rna.get_last_transcribed_pos());
    }
  }

  void GeneticUnit::print_proteins( void ) const
  {
    printf( "  LEADING ( %" PRId32 " )\n", static_cast<int32_t>(_protein_list[LEADING].size()));
    for (const auto& prot: _protein_list[LEADING])
      printf( "    Gene on LEADING at %" PRId32 " (%" PRId32 ") (%f %f %f) (%f) %s\n",
              prot.get_shine_dal_pos(), prot.get_length(),
              prot.get_mean(), prot.get_width(), prot.get_height(), prot.get_concentration(),
              prot.get_is_functional() ? "functional" : "non functional" );


    printf( "  LAGGING ( %" PRId32 " )\n", static_cast<int32_t>(_protein_list[LAGGING].size()));
    for (const auto& prot: _protein_list[LAGGING])
      printf( "    Gene on LAGGING at %" PRId32 " (%" PRId32 ") (%f %f %f) (%f) %s\n",
              prot.get_shine_dal_pos(), prot.get_length(),
              prot.get_mean(), prot.get_width(), prot.get_height(), prot.get_concentration(),
              prot.get_is_functional() ? "functional" : "non functional" );
  }

  bool GeneticUnit::is_promoter( Strand strand, int32_t pos, int8_t& dist ) const
  {
    //~ printf( "=============================================== is_promoter\n" );
    //~ printf( "pos : %" PRId32 "\n", pos );

    const char* genome  = _dna->get_data();
    int32_t  len        = _dna->get_length();

    dist = 0;

    if ( strand == LEADING )
    {
      //~ printf( "LEADING\n" );
      for ( int16_t i = 0 ; i < PROM_SIZE ; i++ )
      {
        //~ printf( "  i : %" PRId32 " dist : %"PRId8"\n", i, dist );
        if ( genome[(pos+i)%len] != PROM_SEQ[i] )
        {
          dist++;
          if ( dist > PROM_MAX_DIFF )
          {
            //~ printf( "=============================================== END is_promoter\n" );
            return false;
          }
        }
      }
    }
    else // ( strand == LAGGING )
    {
      //~ printf( "LAGGING\n" );
      for ( int16_t i = 0 ; i < PROM_SIZE ; i++ )
      {
        //~ printf( "  i : %"PRId32" dist : %"PRId8"\n", i, dist );
        if ( genome[Utils::mod((pos-i),len)] == PROM_SEQ[i] ) // == and not != because we are on the complementary strand...
        {
          dist++;
          if ( dist > PROM_MAX_DIFF )
          {
            //~ printf( "=============================================== END is_promoter\n" );
            return false;
          }
        }
      }
    }


    //~ printf( "=============================================== END is_promoter\n" );
    return true;
  }

  bool GeneticUnit::is_terminator( Strand strand, int32_t pos ) const
  {
    const char* genome  = _dna->get_data();
    int32_t  len        = _dna->get_length();

    if ( strand == LEADING )
    {
      for ( int16_t i = 0 ; i < TERM_STEM_SIZE ; i++ )
      {
        if ( genome[Utils::mod(pos+i,len)] == genome[Utils::mod(pos+(TERM_SIZE-1)-i,len)] ) return false;
      }
    }
    else // ( strand == LAGGING )
    {
      for ( int16_t i = 0 ; i < TERM_STEM_SIZE ; i++ )
      {
        if ( genome[Utils::mod(pos-i,len)] == genome[Utils::mod(pos-(TERM_SIZE-1)+i,len)] ) return false;
      }
    }

    return true;
  }

  bool GeneticUnit::is_shine_dalgarno( Strand strand, int32_t pos ) const
  {
    const char* genome  = _dna->get_data();
    int32_t  len        = _dna->get_length();

    if ( strand == LEADING )
    {
      for ( int8_t i = 0 ; i < SHINE_DAL_SIZE ; i++ )
      {
        if ( genome[Utils::mod((pos+i),len)] != SHINE_DAL_SEQ[i] )
        {
          return false;
        }
      }
    }
    else // ( strand == LAGGING )
    {
      for ( int8_t i = 0 ; i < SHINE_DAL_SIZE ; i++ )
      {
        if ( genome[Utils::mod((pos-i),len)] == SHINE_DAL_SEQ[i] ) // == and not != because we are on the complementary strand...
        {
          return false;
        }
      }
    }

    return true;
  }

  int8_t GeneticUnit::get_codon( Strand strand, int32_t pos ) const
  {
    const char* genome  = _dna->get_data();
    int32_t  len        = _dna->get_length();
    int8_t codon        = 0;

    if ( strand == LEADING )
    {
      for ( int8_t i = 0 ; i < CODON_SIZE ; i++ )
      {
        if ( genome[Utils::mod((pos+i),len)] == '1' )
        {
          codon += 1 << ( CODON_SIZE - i - 1 ); //pow( 2, CODON_SIZE - i - 1 );
        }
      }
    }
    else // ( strand == LAGGING )
    {
      for ( int8_t i = 0 ; i < CODON_SIZE ; i++ )
      {
        if ( genome[Utils::mod((pos-i),len)] != '1' ) // == and not != because we are on the complementary strand...
        {
          codon += 1 << ( CODON_SIZE - i - 1 ); //pow( 2, CODON_SIZE - i - 1 );
        }
      }
    }

    return codon;
  }

  void GeneticUnit::compute_non_coding( void )
  {
    if (_non_coding_computed) return;
    _non_coding_computed = true;

    // Create a table of <genome_length> bools initialized to false (non-coding)
    int32_t genome_length = _dna->get_length();

    // Including Shine-Dalgarno, spacer, START and STOP
    bool* belongs_to_CDS;
    bool* belongs_to_functional_CDS;
    bool* belongs_to_non_functional_CDS; // non-functional CDSs are those that have a null area or that lack a kind of codons (M, W or H)

    // Including Promoters and terminators
    bool* belongs_to_RNA;
    bool* belongs_to_coding_RNA;
    bool* belongs_to_non_coding_RNA;

    // Genes + prom + term (but not UTRs)
    bool* is_essential_DNA;
    bool* is_essential_DNA_including_nf_genes; // Adds non-functional genes + promoters & terminators

    bool* is_not_neutral;            // prom + term + everything in between (as opposed to neutral)


    belongs_to_CDS                      = new bool[genome_length];
    belongs_to_functional_CDS           = new bool[genome_length];
    belongs_to_non_functional_CDS       = new bool[genome_length];
    belongs_to_RNA                      = new bool[genome_length];
    belongs_to_coding_RNA               = new bool[genome_length];
    belongs_to_non_coding_RNA           = new bool[genome_length];
    is_essential_DNA                    = new bool[genome_length];
    is_essential_DNA_including_nf_genes = new bool[genome_length];
    is_not_neutral                      = new bool[genome_length];

    memset( belongs_to_CDS,                      0, genome_length );
    memset( belongs_to_functional_CDS,           0, genome_length );
    memset( belongs_to_non_functional_CDS,       0, genome_length );
    memset( belongs_to_RNA,                      0, genome_length );
    memset( belongs_to_coding_RNA,               0, genome_length );
    memset( belongs_to_non_coding_RNA,           0, genome_length );
    memset( is_essential_DNA,                    0, genome_length );
    memset( is_essential_DNA_including_nf_genes, 0, genome_length );
    memset( is_not_neutral,                      0, genome_length );


    // Parse protein lists and mark the corresponding bases as coding
    for (auto strand: {LEADING, LAGGING}) {
      for (const auto& prot: _protein_list[strand]) {
        int32_t first;
        int32_t last;

        switch (strand) {
          case LEADING:
            first = prot.get_shine_dal_pos();
            last  = prot.get_last_STOP_base_pos();
            break;
          case LAGGING:
            last  = prot.get_shine_dal_pos();
            first = prot.get_last_STOP_base_pos();
            break;
          default:
            assert(false); // error: should never happen
        }

        if ( first <= last )
        {
          for ( int32_t i = first ; i <= last ; i++ )
          {
            belongs_to_CDS[i] = true;
            if ( prot.get_is_functional() ) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
        }
        else
        {
          for ( int32_t i = first ; i < genome_length ; i++ )
          {
            belongs_to_CDS[i] = true;
            if ( prot.get_is_functional() ) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
          for ( int32_t i = 0 ; i <= last ; i++ )
          {
            belongs_to_CDS[i] = true;
            if ( prot.get_is_functional() ) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
        }

        // Include the promoter and terminator to essential DNA
        // Mark everything between promoter and terminator as not neutral
        for (const auto& rna: prot.get_rna_list()) {

          int32_t prom_first;
          int32_t prom_last;
          int32_t term_first;
          int32_t term_last;
          int32_t rna_first;
          int32_t rna_last;

          if ( strand == LEADING )
          {
            prom_first  = rna->get_promoter_pos();
            prom_last   = Utils::mod( prom_first + PROM_SIZE - 1, _dna->get_length() );
            term_last   = rna->get_last_transcribed_pos();
            term_first  = Utils::mod( term_last - TERM_SIZE + 1, _dna->get_length() );
            rna_first   = prom_first;
            rna_last    = term_last;
          }
          else
          {
            prom_last   = rna->get_promoter_pos();
            prom_first  = Utils::mod( prom_last - PROM_SIZE + 1, _dna->get_length() );
            term_first  = rna->get_last_transcribed_pos();
            term_last   = Utils::mod( term_first + TERM_SIZE - 1, _dna->get_length() );
            rna_first   = term_first;
            rna_last    = prom_last;
          }

          // Let us begin with "non-neutral" regions...
          if ( rna_first <= rna_last )
          {
            for ( int32_t i = rna_first ; i <= rna_last ; i++ ) { is_not_neutral[i] = true; }
          }
          else
          {
            for ( int32_t i = rna_first ; i < genome_length ; i++ ) { is_not_neutral[i] = true; }
            for ( int32_t i = 0 ; i <= rna_last ; i++ ) { is_not_neutral[i] = true; }
          }

          // ...and go on with essential DNA
          if ( prom_first <= prom_last )
          {
            for ( int32_t i = prom_first ; i <= prom_last ; i++ )
            {
              //~ printf( "%ld ", i );
              if ( prot.get_is_functional() ) is_essential_DNA[i] = true;
              is_essential_DNA_including_nf_genes[i] = true;
            }
          }
          else
          {
            for ( int32_t i = prom_first ; i < genome_length ; i++ )
            {
              //~ printf( "%ld ", i );
              if ( prot.get_is_functional() ) is_essential_DNA[i] = true;
              is_essential_DNA_including_nf_genes[i] = true;
            }
            for ( int32_t i = 0 ; i <= prom_last ; i++ )
            {
              //~ printf( "%ld ", i );
              if ( prot.get_is_functional() ) is_essential_DNA[i] = true;
              is_essential_DNA_including_nf_genes[i] = true;
            }
          }
          //~ printf( "\n" );

          //~ printf( "term " );
          if ( term_first <= term_last )
          {
            for ( int32_t i = term_first ; i <= term_last ; i++ )
            {
              //~ printf( "%ld ", i );
              if ( prot.get_is_functional() ) is_essential_DNA[i] = true;
              is_essential_DNA_including_nf_genes[i] = true;
            }
          }
          else
          {
            for ( int32_t i = term_first ; i < genome_length ; i++ )
            {
              //~ printf( "%ld ", i );
              if ( prot.get_is_functional() ) is_essential_DNA[i] = true;
              is_essential_DNA_including_nf_genes[i] = true;
            }
            for ( int32_t i = 0 ; i <= term_last ; i++ )
            {
              //~ printf( "%ld ", i );
              if ( prot.get_is_functional() ) is_essential_DNA[i] = true;
              is_essential_DNA_including_nf_genes[i] = true;
            }
          }
          //~ printf( "\n" );
          //~ getchar();
        }




        if ( prot.get_is_functional() )
        {
          if ( first <= last )
          {
            for ( int32_t i = first ; i <= last ; i++ )
            {
              belongs_to_functional_CDS[i] = true;
            }
          }
          else
          {
            for ( int32_t i = first ; i < genome_length ; i++ )
            {
              belongs_to_functional_CDS[i] = true;
            }
            for ( int32_t i = 0 ; i <= last ; i++ )
            {
              belongs_to_functional_CDS[i] = true;
            }
          }
        }
        else // degenerated protein
        {
          if ( first <= last )
          {
            for ( int32_t i = first ; i <= last ; i++ )
            {
              belongs_to_non_functional_CDS[i] = true;
            }
          }
          else
          {
            for ( int32_t i = first ; i < genome_length ; i++ )
            {
              belongs_to_non_functional_CDS[i] = true;
            }
            for ( int32_t i = 0 ; i <= last ; i++ )
            {
              belongs_to_non_functional_CDS[i] = true;
            }
          }
        }
      }
    }


    // Parse RNA lists and mark the corresponding bases as coding (only for the coding RNAs)
    // TODO vld: this block cries for refactoring
    for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ ) {
      for (auto& rna: _rna_list[strand]) {
        int32_t first;
        int32_t last;

        if ( strand == LEADING ) {
          first = rna.get_promoter_pos();
          last  = rna.get_last_transcribed_pos();
        }
        else { // ( strand == LAGGING )
          first = rna.get_last_transcribed_pos();
          last  = rna.get_promoter_pos();
        }

        assert(first < _indiv->get_amount_of_dna());
        assert(last < _indiv->get_amount_of_dna());

        if ( first <= last ) {
          for ( int32_t i = first ; i <= last ; i++ ) {
            belongs_to_RNA[i] = true;
          }
        }
        else {
          for ( int32_t i = first ; i < genome_length ; i++ )
            belongs_to_RNA[i] = true;
          for ( int32_t i = 0 ; i <= last ; i++ )
            belongs_to_RNA[i] = true;
        }

        if (not rna.get_transcribed_proteins().empty()) { // coding RNA
          if ( first <= last ) {
            for ( int32_t i = first ; i <= last ; i++ )
              belongs_to_coding_RNA[i] = true;
          }
          else
          {
            for ( int32_t i = first ; i < genome_length ; i++ )
              belongs_to_coding_RNA[i] = true;
            for ( int32_t i = 0 ; i <= last ; i++ )
              belongs_to_coding_RNA[i] = true;
          }
        }
        else // non coding RNA
        {
          if ( first <= last ) {
            for ( int32_t i = first ; i <= last ; i++ )
              belongs_to_non_coding_RNA[i] = true;
          }
          else {
            for ( int32_t i = first ; i < genome_length ; i++ )
              belongs_to_non_coding_RNA[i] = true;
            for ( int32_t i = 0 ; i <= last ; i++ )
              belongs_to_non_coding_RNA[i] = true;
          }
        }
      }
    }

    // Count non-coding bases
    _nb_bases_in_0_CDS                = 0;
    _nb_bases_in_0_functional_CDS     = 0;
    _nb_bases_in_0_non_functional_CDS = 0;
    _nb_bases_in_0_RNA                = 0;
    _nb_bases_in_0_coding_RNA         = 0;
    _nb_bases_in_0_non_coding_RNA     = 0;
    _nb_bases_non_essential                     = 0;
    _nb_bases_non_essential_including_nf_genes  = 0;
    _nb_bases_in_neutral_regions        = 0;
    _nb_neutral_regions                 = 0;

    // We do not know how many neutral regions there will be, but
    // there should be less than _nb_coding_RNAs + 1
    // As we will see, there may be a shift in values so we take size _nb_coding_RNAs + 2
    int32_t* tmp_beginning_neutral_regions = new int32_t [ _nb_coding_RNAs + 2 ];
    int32_t* tmp_end_neutral_regions = new int32_t [ _nb_coding_RNAs + 2 ];
    memset( tmp_beginning_neutral_regions, -1, _nb_coding_RNAs + 2 );
    memset( tmp_end_neutral_regions, -1, _nb_coding_RNAs + 2 );

    for ( int32_t i = 0 ; i < genome_length ; i++ )
    {
      if ( belongs_to_CDS[i] == false )
      {
        _nb_bases_in_0_CDS++;
      }
      if ( belongs_to_functional_CDS[i] == false )
      {
        _nb_bases_in_0_functional_CDS++;
      }
      if ( belongs_to_non_functional_CDS[i] == false )
      {
        _nb_bases_in_0_non_functional_CDS++;
      }
      if ( belongs_to_RNA[i] == false )
      {
        _nb_bases_in_0_RNA++;
      }
      if ( belongs_to_coding_RNA[i] == false )
      {
        _nb_bases_in_0_coding_RNA++;
      }
      if ( belongs_to_non_coding_RNA[i] == false )
      {
        _nb_bases_in_0_non_coding_RNA++;
      }
      if ( is_essential_DNA[i] == false )
      {
        _nb_bases_non_essential++;
      }
      if ( is_essential_DNA_including_nf_genes[i] == false )
      {
        _nb_bases_non_essential_including_nf_genes++;
      }
      if ( is_not_neutral[i] == false )
      {
        _nb_bases_in_neutral_regions++;
      }
      if ( i != 0 )
      {
        if ( is_not_neutral[i] != is_not_neutral[i-1] )
        {
          if ( is_not_neutral[i-1] == true ) // beginning of a neutral region
          {
            tmp_beginning_neutral_regions [ _nb_neutral_regions ] = i;
          }
          else // end of a neutral region
          {
            tmp_end_neutral_regions [ _nb_neutral_regions ] = i-1;
            _nb_neutral_regions++;
          }
        }
      }
      else // i = 0
      {
        // we arbitrarily set 0 as the beginning of a neutral region (linkage with end of genetic unit
        // will be done later)
        if ( is_not_neutral[0] == false ) { tmp_beginning_neutral_regions [ _nb_neutral_regions ] = 0; }
      }
    }

    // we have to treat specifically the last base of the genetic unit in order to link neutral regions
    // at the end and the beginning of genetic unit
    int32_t shift = 0;
    if ( is_not_neutral[genome_length-1] == false )
    {
      if ( is_not_neutral[0] == true ) // end of a neutral region
      {
        tmp_end_neutral_regions[ _nb_neutral_regions ] = genome_length-1;
        _nb_neutral_regions++;
      }
      else // neutral region goes on after base 0, linkage to be done
      {
        if ( _nb_neutral_regions != 0 )
        {
          tmp_end_neutral_regions[ _nb_neutral_regions ] = tmp_end_neutral_regions[ 0 ];
          // the first neutral region is only a subpart of the last one, it should not be
          // taken into account. When we transfer values to the final array, we introduce a shift
          shift = 1;
          // we do not ++ _nb_neutral_regions as it was already counted
        }
        else // no neutral region detected until now -> all the genetic unit is neutral
        {
          // as all the chromosome is neutral, we indicate 0 as the beginning of the region
          // and genome_length - 1 as its end
          tmp_end_neutral_regions[ 0 ] = genome_length - 1;
          _nb_neutral_regions++;
        }
      }
    }

    // now that we know how many neutral regions there are, we can transfer data to correctly sized arrays
    assert( _nb_neutral_regions <= _nb_coding_RNAs + 1 );
    if ( _beginning_neutral_regions != NULL ) { delete [] _beginning_neutral_regions; }
    if ( _end_neutral_regions != NULL )       { delete [] _end_neutral_regions; }

    if ( _nb_neutral_regions > 0 ) // as unlikely as it seems, there may be no neutral region
    {
      _beginning_neutral_regions = new int32_t [ _nb_neutral_regions ];
      _end_neutral_regions = new int32_t [ _nb_neutral_regions ];
      // transfer from tmp to attributes
      for (int32_t i = 0; i < _nb_neutral_regions; i++)
      {
        _beginning_neutral_regions[i] = tmp_beginning_neutral_regions[i+shift];
        _end_neutral_regions[i]       = tmp_end_neutral_regions[i+shift];
      }
    }
    else // _nb_neutral_regions == 0
    {
      _beginning_neutral_regions = NULL;
      _end_neutral_regions       = NULL;
    }

    delete [] tmp_beginning_neutral_regions;
    delete [] tmp_end_neutral_regions;

    delete [] belongs_to_CDS;
    delete [] belongs_to_functional_CDS;
    delete [] belongs_to_non_functional_CDS;
    delete [] belongs_to_RNA;
    delete [] belongs_to_coding_RNA;
    delete [] belongs_to_non_coding_RNA;
    delete [] is_essential_DNA;
    delete [] is_essential_DNA_including_nf_genes;
    delete [] is_not_neutral;
  }

///
///
/// `duplicated_promoters` is an output parameter, it should be initially empty
  void GeneticUnit::duplicate_promoters_included_in(int32_t pos_1,
                                                    int32_t pos_2,
                                                    Promoters2Strands& duplicated_promoters)
  {
    // 1) Get promoters to be duplicated
    Promoters2Strands retrieved_promoters = {{},{}};
    get_promoters_included_in(pos_1, pos_2, retrieved_promoters);

    // 2) Set RNAs' position as their position on the duplicated segment
    for (auto& strand: {LEADING, LAGGING}) {
      for (auto& rna: retrieved_promoters[strand]) {
        // Make a copy of current RNA inside container
        duplicated_promoters[strand].emplace_back(this, rna);

        // Set RNA's position as it's position on the duplicated segment
        duplicated_promoters[strand].back().shift_position(-pos_1, _dna->get_length());
      }
    }
  }

  void GeneticUnit::get_promoters_included_in(int32_t pos_1,
                                              int32_t pos_2,
                                              Promoters2Strands& promoters)
  {
    assert( pos_1 >= 0 && pos_1 <= _dna->get_length() && pos_2 >= 0 && pos_2 <= _dna->get_length() );

    if ( pos_1 < pos_2 )
    {
      int32_t seg_length = pos_2 - pos_1;

      if ( seg_length >= PROM_SIZE )
      {
        get_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1, promoters[LEADING]);
        get_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1, promoters[LAGGING]);
      }
    }
    else
    {
      int32_t seg_length = _dna->get_length() + pos_2 - pos_1;

      if ( seg_length >= PROM_SIZE )
      {
        bool is_near_end_of_genome        = ( pos_1 + PROM_SIZE > _dna->get_length() );
        bool is_near_beginning_of_genome  = ( pos_2 - PROM_SIZE < 0 );

        if ( !is_near_end_of_genome && !is_near_beginning_of_genome )
        {
          get_promoters(LEADING, AFTER, pos_1, -1, promoters[LEADING] );
          get_promoters(LEADING, BEFORE, -1, pos_2 - PROM_SIZE + 1, promoters[LEADING] );
          get_promoters(LAGGING, AFTER, pos_2, -1, promoters[LAGGING] );
          get_promoters(LAGGING, BEFORE, -1, pos_1 + PROM_SIZE - 1, promoters[LAGGING] );
        }
        else if ( !is_near_end_of_genome ) // => && is_near_beginning_of_genome
        {
          // get_promoters(leading, between, pos_1, pos_2 + _dna->get_length() - PROM_SIZE + 1,
          //                                         promoters[LEADING] );
          get_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1 + _dna->get_length(),
                        promoters[LEADING]);
          get_promoters(LAGGING, AFTER, pos_2, -1, promoters[LAGGING] );
          get_promoters(LAGGING, BEFORE, -1, pos_1 + PROM_SIZE - 1, promoters[LAGGING] );
        }
        else if ( !is_near_beginning_of_genome ) // => && is_near_end_of_genome
        {
          get_promoters(LEADING, AFTER, pos_1, -1, promoters[LEADING] );
          get_promoters(LEADING, BEFORE, -1, pos_2 - PROM_SIZE + 1, promoters[LEADING] );
          get_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1 - _dna->get_length(),
                        promoters[LAGGING] );
        }
        else // is_near_end_of_genome && is_near_beginning_of_genome
        {
          // get_promoters(leading, between, pos_1, pos_2 + _dna->get_length() - PROM_SIZE + 1,
          //                                         promoters[LEADING] );
          get_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1 + _dna->get_length(),
                        promoters[LEADING]);
          get_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1 - _dna->get_length(),
                        promoters[LAGGING] );
        }
      }
    }
  }

/** Get promoters whose starting position are between/before/after
 * pos_1 and pos_2.
 *
 * The promoters will be ordered with regard to the strand's reading direction
 */
  void GeneticUnit::get_promoters(Strand strand_id,
                                  Position before_after_btw, // with regard to the strand's reading direction
                                  int32_t pos1,
                                  int32_t pos2,
                                  Promoters1Strand& promoters) {
    // TODO vld: First try, the parameter list could be cleverer.

    // TODO vld: These find_if puns are not very nice. Could just negate
    // return if LAGGING or something in that spirit.

    assert((before_after_btw == BETWEEN and pos1 >= 0 and pos2 >= 0 and
            pos1 <= _dna->get_length() and pos2 <= _dna->get_length()) or
        (before_after_btw == BEFORE and pos2 >= 0 and pos2 <= _dna->get_length()) or
        (before_after_btw == AFTER and pos1 >= 0 and pos1 <= _dna->get_length()));

    auto strand = _rna_list[strand_id];
    auto it_begin = strand.begin();
    auto it_end   = strand.end();

    if (before_after_btw != BEFORE) {
      it_begin = find_if(strand.begin(),
                         strand.end(),
                         [pos1, strand_id](Rna & p) {
                           if (strand_id == LEADING)
                             return p.get_promoter_pos() >= pos1;
                           else
                             return p.get_promoter_pos() < pos1;
                         });
    }

    if (before_after_btw != AFTER) {
      it_end = find_if(it_begin,
                       strand.end(),
                       [pos2, strand_id](Rna & p) {
                         if (strand_id == LEADING)
                           return p.get_promoter_pos() >= pos2;
                         else
                           return p.get_promoter_pos() < pos2;
                       });
    }

    promoters.insert(promoters.end(), it_begin, it_end);
  }

/// Invert all the promoters of promoter_lists for a sequence of
/// length seq_length.
/*static*/ void GeneticUnit::invert_promoters(Promoters2Strands& promoter_lists,
                                              int32_t seq_length)
  {
    GeneticUnit::invert_promoters(promoter_lists, 0, seq_length);
  }

/// Invert all the promoters of promoter_lists knowing that they
/// represent the promoters of a subsequence beginning at pos_1 and
/// ending at pos_2.
///
/// WARNING : This function is pretty specific, make sure you
/// understand its precise behaviour before using it.
/*static*/ void GeneticUnit::invert_promoters(Promoters2Strands& promoter_lists,
                                              int32_t pos1,
                                              int32_t pos2 )
  {
    assert( pos1 >= 0 && pos1 <= pos2 ); // Could check (pos2 < length) but another parameter would be necessary

    // Exchange LEADING and LAGGING lists
    promoter_lists[LEADING].swap(promoter_lists[LAGGING]);

    // Update the position and strand of each promoter to be inverted...
    for (auto& strand: {LEADING, LAGGING})
      for (auto& rna: promoter_lists[strand]) {
        assert(rna.get_strand() != strand); // strands have just been swapped
        assert(rna.get_promoter_pos() >= pos1);
        assert(rna.get_promoter_pos() < pos2);

        rna.set_promoter_pos(pos1 + pos2 - rna.get_promoter_pos() - 1);
        rna.set_strand(strand);
      }
  }

  void GeneticUnit::invert_promoters_included_in(int32_t pos1,
                                                 int32_t pos2)
  {
    assert(pos1 >= 0);
    assert(pos1 <= pos2);
    assert(pos2 <= _dna->get_length());

    int32_t segment_length = pos2 - pos1;

    if (segment_length < PROM_SIZE)
      return;

    Promoters2Strands inverted_promoters = {{},{}};

    // 1) Extract the promoters completely included on the segment to be inverted
    extract_promoters_included_in(pos1, pos2, inverted_promoters);

    // 2) Invert segment's promoters
    GeneticUnit::invert_promoters(inverted_promoters, pos1, pos2);

    // 3) Reinsert the inverted promoters
    insert_promoters(inverted_promoters);
  }

// TODO vld: should it append extracted promoters to extracted_promoters or replace its content
  void GeneticUnit::extract_leading_promoters_starting_between(int32_t pos_1,
                                                               int32_t pos_2,
                                                               Promoters1Strand& extracted_promoters)
  {
    assert(pos_1 >= 0);
    assert(pos_1 < pos_2);
    assert(pos_2 <= _dna->get_length());

    // Find the first promoters in the interval
    auto& strand = _rna_list[LEADING];
    auto first = find_if(strand.begin(),
                         strand.end(),
                         [pos_1](Rna & p)
                         { return p.get_promoter_pos() >= pos_1; });

    if (first == strand.end() or first->get_promoter_pos() >= pos_2)
      return;

    // Find the last promoters in the interval
    auto end = find_if(first,
                        strand.end(),
                        [pos_2](Rna & p)
                        { return p.get_promoter_pos() >= pos_2; });

    // Extract the promoters (remove them from the individual's list and put them in extracted_promoters)
    extracted_promoters.insert(extracted_promoters.end(), first, end);
    strand.erase(first, end);
  }

  void GeneticUnit::extract_lagging_promoters_starting_between(int32_t pos_1,
                                                               int32_t pos_2,
                                                               Promoters1Strand& extracted_promoters )
  {
    assert(pos_1 >= 0);
    assert(pos_1 < pos_2);
    assert(pos_2 <= _dna->get_length());

    // Find the first promoters in the interval (if any)
    auto& strand = _rna_list[LAGGING];
    auto first = find_if(strand.begin(),
                         strand.end(),
                         [pos_2](Rna & r)
                         { return r.get_promoter_pos() < pos_2; });

    if (first == strand.end() or first->get_promoter_pos() < pos_1)
      return;

    // Find the last promoters in the interval
    auto end = find_if(first,
                        strand.end(),
                        [pos_1](Rna & r)
                        { return r.get_promoter_pos() < pos_1; });

    // Extract the promoters (remove them from the individual's list and put the in extracted_promoters)
    extracted_promoters.insert(extracted_promoters.end(), first, end);
    strand.erase(first, end);
  }




















/// Shift all the promoters in <promoters_to_shift> by <delta_pos> in a sequence of length <seq_length>.
///
/// Every promoter in double stranded list <promoters_to_shift> will
///  be shifted by <delta_pos>, then a modulo <seq_length> will be
///  applied
/*static*/ void GeneticUnit::shift_promoters(Promoters2Strands& promoters_to_shift,
                                             int32_t delta_pos,
                                             int32_t seq_length )
  {
    for (auto& strand: {LEADING, LAGGING})
      for (auto& rna: promoters_to_shift[strand])
        rna.shift_position(delta_pos, seq_length);
  }

/// Insert promoters in double stranded list `promoters_to_insert` into `this->_rna_list`.
///
/// The promoters in <promoters_to_insert> must already be at their
/// rightful position according to <this>
/// and the positions of the promoters from <promoters_to_insert> and
/// <this->_rna_list> must not be interlaced
/// i.e. no promoter in <this->_rna_list> must have a position in
/// [first_prom_to_insert->pos ; last_prom_to_insert->pos]
  void GeneticUnit::insert_promoters(Promoters2Strands& promoters_to_insert)
  {
    // TODO vld: to be merged with insert_promoters_at(...)
    for (auto strand: {LEADING, LAGGING}) {
      if (promoters_to_insert[strand].size() <= 0)
        continue;
      // Get to the right position in individual's list (first promoter after the inserted segment)
      int32_t from_pos = promoters_to_insert[strand].back().get_promoter_pos();

      auto pos = find_if(_rna_list[strand].begin(),
                         _rna_list[strand].end(),
                         [from_pos, strand] (Rna & r) {
                             if (strand == LEADING)
                               return r.get_promoter_pos() >= from_pos;
                             else
                               return r.get_promoter_pos() < from_pos;
                         });

      // Insert the promoters in the individual's RNA list
      for (auto& to_insert: promoters_to_insert[strand])
        // TODO vld: could be compacted in a unique emplace(pos, to_insert) ?
        if (pos != _rna_list[strand].end())
          _rna_list[strand].insert(pos, to_insert);
        else
          _rna_list[strand].push_back(to_insert);
    }
  }

/// Insert promoters in double stranded list `promoters_to_insert`
/// into `this->_rna_list` at position `pos`
///
/// The promoters in `promoters_to_insert` must be at their rightful
/// position according to a stand-alone sequence (i.e. at a RELATIVE
/// position). Their position will be updated automatically.
  void GeneticUnit::insert_promoters_at(Promoters2Strands& promoters_to_insert,
                                        int32_t pos ) {
    for (auto strand: {LEADING, LAGGING}) {
      if (promoters_to_insert[strand].size() <= 0)
        continue;
      // Get to the right position in individual's list (first promoter after the inserted segment)
      auto first = find_if(_rna_list[strand].begin(),
                           _rna_list[strand].end(),
                           [pos, strand](Rna & r) {
                               if (strand == LEADING)
                                 return r.get_promoter_pos() >= pos;
                               else
                                 return r.get_promoter_pos() < pos; });

      // Insert the promoters in the individual's RNA list
      for (auto& to_insert: promoters_to_insert[strand]) {
        // Update promoter position
        to_insert.shift_position( pos, _dna->get_length() );
        // Insert
        if (first != _rna_list[strand].end())
          _rna_list[strand].insert(first, to_insert);
        else
          _rna_list[strand].push_back(to_insert);
      }
    }
  }


/// Remove the RNAs of the LEADING strand whose starting positions lie
/// in [pos_1 ; pos_2[
  void GeneticUnit::remove_leading_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
    assert(pos_1 >= 0);
    assert(pos_1 < _dna->get_length());
    assert(pos_2 >= 0);
    assert(pos_2 <= _dna->get_length());

    if ( pos_1 > pos_2 ) {
      remove_leading_promoters_starting_after(pos_1);
      remove_leading_promoters_starting_before(pos_2);
    } else {
      auto& strand = _rna_list[LEADING];
      // Delete RNAs until we pass pos_2 (or we reach the end of the list)
      // STL Warning: don't erase the current iterator in the for-loop!
      for (auto it = find_if(strand.begin(),
                             strand.end(),
                             [pos_1](Rna & r) { return r.get_promoter_pos() >= pos_1; }),
               nextit = it;
           it != strand.end() and it->get_promoter_pos() < pos_2;
           it = nextit) {
        nextit = next(it);
        strand.erase(it);
      }
    }
  }


/// Remove the RNAs of the LAGGING strand whose starting positions lie
/// in [pos_1 ; pos_2[
///
/// NOTE : A lagging promoter whose starting position is pos spans
/// [pos-PROM_SIZE+1 ; pos], not [pos-PROM_SIZE ; pos[
///
/// Assuming (PROM_SIZE == 4), the LAGGING promoter whose starting
/// position is pos spans the cells filled with X on the following
/// cartoon:
/// \verbatim
///    -------------------------------
///   |   |   | X | X | X | X |   |   |
///    -------------------------------
///                       ^
///                      pos
/// \endverbatim
  void GeneticUnit::remove_lagging_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
    assert( pos_1 >= 0 && pos_1 <= _dna->get_length() && pos_2 >= 0 && pos_2 <= _dna->get_length() );

    if ( pos_1 == _dna->get_length() ) pos_1 = 0;
    if ( pos_2 == 0 )                  pos_2 = _dna->get_length();
    if ( pos_1 > pos_2 ) { // vld: that's a weird case... really do this? used from remove_promoters_around()
      remove_lagging_promoters_starting_after( pos_1 );
      remove_lagging_promoters_starting_before( pos_2 );
    } else {
      auto& strand = _rna_list[LAGGING];
      // Delete RNAs until we pass pos_1 (or we reach the end of the list)
      for (auto it = find_if(strand.begin(),
                             strand.end(),
                             [pos_2](Rna & r) { return r.get_promoter_pos() < pos_2; }),
               nextit = it;
           it != strand.end() and it->get_promoter_pos() >= pos_1;
           it = nextit) {
        nextit = next(it);
        strand.erase(it);
      }
    }
  }


/// Remove the promoters from the LEADING strand whose starting
/// positions are < pos
  void GeneticUnit::remove_leading_promoters_starting_before( int32_t pos )
  {
    assert( pos >= 0 && pos < _dna->get_length() );

    auto& strand = _rna_list[LEADING];
    // Delete RNAs until we reach pos (or we reach the end of the list)
    for (auto it = strand.begin(),
             nextit = it;
         it != strand.end() and it->get_promoter_pos() < pos;
         it = nextit) {
      nextit = next(it);
      strand.erase(it);
    }
  }


/// Remove the promoters from the LAGGING strand whose starting
/// positions are < pos
///
/// NOTE : A lagging promoter whose starting position is pos spans
/// [pos-PROM_SIZE+1 ; pos], not [pos-PROM_SIZE ; pos[
///
/// Assuming (PROM_SIZE == 4), the LAGGING promoter whose starting
/// position is pos spans the cells filled with X on the following
/// cartoon:
/// \verbatim
///    -------------------------------
///   |   |   | X | X | X | X |   |   |
///    -------------------------------
///                       ^
///                      pos
/// \endverbatim
  void GeneticUnit::remove_lagging_promoters_starting_before(int32_t pos) {
    assert(pos >= 0 and pos < _dna->get_length());

    auto& strand = _rna_list[LAGGING];
    // Delete RNAs until we reach pos (or we reach the end of the list)
    // TODO: optimize by starting from the end (with reverse iterators)
    for (auto it = find_if(strand.begin(),
                           strand.end(),
                           [pos](Rna & r)
                           { return r.get_promoter_pos() < pos; }),
             nextit = it;
         it != strand.end();
         it = nextit) {
      nextit = next(it);
      strand.erase(it);
    }
  }


/// Remove the promoters from the LEADING strand whose starting
/// positions are >= pos
  void GeneticUnit::remove_leading_promoters_starting_after( int32_t pos ) {
    assert( pos >= 0 && pos < _dna->get_length() );

    auto& strand = _rna_list[LEADING];
    // TODO: optimize by starting from the end (with reverse iterators)
    for (auto it = find_if(strand.begin(), strand.end(), [pos](Rna & r)
    { return r.get_promoter_pos() >= pos; }),
             nextit = it;
         it != strand.end();
         it = nextit) {
      nextit = next(it);
      strand.erase(it);
    }
  }


/// Remove the promoters from the LAGGING strand whose starting
/// positions are >= pos
///
/// NOTE : A lagging promoter whose starting position is pos spans
/// [pos-PROM_SIZE+1 ; pos], not [pos-PROM_SIZE ; pos[
///
/// Assuming (PROM_SIZE == 4), the LAGGING promoter whose starting
/// position is pos spans the cells filled with X on the following
/// cartoon:
/// \verbatim
///    -------------------------------
///   |   |   | X | X | X | X |   |   |
///    -------------------------------
///                       ^
///                      pos
/// \endverbatim
  void GeneticUnit::remove_lagging_promoters_starting_after(int32_t pos) {
    assert( pos < _dna->get_length() && pos >= 0 );

    auto& strand = _rna_list[LAGGING];
    // Delete RNAs until we pass pos (or we reach the end of the list)
    for (auto it = strand.begin(),
             nextit = it;
         it != strand.end() and it->get_promoter_pos() >= pos;
         it = nextit) {
      nextit = next(it);
      strand.erase(it);
    }
  }


/// Look for new promoters on the LEADING strand whose starting
/// positions would lie in [pos_1 ; pos_2[
  void GeneticUnit::look_for_new_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2 ) {
    assert( pos_1 >= 0 && pos_1 < _dna->get_length() && pos_2 >= 0 && pos_2 < _dna->get_length() );

    // When pos_1 > pos_2, we will perform the search in 2 steps.
    // As positions  0 and _dna->get_length() are equivalent, it's preferable to
    // keep 0 for pos_1 and _dna->get_length() for pos_2.

    if ( pos_1 >= pos_2 ) {
      look_for_new_leading_promoters_starting_after( pos_1 );
      look_for_new_leading_promoters_starting_before( pos_2 );
      return;
    }
    int8_t dist; // Hamming distance of the sequence from the promoter consensus

    for ( int32_t i = pos_1 ; i < pos_2 ; i++ )
    {
      if ( is_promoter( LEADING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
      {
        //~ char tmp[255];
        //~ memcpy( tmp, &_dna->get_data()[i], PROM_SIZE * sizeof(char) );
        //~ printf( "new promoter found on the LEADING strand at position %"PRId32" : %s\n", i, tmp );

        // Look for the right place to insert the new promoter in the list
        auto& strand = _rna_list[LEADING];
        auto first = find_if(strand.begin(),
                             strand.end(),
                             [i](Rna & r)
                             { return r.get_promoter_pos() >= i; });
        if (first == strand.end() or first->get_promoter_pos() != i) {
#ifndef __REGUL
          _rna_list[LEADING].emplace(first, this, LEADING, i, dist);
#else
          _rna_list[LEADING].emplace(first, this, LEADING, i, dist);
#endif
        }
      }
    }
  }


/*!
  \brief Look for new promoters on the LAGGING strand whose starting positions would lie in [pos_1 ; pos_2[

  NOTE : A lagging promoter whose starting position is pos spans [pos-PROM_SIZE+1 ; pos], not [pos-PROM_SIZE ; pos[
  Assuming (PROM_SIZE == 4), the LAGGING promoter whose starting position is pos spans the cells filled with X on the following cartoon:
  \verbatim
     -------------------------------
    |   |   | X | X | X | X |   |   |
     -------------------------------
                        ^
                       pos
  \endverbatim
*/
  void GeneticUnit::look_for_new_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2 ) {
    assert( pos_1 >= 0 && pos_1 < _dna->get_length() && pos_2 >= 0 && pos_2 < _dna->get_length() );

    // When pos_1 > pos_2, we will perform the search in 2 steps.
    // As positions  0 and _dna->get_length() are equivalent, it's preferable to
    // keep 0 for pos_1 and _dna->get_length() for pos_2.

    if ( pos_1 >= pos_2 ) {
      look_for_new_lagging_promoters_starting_after( pos_1 );
      look_for_new_lagging_promoters_starting_before( pos_2 );
      return;
    }

    int8_t dist; // Hamming distance of the sequence from the promoter consensus
    for ( int32_t i = pos_2 - 1 ; i >= pos_1 ; i-- )
    {
      if ( is_promoter( LAGGING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
      {
        assert ( i >= 0 && i < _dna->get_length() );

        // Look for the right place to insert the new promoter in the list
        auto& strand = _rna_list[LAGGING];
        auto first = find_if(strand.begin(),
                             strand.end(),
                             [i](Rna & r)
                             { return r.get_promoter_pos() <= i; });
        if (first == strand.end() or first->get_promoter_pos() != i) {
#ifndef __REGUL
          _rna_list[LAGGING].emplace(first, this, LAGGING, i, dist);
#else
          _rna_list[LAGGING].emplace(first, this, LAGGING, i, dist);
#endif // __REGUL
        }
      }
    }
  }


/// Look for new promoters on the LEADING strand whose starting positions would be >= pos
  void GeneticUnit::look_for_new_leading_promoters_starting_after( int32_t pos ) {
    assert( pos >= 0 && pos < _dna->get_length() );

    // Hamming distance of the sequence from the promoter consensus
    int8_t dist;

    // rna list node used to find the new promoter's place in the list

    for ( int32_t i = pos ; i < _dna->get_length() ; i++ )
    {
      if (is_promoter(LEADING, i, dist)) { // dist takes the hamming distance of the sequence from the consensus
        // Look for the right place to insert the new promoter in the list
        auto& strand = _rna_list[LEADING];
        auto first = find_if(strand.begin(),
                             strand.end(),
                             [i](Rna & r) {
                                 return r.get_promoter_pos() >= i;
                             });
        if (first == strand.end() or first->get_promoter_pos() != i) {
#ifndef __REGUL
          _rna_list[LEADING].emplace(first, this, LEADING, i, dist);
#else
          _rna_list[LEADING].emplace(first, this, LEADING, i, dist);
#endif
        }
      }
    }
  }


/// Look for new promoters on the LAGGING strand whose starting
/// positions would be >= pos
///
/// NOTE : A lagging promoter whose starting position is pos spans
/// [pos-PROM_SIZE+1 ; pos], not [pos-PROM_SIZE ; pos[
///
/// Assuming (PROM_SIZE == 4), the LAGGING promoter whose starting
/// position is pos spans the cells filled with X on the following
/// cartoon:
/// \verbatim
///    -------------------------------
///   |   |   | X | X | X | X |   |   |
///    -------------------------------
///                       ^
///                      pos
/// \endverbatim
  void GeneticUnit::look_for_new_lagging_promoters_starting_after( int32_t pos )
  {
    assert( pos >= 0 && pos < _dna->get_length() );

    // Hamming distance of the sequence from the promoter consensus
    int8_t dist;
    auto& strand = _rna_list[LAGGING];
    auto first = strand.begin();

    for ( int32_t i = _dna->get_length() - 1 ; i >= pos ; i-- ) {
      if ( is_promoter( LAGGING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
      {
        assert ( i >= 0 && i < _dna->get_length() );
        // Look for the right place to insert the new promoter in the list
        first = find_if(first,
                        strand.end(),
                        [i](Rna & r)
                        { return r.get_promoter_pos() <= i; });
        if (first == strand.end() or first->get_promoter_pos() != i) {
#ifndef __REGUL
          _rna_list[LAGGING].emplace(first, this, LAGGING, i, dist);
#else
          _rna_list[LAGGING].emplace(first, this, LAGGING, i, dist);
#endif
        }
      }
    }
  }


/// Look for new promoters on the LEADING strand whose starting
/// positions would be < pos
  void GeneticUnit::look_for_new_leading_promoters_starting_before( int32_t pos )
  {
    assert( pos >= 0 && pos < _dna->get_length() );

    // Hamming distance of the sequence from the promoter consensus
    int8_t dist;

    auto& strand = _rna_list[LEADING];
    auto first = strand.begin(); // TODO vld: should it not be reset at each loop step?

    for ( int32_t i = 0 ; i < pos ; i++ ) {
      if ( is_promoter( LEADING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
      {
        // Look for the right place to insert the new promoter in the list
        first = find_if(first,
                        strand.end(),
                        [i](Rna & r)
                        { return r.get_promoter_pos() >= i; });
        if (first == strand.end() or first->get_promoter_pos() != i) {
#ifndef __REGUL
          _rna_list[LEADING].emplace(first, this, LEADING, i, dist);
#else
          _rna_list[LEADING].emplace(first, this, LEADING, i, dist);
#endif
        }
      }
    }
  }


/// Look for new promoters on the LAGGING strand whose starting positions would be < pos
///
/// NOTE : A lagging promoter whose starting position is pos spans
/// [pos-PROM_SIZE+1 ; pos], not [pos-PROM_SIZE ; pos[
///
/// Assuming (PROM_SIZE == 4), the LAGGING promoter whose starting
/// position is pos spans the cells filled with X on the following
/// cartoon:
/// \verbatim
///    -------------------------------
///   |   |   | X | X | X | X |   |   |
///    -------------------------------
///                       ^
///                      pos
/// \endverbatim
  void GeneticUnit::look_for_new_lagging_promoters_starting_before( int32_t pos )
  {
    assert( pos >= 0 && pos < _dna->get_length() );

    // Hamming distance of the sequence from the promoter consensus
    int8_t dist;

    // rna list node used to find the new promoter's place in the list
    auto& strand = _rna_list[LAGGING];
    auto first = strand.begin();

    for ( int32_t i = pos - 1 ; i >= 0 ; i-- )
    {
      if ( is_promoter( LAGGING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
      {
        assert ( i >= 0 && i < _dna->get_length() );
        // Look for the right place to insert the new promoter in the list
        first = find_if(first,
                        strand.end(),
                        [i](Rna & r) {
                            return r.get_promoter_pos() <= i;
                        });
        if (first == strand.end() or first->get_promoter_pos() != i) {
#ifndef __REGUL
          _rna_list[LAGGING].emplace(first, this, LAGGING, i, dist);
#else
          _rna_list[LAGGING].emplace(first, this, LAGGING, i, dist);
#endif
        }
      }
    }
  }


/// Shift (by delta_post) the positions of the promoters from the
/// LEADING strand whose starting positions are >= pos.
  void GeneticUnit::move_all_leading_promoters_after( int32_t pos, int32_t delta_pos ) {
    auto& strand = _rna_list[LEADING];
    for (auto rna = find_if(strand.begin(), strand.end(), [pos](Rna & r) { return r.get_promoter_pos() >= pos; });
         rna != strand.end();
         ++rna)
      rna->shift_position( delta_pos, _dna->get_length());
  }


/// Shift (by delta_post) the positions of the promoters from the
/// LAGGING strand whose starting positions are >= pos.
///
/// NOTE : A lagging promoter whose starting position is pos spans
/// [pos-PROM_SIZE+1 ; pos], not [pos-PROM_SIZE ; pos[
///
/// Assuming (PROM_SIZE == 4), the LAGGING promoter whose starting
/// position is pos spans the cells filled with X on the following
/// cartoon:
/// \verbatim
///    -------------------------------
///   |   |   | X | X | X | X |   |   |
///    -------------------------------
///                       ^
///                      pos
/// \endverbatim
  void GeneticUnit::move_all_lagging_promoters_after( int32_t pos, int32_t delta_pos ) {
    auto& strand = _rna_list[LAGGING];
    // Update RNAs until we pass pos (or we reach the end of the list)
    for (auto rna = strand.begin();
         rna != strand.end() and rna->get_promoter_pos() >= pos;
         ++rna)
      rna->shift_position(delta_pos, _dna->get_length());
  }

/// Copy (into new_promoter_list) the promoters from the LEADING
/// strand whose starting positions lie in [pos_1 ; pos_2[
  void GeneticUnit::copy_leading_promoters_starting_between(int32_t pos_1,
                                                            int32_t pos_2,
                                                            Promoters1Strand& new_promoter_list ) {
    // 1) Go to first RNA to copy
    auto& strand = _rna_list[LEADING];
    const auto& first = find_if(strand.begin(),
                                strand.end(),
                                [pos_1](Rna & r) { return r.get_promoter_pos() >= pos_1; });

    // 2) Copy RNAs
    if ( pos_1 < pos_2 ) {
      // Copy from pos_1 to pos_2
      for (auto rna = first; rna != strand.end() and rna->get_promoter_pos() < pos_2; ++rna)
        new_promoter_list.emplace_back(this, *rna);
    }
    else {
      // Copy from pos_1 to the end of the list
      for (auto rna = first; rna != strand.end(); ++rna)
        new_promoter_list.emplace_back(this, *rna);

      // Copy from the beginning of the list to pos_2
      for (auto& rna: strand) {
        if (rna.get_promoter_pos() >= pos_2)
          break;
        new_promoter_list.emplace_back(this, rna);
      }
    }
  }

/// Copy (into new_promoter_list) the promoters from the LAGGING
/// strand whose starting positions lie in [pos_1 ; pos_2[
///
/// NOTE : A lagging promoter whose starting position is pos spans
/// [pos-PROM_SIZE+1 ; pos], not [pos-PROM_SIZE ; pos[
///
/// Assuming (PROM_SIZE == 4), the LAGGING promoter whose starting
/// position is pos spans the cells filled with X on the following
/// cartoon:
/// \verbatim
///    -------------------------------
///   |   |   | X | X | X | X |   |   |
///    -------------------------------
///                       ^
///                      pos
/// \endverbatim
void GeneticUnit::copy_lagging_promoters_starting_between(int32_t pos_1,
                                                          int32_t pos_2,
                                                          Promoters1Strand& new_promoter_list) {
  // Go to first RNA to copy
  auto& strand = _rna_list[LAGGING];
  const auto& first = find_if(strand.rbegin(),
                              strand.rend(),
                              [pos_1](Rna & r) { return r.get_promoter_pos() >= pos_1; });

  // Copy RNAs
  if (pos_1 < pos_2) {
    // Copy from pos_1 to pos_2
    for (auto rna = first; rna != strand.rend() and rna->get_promoter_pos() < pos_2 ; ++rna)
      new_promoter_list.emplace_front(this, *rna);
  }
  else
  {
    // Copy from pos_1 to the beginning of the list (we are going backwards)
    for (auto rna = first; rna != strand.rend(); ++rna)
      new_promoter_list.emplace_front(this, *rna);

    // Copy from the end of the list to pos_2 (we are going backwards)
    for (auto rna = strand.rbegin(); rna != strand.rend() and rna->get_promoter_pos() < pos_2; ++rna)
      new_promoter_list.emplace_front(this, *rna);
  }
}

void GeneticUnit::save(gzFile backup_file) const
{
  _dna->save(backup_file);
  gzwrite(backup_file, &_min_gu_length, sizeof(_min_gu_length));
  gzwrite(backup_file, &_max_gu_length, sizeof(_max_gu_length));
}

int32_t GeneticUnit::get_nb_terminators() {
  int32_t nb_term = 0;
  if (_dna->get_length() >= TERM_SIZE)
    for (int32_t i = 0 ; i < _dna->get_length() ; i++)
      if (is_terminator(LEADING, i))  // No need to count on both the
        // LEADING and the LAGGING
        // strand as terminators are
        // "shared"
        nb_term++;
  return nb_term;
}

#ifdef DEBUG
/// Compare current _rna_list with locate_promoters-generated _rna_list
void GeneticUnit::assert_promoters() {
  // Check that the lists are ordered correctly
  assert_promoters_order();

  // Make a backup of the genetic unit's lists of RNAs
  auto backup = _rna_list;
  _rna_list = {{},{}};
  locate_promoters(); // regenerate _rna_list

  // Compare lists
  for (auto strand: {LEADING, LAGGING}) {
    if(backup[strand].size() != _rna_list[strand].size()) {
      printf("Individual %" PRId32 "\n", _indiv->get_id());
      printf("***************** FOUND *******************");
      print_rnas(backup[strand], strand);
      printf("***************** EXPECTED *******************");
      print_rnas(_rna_list[strand], strand);
      printf("******************************************");
      assert(false);
    }

    auto node_old = backup[strand].begin();
    auto node_new = node_old; // just for the type
    for (node_old = backup[strand].begin(), node_new = _rna_list[strand].begin();
         node_old != backup[strand].end(); // _rna_list is the same size as backup
         ++node_old, ++node_new) {
      // TODO vld: to factor
      if ( node_old->get_strand() != node_new->get_strand() ) {
        printf("Individual %" PRId32 "\n", _indiv->get_id());
        printf( "****************************** STRAND problem ******************************\n" );
        printf( "should be : \n" );
        print_rnas( _rna_list );
        printf( "is : \n" );
        print_rnas( backup );
        printf( "****************************************************************************\n" );
        printf( "  %" PRId32 " (%s) : %f    vs    %" PRId32 " (%s) : %f\n",
                node_old->get_promoter_pos(), StrandName[strand], node_old->get_basal_level(),
                node_new->get_promoter_pos(), StrandName[strand], node_new->get_basal_level() );
        printf( "  genome length : %" PRId32 "\n", _dna->get_length() );
        assert( node_old->get_strand() == node_new->get_strand() );
      }

      if ( node_old->get_promoter_pos() != node_new->get_promoter_pos() ) {
        printf("Individual %" PRId32 "\n", _indiv->get_id());
        printf( "***************************** POSITION problem *****************************\n" );
        printf( "should be : \n" );
        print_rnas( _rna_list );
        printf( "is : \n" );
        print_rnas( backup );
        printf( "****************************************************************************\n" );
        printf( "  %" PRId32 " (%s) : %f    vs    %" PRId32 " (%s) : %f\n",
                node_old->get_promoter_pos(), StrandName[strand], node_old->get_basal_level(),
                node_new->get_promoter_pos(), StrandName[strand], node_new->get_basal_level() );
        printf( "  genome length : %" PRId32 "\n", _dna->get_length() );
        assert( node_old->get_promoter_pos() == node_new->get_promoter_pos()  );
      }

      if ( node_old->get_basal_level() != node_new->get_basal_level() ) {
        printf("Individual %" PRId32 "\n", _indiv->get_id());
        printf( "*************************** BASAL LEVEL problem ****************************\n" );
        printf( "should be : \n" );
        print_rnas( _rna_list );
        printf( "is : \n" );
        print_rnas( backup );
        printf( "****************************************************************************\n" );
        printf( "  %" PRId32 " (%s) : %f    vs    %" PRId32 " (%s) : %f\n",
                node_old->get_promoter_pos(), StrandName[strand], node_old->get_basal_level(),
                node_new->get_promoter_pos(), StrandName[strand], node_new->get_basal_level() );
        printf( "  genome length : %" PRId32 "\n", _dna->get_length() );
        assert( node_old->get_basal_level() == node_new->get_basal_level() );
      }
    }
  }

  _rna_list[LEADING].clear();
  _rna_list[LAGGING].clear();
  _rna_list = backup;
}

void GeneticUnit::assert_promoters_order() {
  for (auto strand: {LEADING, LAGGING}) {
    if (_rna_list[strand].size() < 2)
      continue;

    for (auto it =_rna_list[strand].begin(); next(it) != _rna_list[strand].end(); ++it) {
      if ((strand == LEADING and it->get_promoter_pos() >= next(it)->get_promoter_pos()) or
          (strand == LAGGING and it->get_promoter_pos() <= next(it)->get_promoter_pos())) {
        printf("Individual %" PRId32 "\n", _indiv->get_id());
        printf("********************** ORDER problem (%s) ***********************\n", StrandName[strand]);
        print_rnas();
        printf("****************************************************************************\n");
        assert(false);
      }
    }
  }
}
#endif

/// Retrieve for each base if it belongs or not to coding RNA.
///
/// \return Boolean table of sequence length size with for each base if it belongs or not to coding RNA
bool* GeneticUnit::is_belonging_to_coding_RNA( void ) {
  int32_t genome_length = _dna->get_length();
  bool* belongs_to_coding_RNA = new bool[genome_length];
  memset(belongs_to_coding_RNA,0, genome_length);

  // Parse RNA lists and mark the corresponding bases as coding (only for the coding RNAs)
  for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ ) {
    for (auto& rna: _rna_list[strand]) {
      int32_t first;
      int32_t last;

      if ( strand == LEADING )
      {
        first = rna.get_promoter_pos();
        last  = rna.get_last_transcribed_pos();
      }
      else // ( strand == LAGGING )
      {
        last  = rna.get_promoter_pos();
        first = rna.get_last_transcribed_pos();
      }

      if (not rna.get_transcribed_proteins().empty()) // coding RNA
      {
        if ( first <= last )
        {
          for ( int32_t i = first ; i <= last ; i++ )
          {
            belongs_to_coding_RNA[i] = true;
          }
        }
        else
        {
          for ( int32_t i = first ; i < genome_length ; i++ )
          {
            belongs_to_coding_RNA[i] = true;
          }
          for ( int32_t i = 0 ; i <= last ; i++ )
          {
            belongs_to_coding_RNA[i] = true;
          }
        }
      }
    }
  }
  return belongs_to_coding_RNA;
}

/**
 * Remove the bases that are not in coding RNA
 *
 * Remove the bases that are not in coding RNA and test at each loss that fitness is not changed
 */
void GeneticUnit::remove_non_coding_bases( void)
{
  // TODO <david.parsons@inria.fr> Restore method (deal with checking that the fitness remains untouched)
//  Environment* env = _exp_m->get_env() ;
//
//  reset_expression();
//  locate_promoters();
//  _distance_to_target_computed        = false;
//  _fitness_computed                   = false;
//  compute_phenotypic_contribution();
//  if(_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
//  {
//    compute_distance_to_target(env);
//  }
//  compute_fitness(env);
//  double initial_fitness = get_fitness();
//
//  int32_t genome_length = _dna->get_length();
//  bool* belongs_to_coding_RNA = is_belonging_to_coding_RNA();
//
//  int32_t non_coding_bases_nb = get_nb_bases_in_0_coding_RNA();
//  int32_t start = 0;
//  int32_t end = 0;
//
//
//  for ( int32_t i = genome_length-1 ; i >= 0 ; i-- )
//  {
//    if ( belongs_to_coding_RNA[i] == false )
//    {
//      end = i+1;
//      while((i-1) >=0 && belongs_to_coding_RNA[(i-1)] == false )
//      {
//        i--;
//      }
//      start = i;
//      _dna->remove(start,end);
//
//      locate_promoters();
//      reset_expression();
//      _distance_to_target_computed        = false;
//      _fitness_computed                   = false;
//      compute_phenotypic_contribution();
//      if(_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
//      {
//        compute_distance_to_target(env);
//      }
//      compute_fitness(env);
//      assert(get_fitness()==initial_fitness);
//
//      _non_coding_computed = false;
//      non_coding_bases_nb = get_nb_bases_in_0_coding_RNA();
//    }
//  }
//
//  locate_promoters();
//  reset_expression();
//  _distance_to_target_computed        = false;
//  _fitness_computed                   = false;
//  compute_phenotypic_contribution();
//  if(_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
//  {
//    compute_distance_to_target(env);
//  }
//  compute_fitness(env);
//  assert(get_fitness()==initial_fitness);
//
//  _non_coding_computed = false;
//  non_coding_bases_nb = get_nb_bases_in_0_coding_RNA();
//  assert(non_coding_bases_nb==0);
//
//  delete [] belongs_to_coding_RNA;
}

/**
 * Double the bases that are not in coding RNA
 *
 * Double the bases that are not in coding RNA by addition of random bases and test at each addition that fitness is not changed
 */
void GeneticUnit::double_non_coding_bases(void)
{
  // TODO <david.parsons@inria.fr> Restore method (deal with checking that the fitness remains untouched)
//  Environment* env = _exp_m->get_env() ;
//
//  reset_expression();
//  locate_promoters();
//  _distance_to_target_computed        = false;
//  _fitness_computed                   = false;
//  compute_phenotypic_contribution();
//  if(_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
//  {
//    compute_distance_to_target(env);
//  }
//  compute_fitness(env);
//  double initial_fitness = get_fitness();
//
//  int32_t genome_length = _dna->get_length();
//  bool* belongs_to_coding_RNA = is_belonging_to_coding_RNA();
//
//  _non_coding_computed = false;
//  int32_t inital_non_coding_bases_nb = get_nb_bases_in_0_coding_RNA();
//  int32_t start = 0;
//  int32_t end = 0;
//  int32_t length = 0;
//  int32_t pos = 0;
//  char * random_portion = NULL;
//  bool insertion_ok = false;
//  int32_t non_coding_bases_nb_before_fitness = get_nb_bases_in_0_coding_RNA();
//  int32_t non_coding_bases_nb_after_fitness = get_nb_bases_in_0_coding_RNA();
//
//  int32_t non_coding_bases_nb = get_nb_bases_in_0_coding_RNA();
//
//  for ( int32_t i = genome_length-1 ; i >= 0 ; i-- )
//  {
//    if ( belongs_to_coding_RNA[i] == false )
//    {
//      end = i+1;
//      while((i-1) >=0 && belongs_to_coding_RNA[(i-1)] == false )
//      {
//        i--;
//      }
//      start = i;
//      length = end-start;
//
//      insertion_ok = false;
//      while(not insertion_ok)
//      {
//        random_portion = new char [length+1];
//        for ( int32_t j = 0 ; j < length ; j++ )
//        {
//          random_portion[j] = '0' + _indiv->get_mut_prng()->random( NB_BASE );
//        }
//        random_portion[length] = 0;
//        pos = _indiv->get_mut_prng()->random( length )+start;
//        _dna->insert(pos, random_portion);
//
//        _non_coding_computed = false;
//        non_coding_bases_nb_before_fitness = get_nb_bases_in_0_coding_RNA();
//
//        locate_promoters();
//        reset_expression();
//        _distance_to_target_computed        = false;
//        _fitness_computed                   = false;
//        compute_phenotypic_contribution();
//        if(_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
//        {
//          compute_distance_to_target(env);
//        }
//        compute_fitness(env);
//        assert(get_fitness()==initial_fitness);
//
//        _non_coding_computed = false;
//        non_coding_bases_nb_after_fitness = get_nb_bases_in_0_coding_RNA();
//
//        if (non_coding_bases_nb_before_fitness != non_coding_bases_nb_after_fitness)
//        {
//          _dna->remove(pos, pos + length);
//        }
//        else
//        {
//          insertion_ok = true;
//        }
//
//      }
//      _non_coding_computed = false;
//
//      delete [] random_portion;
//      random_portion = NULL;
//    }
//  }
//
//  locate_promoters();
//  reset_expression();
//  _distance_to_target_computed        = false;
//  _fitness_computed                   = false;
//  compute_phenotypic_contribution();
//  if(_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
//  {
//    compute_distance_to_target(env);
//  }
//  compute_fitness(env);
//  assert(get_fitness()==initial_fitness);
//
//  _non_coding_computed = false;
//  non_coding_bases_nb = get_nb_bases_in_0_coding_RNA();
//  assert(non_coding_bases_nb == 2*inital_non_coding_bases_nb);
//
//  delete [] belongs_to_coding_RNA;
}


// This is an auxiliary function for the method ae_genetic_unit::compute_nb_of_affected_genes.
// The gene should go from 'first' to 'last' in the clockwise sense, implying that 'first' should be:
//  - the beginning of the promoter if the gene is on the leading strand
//  - or the end of the terminator if the gene is on the lagging strand.
  static bool breakpoint_inside_gene(int32_t pos_brkpt, int32_t first, int32_t last)
  {
    if (first < last) // most frequent case
    {
      if( (first <= pos_brkpt) && (pos_brkpt <= last)) {return true;}
      else {return false;}
    }
    else // special case where the RNA overlaps ori
    {
      if( (first <= pos_brkpt) || (pos_brkpt <= last) ) {return true;}
      else {return false;}
    }
  }


// This is an auxiliary function for the method ae_genetic_unit::compute_nb_of_affected_genes.
// It return true if the gene is totally included in the segment [pos1, pos2].
// The gene should go from 'first' to 'last' in the clockwise sense, implying that 'first' should be:
//  - the beginning of the promoter if the gene is on the leading strand
//  - or the end of the terminator if the gene is on the lagging strand.
  static bool gene_totally_in_segment(int32_t pos1, int32_t pos2, int32_t first, int32_t last)
  {
    if ( (first < last)  && (pos1 <= pos2) )
    {
      if ( ((first >= pos1) && (first <= pos2)) && ((last >= pos1) && (last <= pos2)) ) {return true; }
      else {return false;}
    }
    else if ( (first < last) && (pos1 > pos2) )  // mut seg in 2 pieces but not the gene
    {
      if ( (first >= pos1) || (last <= pos2) )  // the gene is either completely in [pos1, genlen-1] or completely in [0, pos2]
      {
        return true;
      }
      else return false;
    }
    else if ( (first > last) && (pos1 <= pos2) )  // gene in two pieces but not mut seg, the gene cannot be totally included
    {
      return false;
    }
    else // both mut seg and the gene are in 2 pieces
    {
      if ((first >= pos1) && (last <= pos2)) {return true;}
      else {return false;}
    }
  }



// WARNING: This method works properly only in the case of a single genetic unit (no plasmid).
// Translocations between different genetic units is not managed.
  void GeneticUnit::compute_nb_of_affected_genes(const Mutation * mut, int & nb_genes_at_breakpoints, int & nb_genes_in_segment, int & nb_genes_in_replaced_segment)
  {
    nb_genes_at_breakpoints = 0;
    nb_genes_in_segment = 0;
    nb_genes_in_replaced_segment = 0;
    int32_t genlen = get_seq_length(); // length of the genetic unit, in bp

    int32_t pos0 = -1, pos1 = -1, pos2 = -1, pos3 = -1, mutlength = -1;
    int32_t pos3donor = -1;
    char * seq = NULL;
    int32_t first, last;
    MutationType type = mut->get_mut_type();
    switch(type) {
      case SWITCH:
        pos0 = dynamic_cast<const PointMutation*>(mut)->pos();
        mutlength = 1;
        break;
      case S_INS : {
        const auto* s_ins = dynamic_cast<const SmallInsertion*>(mut);
        pos0 = s_ins->pos();
        mutlength = s_ins->length();
        break;
      }
      case S_DEL : {
        const auto* s_del = dynamic_cast<const SmallDeletion*>(mut);
        pos0 = s_del->pos();
        mutlength = s_del->length();
        break;
      }
      case DUPL : {
        const auto& dupl = dynamic_cast<const Duplication*>(mut);
        pos1 = dupl->pos1();
        pos2 = Utils::mod(dupl->pos2() - 1, genlen);
        pos0 = dupl->pos3();
        mutlength = dupl->length();
      }
      case DEL : {
        const auto& del = dynamic_cast<const Deletion*>(mut);
        pos1 = del->pos1();
        pos2 = Utils::mod(del->pos2() - 1, genlen);
        mutlength = del->length();
        break;
      }
      case TRANS : {
        const auto& trans = dynamic_cast<const Translocation*>(mut);
        pos1 = trans->pos1();
        pos2 = Utils::mod(trans->pos2() - 1, genlen);
        pos3 = trans->pos3();
        pos0 = trans->pos4();
        mutlength = trans->length();
        break;
      }
      case INV : {
        const auto& inv = dynamic_cast<const Inversion*>(mut);
        pos1 = inv->pos1();
        pos2 = Utils::mod(inv->pos2() - 1, genlen);
        mutlength = inv->length();
        break;
      }
      case INS_HT: {
        const auto& ins_ht = dynamic_cast<const InsertionHT*>(mut);
        pos0 = ins_ht->exogenote_pos();// TODO <david.parsons@inria.fr> weird !
        pos3donor = ins_ht->receiver_pos();// TODO <david.parsons@inria.fr> weird !
        mutlength = ins_ht->length();
        seq = ins_ht->seq();

        // Make a temporary genetic unit and translate it to count how many genes were on the exogenote
        GeneticUnit * tmpunit = new GeneticUnit( _indiv, seq, mutlength);
        tmpunit->do_transcription();
        tmpunit->do_translation();
        nb_genes_in_segment = tmpunit->get_nb_coding_RNAs();

        // Check whether the pos3donor breakpoint killed one or several of those genes, in that case decrement nb_genes_in_segment
        for (auto& strand: tmpunit->_rna_list)
          for (auto& rna: strand) {
            if (not rna.is_coding())
              continue;
            first = rna.get_promoter_pos();
            last = rna.get_last_transcribed_pos();
            if (breakpoint_inside_gene(pos3donor, first, last))
              nb_genes_in_segment --;
          }

        delete tmpunit;
        seq = NULL;
        break;
      }
      case REPL_HT: {
        const auto& repl_ht = dynamic_cast<const ReplacementHT*>(mut);
        pos1 = repl_ht->receiver_pos1();
        pos2 = Utils::mod(repl_ht->receiver_pos2() - 1, genlen);
        mutlength = repl_ht->length();
        seq = repl_ht->seq();

        // Make a temporary genetic unit and translate it to count how many genes were on the donor segment
        GeneticUnit * tmpunit = new GeneticUnit( _indiv, seq, mutlength);
        tmpunit->do_transcription();
        tmpunit->do_translation();
        nb_genes_in_segment = tmpunit->get_nb_coding_RNAs();

        // Remove the genes that overlap ori in this temporary unit, as the transferred segment was actually linear
        for (auto& strand: tmpunit->_rna_list)
          for (auto& rna: strand) {
            if (not rna.is_coding())
              continue;
            first = rna.get_promoter_pos();
            last = rna.get_last_transcribed_pos();
            if (first > last)
              nb_genes_in_segment --;
          }

        delete tmpunit;
        seq = NULL;
        break;
      }
      default:
        fprintf(stderr, "Error: unknown mutation type in ae_genetic_unit::compute_nb_of_affected_genes.\n");
    }

    for (auto strand: {LEADING, LAGGING})
      for (auto& rna: _rna_list[strand]) {
        if (rna.is_coding() == false) continue;

        switch (strand) {
          case LEADING:
            first = rna.get_promoter_pos();
            last = rna.get_last_transcribed_pos();
            break;
          case LAGGING:
            first = rna.get_last_transcribed_pos();
            last = rna.get_promoter_pos();
        };

        // TODO vld: reoder lines (if invariant) in cases DUPL and
        // DEL/INV and merge with S_DEL thanks to fallthrough (etc)
        switch(type) {
          case SWITCH: // fall through
          case S_INS:  // fall through
          case INSERT: // fall through
          case INS_HT: // fall through
          case S_DEL:
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          case DUPL:
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_segment ++;
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          case DEL:
          case INV:
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_segment ++;
            if (breakpoint_inside_gene(pos1, first, last)) nb_genes_at_breakpoints ++;
            else if (breakpoint_inside_gene(pos2, first, last)) nb_genes_at_breakpoints ++;     // else because the gene must not be counted twice if both p1 and p2 are in the same gene
            break;
          case TRANS:
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_segment ++;
            if (breakpoint_inside_gene(pos1, first, last)) nb_genes_at_breakpoints ++;   // beginning of the excised segment
            else if (breakpoint_inside_gene(pos2, first, last)) nb_genes_at_breakpoints ++;   // end of the excised segment
            else if (breakpoint_inside_gene(pos3, first, last)) nb_genes_at_breakpoints ++;   // breakpoint inside the segment for the reinsertion
            else if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;   // reinsertion point in the genetic unit
            break;
          case REPL_HT:
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_replaced_segment ++; // the whole gene sequence was replaced by the donor DNA
            if (breakpoint_inside_gene(pos1, first, last)) nb_genes_at_breakpoints ++;  // the gene was disrupted by the breakpoint p1
            else if (breakpoint_inside_gene(pos2, first, last)) nb_genes_at_breakpoints ++;  // the gene was disrupted by the breakpoint p2
            break;
          default:
            // Only simple mutation types are considered.
            break;
        }
      }

    //  if (type == REPL_HT) printf("%d genes in the replaced segment, %d in the donor\n", nb_genes_in_replaced_segment, nb_genes_in_segment);

    if (seq != NULL) delete [] seq;
  }



// =================================================================
//                           Protected Methods
// =================================================================
  void GeneticUnit::init_statistical_data( void ) // TODO : integrate into compute_statistical_data
  {
    //~ _nb_promoters[LEADING]        = 0;
    //~ _nb_promoters[LAGGING]        = 0;
    //~ _nb_genes[LEADING]            = 0;
    //~ _nb_genes[LAGGING]            = 0;
    //~ _average_gene_size            = 0;
    //~ _average_functional_gene_size = 0;
    //~ _nb_coding_bp                 = 0;
    //~ _clustering                   = 0;

    _nb_coding_RNAs               = 0;
    _nb_non_coding_RNAs           = 0;
    _overall_size_coding_RNAs     = 0;
    _overall_size_non_coding_RNAs = 0;
    _nb_genes_activ               = 0;
    _nb_genes_inhib               = 0;
    _nb_fun_genes                 = 0;
    _nb_non_fun_genes             = 0;
    _overall_size_fun_genes       = 0;
    _overall_size_non_fun_genes   = 0;

    _nb_bases_in_0_CDS                = -1;
    _nb_bases_in_0_functional_CDS     = -1;
    _nb_bases_in_0_non_functional_CDS = -1;
    _nb_bases_in_0_RNA                = -1;
    _nb_bases_in_0_coding_RNA         = -1;
    _nb_bases_in_0_non_coding_RNA     = -1;
    _nb_bases_non_essential                     = -1;
    _nb_bases_non_essential_including_nf_genes  = -1;

    _modularity = -1;

    _beginning_neutral_regions = NULL;
    _end_neutral_regions = NULL;
  }

} // namespace aevol
