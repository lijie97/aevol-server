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
  return exp_m_;
}

Individual* GeneticUnit::get_indiv() const
{
  return indiv_;
}

Dna* GeneticUnit::get_dna() const
{
  assert(dna_->length() != 0);
  return dna_;
}

const Promoters2Strands& GeneticUnit::get_rna_list() const {
  return rna_list_;
}

std::list<Protein>& GeneticUnit::get_protein_list(Strand strand) {
  return protein_list_[strand];
}

Fuzzy* GeneticUnit::get_activ_contribution() const
{
  return activ_contribution_;
}

Fuzzy* GeneticUnit::get_inhib_contribution() const
{
  return inhib_contribution_;
}

Fuzzy* GeneticUnit::get_phenotypic_contribution() const
{
  assert(phenotypic_contribution_ != NULL);
  return phenotypic_contribution_;
}

/*!
  Returns the DNA sequence
*/
const char* GeneticUnit::get_sequence( void ) const
{
  return dna_->data();
}

/*!
  Returns the DNA sequence length
*/
int32_t GeneticUnit::get_seq_length( void ) const
{
  return dna_->length();
}

int32_t GeneticUnit::get_nb_coding_RNAs( void ) const
{
  return nb_coding_RNAs_;
}

int32_t GeneticUnit::get_nb_non_coding_RNAs( void ) const
{
  return nb_non_coding_RNAs_;
}

double GeneticUnit::get_overall_size_coding_RNAs( void ) const
{
  return overall_size_coding_RNAs_;
}

double GeneticUnit::get_av_size_coding_RNAs( void ) const
{
  if ( nb_coding_RNAs_ != 0 )
  {
    return overall_size_coding_RNAs_ / nb_coding_RNAs_;
  }
  else return 0.0;
}

double GeneticUnit::get_overall_size_non_coding_RNAs( void ) const
{
  return overall_size_non_coding_RNAs_;
}

double GeneticUnit::get_av_size_non_coding_RNAs( void ) const
{
  if ( nb_non_coding_RNAs_ != 0 )
  {
    return overall_size_non_coding_RNAs_ / nb_non_coding_RNAs_;
  }
  else return 0.0;
}

int32_t GeneticUnit::get_nb_genes_activ( void ) const
{
  return nb_genes_activ_;
}

int32_t GeneticUnit::get_nb_genes_inhib( void ) const
{
  return nb_genes_inhib_;
}

int32_t GeneticUnit::get_nb_functional_genes( void ) const
{
  return nb_fun_genes_;
}

int32_t GeneticUnit::get_nb_non_functional_genes( void ) const
{
  return nb_non_fun_genes_;
}

double GeneticUnit::get_overall_size_functional_genes( void ) const
{
  return overall_size_fun_genes_;
}

double GeneticUnit::get_av_size_functional_genes( void ) const
{
  if ( nb_fun_genes_ != 0 )
  {
    return overall_size_fun_genes_ / nb_fun_genes_;
  }
  else return 0.0;
}

double GeneticUnit::get_overall_size_non_functional_genes( void ) const
{
  return overall_size_non_fun_genes_;
}

double GeneticUnit::get_av_size_non_functional_genes( void ) const
{
  if ( nb_non_fun_genes_ != 0 )
  {
    return overall_size_non_fun_genes_ / nb_non_fun_genes_;
  }
  else return 0.0;
}

int32_t GeneticUnit::get_nb_bases_in_0_CDS() const
{
  assert (non_coding_computed_);
  return nb_bases_in_0_CDS_;
}

int32_t GeneticUnit::get_nb_bases_in_0_functional_CDS() const
{
  assert (non_coding_computed_);
  return nb_bases_in_0_functional_CDS_;
}

int32_t GeneticUnit::get_nb_bases_in_0_non_functional_CDS() const
{
  assert (non_coding_computed_);
  return nb_bases_in_0_non_functional_CDS_;
}

int32_t GeneticUnit::get_nb_bases_in_0_RNA() const
{
  assert (non_coding_computed_);
  return nb_bases_in_0_RNA_;
}

int32_t GeneticUnit::get_nb_bases_in_0_coding_RNA() const
{
  assert (non_coding_computed_);
  return nb_bases_in_0_coding_RNA_;
}

int32_t GeneticUnit::get_nb_bases_in_0_non_coding_RNA() const
{
  assert (non_coding_computed_);
  return nb_bases_in_0_non_coding_RNA_;
}

int32_t GeneticUnit::get_nb_bases_non_essential() const
{
  assert (non_coding_computed_);
  return nb_bases_non_essential_;
}

int32_t GeneticUnit::get_nb_bases_non_essential_including_nf_genes() const
{
  assert (non_coding_computed_);
  return nb_bases_non_essential_including_nf_genes_;
}

int32_t GeneticUnit::get_nb_bases_in_neutral_regions() const
{
  assert (non_coding_computed_);
  return nb_bases_in_neutral_regions_;
}

int32_t GeneticUnit::get_nb_neutral_regions() const
{
  assert (non_coding_computed_);
  return nb_neutral_regions_;
}

int32_t* GeneticUnit::get_beginning_neutral_regions() const
{
  assert (non_coding_computed_);
  return beginning_neutral_regions_;
}

int32_t* GeneticUnit::get_end_neutral_regions() const
{
  assert (non_coding_computed_);
  return end_neutral_regions_;
}

double GeneticUnit::get_modularity() const
{
  return modularity_;
}

double GeneticUnit::get_dist_to_target_by_feature( PhenotypicFeature feature ) const
{
  assert(distance_to_target_computed_);

  return dist_to_target_by_feature_[feature];
}

double GeneticUnit::get_fitness() const
{
  assert(fitness_computed_);

  return fitness_;
}

double GeneticUnit::get_fitness_by_feature( PhenotypicFeature feature ) const
{
  assert(fitness_computed_);

  return fitness_by_feature_[feature];
}

int32_t GeneticUnit::get_min_gu_length( void ) const
{
  return min_gu_length_;
}

int32_t GeneticUnit::get_max_gu_length( void ) const
{
  return max_gu_length_;
}

void GeneticUnit::set_min_gu_length( int32_t min_gu_length )
{
  min_gu_length_ = min_gu_length;
}

void GeneticUnit::set_max_gu_length( int32_t max_gu_length )
{
  max_gu_length_ = max_gu_length;
}

void GeneticUnit::set_exp_m( ExpManager * exp_m )
{
  exp_m_ = exp_m;
}

// =====================================================================
//                       functions' definition
// =====================================================================
void GeneticUnit::print_rnas() const {
  print_rnas(rna_list_);
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
  rna_list_[LEADING].clear();
  rna_list_[LAGGING].clear();
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
  assert(pos_2 <= dna_->length());

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
  if (dna_->length() >= PROM_SIZE )
  {
    remove_leading_promoters_starting_between(Utils::mod(pos - PROM_SIZE + 1,
                                                         dna_->length()),
                                              pos );
    remove_lagging_promoters_starting_between(pos,
                                              Utils::mod(pos + PROM_SIZE - 1,
                                                         dna_->length()));
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
  if ( Utils::mod(pos_1 - pos_2, dna_->length()) >= PROM_SIZE )
  {
    remove_leading_promoters_starting_between( Utils::mod(pos_1 - PROM_SIZE + 1,
                                                          dna_->length()), pos_2 );
    remove_lagging_promoters_starting_between( pos_1, Utils::mod(pos_2 + PROM_SIZE - 1,
                                                                 dna_->length()) );
  }
  else
  {
    remove_all_promoters();
  }
}


/// Look for promoters that are astride pos and add them to the list
/// of promoters (rna_list_).
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
  assert(pos >= 0 && pos <= dna_->length());

  if (dna_->length() >= PROM_SIZE)
  {
    look_for_new_leading_promoters_starting_between(
        Utils::mod(pos - PROM_SIZE + 1, dna_->length()),
        pos);
    look_for_new_lagging_promoters_starting_between(
        pos,
        Utils::mod(pos + PROM_SIZE - 1, dna_->length()));
  }
}


/// Look for promoters that contain at least 1 base lying in [pos_1 ;
/// pos_2[ and add them to the list of promoters (rna_list_).
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
  //~ if ( Utils::mod( pos_1 - pos_2, dna_->length()) == PROM_SIZE - 1 )
  //~ {
  //~ // We have to look at every possible position on the genome.
  //~ locate_promoters();
  //~ }
  /*else*/ if (dna_->length() >= PROM_SIZE )
  {
    look_for_new_leading_promoters_starting_between( Utils::mod(pos_1 - PROM_SIZE + 1,
                                                                dna_->length()), pos_2 );
    look_for_new_lagging_promoters_starting_between( pos_1, Utils::mod(pos_2 + PROM_SIZE - 1,
                                                                       dna_->length()) );
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
  if ( Utils::mod( pos_2 - pos_1 - 1, dna_->length() ) + 1 >= PROM_SIZE )
  {
    copy_leading_promoters_starting_between( pos_1, Utils::mod( pos_2 - PROM_SIZE + 1,
                                                                dna_->length() ), new_promoter_lists[LEADING] );
    copy_lagging_promoters_starting_between( Utils::mod( pos_1 + PROM_SIZE - 1,
                                                         dna_->length() ), pos_2, new_promoter_lists[LAGGING] );
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
  indiv_ = indiv;
  exp_m_ = indiv->get_exp_m();

  transcribed_ = false;
  translated_ = false;
  phenotypic_contributions_computed_ = false;
  non_coding_computed_ = false;
  distance_to_target_computed_ = false;
  fitness_computed_ = false;

  min_gu_length_ = -1;
  max_gu_length_ = -1;

  dna_ = new Dna( this, length, prng );

  // Create empty fuzzy sets for the phenotypic contributions
  activ_contribution_ = new Fuzzy();
  inhib_contribution_ = new Fuzzy();
  phenotypic_contribution_ = NULL;
  // NB : phenotypic_contribution_ is only an indicative value,
  //      it is not used for the whole phenotype computation

  // dist_to_target_per_segment_ depends on the segmentation of the environment
  // and will hence be newed at evaluation time
  dist_to_target_per_segment_ = NULL;

  dist_to_target_by_feature_ = new double [NB_FEATURES];
  fitness_by_feature_ = new double [NB_FEATURES];

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    dist_to_target_by_feature_[i] = 0.0;
    fitness_by_feature_[i]        = 0.0;
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
  exp_m_ = indiv->get_exp_m();
  indiv_ = indiv;

  transcribed_ = false;
  translated_ = false;
  phenotypic_contributions_computed_ = false;
  non_coding_computed_ = false;
  distance_to_target_computed_ = false;
  fitness_computed_ = false;

  min_gu_length_ = -1;
  max_gu_length_ = -1;

  dna_ = new Dna( this, seq, length );

  if (not prom_list[LEADING].empty() and not prom_list[LAGGING].empty()) { // if not default `prom_list`
    // Copy rna lists
    rna_list_ = prom_list;
    Dna::set_GU(rna_list_, this);
  } else {
    for (auto& strand: {LEADING, LAGGING})
      assert(rna_list_[strand].empty());

    // Look for promoters
    locate_promoters();
  }

  // Create empty fuzzy sets for the phenotypic contributions
  activ_contribution_ = new Fuzzy();
  inhib_contribution_ = new Fuzzy();
  phenotypic_contribution_ = NULL;
  // NB : phenotypic_contribution_ is only an indicative value,
  //      it is not used for the whole phenotype computation

  // Initialize all the fitness-related stuff
  dist_to_target_per_segment_ = NULL;

  dist_to_target_by_feature_ = new double [NB_FEATURES];
  fitness_by_feature_ = new double [NB_FEATURES];

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    dist_to_target_by_feature_[i] = 0.0;
    fitness_by_feature_[i]        = 0.0;
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
  exp_m_ = indiv->get_exp_m();
  indiv_ = indiv;

  transcribed_ = false;
  translated_ = false;
  phenotypic_contributions_computed_ = false;
  non_coding_computed_ = false;
  distance_to_target_computed_ = model.distance_to_target_computed_;
  fitness_computed_ = model.fitness_computed_;

  min_gu_length_ = model.min_gu_length_;
  max_gu_length_ = model.max_gu_length_;

  // Copy DNA
  dna_ = new Dna( this, *(model.dna_) );

  // Create empty fuzzy sets for the phenotypic contributions
  activ_contribution_ = new Fuzzy();
  inhib_contribution_ = new Fuzzy();
  phenotypic_contribution_ = NULL;
  // NB : phenotypic_contribution_ is only an indicative value, not used for the whole phenotype computation
  dist_to_target_per_segment_ = NULL;
  dist_to_target_by_feature_ = new double [NB_FEATURES];
  fitness_by_feature_ = new double [NB_FEATURES];

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    dist_to_target_by_feature_[i] = 0.0;
    fitness_by_feature_[i]        = 0.0;
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
  exp_m_ = indiv->get_exp_m();
  indiv_ = indiv;

  transcribed_ = false;
  translated_ = false;
  phenotypic_contributions_computed_ = false;
  non_coding_computed_ = false;
  distance_to_target_computed_ = false;
  fitness_computed_ = false;

  min_gu_length_ = parent->min_gu_length_;
  max_gu_length_ = parent->max_gu_length_;

  // Copy DNA
  dna_ = new Dna( this, parent->dna_);

  // Copy promoter list (rna_list_)
  // Note that the length of the RNA will have to be recomputed (do_transcription)
  for (auto& strand: {LEADING, LAGGING}) {
    for (auto& rna: parent->rna_list_[strand]) {
#ifndef __REGUL
      rna_list_[strand].emplace_back(this, rna);
#else
      rna_list_[strand].emplace_back(this, (dynamic_cast<ae_rna_R>(rna)));
#endif
    }
  }

  // Create empty fuzzy sets for the phenotypic contributions
  activ_contribution_ = new Fuzzy();
  inhib_contribution_ = new Fuzzy();
  phenotypic_contribution_ = NULL;
  // NB : phenotypic_contribution_ is only an indicative value, not used for the whole phenotype computation

  // Initialize all the fitness-related stuff
  dist_to_target_per_segment_ = NULL;

  dist_to_target_by_feature_ = new double [NB_FEATURES];
  fitness_by_feature_ = new double [NB_FEATURES];

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    dist_to_target_by_feature_[i] = 0.0;
    fitness_by_feature_[i]        = 0.0;
  }

  init_statistical_data();
}

GeneticUnit::GeneticUnit( Individual * indiv, gzFile backup_file )
{
  exp_m_ = indiv->get_exp_m();
  indiv_ = indiv;

  transcribed_ = false;
  translated_ = false;
  phenotypic_contributions_computed_ = false;
  non_coding_computed_ = false;
  distance_to_target_computed_ = false;
  fitness_computed_ = false;

  dna_ = new Dna( this, backup_file );

  gzread( backup_file, &min_gu_length_, sizeof(min_gu_length_) );
  gzread( backup_file, &max_gu_length_, sizeof(max_gu_length_) );

  // Create empty fuzzy sets for the phenotypic contributions
  activ_contribution_ = new Fuzzy();
  inhib_contribution_ = new Fuzzy();
  phenotypic_contribution_ = NULL;
  // NB : phenotypic_contribution_ is only an indicative value, not used for the whole phenotype computation

  // Initialize all the fitness-related stuff
  dist_to_target_per_segment_ = NULL;

  dist_to_target_by_feature_ = new double [NB_FEATURES];
  fitness_by_feature_ = new double [NB_FEATURES];

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    dist_to_target_by_feature_[i] = 0.0;
    fitness_by_feature_[i]        = 0.0;
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
  exp_m_ = indiv->get_exp_m();
  indiv_ = indiv;

  transcribed_ = false;
  translated_ = false;
  phenotypic_contributions_computed_ = false;
  non_coding_computed_ = false;
  distance_to_target_computed_ = false;
  fitness_computed_ = false;

  dna_ = new Dna( this, organism_file_name );

  // Create empty fuzzy sets for the phenotypic contributions
  activ_contribution_ = new Fuzzy();
  inhib_contribution_ = new Fuzzy();
  phenotypic_contribution_ = NULL;
  // NB : phenotypic_contribution_ is only an indicative value,
  //      it is not used for the whole phenotype computation

  // Initialize all the fitness-related stuff
  dist_to_target_per_segment_ = NULL;

  dist_to_target_by_feature_ = new double [NB_FEATURES];
  fitness_by_feature_ = new double [NB_FEATURES];

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    dist_to_target_by_feature_[i] = 0.0;
    fitness_by_feature_[i]        = 0.0;
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
  delete dna_;
  delete activ_contribution_;
  delete inhib_contribution_;
  if ( phenotypic_contribution_ != NULL ) delete phenotypic_contribution_;

  if ( dist_to_target_per_segment_ != NULL ) delete [] dist_to_target_per_segment_;

  assert( dist_to_target_by_feature_ != NULL );
  delete [] dist_to_target_by_feature_;
  assert( fitness_by_feature_ != NULL );
  delete [] fitness_by_feature_;

  delete [] beginning_neutral_regions_;
  delete [] end_neutral_regions_;
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
  for (auto& strand: rna_list_) {
    strand.clear();
  }

  if (dna_->length() < PROM_SIZE )
    return;
  for ( int32_t i = 0 ; i < dna_->length() ; i++ ) {
#ifndef __REGUL
    if (is_promoter(LEADING, i, dist)) // dist takes the hamming distance of the sequence from the consensus
      rna_list_[LEADING].emplace_back(this, LEADING, i, dist);
    if (is_promoter(LAGGING, dna_->length() - i - 1, dist))
      rna_list_[LAGGING].emplace_back(this, LAGGING, dna_->length() - i - 1, dist);
#else
    if (is_promoter(LEADING, i, dist))
      _rna_list[LEADING].emplace_back(this, LEADING, i, dist);
    if (is_promoter(LAGGING, _dna->length() - i - 1, dist))
      _rna_list[LAGGING].emplace_back(this, LAGGING, _dna->length() - i - 1, dist);
#endif
  }
}

void GeneticUnit::do_transcription( void )
{
  if (transcribed_) return;
  transcribed_ = true;

  int32_t transcript_start  = -1;
  int32_t genome_length     = dna_->length();

  // If the genome is not long enough to bear a promoter and a terminator,
  // we set all its RNAs to a length of -1
  // (NB but a terminator can share code with the promoter, making it
  // possible for the genome to be no longer than the promoter)
  if (genome_length < PROM_SIZE) {
    for (auto& strand: rna_list_)
      for (auto& rna: strand)
        rna.set_transcript_length(-1);
    return;
  }

  for (auto& strand_id: {LEADING, LAGGING}) {
    auto& strand = rna_list_[strand_id];
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
  if (translated_)
    return;
  translated_ = true;
  if (not transcribed_)
    do_transcription();

  int32_t transcript_start  = -1;
  int32_t transcript_length = -1;
  int32_t genome_length     = dna_->length();

  for (auto strand: {LEADING, LAGGING}) {
    for (auto& rna: rna_list_[strand]) {
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
          auto& protein_strand = protein_list_[strand];
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
              auto codon = new Codon(dna_,
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
                    nb_fun_genes_++;
                    overall_size_fun_genes_ += protein.get_length() * CODON_SIZE;

                    if (protein.get_height() > 0) nb_genes_activ_++;
                    else                           nb_genes_inhib_++;
                  }
                  else {
                    nb_non_fun_genes_++;
                    overall_size_non_fun_genes_ += (strand == LEADING) ?
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
        nb_coding_RNAs_++;
        overall_size_coding_RNAs_ += rna.get_transcript_length();
      }
      else { // non-coding RNA
        nb_non_coding_RNAs_++;
        overall_size_non_coding_RNAs_ += rna.get_transcript_length();
      }
    }
  }
}

void GeneticUnit::compute_phenotypic_contribution( void )
{
  if (phenotypic_contributions_computed_) return;
  phenotypic_contributions_computed_ = true;
  if ( !translated_) do_translation();

  for (const auto& strand: protein_list_) // two strands: LEADING & LAGGING
    for (const auto& prot: strand)
      if (prot.get_is_functional())
        ((prot.get_height() > 0) ? activ_contribution_ : inhib_contribution_)
            ->add_triangle(prot.get_mean(),
                           prot.get_width(),
                           prot.get_height() * prot.get_concentration());
  // if ( prot->get_height() > 0 )
  //   activ_contribution_->add_triangle(prot->get_mean(),
  //                                     prot->get_width(),
  //                                     prot->get_height() * prot->get_concentration() );
  // else
  //   inhib_contribution_->add_triangle(prot->get_mean(),
  //                                     prot->get_width(),
  //                                     prot->get_height() * prot->get_concentration() );

  // It is not necessary to add a lower bound to activ_contribution_ as there can be no negative y
  // The same goes for the upper bound for inhib_contribution_
  activ_contribution_->clip(Fuzzy::max,   Y_MAX );
  inhib_contribution_->clip(Fuzzy::min, - Y_MAX );
  activ_contribution_->simplify();
  inhib_contribution_->simplify();

  if ( exp_m_->get_output_m()->get_compute_phen_contrib_by_GU() )
  {
    phenotypic_contribution_ = new Phenotype();
    phenotypic_contribution_->add( *activ_contribution_);
    phenotypic_contribution_->add( *inhib_contribution_);
    phenotypic_contribution_->simplify();
  }
}

/*!
  \brief Compute the areas between the phenotype and the environment for each environmental segment.

  If the environment is not segmented, the total area is computed
*/
void GeneticUnit::compute_distance_to_target(const PhenotypicTarget& target) {
  if (distance_to_target_computed_) return; // _distance_to_target has already been computed, nothing to do.
  distance_to_target_computed_ = true;

  compute_phenotypic_contribution();

  // Compute the difference between the (whole) phenotype and the environment
  Fuzzy* delta = new Fuzzy( *phenotypic_contribution_);
  delta->sub( target );

  PhenotypicSegment ** segments = target.segments();

  // TODO <david.parsons@inria.fr> We should take into account that we compute the areas in order (from the leftmost segment, rightwards)
  //   => We shouldn't parse the whole list of points on the left of the segment we are considering (we have
  //      already been through them!)

  if ( dist_to_target_per_segment_ == NULL )
  {
    dist_to_target_per_segment_ = new double [target.nb_segments()]; // Can not be allocated in constructor because number of segments is then unknow
  }
  for (int8_t i = 0 ; i < target.nb_segments() ; i++)
  {
    dist_to_target_per_segment_[i] = delta->get_geometric_area( segments[i]->start, segments[i]->stop );
    dist_to_target_by_feature_[segments[i]->feature] += dist_to_target_per_segment_[i];
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
  if (fitness_computed_) return; // Fitness has already been computed, nothing to do.
  fitness_computed_ = true;

#ifdef NORMALIZED_FITNESS

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    if (target.get_area_by_feature(i)==0.)
    {
      fitness_by_feature_[i] = 0.;
    }
    else
    {
      fitness_by_feature_[i] =  ( target.get_area_by_feature(i) - dist_to_target_by_feature_[i] ) / target.get_area_by_feature(i);
      if ( (fitness_by_feature_[i] < 0.) && (i != METABOLISM) ) // non-metabolic fitness can NOT be lower than zero (we do not want individual to secrete a negative quantity of public good)
      {
        fitness_by_feature_[i] = 0.;
      }
    }
  }

  if ((! indiv_->get_placed_in_population()) || (! exp_m_->get_with_secretion() ))
  {
    fitness_ = fitness_by_feature_[METABOLISM];
  }
  else
  {
    fitness_ =  fitness_by_feature_[METABOLISM] * ( 1 + exp_m_->get_secretion_contrib_to_fitness() * ( indiv_->get_grid_cell()->get_compound_amount() - exp_m_->get_secretion_cost() * fitness_by_feature_[SECRETION] ) );
  }

  if ( exp_m_->get_selection_scheme() == FITNESS_PROPORTIONATE ) // Then the exponential selection is integrated inside the fitness value
  {
    fitness_ = exp( -exp_m_->get_selection_pressure() * (1 - fitness_) );
  }

#else

    for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
    {
      if ( i == SECRETION )
      {
        fitness_by_feature_[SECRETION] =  exp( -exp_m_->get_selection_pressure() * dist_to_target_by_feature_[SECRETION] )
                                          - exp( -exp_m_->get_selection_pressure() * target.area_by_feature(SECRETION) );

        if ( fitness_by_feature_[i] < 0 )
        {
          fitness_by_feature_[i] = 0;
        }
      }
      else
      {
        fitness_by_feature_[i] = exp( -exp_m_->get_selection_pressure() * dist_to_target_by_feature_[i] );
      }
    }

    // Calculate combined, total fitness here!
    // Multiply the contribution of metabolism and the amount of compound in the environment
    if ((! indiv_->get_placed_in_population()) || (! exp_m_->get_with_secretion() ))
    {
      fitness_ = fitness_by_feature_[METABOLISM] ;
    }
    else
    {
      fitness_ =  fitness_by_feature_[METABOLISM] *
                  ( 1 + exp_m_->get_secretion_contrib_to_fitness() * indiv_->get_grid_cell()->compound_amount()
                    - exp_m_->get_secretion_cost() * fitness_by_feature_[SECRETION] );
    }

#endif

  }


  void GeneticUnit::reset_expression( void )
  {
    // useful if the DNA sequence has changed (cf post-treatment programs
    // which replay mutations)

    transcribed_ = false;
    translated_ = false;
    phenotypic_contributions_computed_ = false;
    non_coding_computed_ = false;
    distance_to_target_computed_ = false;
    fitness_computed_ = false;


    // I do not erase the RNA lists, because they were updated
    // during the mutations (cf Dna::undergo_this_mutation)
    // TODO : Reinitialize _transcribed proteins ?

    if ( activ_contribution_ != NULL )
    {
      delete activ_contribution_;
      activ_contribution_ = new Fuzzy();
    }

    if ( inhib_contribution_ != NULL )
    {
      delete inhib_contribution_;
      inhib_contribution_ = new Fuzzy();
    }

    if ( phenotypic_contribution_ != NULL )
    {
      delete phenotypic_contribution_; // Not re-created now, will be conditionally allocated in compute_phenotypic_contribution
      phenotypic_contribution_ = NULL;
    }

    init_statistical_data();
  }


  void GeneticUnit::print_coding_rnas() {
    for (int strand_id = LEADING; strand_id <= LAGGING; ++strand_id) {
      auto& strand = rna_list_[strand_id];
      printf("  %s \n", StrandName[strand_id]);
      for (const auto& rna: strand)
        if (rna.is_coding())
          printf("Promoter at %" PRId32 ", last transcribed position at %" PRId32 "\n", rna.get_promoter_pos(), rna.get_last_transcribed_pos());
    }
  }

  void GeneticUnit::print_proteins( void ) const
  {
    printf( "  LEADING ( %" PRId32 " )\n", static_cast<int32_t>(protein_list_[LEADING].size()));
    for (const auto& prot: protein_list_[LEADING])
      printf( "    Gene on LEADING at %" PRId32 " (%" PRId32 ") (%f %f %f) (%f) %s\n",
              prot.get_shine_dal_pos(), prot.get_length(),
              prot.get_mean(), prot.get_width(), prot.get_height(), prot.get_concentration(),
              prot.get_is_functional() ? "functional" : "non functional" );


    printf( "  LAGGING ( %" PRId32 " )\n", static_cast<int32_t>(protein_list_[LAGGING].size()));
    for (const auto& prot: protein_list_[LAGGING])
      printf( "    Gene on LAGGING at %" PRId32 " (%" PRId32 ") (%f %f %f) (%f) %s\n",
              prot.get_shine_dal_pos(), prot.get_length(),
              prot.get_mean(), prot.get_width(), prot.get_height(), prot.get_concentration(),
              prot.get_is_functional() ? "functional" : "non functional" );
  }

  bool GeneticUnit::is_promoter( Strand strand, int32_t pos, int8_t& dist ) const
  {
    //~ printf( "=============================================== is_promoter\n" );
    //~ printf( "pos : %" PRId32 "\n", pos );

    const char* genome  = dna_->data();
    int32_t  len        = dna_->length();

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
    const char* genome  = dna_->data();
    int32_t  len        = dna_->length();

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
    const char* genome  = dna_->data();
    int32_t  len        = dna_->length();

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
    const char* genome  = dna_->data();
    int32_t  len        = dna_->length();
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
    if (non_coding_computed_) return;
    non_coding_computed_ = true;

    // Create a table of <genome_length> bools initialized to false (non-coding)
    int32_t genome_length = dna_->length();

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
      for (const auto& prot: protein_list_[strand]) {
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
            prom_last   = Utils::mod( prom_first + PROM_SIZE - 1, dna_->length() );
            term_last   = rna->get_last_transcribed_pos();
            term_first  = Utils::mod( term_last - TERM_SIZE + 1, dna_->length() );
            rna_first   = prom_first;
            rna_last    = term_last;
          }
          else
          {
            prom_last   = rna->get_promoter_pos();
            prom_first  = Utils::mod( prom_last - PROM_SIZE + 1, dna_->length() );
            term_first  = rna->get_last_transcribed_pos();
            term_last   = Utils::mod( term_first + TERM_SIZE - 1, dna_->length() );
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
      for (auto& rna: rna_list_[strand]) {
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

        assert(first < indiv_->get_amount_of_dna());
        assert(last < indiv_->get_amount_of_dna());

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
    nb_bases_in_0_CDS_ = 0;
    nb_bases_in_0_functional_CDS_ = 0;
    nb_bases_in_0_non_functional_CDS_ = 0;
    nb_bases_in_0_RNA_ = 0;
    nb_bases_in_0_coding_RNA_ = 0;
    nb_bases_in_0_non_coding_RNA_ = 0;
    nb_bases_non_essential_ = 0;
    nb_bases_non_essential_including_nf_genes_ = 0;
    nb_bases_in_neutral_regions_ = 0;
    nb_neutral_regions_ = 0;

    // We do not know how many neutral regions there will be, but
    // there should be less than nb_coding_RNAs_ + 1
    // As we will see, there may be a shift in values so we take size nb_coding_RNAs_ + 2
    int32_t* tmp_beginning_neutral_regions = new int32_t [ nb_coding_RNAs_ + 2 ];
    int32_t* tmp_end_neutral_regions = new int32_t [ nb_coding_RNAs_ + 2 ];
    memset( tmp_beginning_neutral_regions, -1, nb_coding_RNAs_ + 2 );
    memset( tmp_end_neutral_regions, -1, nb_coding_RNAs_ + 2 );

    for ( int32_t i = 0 ; i < genome_length ; i++ )
    {
      if ( belongs_to_CDS[i] == false )
      {
        nb_bases_in_0_CDS_++;
      }
      if ( belongs_to_functional_CDS[i] == false )
      {
        nb_bases_in_0_functional_CDS_++;
      }
      if ( belongs_to_non_functional_CDS[i] == false )
      {
        nb_bases_in_0_non_functional_CDS_++;
      }
      if ( belongs_to_RNA[i] == false )
      {
        nb_bases_in_0_RNA_++;
      }
      if ( belongs_to_coding_RNA[i] == false )
      {
        nb_bases_in_0_coding_RNA_++;
      }
      if ( belongs_to_non_coding_RNA[i] == false )
      {
        nb_bases_in_0_non_coding_RNA_++;
      }
      if ( is_essential_DNA[i] == false )
      {
        nb_bases_non_essential_++;
      }
      if ( is_essential_DNA_including_nf_genes[i] == false )
      {
        nb_bases_non_essential_including_nf_genes_++;
      }
      if ( is_not_neutral[i] == false )
      {
        nb_bases_in_neutral_regions_++;
      }
      if ( i != 0 )
      {
        if ( is_not_neutral[i] != is_not_neutral[i-1] )
        {
          if ( is_not_neutral[i-1] == true ) // beginning of a neutral region
          {
            tmp_beginning_neutral_regions [nb_neutral_regions_] = i;
          }
          else // end of a neutral region
          {
            tmp_end_neutral_regions [nb_neutral_regions_] = i-1;
            nb_neutral_regions_++;
          }
        }
      }
      else // i = 0
      {
        // we arbitrarily set 0 as the beginning of a neutral region (linkage with end of genetic unit
        // will be done later)
        if ( is_not_neutral[0] == false ) { tmp_beginning_neutral_regions [nb_neutral_regions_] = 0; }
      }
    }

    // we have to treat specifically the last base of the genetic unit in order to link neutral regions
    // at the end and the beginning of genetic unit
    int32_t shift = 0;
    if ( is_not_neutral[genome_length-1] == false )
    {
      if ( is_not_neutral[0] == true ) // end of a neutral region
      {
        tmp_end_neutral_regions[nb_neutral_regions_] = genome_length-1;
        nb_neutral_regions_++;
      }
      else // neutral region goes on after base 0, linkage to be done
      {
        if ( nb_neutral_regions_ != 0 )
        {
          tmp_end_neutral_regions[nb_neutral_regions_] = tmp_end_neutral_regions[ 0 ];
          // the first neutral region is only a subpart of the last one, it should not be
          // taken into account. When we transfer values to the final array, we introduce a shift
          shift = 1;
          // we do not ++ nb_neutral_regions_ as it was already counted
        }
        else // no neutral region detected until now -> all the genetic unit is neutral
        {
          // as all the chromosome is neutral, we indicate 0 as the beginning of the region
          // and genome_length - 1 as its end
          tmp_end_neutral_regions[ 0 ] = genome_length - 1;
          nb_neutral_regions_++;
        }
      }
    }

    // now that we know how many neutral regions there are, we can transfer data to correctly sized arrays
    assert( nb_neutral_regions_ <= nb_coding_RNAs_ + 1 );
    if ( beginning_neutral_regions_ != NULL ) { delete [] beginning_neutral_regions_; }
    if ( end_neutral_regions_ != NULL )       { delete [] end_neutral_regions_; }

    if ( nb_neutral_regions_ > 0 ) // as unlikely as it seems, there may be no neutral region
    {
      beginning_neutral_regions_ = new int32_t [nb_neutral_regions_];
      end_neutral_regions_ = new int32_t [nb_neutral_regions_];
      // transfer from tmp to attributes
      for (int32_t i = 0; i < nb_neutral_regions_; i++)
      {
        beginning_neutral_regions_[i] = tmp_beginning_neutral_regions[i+shift];
        end_neutral_regions_[i]       = tmp_end_neutral_regions[i+shift];
      }
    }
    else // nb_neutral_regions_ == 0
    {
      beginning_neutral_regions_ = NULL;
      end_neutral_regions_ = NULL;
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
        duplicated_promoters[strand].back().shift_position(-pos_1,
                                                           dna_->length());
      }
    }
  }

  void GeneticUnit::get_promoters_included_in(int32_t pos_1,
                                              int32_t pos_2,
                                              Promoters2Strands& promoters)
  {
    assert(pos_1 >= 0 && pos_1 <= dna_->length() && pos_2 >= 0 &&
           pos_2 <= dna_->length() );

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
      int32_t seg_length = dna_->length() + pos_2 - pos_1;

      if ( seg_length >= PROM_SIZE )
      {
        bool is_near_end_of_genome = ( pos_1 + PROM_SIZE > dna_->length() );
        bool is_near_beginning_of_genome = ( pos_2 - PROM_SIZE < 0 );

        if ( !is_near_end_of_genome && !is_near_beginning_of_genome )
        {
          get_promoters(LEADING, AFTER, pos_1, -1, promoters[LEADING] );
          get_promoters(LEADING, BEFORE, -1, pos_2 - PROM_SIZE + 1, promoters[LEADING] );
          get_promoters(LAGGING, AFTER, pos_2, -1, promoters[LAGGING] );
          get_promoters(LAGGING, BEFORE, -1, pos_1 + PROM_SIZE - 1, promoters[LAGGING] );
        }
        else if ( !is_near_end_of_genome ) // => && is_near_beginning_of_genome
        {
          // get_promoters(leading, between, pos_1, pos_2 + dna_->length() - PROM_SIZE + 1,
          //                                         promoters[LEADING] );
          get_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1 +
                            dna_->length(),
                        promoters[LEADING]);
          get_promoters(LAGGING, AFTER, pos_2, -1, promoters[LAGGING] );
          get_promoters(LAGGING, BEFORE, -1, pos_1 + PROM_SIZE - 1, promoters[LAGGING] );
        }
        else if ( !is_near_beginning_of_genome ) // => && is_near_end_of_genome
        {
          get_promoters(LEADING, AFTER, pos_1, -1, promoters[LEADING] );
          get_promoters(LEADING, BEFORE, -1, pos_2 - PROM_SIZE + 1, promoters[LEADING] );
          get_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1 -
                            dna_->length(),
                        promoters[LAGGING] );
        }
        else // is_near_end_of_genome && is_near_beginning_of_genome
        {
          // get_promoters(leading, between, pos_1, pos_2 + dna_->length() - PROM_SIZE + 1,
          //                                         promoters[LEADING] );
          get_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1 +
                            dna_->length(),
                        promoters[LEADING]);
          get_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1 -
                            dna_->length(),
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
            pos1 <= dna_->length() and pos2 <= dna_->length()) or
        (before_after_btw == BEFORE and pos2 >= 0 and pos2 <= dna_->length()) or
        (before_after_btw == AFTER and pos1 >= 0 and pos1 <= dna_->length()));

    auto strand = rna_list_[strand_id];
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
    assert(pos2 <= dna_->length());

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
    assert(pos_2 <= dna_->length());

    // Find the first promoters in the interval
    auto& strand = rna_list_[LEADING];
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
    assert(pos_2 <= dna_->length());

    // Find the first promoters in the interval (if any)
    auto& strand = rna_list_[LAGGING];
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

/// Insert promoters in double stranded list `promoters_to_insert` into `this->rna_list_`.
///
/// The promoters in <promoters_to_insert> must already be at their
/// rightful position according to <this>
/// and the positions of the promoters from <promoters_to_insert> and
/// <this->rna_list_> must not be interlaced
/// i.e. no promoter in <this->rna_list_> must have a position in
/// [first_prom_to_insert->pos ; last_prom_to_insert->pos]
  void GeneticUnit::insert_promoters(Promoters2Strands& promoters_to_insert)
  {
    // TODO vld: to be merged with insert_promoters_at(...)
    for (auto strand: {LEADING, LAGGING}) {
      if (promoters_to_insert[strand].size() <= 0)
        continue;
      // Get to the right position in individual's list (first promoter after the inserted segment)
      int32_t from_pos = promoters_to_insert[strand].back().get_promoter_pos();

      auto pos = find_if(rna_list_[strand].begin(),
                         rna_list_[strand].end(),
                         [from_pos, strand] (Rna & r) {
                             if (strand == LEADING)
                               return r.get_promoter_pos() >= from_pos;
                             else
                               return r.get_promoter_pos() < from_pos;
                         });

      // Insert the promoters in the individual's RNA list
      for (auto& to_insert: promoters_to_insert[strand])
        // TODO vld: could be compacted in a unique emplace(pos, to_insert) ?
        if (pos != rna_list_[strand].end())
          rna_list_[strand].insert(pos, to_insert);
        else
          rna_list_[strand].push_back(to_insert);
    }
  }

/// Insert promoters in double stranded list `promoters_to_insert`
/// into `this->rna_list_` at position `pos`
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
      auto first = find_if(rna_list_[strand].begin(),
                           rna_list_[strand].end(),
                           [pos, strand](Rna & r) {
                               if (strand == LEADING)
                                 return r.get_promoter_pos() >= pos;
                               else
                                 return r.get_promoter_pos() < pos; });

      // Insert the promoters in the individual's RNA list
      for (auto& to_insert: promoters_to_insert[strand]) {
        // Update promoter position
        to_insert.shift_position( pos, dna_->length() );
        // Insert
        if (first != rna_list_[strand].end())
          rna_list_[strand].insert(first, to_insert);
        else
          rna_list_[strand].push_back(to_insert);
      }
    }
  }


/// Remove the RNAs of the LEADING strand whose starting positions lie
/// in [pos_1 ; pos_2[
  void GeneticUnit::remove_leading_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
    assert(pos_1 >= 0);
    assert(pos_1 < dna_->length());
    assert(pos_2 >= 0);
    assert(pos_2 <= dna_->length());

    if ( pos_1 > pos_2 ) {
      remove_leading_promoters_starting_after(pos_1);
      remove_leading_promoters_starting_before(pos_2);
    } else {
      auto& strand = rna_list_[LEADING];
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
    assert(pos_1 >= 0 && pos_1 <= dna_->length() && pos_2 >= 0 &&
           pos_2 <= dna_->length() );

    if ( pos_1 == dna_->length() ) pos_1 = 0;
    if ( pos_2 == 0 ) pos_2 = dna_->length();
    if ( pos_1 > pos_2 ) { // vld: that's a weird case... really do this? used from remove_promoters_around()
      remove_lagging_promoters_starting_after( pos_1 );
      remove_lagging_promoters_starting_before( pos_2 );
    } else {
      auto& strand = rna_list_[LAGGING];
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
    assert( pos >= 0 && pos < dna_->length() );

    auto& strand = rna_list_[LEADING];
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
    assert(pos >= 0 and pos < dna_->length());

    auto& strand = rna_list_[LAGGING];
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
    assert( pos >= 0 && pos < dna_->length() );

    auto& strand = rna_list_[LEADING];
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
    assert( pos < dna_->length() && pos >= 0 );

    auto& strand = rna_list_[LAGGING];
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
    assert(pos_1 >= 0 && pos_1 < dna_->length() && pos_2 >= 0 &&
           pos_2 < dna_->length() );

    // When pos_1 > pos_2, we will perform the search in 2 steps.
    // As positions  0 and _dna->length() are equivalent, it's preferable to
    // keep 0 for pos_1 and _dna->length() for pos_2.

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
        //~ memcpy( tmp, &dna_->data()[i], PROM_SIZE * sizeof(char) );
        //~ printf( "new promoter found on the LEADING strand at position %"PRId32" : %s\n", i, tmp );

        // Look for the right place to insert the new promoter in the list
        auto& strand = rna_list_[LEADING];
        auto first = find_if(strand.begin(),
                             strand.end(),
                             [i](Rna & r)
                             { return r.get_promoter_pos() >= i; });
        if (first == strand.end() or first->get_promoter_pos() != i) {
#ifndef __REGUL
          rna_list_[LEADING].emplace(first, this, LEADING, i, dist);
#else
          rna_list_[LEADING].emplace(first, this, LEADING, i, dist);
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
    assert(pos_1 >= 0 && pos_1 < dna_->length() && pos_2 >= 0 &&
           pos_2 < dna_->length());

    // When pos_1 > pos_2, we will perform the search in 2 steps.
    // As positions  0 and dna_->length() are equivalent, it's preferable to
    // keep 0 for pos_1 and dna_->length() for pos_2.

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
        assert ( i >= 0 && i < dna_->length() );

        // Look for the right place to insert the new promoter in the list
        auto& strand = rna_list_[LAGGING];
        auto first = find_if(strand.begin(),
                             strand.end(),
                             [i](Rna & r)
                             { return r.get_promoter_pos() <= i; });
        if (first == strand.end() or first->get_promoter_pos() != i) {
#ifndef __REGUL
          rna_list_[LAGGING].emplace(first, this, LAGGING, i, dist);
#else
          rna_list_[LAGGING].emplace(first, this, LAGGING, i, dist);
#endif // __REGUL
        }
      }
    }
  }


/// Look for new promoters on the LEADING strand whose starting positions would be >= pos
  void GeneticUnit::look_for_new_leading_promoters_starting_after( int32_t pos ) {
    assert( pos >= 0 && pos < dna_->length() );

    // Hamming distance of the sequence from the promoter consensus
    int8_t dist;

    // rna list node used to find the new promoter's place in the list

    for ( int32_t i = pos ; i < dna_->length() ; i++ )
    {
      if (is_promoter(LEADING, i, dist)) { // dist takes the hamming distance of the sequence from the consensus
        // Look for the right place to insert the new promoter in the list
        auto& strand = rna_list_[LEADING];
        auto first = find_if(strand.begin(),
                             strand.end(),
                             [i](Rna & r) {
                                 return r.get_promoter_pos() >= i;
                             });
        if (first == strand.end() or first->get_promoter_pos() != i) {
#ifndef __REGUL
          rna_list_[LEADING].emplace(first, this, LEADING, i, dist);
#else
          rna_list_[LEADING].emplace(first, this, LEADING, i, dist);
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
    assert( pos >= 0 && pos < dna_->length() );

    // Hamming distance of the sequence from the promoter consensus
    int8_t dist;
    auto& strand = rna_list_[LAGGING];
    auto first = strand.begin();

    for ( int32_t i = dna_->length() - 1 ; i >= pos ; i-- ) {
      if ( is_promoter( LAGGING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
      {
        assert ( i >= 0 && i < dna_->length() );
        // Look for the right place to insert the new promoter in the list
        first = find_if(first,
                        strand.end(),
                        [i](Rna & r)
                        { return r.get_promoter_pos() <= i; });
        if (first == strand.end() or first->get_promoter_pos() != i) {
#ifndef __REGUL
          rna_list_[LAGGING].emplace(first, this, LAGGING, i, dist);
#else
          rna_list_[LAGGING].emplace(first, this, LAGGING, i, dist);
#endif
        }
      }
    }
  }


/// Look for new promoters on the LEADING strand whose starting
/// positions would be < pos
  void GeneticUnit::look_for_new_leading_promoters_starting_before( int32_t pos )
  {
    assert( pos >= 0 && pos < dna_->length() );

    // Hamming distance of the sequence from the promoter consensus
    int8_t dist;

    auto& strand = rna_list_[LEADING];
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
          rna_list_[LEADING].emplace(first, this, LEADING, i, dist);
#else
          rna_list_[LEADING].emplace(first, this, LEADING, i, dist);
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
    assert( pos >= 0 && pos < dna_->length() );

    // Hamming distance of the sequence from the promoter consensus
    int8_t dist;

    // rna list node used to find the new promoter's place in the list
    auto& strand = rna_list_[LAGGING];
    auto first = strand.begin();

    for ( int32_t i = pos - 1 ; i >= 0 ; i-- )
    {
      if ( is_promoter( LAGGING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
      {
        assert ( i >= 0 && i < dna_->length() );
        // Look for the right place to insert the new promoter in the list
        first = find_if(first,
                        strand.end(),
                        [i](Rna & r) {
                            return r.get_promoter_pos() <= i;
                        });
        if (first == strand.end() or first->get_promoter_pos() != i) {
#ifndef __REGUL
          rna_list_[LAGGING].emplace(first, this, LAGGING, i, dist);
#else
          rna_list_[LAGGING].emplace(first, this, LAGGING, i, dist);
#endif
        }
      }
    }
  }


/// Shift (by delta_post) the positions of the promoters from the
/// LEADING strand whose starting positions are >= pos.
  void GeneticUnit::move_all_leading_promoters_after( int32_t pos, int32_t delta_pos ) {
    auto& strand = rna_list_[LEADING];
    for (auto rna = find_if(strand.begin(), strand.end(), [pos](Rna & r) { return r.get_promoter_pos() >= pos; });
         rna != strand.end();
         ++rna)
      rna->shift_position( delta_pos, dna_->length());
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
    auto& strand = rna_list_[LAGGING];
    // Update RNAs until we pass pos (or we reach the end of the list)
    for (auto rna = strand.begin();
         rna != strand.end() and rna->get_promoter_pos() >= pos;
         ++rna)
      rna->shift_position(delta_pos, dna_->length());
  }

/// Copy (into new_promoter_list) the promoters from the LEADING
/// strand whose starting positions lie in [pos_1 ; pos_2[
  void GeneticUnit::copy_leading_promoters_starting_between(int32_t pos_1,
                                                            int32_t pos_2,
                                                            Promoters1Strand& new_promoter_list ) {
    // 1) Go to first RNA to copy
    auto& strand = rna_list_[LEADING];
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
  auto& strand = rna_list_[LAGGING];
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
  dna_->save(backup_file);
  gzwrite(backup_file, &min_gu_length_, sizeof(min_gu_length_));
  gzwrite(backup_file, &max_gu_length_, sizeof(max_gu_length_));
}

int32_t GeneticUnit::get_nb_terminators() {
  int32_t nb_term = 0;
  if (dna_->length() >= TERM_SIZE)
    for (int32_t i = 0 ; i < dna_->length() ; i++)
      if (is_terminator(LEADING, i))  // No need to count on both the
        // LEADING and the LAGGING
        // strand as terminators are
        // "shared"
        nb_term++;
  return nb_term;
}

#ifdef DEBUG
/// Compare current rna_list_ with locate_promoters-generated rna_list_
void GeneticUnit::assert_promoters() {
  // Check that the lists are ordered correctly
  assert_promoters_order();

  // Make a backup of the genetic unit's lists of RNAs
  auto backup = rna_list_;
  rna_list_ = {{},{}};
  locate_promoters(); // regenerate rna_list_

  // Compare lists
  for (auto strand: {LEADING, LAGGING}) {
    if(backup[strand].size() != rna_list_[strand].size()) {
      printf("Individual %" PRId32 "\n", indiv_->get_id());
      printf("***************** FOUND *******************");
      print_rnas(backup[strand], strand);
      printf("***************** EXPECTED *******************");
      print_rnas(rna_list_[strand], strand);
      printf("******************************************");
      assert(false);
    }

    auto node_old = backup[strand].begin();
    auto node_new = node_old; // just for the type
    for (node_old = backup[strand].begin(), node_new = rna_list_[strand].begin();
         node_old != backup[strand].end(); // rna_list_ is the same size as backup
         ++node_old, ++node_new) {
      // TODO vld: to factor
      if ( node_old->get_strand() != node_new->get_strand() ) {
        printf("Individual %" PRId32 "\n", indiv_->get_id());
        printf( "****************************** STRAND problem ******************************\n" );
        printf( "should be : \n" );
        print_rnas(rna_list_);
        printf( "is : \n" );
        print_rnas( backup );
        printf( "****************************************************************************\n" );
        printf( "  %" PRId32 " (%s) : %f    vs    %" PRId32 " (%s) : %f\n",
                node_old->get_promoter_pos(), StrandName[strand], node_old->get_basal_level(),
                node_new->get_promoter_pos(), StrandName[strand], node_new->get_basal_level() );
        printf( "  genome length : %" PRId32 "\n", dna_->length() );
        assert( node_old->get_strand() == node_new->get_strand() );
      }

      if ( node_old->get_promoter_pos() != node_new->get_promoter_pos() ) {
        printf("Individual %" PRId32 "\n", indiv_->get_id());
        printf( "***************************** POSITION problem *****************************\n" );
        printf( "should be : \n" );
        print_rnas(rna_list_);
        printf( "is : \n" );
        print_rnas( backup );
        printf( "****************************************************************************\n" );
        printf( "  %" PRId32 " (%s) : %f    vs    %" PRId32 " (%s) : %f\n",
                node_old->get_promoter_pos(), StrandName[strand], node_old->get_basal_level(),
                node_new->get_promoter_pos(), StrandName[strand], node_new->get_basal_level() );
        printf( "  genome length : %" PRId32 "\n", dna_->length() );
        assert( node_old->get_promoter_pos() == node_new->get_promoter_pos()  );
      }

      if ( node_old->get_basal_level() != node_new->get_basal_level() ) {
        printf("Individual %" PRId32 "\n", indiv_->get_id());
        printf( "*************************** BASAL LEVEL problem ****************************\n" );
        printf( "should be : \n" );
        print_rnas(rna_list_);
        printf( "is : \n" );
        print_rnas( backup );
        printf( "****************************************************************************\n" );
        printf( "  %" PRId32 " (%s) : %f    vs    %" PRId32 " (%s) : %f\n",
                node_old->get_promoter_pos(), StrandName[strand], node_old->get_basal_level(),
                node_new->get_promoter_pos(), StrandName[strand], node_new->get_basal_level() );
        printf( "  genome length : %" PRId32 "\n", dna_->length() );
        assert( node_old->get_basal_level() == node_new->get_basal_level() );
      }
    }
  }

  rna_list_[LEADING].clear();
  rna_list_[LAGGING].clear();
  rna_list_ = backup;
}

void GeneticUnit::assert_promoters_order() {
  for (auto strand: {LEADING, LAGGING}) {
    if (rna_list_[strand].size() < 2)
      continue;

    for (auto it = rna_list_[strand].begin(); next(it) != rna_list_[strand].end(); ++it) {
      if ((strand == LEADING and it->get_promoter_pos() >= next(it)->get_promoter_pos()) or
          (strand == LAGGING and it->get_promoter_pos() <= next(it)->get_promoter_pos())) {
        printf("Individual %" PRId32 "\n", indiv_->get_id());
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
  int32_t genome_length = dna_->length();
  bool* belongs_to_coding_RNA = new bool[genome_length];
  memset(belongs_to_coding_RNA,0, genome_length);

  // Parse RNA lists and mark the corresponding bases as coding (only for the coding RNAs)
  for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ ) {
    for (auto& rna: rna_list_[strand]) {
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
//  Environment* env = exp_m_->get_env() ;
//
//  reset_expression();
//  locate_promoters();
//  distance_to_target_computed_        = false;
//  fitness_computed_                   = false;
//  compute_phenotypic_contribution();
//  if(exp_m_->get_output_m()->get_compute_phen_contrib_by_GU())
//  {
//    compute_distance_to_target(env);
//  }
//  compute_fitness(env);
//  double initial_fitness = get_fitness();
//
//  int32_t genome_length = dna_->length();
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
//      dna_->remove(start,end);
//
//      locate_promoters();
//      reset_expression();
//      distance_to_target_computed_        = false;
//      fitness_computed_                   = false;
//      compute_phenotypic_contribution();
//      if(exp_m_->get_output_m()->get_compute_phen_contrib_by_GU())
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
//  distance_to_target_computed_        = false;
//  fitness_computed_                   = false;
//  compute_phenotypic_contribution();
//  if(exp_m_->get_output_m()->get_compute_phen_contrib_by_GU())
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
//  Environment* env = exp_m_->get_env() ;
//
//  reset_expression();
//  locate_promoters();
//  distance_to_target_computed_        = false;
//  fitness_computed_                   = false;
//  compute_phenotypic_contribution();
//  if(exp_m_->get_output_m()->get_compute_phen_contrib_by_GU())
//  {
//    compute_distance_to_target(env);
//  }
//  compute_fitness(env);
//  double initial_fitness = get_fitness();
//
//  int32_t genome_length = dna_->length();
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
//          random_portion[j] = '0' + indiv_->get_mut_prng()->random( NB_BASE );
//        }
//        random_portion[length] = 0;
//        pos = indiv_->get_mut_prng()->random( length )+start;
//        dna_->insert(pos, random_portion);
//
//        _non_coding_computed = false;
//        non_coding_bases_nb_before_fitness = get_nb_bases_in_0_coding_RNA();
//
//        locate_promoters();
//        reset_expression();
//        distance_to_target_computed_        = false;
//        fitness_computed_                   = false;
//        compute_phenotypic_contribution();
//        if(exp_m_->get_output_m()->get_compute_phen_contrib_by_GU())
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
//          dna_->remove(pos, pos + length);
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
//  distance_to_target_computed_        = false;
//  fitness_computed_                   = false;
//  compute_phenotypic_contribution();
//  if(exp_m_->get_output_m()->get_compute_phen_contrib_by_GU())
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
        GeneticUnit * tmpunit = new GeneticUnit(indiv_, seq, mutlength);
        tmpunit->do_transcription();
        tmpunit->do_translation();
        nb_genes_in_segment = tmpunit->get_nb_coding_RNAs();

        // Check whether the pos3donor breakpoint killed one or several of those genes, in that case decrement nb_genes_in_segment
        for (auto& strand: tmpunit->rna_list_)
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
        GeneticUnit * tmpunit = new GeneticUnit(indiv_, seq, mutlength);
        tmpunit->do_transcription();
        tmpunit->do_translation();
        nb_genes_in_segment = tmpunit->get_nb_coding_RNAs();

        // Remove the genes that overlap ori in this temporary unit, as the transferred segment was actually linear
        for (auto& strand: tmpunit->rna_list_)
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
      for (auto& rna: rna_list_[strand]) {
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

    nb_coding_RNAs_ = 0;
    nb_non_coding_RNAs_ = 0;
    overall_size_coding_RNAs_ = 0;
    overall_size_non_coding_RNAs_ = 0;
    nb_genes_activ_ = 0;
    nb_genes_inhib_ = 0;
    nb_fun_genes_ = 0;
    nb_non_fun_genes_ = 0;
    overall_size_fun_genes_ = 0;
    overall_size_non_fun_genes_ = 0;

    nb_bases_in_0_CDS_ = -1;
    nb_bases_in_0_functional_CDS_ = -1;
    nb_bases_in_0_non_functional_CDS_ = -1;
    nb_bases_in_0_RNA_ = -1;
    nb_bases_in_0_coding_RNA_ = -1;
    nb_bases_in_0_non_coding_RNA_ = -1;
    nb_bases_non_essential_ = -1;
    nb_bases_non_essential_including_nf_genes_ = -1;

    modularity_ = -1;

    beginning_neutral_regions_ = NULL;
    end_neutral_regions_ = NULL;
  }

void GeneticUnit::clear_transcribed_proteins() {
  for (auto strand: {LEADING, LAGGING}) {
    for (auto& rna: rna_list_[strand])
      rna.clear_transcribed_proteins();
    protein_list_[strand].clear();
  }
}
} // namespace aevol
