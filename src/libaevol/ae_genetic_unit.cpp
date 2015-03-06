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
//                              Libraries
// =================================================================
#include <assert.h>


// =================================================================
//                            Project Files
// =================================================================
#include "ae_genetic_unit.h"


#include "ae_exp_manager.h"
#include "ae_exp_setup.h"
#include "ae_codon.h"
#include "ae_mutation.h"
#include "ae_enums.h"

#ifdef __REGUL
  #include "ae_individual_R.h"
#else
  #include "ae_individual.h"
#endif

#include "fuzzy.h"

namespace aevol {


// =================================================================
//                       Miscellaneous Functions
// =================================================================
int compare_prot_pos( const void* pos, const void* prot ) // This function has to be a plain int
                                                            // to comply with the definition of bsearch()
{
  if ( ((ae_protein*)prot)->get_shine_dal_pos() == *(int32_t*)pos ) return 0;
  else return 1;
}

//##############################################################################
//                                                                             #
//                            Class ae_genetic_unit                            #
//                                                                             #
//##############################################################################
// =====================================================================
//                          Accessors' definitions
// =====================================================================
ae_exp_manager* ae_genetic_unit::get_exp_m( void ) const
{
  return _exp_m;
}

ae_individual* ae_genetic_unit::get_indiv( void ) const
{
  return _indiv;
}

ae_dna* ae_genetic_unit::get_dna( void ) const
{
  return _dna;
}

/*const*/ vector<list<ae_rna*>> ae_genetic_unit::get_rna_list() const
{
  vector<list<ae_rna*>> r(LAGGING - LEADING + 1);
  for (int8_t strand = LEADING ; strand <= LAGGING ; strand++)
    for (ae_list_node<ae_rna*>* rna_node = _rna_list[strand]->get_first();
         rna_node != NULL;
         rna_node = rna_node->get_next())
      r[strand].push_back(rna_node->get_obj());
  return r;//std::move(r);
}

const list<ae_protein*> ae_genetic_unit::get_protein_list(ae_strand strand) const
{
  list<ae_protein*> r;
    for (ae_list_node<ae_protein*>* protein_node = _protein_list[strand]->get_first() ;
         protein_node != NULL ;
         protein_node = protein_node->get_next())
      r.push_back(protein_node->get_obj());
  return r;
}

Fuzzy* ae_genetic_unit::get_activ_contribution( void ) const
{
  return _activ_contribution;
}

Fuzzy* ae_genetic_unit::get_inhib_contribution( void ) const
{
  return _inhib_contribution;
}

Fuzzy* ae_genetic_unit::get_phenotypic_contribution( void ) const
{
  assert(_phenotypic_contribution != NULL);
  return _phenotypic_contribution;
}

/*!
  Returns the DNA sequence
*/
const char* ae_genetic_unit::get_sequence( void ) const
{
  return _dna->get_data();
}

/*!
  Returns the DNA sequence length
*/
int32_t ae_genetic_unit::get_seq_length( void ) const
{
  return _dna->get_length();
}

int32_t ae_genetic_unit::get_nb_coding_RNAs( void ) const
{
  return _nb_coding_RNAs;
}

int32_t ae_genetic_unit::get_nb_non_coding_RNAs( void ) const
{
  return _nb_non_coding_RNAs;
}

double ae_genetic_unit::get_overall_size_coding_RNAs( void ) const
{
  return _overall_size_coding_RNAs;
}

double ae_genetic_unit::get_av_size_coding_RNAs( void ) const
{
  if ( _nb_coding_RNAs != 0 )
  {
    return _overall_size_coding_RNAs / _nb_coding_RNAs;
  }
  else return 0.0;
}

double ae_genetic_unit::get_overall_size_non_coding_RNAs( void ) const
{
  return _overall_size_non_coding_RNAs;
}

double ae_genetic_unit::get_av_size_non_coding_RNAs( void ) const
{
  if ( _nb_non_coding_RNAs != 0 )
  {
    return _overall_size_non_coding_RNAs / _nb_non_coding_RNAs;
  }
  else return 0.0;
}

int32_t ae_genetic_unit::get_nb_genes_activ( void ) const
{
  return _nb_genes_activ;
}

int32_t ae_genetic_unit::get_nb_genes_inhib( void ) const
{
  return _nb_genes_inhib;
}

int32_t ae_genetic_unit::get_nb_functional_genes( void ) const
{
  return _nb_fun_genes;
}

int32_t ae_genetic_unit::get_nb_non_functional_genes( void ) const
{
  return _nb_non_fun_genes;
}

double ae_genetic_unit::get_overall_size_functional_genes( void ) const
{
  return _overall_size_fun_genes;
}

double ae_genetic_unit::get_av_size_functional_genes( void ) const
{
  if ( _nb_fun_genes != 0 )
  {
    return _overall_size_fun_genes / _nb_fun_genes;
  }
  else return 0.0;
}

double ae_genetic_unit::get_overall_size_non_functional_genes( void ) const
{
  return _overall_size_non_fun_genes;
}

double ae_genetic_unit::get_av_size_non_functional_genes( void ) const
{
  if ( _nb_non_fun_genes != 0 )
  {
    return _overall_size_non_fun_genes / _nb_non_fun_genes;
  }
  else return 0.0;
}

int32_t ae_genetic_unit::get_nb_bases_in_0_CDS( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_CDS;
}

int32_t ae_genetic_unit::get_nb_bases_in_0_functional_CDS( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_functional_CDS;
}

int32_t ae_genetic_unit::get_nb_bases_in_0_non_functional_CDS( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_non_functional_CDS;
}

int32_t ae_genetic_unit::get_nb_bases_in_0_RNA( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_RNA;
}

int32_t ae_genetic_unit::get_nb_bases_in_0_coding_RNA( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_coding_RNA;
}

int32_t ae_genetic_unit::get_nb_bases_in_0_non_coding_RNA( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_non_coding_RNA;
}

int32_t ae_genetic_unit::get_nb_bases_non_essential( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_non_essential;
}

int32_t ae_genetic_unit::get_nb_bases_non_essential_including_nf_genes( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_non_essential_including_nf_genes;
}

int32_t ae_genetic_unit::get_nb_bases_in_neutral_regions( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_neutral_regions;
}

int32_t ae_genetic_unit::get_nb_neutral_regions( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_neutral_regions;
}

int32_t* ae_genetic_unit::get_beginning_neutral_regions( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _beginning_neutral_regions;
}

int32_t* ae_genetic_unit::get_end_neutral_regions( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _end_neutral_regions;
}

double ae_genetic_unit::get_modularity( void ) const
{
  return _modularity;
}

double ae_genetic_unit::get_dist_to_target_by_feature( ae_env_axis_feature feature ) const
{
   assert( _distance_to_target_computed );

  return _dist_to_target_by_feature[feature];
}

double ae_genetic_unit::get_fitness( void ) const
{
  assert( _fitness_computed );

  return _fitness;
}

double ae_genetic_unit::get_fitness_by_feature( ae_env_axis_feature feature ) const
{
  assert( _fitness_computed );

  return _fitness_by_feature[feature];
}

int32_t ae_genetic_unit::get_min_gu_length( void ) const
{
  return _min_gu_length;
}

int32_t ae_genetic_unit::get_max_gu_length( void ) const
{
  return _max_gu_length;
}

void ae_genetic_unit::set_min_gu_length( int32_t min_gu_length )
{
  _min_gu_length = min_gu_length;
}

void ae_genetic_unit::set_max_gu_length( int32_t max_gu_length )
{
  _max_gu_length = max_gu_length;
}

void ae_genetic_unit::set_exp_m( ae_exp_manager* exp_m )
{
  _exp_m = exp_m;
}

// =====================================================================
//                       functions' definition
// =====================================================================
void ae_genetic_unit::print_rnas( void ) const
{
  print_rnas( _rna_list );
}

/* static */ void ae_genetic_unit::print_rnas( ae_list<ae_rna*>** rnas )
{
  print_rnas( rnas[LEADING], LEADING );
  print_rnas( rnas[LAGGING], LAGGING );
}

bool ae_genetic_unit::is_start( ae_strand strand, int32_t index ) const
{
  return ( get_codon( strand, index ) == CODON_START );
}

bool ae_genetic_unit::is_stop( ae_strand strand, int32_t index ) const
{
  return ( get_codon( strand, index ) == CODON_STOP );
}

void ae_genetic_unit::remove_all_promoters( void )
{
  _rna_list[LEADING]->erase( true );
  _rna_list[LAGGING]->erase( true );
}

void ae_genetic_unit::move_all_promoters_after( int32_t pos, int32_t delta_pos )
{
  move_all_leading_promoters_after( pos, delta_pos );
  move_all_lagging_promoters_after( pos, delta_pos );
}

void ae_genetic_unit::extract_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>** extracted_promoters )
{
  assert( pos_1 >= 0 && pos_1 < pos_2 && pos_2 <= _dna->get_length() );
  if ( pos_2 - pos_1 >= PROM_SIZE )
  {
    extract_leading_promoters_starting_between( pos_1, pos_2 - PROM_SIZE + 1, extracted_promoters[LEADING] );
    extract_lagging_promoters_starting_between( pos_1 + PROM_SIZE - 1, pos_2, extracted_promoters[LAGGING] );
  }
}

void ae_genetic_unit::extract_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>** extracted_promoters )
{
  extract_leading_promoters_starting_between( pos_1, pos_2, extracted_promoters[LEADING] );
  extract_lagging_promoters_starting_between( pos_1, pos_2, extracted_promoters[LAGGING] );
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
void ae_genetic_unit::remove_promoters_around( int32_t pos )
{
  if ( _dna->get_length() >= PROM_SIZE )
  {
    remove_leading_promoters_starting_between( ae_utils::mod(pos - PROM_SIZE + 1, _dna->get_length()), pos );
    remove_lagging_promoters_starting_between( pos, ae_utils::mod(pos + PROM_SIZE - 1, _dna->get_length()) );
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
void ae_genetic_unit::remove_promoters_around( int32_t pos_1, int32_t pos_2 )
{
  if ( ae_utils::mod(pos_1 - pos_2, _dna->get_length()) >= PROM_SIZE )
  {
    remove_leading_promoters_starting_between( ae_utils::mod(pos_1 - PROM_SIZE + 1, _dna->get_length()), pos_2 );
    remove_lagging_promoters_starting_between( pos_1, ae_utils::mod(pos_2 + PROM_SIZE - 1, _dna->get_length()) );
  }
  else
  {
    remove_all_promoters();
  }
}


/*!
  \brief  Look for promoters that are astride pos and add them to the list of promoters (_rna_list).

  Look for promoters that include BOTH the base before AND after pos (marked X in the cartoon below).
  If the genome is smaller than the size of a promoter, no search is performed.

  \verbatim
     -------------------------------------------------------
    |   |   |   |   | X | X |   |   |   |   |   |   |   |   |
     -------------------------------------------------------
    ^                   ^
    0                  pos
  \endverbatim
*/
void ae_genetic_unit::look_for_new_promoters_around( int32_t pos )
{
  assert( pos >= 0 && pos <= _dna->get_length() );

  if ( _dna->get_length() >= PROM_SIZE )
  {
    look_for_new_leading_promoters_starting_between(  ae_utils::mod(pos - PROM_SIZE + 1, _dna->get_length()),
                                                      ae_utils::mod(pos                , _dna->get_length()) );
    look_for_new_lagging_promoters_starting_between(  ae_utils::mod(pos                , _dna->get_length()),
                                                      ae_utils::mod(pos + PROM_SIZE - 1, _dna->get_length()) );
  }
}


/*!
  \brief  Look for promoters that contain at least 1 base lying in [pos_1 ; pos_2[ and add them to the list of promoters (_rna_list).

  Look for promoters that   * include BOTH the base before AND after pos_1 (marked X in the cartoon below).
                            * include BOTH the base before AND after pos_2 (marked Y in the cartoon below).
                            * are completely contained between pos_1 and pos_2.
  If the genome is smaller than the size of a promoter, no search is performed.

  \verbatim
     -------------------------------------------------------
    |   |   |   |   | X | X |   |   |   | Y | Y |   |   |   |
     -------------------------------------------------------
    ^                   ^                   ^
    0                 pos_1               pos_2
  \endverbatim
*/
void ae_genetic_unit::look_for_new_promoters_around( int32_t pos_1, int32_t pos_2 )
{
  //~ if ( ae_utils::mod( pos_1 - pos_2, _dna->get_length()) == PROM_SIZE - 1 )
  //~ {
    //~ // We have to look at every possible position on the genome.
    //~ locate_promoters();
  //~ }
  /*else*/ if ( _dna->get_length() >= PROM_SIZE )
  {
    look_for_new_leading_promoters_starting_between( ae_utils::mod(pos_1 - PROM_SIZE + 1, _dna->get_length()), pos_2 );
    look_for_new_lagging_promoters_starting_between( pos_1, ae_utils::mod(pos_2 + PROM_SIZE - 1, _dna->get_length()) );
  }
}

//~ void ae_genetic_unit::duplicate_promoters_starting_between( int32_t pos_1, int32_t pos_2, int32_t delta_pos )
//~ {
  //~ duplicate_leading_promoters_starting_between( pos_1, pos_2, delta_pos );
  //~ duplicate_lagging_promoters_starting_between( pos_1, pos_2, delta_pos );
//~ }

void ae_genetic_unit::copy_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>** new_promoter_lists )
{
  if ( ae_utils::mod( pos_2 - pos_1 - 1, _dna->get_length() ) + 1 >= PROM_SIZE )
  {
    copy_leading_promoters_starting_between( pos_1, ae_utils::mod( pos_2 - PROM_SIZE + 1, _dna->get_length() ), new_promoter_lists[LEADING] );
    copy_lagging_promoters_starting_between( ae_utils::mod( pos_1 + PROM_SIZE - 1, _dna->get_length() ), pos_2, new_promoter_lists[LAGGING] );
  }
}

void ae_genetic_unit::copy_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>** new_promoter_lists )
{
  copy_leading_promoters_starting_between( pos_1, pos_2, new_promoter_lists[LEADING] );
  copy_lagging_promoters_starting_between( pos_1, pos_2, new_promoter_lists[LAGGING] );
}

//~ void ae_genetic_unit::copy_all_promoters( ae_list<ae_rna*>** new_promoter_lists )
//~ {
  //~ ae_list_node<ae_rna*>* rna_node = NULL;

  //~ for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ )
  //~ {
    //~ rna_node = _rna_list[strand]->get_first();

    //~ while ( rna_node != NULL )
    //~ {
      //~ new_promoter_lists[strand]->add( new ae_rna( this, *(rna_node->get_obj()) ) );

      //~ rna_node = rna_node->get_next();
    //~ }
  //~ }
//~ }

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
ae_genetic_unit::ae_genetic_unit( ae_individual* indiv, int32_t length, ae_jumping_mt * prng )
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

  _dna = new ae_dna( this, length, prng );

  // Create empty rna and protein lists
  _rna_list           = new ae_list<ae_rna*>*[2];
  _rna_list[LEADING]  = new ae_list<ae_rna*>();
  _rna_list[LAGGING]  = new ae_list<ae_rna*>();

  _protein_list           = new ae_list<ae_protein*>*[2];
  _protein_list[LEADING]  = new ae_list<ae_protein*>();
  _protein_list[LAGGING]  = new ae_list<ae_protein*>();

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

/*!
  \brief Create a new genetic unit for indiv with sequence seq [of size length] [and containing promoters prom_list]

  Promoters will be looked for if prom_list is not provided (this may take some time).

  WARNING :
    seq will be used directly which means the caller must not delete it
    The same goes for prom_list if it is provided.
*/
ae_genetic_unit::ae_genetic_unit( ae_individual* indiv, char* seq, int32_t length, ae_list<ae_rna*>** prom_list /*= NULL*/ )
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

  _dna = new ae_dna( this, seq, length );

  if ( prom_list != NULL )
  {
    // Copy rna lists
    _rna_list = prom_list;

    // TODO vld: remove conversion that should not happen here
    // (hack in the process of removing ae_lists from ae_dna)
    ae_dna::set_GU(get_rna_list(), this);
  }
  else
  {
    // Create empty rna lists
    _rna_list           = new ae_list<ae_rna*>*[2];
    _rna_list[LEADING]  = new ae_list<ae_rna*>();
    _rna_list[LAGGING]  = new ae_list<ae_rna*>();

    // Look for promoters
    locate_promoters();
  }

  // Create empty protein lists
  _protein_list           = new ae_list<ae_protein*>*[2];
  _protein_list[LEADING]  = new ae_list<ae_protein*>();
  _protein_list[LAGGING]  = new ae_list<ae_protein*>();

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
ae_genetic_unit::ae_genetic_unit(ae_individual* indiv, const ae_genetic_unit& model)
{
  _exp_m = indiv->get_exp_m();
  _indiv = indiv;

  _transcribed                        = false;
  _translated                         = false;
  _phenotypic_contributions_computed  = false;
  _non_coding_computed                = false;
  _distance_to_target_computed        = false;
  _fitness_computed                   = false;

  _min_gu_length = model._min_gu_length;
  _max_gu_length = model._max_gu_length;

  // Copy DNA
  _dna = new ae_dna( this, *(model._dna) );

  // Create empty rna and protein lists
  _rna_list           = new ae_list<ae_rna*>*[2];
  _rna_list[LEADING]  = new ae_list<ae_rna*>();
  _rna_list[LAGGING]  = new ae_list<ae_rna*>();

  _protein_list           = new ae_list<ae_protein*>*[2];
  _protein_list[LEADING]  = new ae_list<ae_protein*>();
  _protein_list[LAGGING]  = new ae_list<ae_protein*>();

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

ae_genetic_unit::ae_genetic_unit(ae_individual* indiv, const ae_genetic_unit* parent)
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
  _dna = new ae_dna( this, parent->_dna );

  // Copy promoter list (_rna_list)
  // Note that the length of the RNA will have to be recomputed (do_transcription)
  _rna_list     = new ae_list<ae_rna*>*[2];

  for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ )
  {
    _rna_list[strand] = new ae_list<ae_rna*>();

    ae_list_node<ae_rna*>* rna_node = parent->_rna_list[strand]->get_first();
    ae_rna* rna;

    while ( rna_node != NULL )
    {
      rna = rna_node->get_obj();

      #ifndef __REGUL
        _rna_list[strand]->add( new ae_rna( this, *rna ) );
      #else
        _rna_list[strand]->add( new ae_rna_R( this, *(dynamic_cast<ae_rna_R*>(rna)) ) );
      #endif

      rna_node = rna_node->get_next();
    }
  }

  // Create an empty protein list
  _protein_list           = new ae_list<ae_protein*>*[2];
  _protein_list[LEADING]  = new ae_list<ae_protein*>();
  _protein_list[LAGGING]  = new ae_list<ae_protein*>();

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

ae_genetic_unit::ae_genetic_unit( ae_individual* indiv, gzFile backup_file )
{
  _exp_m = indiv->get_exp_m();
  _indiv = indiv;

  _transcribed                        = false;
  _translated                         = false;
  _phenotypic_contributions_computed  = false;
  _non_coding_computed                = false;
  _distance_to_target_computed        = false;
  _fitness_computed                   = false;

  _dna = new ae_dna( this, backup_file );

  gzread( backup_file, &_min_gu_length, sizeof(_min_gu_length) );
  gzread( backup_file, &_max_gu_length, sizeof(_max_gu_length) );

  _rna_list           = new ae_list<ae_rna*>*[2];
  _rna_list[LEADING]  = new ae_list<ae_rna*>();
  _rna_list[LAGGING]  = new ae_list<ae_rna*>();

  _protein_list           = new ae_list<ae_protein*>*[2];
  _protein_list[LEADING]  = new ae_list<ae_protein*>();
  _protein_list[LAGGING]  = new ae_list<ae_protein*>();

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
ae_genetic_unit::ae_genetic_unit( ae_individual* indiv, char* organism_file_name )
{
  _exp_m = indiv->get_exp_m();
  _indiv = indiv;

  _transcribed                        = false;
  _translated                         = false;
  _phenotypic_contributions_computed  = false;
  _non_coding_computed                = false;
  _distance_to_target_computed        = false;
  _fitness_computed                   = false;

  _dna = new ae_dna( this, organism_file_name );

  // Create empty rna and protein lists
  _rna_list           = new ae_list<ae_rna*>*[2];
  _rna_list[LEADING]  = new ae_list<ae_rna*>();
  _rna_list[LAGGING]  = new ae_list<ae_rna*>();

  _protein_list           = new ae_list<ae_protein*>*[2];
  _protein_list[LEADING]  = new ae_list<ae_protein*>();
  _protein_list[LAGGING]  = new ae_list<ae_protein*>();

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
ae_genetic_unit::~ae_genetic_unit( void )
{
  assert( _protein_list           != NULL );
  assert( _protein_list[LEADING]  != NULL );
  assert( _protein_list[LAGGING]  != NULL );
  _protein_list[LEADING]->erase( true );
  _protein_list[LAGGING]->erase( true );
  delete _protein_list[LEADING];
  delete _protein_list[LAGGING];
  delete [] _protein_list;

  assert( _rna_list           != NULL );
  assert( _rna_list[LEADING]  != NULL );
  assert( _rna_list[LAGGING]  != NULL );
  _rna_list[LEADING]->erase( true );
  _rna_list[LAGGING]->erase( true );
  delete _rna_list[LEADING];
  delete _rna_list[LAGGING];
  delete [] _rna_list;

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
void ae_genetic_unit::locate_promoters( void )
{
  // Look for promoters in the genome and create a new ae_rna
  // in the corresponding strand's RNA list
  int8_t dist; // Hamming distance of the sequence from the promoter consensus

  // Empty RNA list
  _rna_list[LEADING]->erase( true );
  _rna_list[LAGGING]->erase( true );

  if ( _dna->get_length() >= PROM_SIZE )
  {
    for ( int32_t i = 0 ; i < _dna->get_length() ; i++ )
    {
      if ( is_promoter( LEADING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
      {
        #ifndef __REGUL
          _rna_list[LEADING]->add( new ae_rna( this, LEADING, i, dist ) );
        #else
          _rna_list[LEADING]->add( new ae_rna_R( this, LEADING, i, dist ) );
        #endif
      }
      if ( is_promoter( LAGGING, _dna->get_length() - i - 1, dist ) )
      {
        #ifndef __REGUL
          _rna_list[LAGGING]->add( new ae_rna( this, LAGGING, _dna->get_length() - i - 1, dist ) );
        #else
          _rna_list[LAGGING]->add( new ae_rna_R( this, LAGGING, _dna->get_length() - i - 1, dist ) );
        #endif
      }
    }
  }
}

void ae_genetic_unit::do_transcription( void )
{
  if ( _transcribed ) return;
  _transcribed = true;

  ae_list_node<ae_rna*>* rna_node    = NULL;
  ae_rna* rna = NULL;
  int32_t transcript_start  = -1;
  int32_t genome_length     = _dna->get_length();

  // If the genome is not long enough to bear a promoter and a terminator,
  // we set all its RNAs to a length of -1
  // (NB but a terminator can share code with the promoter, making it
  // possible for the genome to be no longer than the promoter)
  if ( genome_length < PROM_SIZE )
  {
    rna_node = _rna_list[LEADING]->get_first();
    while ( rna_node != NULL )
    {
      rna_node->get_obj()->set_transcript_length( -1 );
      rna_node = rna_node->get_next();
    }

    rna_node = _rna_list[LAGGING]->get_first();
    while ( rna_node != NULL )
    {
      rna_node->get_obj()->set_transcript_length( -1 );
      rna_node = rna_node->get_next();
    }

    return;
  }

  // ----------------
  //  LEADING strand
  // ----------------
  rna_node = _rna_list[LEADING]->get_first();

  while ( rna_node != NULL )
  {
    rna = rna_node->get_obj();
    transcript_start = rna->get_first_transcribed_pos();
    rna->set_transcript_length( -1 );

    int32_t i;
    for ( i = 0 ; i < genome_length ; i++ )
    {
      if ( is_terminator( LEADING, transcript_start + i ) )
      {
        // Found terminator => set transcript's length
        rna->set_transcript_length( i + TERM_SIZE );

        // Deduce the length of all the RNAs that share the same terminator
        // These are the RNAs whose promoter is entirely (and strictly) included
        // between the promoter and the terminator of the RNA we have just treated.
        // They are hence the RNAs whose promoter starts at most i bases after the
        // current rna's promoter
        ae_list_node<ae_rna*>* rna_node_2  = rna_node->get_next();
        ae_rna* rna_2 = NULL;
        while ( rna_node_2 != NULL )
        {
          rna_2 = rna_node_2->get_obj();

          // We know rna_2 is after rna => rna_2->pos > rna->pos (LEADING strand) because the list is sorted
          if ( rna_2->get_promoter_pos() - rna->get_promoter_pos() <= i )
          {
            rna_2->set_transcript_length( i - (rna_2->get_promoter_pos() - rna->get_promoter_pos()) + TERM_SIZE );

            // Step forward in RNA list
            rna_node = rna_node_2;
          }
          else
          {
            // The promoter of rna_2 is after (or contains a part of) the terminator of rna,
            // we will need to search its own terminator
            break;
          }

          rna_node_2 = rna_node_2->get_next();
        }

        // Terminator found for this RNA, nothing else to do (for this RNA)
        break;
      }
    }

    if (i == genome_length)
      {
        // We have searched the whole genome and found no terminator for this promoter.
        // We consider that no RNA can actually be produced, hence we set the transcript
        // length to -1. This will prevent the search for coding sequences downstream of this promoter.
        // However, we do not destroy the ae_rna object, it must still be kept in memory and
        // transmitted to the offspring in case a mutation recreates a terminator.
        rna->set_transcript_length( -1 );
      }

    rna_node = rna_node->get_next();
  }

  // ----------------
  //  LAGGING strand
  // ----------------
  rna_node = _rna_list[LAGGING]->get_first();

  while ( rna_node != NULL )
  {
    rna = rna_node->get_obj();
    transcript_start = rna->get_first_transcribed_pos();
    rna->set_transcript_length( -1 );

    int32_t i;
    for ( i = 0 ; i < genome_length ; i++ )
    {
      if ( is_terminator( LAGGING, transcript_start - i ) )
      {
        // Found terminator => set transcript's length
        rna->set_transcript_length( i + TERM_SIZE );

        // Deduce the length of all the RNAs that share the same terminator
        // These are the RNAs whose promoter is entirely (and strictly) included
        // between the promoter and the terminator of the RNA we have just treated.
        // They are hence the RNAs whose promoter starts at most i bases after the
        // current rna's promoter
        ae_list_node<ae_rna*>* rna_node_2  = rna_node->get_next();
        ae_rna* rna_2 = NULL;
        while ( rna_node_2 != NULL )
        {
          rna_2 = rna_node_2->get_obj();

          // We know rna_2 is after rna => rna_2->pos < rna->pos (LAGGING strand) because the list is sorted
          if ( rna->get_promoter_pos() - rna_2->get_promoter_pos() <= i )
          {
            rna_2->set_transcript_length( i - (rna->get_promoter_pos() - rna_2->get_promoter_pos()) + TERM_SIZE );

            // Step forward in RNA list
            rna_node = rna_node_2;
          }
          else
          {
            // The promoter of rna_2 is after (or contains a part of) the terminator of rna,
            // we will need to search its own terminator
            break;
          }

          rna_node_2 = rna_node_2->get_next();
        }

        // Terminator found for this RNA, nothing else to do (for this RNA)
        break;
      }
    }

    if (i == genome_length)
      {
        // We have searched the whole genome and found no terminator for this promoter.
        // We consider that no RNA can actually be produced, hence we set the transcript
        // length to -1. This will prevent the search for coding sequences downstream of this promoter.
        // However, we do not destroy the ae_rna object, it must still be kept in memory and
        // transmitted to the offspring in case a mutation recreates a terminator.
        rna->set_transcript_length( -1 );
      }

    rna_node = rna_node->get_next();
  }


  /******************** DEBUG (print rnas' sequences, positions and strands ********************/
  /*for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ )
  {
    rna_node = _rna_list[strand]->get_first();

    while ( rna_node != NULL )
    {
      rna = rna_node->get_obj();

      char* seq = new char[PROM_SIZE + rna->get_transcript_length() + 1];

      //~ if ( rna->get_promoter_pos() + PROM_SIZE + rna->get_transcript_length() < _dna->get_length() )
      //~ {
        //~ memcpy( seq, &_dna->get_data()[rna->get_promoter_pos()], PROM_SIZE + rna->get_transcript_length()  );
      //~ }
      //~ else
      //~ {
        //~ memcpy( seq, &_dna->get_data()[rna->get_promoter_pos()], _dna->get_length() - rna->get_promoter_pos()  );
        //~ memcpy( &seq[_dna->get_length() - rna->get_promoter_pos()], _dna->get_data(),
        //~ PROM_SIZE + rna->get_transcript_length() - (_dna->get_length() - rna->get_promoter_pos())  );
      //~ }

      //~ seq[PROM_SIZE + rna->get_transcript_length() ] = '\0';

      //~ printf( "rna seq : %s\n", seq );
      printf( "RNA at pos : %"PRId32"      length : %"PRId32"\n", rna->get_promoter_pos(), rna->get_transcript_length() );
      printf( "  strand : %s    basal_level : %f\n", (rna->get_strand() == LEADING)?"LEADING":"LAGGING", rna->get_basal_level() );
      //~ getchar();

      rna_node = rna_node->get_next();
    }
  }*/
  /******************** END DEBUG ********************/
}

void ae_genetic_unit::do_translation( void )
{
  if ( _translated ) return;
  _translated = true;
  if ( ! _transcribed ) do_transcription();

  ae_list_node<ae_rna*>* rna_node    = NULL;
  ae_rna* rna = NULL;
  int32_t transcript_start  = -1;
  int32_t transcript_length = -1;
  int32_t genome_length     = _dna->get_length();

  // ----------------
  //  LEADING strand
  // ----------------
  rna_node = _rna_list[LEADING]->get_first();

  while ( rna_node != NULL )
  {
    rna = rna_node->get_obj();
    transcript_start  = rna->get_first_transcribed_pos();
    transcript_length = rna->get_transcript_length();

    // Try every position where a translation process could occur
    // Minimum number of bases needed is SHINE_DAL_SIZE + SHINE_START_SPACER + 3 * CODON_SIZE
    // (3 codons for START + STOP + at least one amino-acid)
    for ( int32_t i = 0 ; transcript_length - i >= SHINE_DAL_SIZE + SHINE_START_SPACER + 3 * CODON_SIZE ; i++ )
    {
      if (  ( is_shine_dalgarno( LEADING, ae_utils::mod(transcript_start + i, genome_length) ) ) &&
            ( is_start( LEADING, ae_utils::mod(transcript_start + i + SHINE_DAL_SIZE + SHINE_START_SPACER, genome_length) ) ) )
      {
        // We found a translation initiation, we can now build the protein until we find a STOP codon or until we reach the end
        // of the transcript (in which case the protein is not valid)


        // First of all, we will check whether this CDS has already been translated (because it is present on another RNA
        // In that case, we don't need to tranlate it again, we only need to increase the protein's concentration according to
        // the promoter transcription level
        int32_t shine_dal_pos = ae_utils::mod(transcript_start + i, genome_length);
        ae_list_node<ae_protein*>* protein_node = _protein_list[LEADING]->bsearch( &shine_dal_pos, compare_prot_pos );

        if ( protein_node != NULL )
        {
          ae_protein* protein = protein_node->get_obj();
          protein->add_RNA( rna );
          rna->add_transcribed_protein( protein );
        }
        else
        {
          // Build codon list and make new protein when stop found
          int32_t j = i + SHINE_DAL_SIZE + SHINE_START_SPACER + CODON_SIZE; // next codon to examine
          ae_codon* codon;
          std::list<ae_codon*> codon_list;
          //~ int32_t nb_m, nb_w, nb_h;  // Number of M, W, H-codons found in the gene
          //~ nb_m = nb_w = nb_h = 0; // TODO : usefull?

          while ( (transcript_length - j >= CODON_SIZE) )
          {
            codon = new ae_codon( _dna, LEADING, ae_utils::mod(transcript_start + j, genome_length) );

            if ( codon->is_stop() )
            {
              if (not codon_list.empty()) // at least one amino-acid
              {
                // The protein is valid, create the corresponding object
                ae_protein* protein;
                #ifndef __REGUL
                  protein = new ae_protein( this, codon_list, LEADING, shine_dal_pos, rna );
                #else
                  protein = new ae_protein_R( this, new ae_list<ae_codon*>(codon_list), LEADING, shine_dal_pos, rna );
                #endif

                _protein_list[LEADING]->add( protein );
                rna->add_transcribed_protein( protein );

                if ( protein->get_is_functional() )
                {
                  _nb_fun_genes++;
                  //~ _overall_size_fun_genes += ( protein->get_length() + 2 ) * CODON_SIZE;
                  _overall_size_fun_genes += protein->get_length() * CODON_SIZE;

                  if ( protein->get_height() > 0 )  _nb_genes_activ++;
                  else                              _nb_genes_inhib++;
                }
                else
                {
                  _nb_non_fun_genes++;
                  //~ _overall_size_non_fun_genes += ( protein->get_length() + 2 ) * CODON_SIZE;
                  _overall_size_non_fun_genes += protein->get_length() * CODON_SIZE;
                }
              }

              delete codon;
              break;
            }
            else
            {
              codon_list.push_back(codon);
            }

            j += CODON_SIZE;
          }
        }
      }
    }

    // Statistics
    if (not rna->get_transcribed_proteins().empty()) // coding RNA
    {
      _nb_coding_RNAs++;
      _overall_size_coding_RNAs += rna->get_transcript_length();
    }
    else // non-coding RNA
    {
      _nb_non_coding_RNAs++;
      _overall_size_non_coding_RNAs += rna->get_transcript_length();
    }

    rna_node = rna_node->get_next();
  }

  // ----------------
  //  LAGGING strand
  // ----------------
  rna_node = _rna_list[LAGGING]->get_first();

  while ( rna_node != NULL )
  {
    rna = rna_node->get_obj();
    transcript_start  = rna->get_first_transcribed_pos();
    transcript_length = rna->get_transcript_length();

    // Try every position where a translation process could occur
    // Minimum number of bases needed is SHINE_DAL_SIZE + SHINE_START_SPACER + 3 * CODON_SIZE
    // (3 codons for START + STOP + at least one amino-acid)
    for ( int32_t i = 0 ; transcript_length - i >= SHINE_DAL_SIZE + SHINE_START_SPACER + 3 * CODON_SIZE ; i++ )
    {
      if (  ( is_shine_dalgarno( LAGGING, ae_utils::mod(transcript_start - i, genome_length) ) ) &&
            ( is_start( LAGGING, ae_utils::mod(transcript_start - i - SHINE_DAL_SIZE - SHINE_START_SPACER, genome_length) ) ) )
      {
        // We found a translation initiation, we can now build the protein until we find a STOP codon or until we reach the end
        // of the transcript (in which case the protein is not valid)

        // First of all, we will check whether this CDS has already been translated (because it is present on another RNA
        // In that case, we don't need to tranlate it again, we only need to increase the protein's concentration according to
        // the promoter strength
        int32_t shine_dal_pos = ae_utils::mod(transcript_start - i, genome_length);
        ae_list_node<ae_protein*>* protein_node = _protein_list[LAGGING]->bsearch( &shine_dal_pos, compare_prot_pos );

        if ( protein_node != NULL )
        {
          ae_protein* protein = protein_node->get_obj();
          protein->add_RNA( rna );
          rna->add_transcribed_protein( protein );
        }
        else
        {
          // Build codon list and make new protein when stop found
          int32_t j = i + SHINE_DAL_SIZE + SHINE_START_SPACER + CODON_SIZE; // next codon to examine
          ae_codon* codon;
          ae_list<ae_codon*>* codon_list = new ae_list<ae_codon*>();
          //~ int32_t nb_m, nb_w, nb_h; // Number of M, W, H-codons found in the gene
          //~ nb_m = nb_w = nb_h = 0;

          while ( (transcript_length - j >= CODON_SIZE) )
          {
            codon = new ae_codon( _dna, LAGGING, ae_utils::mod(transcript_start - j, genome_length) );

            if ( codon->is_stop() )
            {
              if ( codon_list->is_empty() == false ) // at least one amino-acid
              {
                // The protein is valid, create the corresponding object
                ae_protein* protein;
                #ifndef __REGUL
                  protein = new ae_protein(this, aelist_to_stdlist(codon_list), LAGGING, shine_dal_pos, rna);
                #else
                  protein = new ae_protein_R(this, aelist_to_stdlist(codon_list), LAGGING, shine_dal_pos, rna);
                #endif

                // The codon list will be kept in the protein
                codon_list = NULL;

                _protein_list[LAGGING]->add( protein );
                rna->add_transcribed_protein( protein );

                if ( protein->get_is_functional() )
                {
                  _nb_fun_genes++;
                  //~ _overall_size_fun_genes += ( protein->get_length() + 2 ) * CODON_SIZE;
                  _overall_size_fun_genes += protein->get_length() * CODON_SIZE;

                  if ( protein->get_height() > 0 )  _nb_genes_activ++;
                  else                              _nb_genes_inhib++;
                }
                else
                {
                  _nb_non_fun_genes++;
                  _overall_size_non_fun_genes += ( protein->get_length() + 2 ) * CODON_SIZE;
                }
              }

              delete codon;
              break;
            }
            else
            {
              codon_list->add( codon );
            }

            j += CODON_SIZE;
          }

          // The codon list is no longer useful, delete it with all its items
          // TODO : memory leek in RAEVOL?
//          #ifndef __REGUL
            if ( codon_list != NULL )
            {
              codon_list->erase( true );
              delete codon_list;
            }
//          #endif
        }
      }
    }

    // Statistics
    if (not rna->get_transcribed_proteins().empty()) // coding RNA
    {
      _nb_coding_RNAs++;
      _overall_size_coding_RNAs += rna->get_transcript_length();
    }
    else // non-coding RNA
    {
      _nb_non_coding_RNAs++;
      _overall_size_non_coding_RNAs += rna->get_transcript_length();
    }

    rna_node = rna_node->get_next();
  }
}

void ae_genetic_unit::compute_phenotypic_contribution( void )
{
  if ( _phenotypic_contributions_computed ) return;
  _phenotypic_contributions_computed = true;
  if ( ! _translated ) do_translation();

  ae_list_node<ae_protein*>* prot_node;
  ae_protein*   prot;

  // LEADING strand
  prot_node = _protein_list[LEADING]->get_first();
  while ( prot_node != NULL )
  {
    prot = prot_node->get_obj();

    if ( prot->get_is_functional() )
    {
      if ( prot->get_height() > 0 )
      {
        _activ_contribution->add_triangle(  prot->get_mean(),
                                            prot->get_width(),
                                            prot->get_height() * prot->get_concentration() );
      }
      else
      {
        _inhib_contribution->add_triangle(  prot->get_mean(),
                                            prot->get_width(),
                                            prot->get_height() * prot->get_concentration() );
      }
    }

    prot_node = prot_node->get_next();
  }

  // LAGGING strand
  prot_node = _protein_list[LAGGING]->get_first();
  while ( prot_node != NULL )
  {
    prot = prot_node->get_obj();

    if ( prot->get_is_functional() )
    {
      if ( prot->get_height() > 0 )
      {
        _activ_contribution->add_triangle(  prot->get_mean(),
                                            prot->get_width(),
                                            prot->get_height() * prot->get_concentration() );
      }
      else
      {
        _inhib_contribution->add_triangle(  prot->get_mean(),
                                            prot->get_width(),
                                            prot->get_height() * prot->get_concentration() );
      }
    }

    prot_node = prot_node->get_next();
  }


  // It is not necessary to add a lower bound to _activ_contribution as there can be no negative y
  // The same goes for the upper bound for _inhib_contribution
  _activ_contribution->clip(Fuzzy::max,   Y_MAX );
  _inhib_contribution->clip(Fuzzy::min, - Y_MAX );
  _activ_contribution->simplify();
  _inhib_contribution->simplify();

  if ( _exp_m->get_output_m()->get_compute_phen_contrib_by_GU() )
  {
    _phenotypic_contribution = new ae_phenotype();
    _phenotypic_contribution->add( *_activ_contribution );
    _phenotypic_contribution->add( *_inhib_contribution );
    _phenotypic_contribution->simplify();
  }
}

/*!
  \brief Compute the areas between the phenotype and the environment for each environmental segment.

  If the environment is not segmented, the total area is computed
*/
void ae_genetic_unit::compute_distance_to_target( Environment* env )
{
  if ( _distance_to_target_computed ) return; // _distance_to_target has already been computed, nothing to do.
  _distance_to_target_computed = true;

  compute_phenotypic_contribution();

  // Compute the difference between the (whole) phenotype and the environment
  Fuzzy* delta = new Fuzzy( *_phenotypic_contribution );
  delta->sub( *env );

  ae_env_segment** segments = env->get_segments();

  // TODO : We should take into account that we compute the areas in order (from the leftmost segment, rightwards)
  //   => We shouldn't parse the whole list of points on the left of the segment we are considering (we have
  //      already been through them!)

  if ( _dist_to_target_per_segment == NULL )
  {
    _dist_to_target_per_segment = new double [env->get_nb_segments()]; // Can not be allocated in constructor because number of segments is then unknow
  }
  for ( size_t i = 0 ; i < env->get_nb_segments() ; i++ )
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
void ae_genetic_unit::compute_fitness( Environment* env )
{
  if ( _fitness_computed ) return; // Fitness has already been computed, nothing to do.
  _fitness_computed = true;

#ifdef NORMALIZED_FITNESS

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    if (env->get_area_by_feature(i)==0.)
    {
      _fitness_by_feature[i] = 0.;
    }
    else
    {
      _fitness_by_feature[i] =  ( env->get_area_by_feature(i) - _dist_to_target_by_feature[i] ) / env->get_area_by_feature(i);
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
                                      - exp( -_exp_m->get_selection_pressure() * env->get_area_by_feature(SECRETION) );

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
                ( 1 + _exp_m->get_secretion_contrib_to_fitness() * _indiv->get_grid_cell()->get_compound_amount()
                    - _exp_m->get_secretion_cost() * _fitness_by_feature[SECRETION] );
  }

#endif

}


void ae_genetic_unit::reset_expression( void )
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
  // during the mutations (cf ae_dna::undergo_this_mutation)
  // TODO : Reinitialize _transcribed proteins ?

  _protein_list[LEADING]->erase(true);
  _protein_list[LAGGING]->erase(true);

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


void ae_genetic_unit::print_coding_rnas()
{
  printf( "  LEADING \n" );
  ae_list_node<ae_rna*>* a_node = _rna_list[LEADING]->get_first();
  ae_rna * rna = NULL;
  while (a_node != NULL)
    {
      rna = (ae_rna *) a_node->get_obj();
      a_node = a_node->get_next();
      if (rna->is_coding() == false) continue;
      printf("Promoter at %" PRId32 ", last transcribed position at %" PRId32 "\n", rna->get_promoter_pos(), rna->get_last_transcribed_pos());
    }

  printf( "  LAGGING \n" );
  a_node = _rna_list[LAGGING]->get_first();
  rna = NULL;
  while (a_node != NULL)
    {
      rna = (ae_rna *) a_node->get_obj();
      a_node = a_node->get_next();
      if (rna->is_coding() == false) continue;
      printf("Promoter at %" PRId32 ", last transcribed position at %" PRId32 "\n", rna->get_promoter_pos(), rna->get_last_transcribed_pos());
    }
}


void ae_genetic_unit::print_rnas( ae_list<ae_rna*>* rnas, ae_strand strand )
{
  ae_list_node<ae_rna*>* rna_node = NULL;
  ae_rna* rna = NULL;

  printf( "  %s ( %" PRId32 " )\n", strand == LEADING ? " LEADING " : "LAGGING", rnas->get_nb_elts() );
  rna_node = rnas->get_first();

  while ( rna_node != NULL )
  {
    rna = rna_node->get_obj();

    assert( rna->get_strand() == strand );

    printf( "    Promoter on %s at %" PRId32 "\n", strand == LEADING ? " LEADING " : "LAGGING", rna->get_promoter_pos() );
    //~ printf( "      length %" PRId32 "  basal_level %f\n", rna->get_transcript_length(), rna->get_basal_level() );

    rna_node = rna_node->get_next();
  }
}

void ae_genetic_unit::print_proteins( void ) const
{
  ae_list_node<ae_protein*>* prot_node = NULL;
  ae_protein*   prot      = NULL;

  printf( "  LEADING ( %" PRId32 " )\n", _protein_list[LEADING]->get_nb_elts() );
  prot_node = _protein_list[LEADING]->get_first();

  while ( prot_node != NULL )
  {
    prot = prot_node->get_obj();

    printf( "    Gene on LEADING at %" PRId32 " (%" PRId32 ") (%f %f %f) (%f) %s\n",
            prot->get_shine_dal_pos(), prot->get_length(),
            prot->get_mean(), prot->get_width(), prot->get_height(), prot->get_concentration(),
            prot->get_is_functional() ? "functional" : "non functional" );

    prot_node = prot_node->get_next();
  }


  printf( "  LAGGING ( %" PRId32 " )\n", _protein_list[LAGGING]->get_nb_elts() );
  prot_node = _protein_list[LAGGING]->get_first();

  while ( prot_node != NULL )
  {
    prot = prot_node->get_obj();

    printf( "    Gene on LAGGING at %" PRId32 " (%" PRId32 ") (%f %f %f) (%f) %s\n",
            prot->get_shine_dal_pos(), prot->get_length(),
            prot->get_mean(), prot->get_width(), prot->get_height(), prot->get_concentration(),
            prot->get_is_functional() ? "functional" : "non functional" );

    prot_node = prot_node->get_next();
  }
}

bool ae_genetic_unit::is_promoter( ae_strand strand, int32_t pos, int8_t& dist ) const
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
      if ( genome[ae_utils::mod((pos-i),len)] == PROM_SEQ[i] ) // == and not != because we are on the complementary strand...
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

bool ae_genetic_unit::is_terminator( ae_strand strand, int32_t pos ) const
{
  const char* genome  = _dna->get_data();
  int32_t  len        = _dna->get_length();

  if ( strand == LEADING )
  {
    for ( int16_t i = 0 ; i < TERM_STEM_SIZE ; i++ )
    {
      if ( genome[ae_utils::mod(pos+i,len)] == genome[ae_utils::mod(pos+(TERM_SIZE-1)-i,len)] ) return false;
    }
  }
  else // ( strand == LAGGING )
  {
    for ( int16_t i = 0 ; i < TERM_STEM_SIZE ; i++ )
    {
      if ( genome[ae_utils::mod(pos-i,len)] == genome[ae_utils::mod(pos-(TERM_SIZE-1)+i,len)] ) return false;
    }
  }

  return true;
}

bool ae_genetic_unit::is_shine_dalgarno( ae_strand strand, int32_t pos ) const
{
  const char* genome  = _dna->get_data();
  int32_t  len        = _dna->get_length();

  if ( strand == LEADING )
  {
    for ( int8_t i = 0 ; i < SHINE_DAL_SIZE ; i++ )
    {
      if ( genome[ae_utils::mod((pos+i),len)] != SHINE_DAL_SEQ[i] )
      {
        return false;
      }
    }
  }
  else // ( strand == LAGGING )
  {
    for ( int8_t i = 0 ; i < SHINE_DAL_SIZE ; i++ )
    {
      if ( genome[ae_utils::mod((pos-i),len)] == SHINE_DAL_SEQ[i] ) // == and not != because we are on the complementary strand...
      {
        return false;
      }
    }
  }

  return true;
}

int8_t ae_genetic_unit::get_codon( ae_strand strand, int32_t pos ) const
{
  const char* genome  = _dna->get_data();
  int32_t  len        = _dna->get_length();
  int8_t codon        = 0;

  if ( strand == LEADING )
  {
    for ( int8_t i = 0 ; i < CODON_SIZE ; i++ )
    {
      if ( genome[ae_utils::mod((pos+i),len)] == '1' )
      {
        codon += 1 << ( CODON_SIZE - i - 1 ); //pow( 2, CODON_SIZE - i - 1 );
      }
    }
  }
  else // ( strand == LAGGING )
  {
    for ( int8_t i = 0 ; i < CODON_SIZE ; i++ )
    {
      if ( genome[ae_utils::mod((pos-i),len)] != '1' ) // == and not != because we are on the complementary strand...
      {
        codon += 1 << ( CODON_SIZE - i - 1 ); //pow( 2, CODON_SIZE - i - 1 );
      }
    }
  }

  return codon;
}

void ae_genetic_unit::compute_non_coding( void )
{
  if ( _non_coding_computed ) return;
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
  for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ )
  {
    ae_list_node<ae_protein*>* prot_node = _protein_list[strand]->get_first();
    ae_protein* prot = NULL;

    while ( prot_node != NULL )
    {
      prot = prot_node->get_obj();
      int32_t first;
      int32_t last;

      if ( strand == LEADING )
      {
        first = prot->get_shine_dal_pos();
        last  = prot->get_last_STOP_base_pos();
      }
      else // ( strand == LAGGING )
      {
        last  = prot->get_shine_dal_pos();
        first = prot->get_last_STOP_base_pos();
      }

      if ( first <= last )
      {
        for ( int32_t i = first ; i <= last ; i++ )
        {
          belongs_to_CDS[i] = true;
          if ( prot->get_is_functional() ) is_essential_DNA[i] = true;
          is_essential_DNA_including_nf_genes[i] = true;
        }
      }
      else
      {
        for ( int32_t i = first ; i < genome_length ; i++ )
        {
          belongs_to_CDS[i] = true;
          if ( prot->get_is_functional() ) is_essential_DNA[i] = true;
          is_essential_DNA_including_nf_genes[i] = true;
        }
        for ( int32_t i = 0 ; i <= last ; i++ )
        {
          belongs_to_CDS[i] = true;
          if ( prot->get_is_functional() ) is_essential_DNA[i] = true;
          is_essential_DNA_including_nf_genes[i] = true;
        }
      }

      // Include the promoter and terminator to essential DNA
      // Mark everything between promoter and terminator as not neutral
      for (const auto& rna: prot->get_rna_list_std()) {

        int32_t prom_first;
        int32_t prom_last;
        int32_t term_first;
        int32_t term_last;
        int32_t rna_first;
        int32_t rna_last;

        if ( strand == LEADING )
        {
          prom_first  = rna->get_promoter_pos();
          prom_last   = ae_utils::mod( prom_first + PROM_SIZE - 1, _dna->get_length() );
          term_last   = rna->get_last_transcribed_pos();
          term_first  = ae_utils::mod( term_last - TERM_SIZE + 1, _dna->get_length() );
          rna_first   = prom_first;
          rna_last    = term_last;
        }
        else
        {
          prom_last   = rna->get_promoter_pos();
          prom_first  = ae_utils::mod( prom_last - PROM_SIZE + 1, _dna->get_length() );
          term_first  = rna->get_last_transcribed_pos();
          term_last   = ae_utils::mod( term_first + TERM_SIZE - 1, _dna->get_length() );
          rna_first   = term_first;
          rna_last    = prom_last;
        }
        //~ printf( "\n" );
        //~ if ( strand == LEADING ) printf( "LEADING\n" );
        //~ else printf( "LAGGING\n" );
        //~ printf( "prom_first : %ld    prom_last %ld   (size %ld)\n", prom_first, prom_last, _dna->get_length() );
        //~ printf( "term_first : %ld    term_last %ld   (size %ld)\n", term_first, term_last, _dna->get_length() );
        //~ getchar();

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
        //~ printf( "prom " );
        if ( prom_first <= prom_last )
        {
          for ( int32_t i = prom_first ; i <= prom_last ; i++ )
          {
            //~ printf( "%ld ", i );
            if ( prot->get_is_functional() ) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
        }
        else
        {
          for ( int32_t i = prom_first ; i < genome_length ; i++ )
          {
            //~ printf( "%ld ", i );
            if ( prot->get_is_functional() ) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
          for ( int32_t i = 0 ; i <= prom_last ; i++ )
          {
            //~ printf( "%ld ", i );
            if ( prot->get_is_functional() ) is_essential_DNA[i] = true;
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
            if ( prot->get_is_functional() ) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
        }
        else
        {
          for ( int32_t i = term_first ; i < genome_length ; i++ )
          {
            //~ printf( "%ld ", i );
            if ( prot->get_is_functional() ) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
          for ( int32_t i = 0 ; i <= term_last ; i++ )
          {
            //~ printf( "%ld ", i );
            if ( prot->get_is_functional() ) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
        }
        //~ printf( "\n" );
        //~ getchar();
      }




      if ( prot->get_is_functional() )
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

      prot_node = prot_node->get_next();
    }
  }


  // Parse RNA lists and mark the corresponding bases as coding (only for the coding RNAs)
  for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ )
  {
    ae_list_node<ae_rna*>* rna_node  = _rna_list[strand]->get_first();
    ae_rna* rna = NULL;

    while ( rna_node != NULL )
    {
      rna = rna_node->get_obj();

      int32_t first;
      int32_t last;

      if ( strand == LEADING )
      {

        first = rna->get_promoter_pos();
        last  = rna->get_last_transcribed_pos();

      }
      else // ( strand == LAGGING )
      {
        last  = rna->get_promoter_pos();

        first = rna->get_last_transcribed_pos();
      }

      if ( first <= last )
      {
        for ( int32_t i = first ; i <= last ; i++ )
        {
          belongs_to_RNA[i] = true;
        }
      }
      else
      {
        for ( int32_t i = first ; i < genome_length ; i++ )
        {
          belongs_to_RNA[i] = true;
        }
        for ( int32_t i = 0 ; i <= last ; i++ )
        {
          belongs_to_RNA[i] = true;
        }
      }

      if (not rna->get_transcribed_proteins().empty()) // coding RNA
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
      else // non coding RNA
      {
        if ( first <= last )
        {
          for ( int32_t i = first ; i <= last ; i++ )
          {
            belongs_to_non_coding_RNA[i] = true;
          }
        }
        else
        {
          for ( int32_t i = first ; i < genome_length ; i++ )
          {
            belongs_to_non_coding_RNA[i] = true;
          }
          for ( int32_t i = 0 ; i <= last ; i++ )
          {
            belongs_to_non_coding_RNA[i] = true;
          }
        }
      }

      rna_node = rna_node->get_next();
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

void ae_genetic_unit::duplicate_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>** duplicated_promoters )
{
  // 1) Get promoters to be duplicated
  get_promoters_included_in( pos_1, pos_2, duplicated_promoters );

  // 2) Set RNAs' position as their position on the duplicated segment
  ae_list_node<ae_rna*>* rna_node  = NULL;

  // -- LEADING --
  rna_node = duplicated_promoters[LEADING]->get_first();
  while ( rna_node != NULL )
  {

#ifndef __REGUL
    // Make a copy of current RNA
    ae_rna* copy = new ae_rna( this, *(rna_node->get_obj()) );
#else
    // Make a copy of current RNA
    ae_rna_R* copy = new ae_rna_R( this, *(rna_node->get_obj()) );
#endif

    // Set RNA's position as it's position on the duplicated segment
    copy->shift_position( -pos_1, _dna->get_length() );

    // Replace current object by the copy
    // Do not delete the replaced object as it is still a valid promoter in the individual's promoter list
    rna_node->set_obj( copy );

    rna_node = rna_node->get_next();
  }

  // -- LAGGING --
  rna_node = duplicated_promoters[LAGGING]->get_first();
  while ( rna_node != NULL )
  {

#ifndef __REGUL
    // Make a copy of current RNA
    ae_rna* copy = new ae_rna( this, *(rna_node->get_obj()) );
#else
    // Make a copy of current RNA
    ae_rna_R* copy = new ae_rna_R( this, *(rna_node->get_obj()) );
#endif
    // Set RNA's position as it's position on the duplicated segment
    copy->shift_position( -pos_1, _dna->get_length() );

    // Replace current object by the copy (don't delete the replaced object as it is still a valid promoter)
    rna_node->set_obj( copy );

    rna_node = rna_node->get_next();
  }
}

void ae_genetic_unit::get_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>** promoters )
{
  assert( pos_1 >= 0 && pos_1 <= _dna->get_length() && pos_2 >= 0 && pos_2 <= _dna->get_length() );

  if ( pos_1 < pos_2 )
  {
    int32_t seg_length = pos_2 - pos_1;

    if ( seg_length >= PROM_SIZE )
    {
      get_leading_promoters_starting_between( pos_1, pos_2 - PROM_SIZE + 1, promoters[LEADING] );
      get_lagging_promoters_starting_between( pos_1 + PROM_SIZE - 1, pos_2, promoters[LAGGING] );
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
        get_leading_promoters_starting_after( pos_1, promoters[LEADING] );
        get_leading_promoters_starting_before( pos_2 - PROM_SIZE + 1, promoters[LEADING] );
        get_lagging_promoters_starting_before( pos_2, promoters[LAGGING] );
        get_lagging_promoters_starting_after( pos_1 + PROM_SIZE - 1, promoters[LAGGING] );
      }
      else if ( !is_near_end_of_genome ) // => && is_near_beginning_of_genome
      {
        get_leading_promoters_starting_between( pos_1, pos_2 + _dna->get_length() - PROM_SIZE + 1,
                                                promoters[LEADING] );
        get_lagging_promoters_starting_before( pos_2, promoters[LAGGING] );
        get_lagging_promoters_starting_after( pos_1 + PROM_SIZE - 1, promoters[LAGGING] );
      }
      else if ( !is_near_beginning_of_genome ) // => && is_near_end_of_genome
      {
        get_leading_promoters_starting_after( pos_1, promoters[LEADING] );
        get_leading_promoters_starting_before( pos_2 - PROM_SIZE + 1, promoters[LEADING] );
        get_lagging_promoters_starting_between( pos_1 - _dna->get_length() + PROM_SIZE - 1, pos_2,
                                                promoters[LAGGING] );
      }
      else // is_near_end_of_genome && is_near_beginning_of_genome
      {
        get_leading_promoters_starting_between( pos_1, pos_2 + _dna->get_length() - PROM_SIZE + 1,
                                                promoters[LEADING] );
        get_lagging_promoters_starting_between( pos_1 - _dna->get_length() + PROM_SIZE - 1, pos_2,
                                                promoters[LAGGING] );
      }
    }
  }
}

void ae_genetic_unit::get_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>* leading_promoters )
{
  assert( pos_1 >= 0 && pos_1 < pos_2 && pos_2 <= _dna->get_length() );

  // Go to first RNA after pos_1
  ae_list_node<ae_rna*>* rna_node  = _rna_list[LEADING]->get_first();
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_1 )
  {
    rna_node  = rna_node->get_next();
  }

  // Add RNAs to new_list until we pass pos_2 (or we reach the end of the list)
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_2 )
  {
    leading_promoters->add( rna_node->get_obj() );

    rna_node = rna_node->get_next();
  }
}

void ae_genetic_unit::get_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>* lagging_promoters )
{
  assert( pos_1 >= 0 && pos_1 < pos_2 && pos_2 <= _dna->get_length() );

  // Go to first RNA before pos_2
  ae_list_node<ae_rna*>* rna_node  = _rna_list[LAGGING]->get_first();
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos_2 )
  {
    rna_node = rna_node->get_next();
  }

  // Add RNAs to new_list until we pass pos_1 (or we reach the end of the list)
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos_1 )
  {
    lagging_promoters->add( rna_node->get_obj() );

    rna_node = rna_node->get_next();
  }
}

void ae_genetic_unit::get_leading_promoters_starting_after( int32_t pos, ae_list<ae_rna*>* leading_promoters )
{
  assert( pos >= 0 && pos < _dna->get_length() );

  // Go to first RNA after pos
  ae_list_node<ae_rna*>* rna_node  = _rna_list[LEADING]->get_first();
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos )
  {
    rna_node  = rna_node->get_next();
  }

  // Add RNAs to new_list until we reach the end of the list
  while ( rna_node != NULL )
  {
    leading_promoters->add( rna_node->get_obj() );

    rna_node = rna_node->get_next();
  }
}

void ae_genetic_unit::get_leading_promoters_starting_before( int32_t pos, ae_list<ae_rna*>* leading_promoters )
{
  assert( pos >= 0 && pos < _dna->get_length() );

  ae_list_node<ae_rna*>* rna_node  = _rna_list[LEADING]->get_first();

  // Add RNAs to new_list until we pass pos (or we reach the end of the list)
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos )
  {
    leading_promoters->add( rna_node->get_obj() );

    rna_node = rna_node->get_next();
  }
}

void ae_genetic_unit::get_lagging_promoters_starting_before( int32_t pos, ae_list<ae_rna*>* lagging_promoters )
{
  assert( pos >= 0 && pos < _dna->get_length() );

  // Go to first RNA before pos
  ae_list_node<ae_rna*>* rna_node  = _rna_list[LAGGING]->get_first();
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos )
  {
    rna_node  = rna_node->get_next();
  }

  // Add RNAs to new_list until we reach the end of the list
  while ( rna_node != NULL )
  {
    lagging_promoters->add( rna_node->get_obj() );

    rna_node = rna_node->get_next();
  }
}

void ae_genetic_unit::get_lagging_promoters_starting_after( int32_t pos, ae_list<ae_rna*>* lagging_promoters )
{
  assert( pos >= 0 && pos < _dna->get_length() );

  ae_list_node<ae_rna*>* rna_node  = _rna_list[LAGGING]->get_first();

  // Add RNAs to new_list until we pass pos (or we reach the end of the list)
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos )
  {
    lagging_promoters->add( rna_node->get_obj() );

    rna_node = rna_node->get_next();
  }
}

void ae_genetic_unit::invert_promoters_included_in( int32_t pos_1, int32_t pos_2 )
{
  assert( pos_1 >= 0 && pos_1 <= pos_2 && pos_2 <= _dna->get_length() );

  int32_t segment_length = pos_2 - pos_1;

  if ( segment_length >= PROM_SIZE )
  {
    ae_list<ae_rna*>** inverted_promoters = new ae_list<ae_rna*>*[2];
    inverted_promoters[LEADING] = new ae_list<ae_rna*>();
    inverted_promoters[LAGGING] = new ae_list<ae_rna*>();

    //~ print_rnas();

    // 1) Extract the promoters completely included on the segment to be inverted
    extract_promoters_included_in( pos_1, pos_2, inverted_promoters );

    //~ print_rnas( inverted_promoters );

    // 2) Invert segment's promoters
    ae_genetic_unit::invert_promoters( inverted_promoters, pos_1, pos_2 );

    // 3) Reinsert the inverted promoters
    insert_promoters( inverted_promoters );

    delete inverted_promoters[LEADING];
    delete inverted_promoters[LAGGING];
    delete [] inverted_promoters;
  }
}

/*!
  \brief Invert all the promoters of promoter_lists for a sequence of length seq_length.
*/
/*static*/ void ae_genetic_unit::invert_promoters( ae_list<ae_rna*>** promoter_lists, int32_t seq_length )
{
  ae_genetic_unit::invert_promoters( promoter_lists, 0, seq_length );
}

/*!
  \brief Invert all the promoters of promoter_lists knowing that they represent the promoters of a subsequence beginning at pos_1 and ending at pos_2.

  WARNING : This function is pretty specific, make sure you understand its precise behaviour before using it.
*/
/*static*/ void ae_genetic_unit::invert_promoters( ae_list<ae_rna*>** promoter_lists, int32_t pos_1, int32_t pos_2 )
{
  assert( pos_1 >= 0 && pos_1 <= pos_2 ); // Could check (pos_2 < length) but another parameter would be necessary

  // Exchange LEADING and LAGGING lists
  ae_list<ae_rna*>* tmp   = promoter_lists[LEADING];
  promoter_lists[LEADING] = promoter_lists[LAGGING];
  promoter_lists[LAGGING] = tmp;

  // Update the position and strand of each promoter to be inverted...
  ae_list_node<ae_rna*>* rna_node  = NULL;
  ae_rna* rna             = NULL;

  // ...on the former LAGGING strand (becoming the LEADING strand)
  rna_node = promoter_lists[LEADING]->get_first();
  int i = 0;
  while ( rna_node != NULL )
  {
    i++;
    rna = rna_node->get_obj();
    assert( rna->get_strand() == LAGGING );
    assert( rna->get_promoter_pos() >= pos_1 && rna->get_promoter_pos() < pos_2 );
    rna->set_promoter_pos( pos_1 + pos_2 - rna->get_promoter_pos() - 1 );
    rna->set_strand( LEADING );

    rna_node = rna_node->get_next();
  }


  // ... and on the former LEADING strand (becoming the LAGGING strand)
  rna_node = promoter_lists[LAGGING]->get_first();
  i = 0;
  while ( rna_node != NULL )
  {
    i++;
    rna = rna_node->get_obj();
    assert( rna->get_strand() == LEADING );
    assert( rna->get_promoter_pos() >= pos_1 && rna->get_promoter_pos() < pos_2 );
    rna->set_promoter_pos( pos_1 + pos_2 - rna->get_promoter_pos() - 1);
    rna->set_strand( LAGGING );

    rna_node = rna_node->get_next();
  }
}

void ae_genetic_unit::extract_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>* extracted_promoters )
{
  assert( pos_1 >= 0 && pos_1 < pos_2 && pos_2 <= _dna->get_length() );

  ae_list_node<ae_rna*>* rna_node = NULL;

  // Find the first promoters in the interval
  ae_list_node<ae_rna*>* node_first = NULL;
  rna_node = _rna_list[LEADING]->get_first();
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_1 )
  {
    rna_node  = rna_node->get_next();
  }
  if ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_2 )
  {
    node_first = rna_node;
  }

  if ( node_first != NULL )
  {
    // Find the last promoters in the interval
    ae_list_node<ae_rna*>* node_last = node_first;
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_2 )
    {
      node_last = rna_node;
      rna_node  = rna_node->get_next();
    }

    // Extract the promoters (remove them from the individual's list and put them in extracted_promoters)
    //~ print_rnas();
    ae_list<ae_rna*>* tmp_list = _rna_list[LEADING]->extract_sublist( node_first, node_last );
    //~ printf( "----------------------------------------------------\n" );
    //~ print_rnas( tmp_list, LEADING );
    //~ printf( "-------------------------==========-----------------\n" );
    extracted_promoters->merge( tmp_list );
    delete tmp_list;
  }
}

void ae_genetic_unit::extract_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>* extracted_promoters )
{
  assert( pos_1 >= 0 && pos_1 < pos_2 && pos_2 <= _dna->get_length() );

  ae_list_node<ae_rna*>* rna_node = NULL;

  // Find the first promoters in the interval (if any)
  ae_list_node<ae_rna*>* node_first = NULL;
  rna_node = _rna_list[LAGGING]->get_first();
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos_2 )
  {
    rna_node  = rna_node->get_next();
  }
  if ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos_1 )
  {
    node_first = rna_node;
  }

  // BREAKPOINT (contenu de extracted_promoters)

  if ( node_first != NULL )
  {
    // Find the last promoters in the interval
    ae_list_node<ae_rna*>* node_last = NULL;
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos_1 )
    {
      node_last = rna_node;
      rna_node  = rna_node->get_next();
    }

    //~ // <DEBUG>
    //~ printf( "_rna_list[LAGGING]  (first : %p last : %p)\n", _rna_list[LAGGING]->get_first(), _rna_list[LAGGING]->get_last() );
    //~ ae_list_node* anode = _rna_list[LAGGING]->get_first();
    //~ while ( anode )
    //~ {
      //~ printf( "  node : %p\n", anode );
      //~ anode = anode->get_next();
    //~ }
    //~ // </DEBUG>

    // Extract the promoters (remove them from the individual's list and put the in extracted_promoters)
    //~ printf( "tmp_list = _rna_list[LAGGING]->extract_sublist( %p, %p);\n", node_first, node_last );
    ae_list<ae_rna*>* tmp_list = _rna_list[LAGGING]->extract_sublist( node_first, node_last );

    //~ // <DEBUG>
    //~ printf( "tmp_list  (first : %p last : %p)\n", tmp_list->get_first(), tmp_list->get_last() );
    //~ anode = tmp_list->get_first();
    //~ while ( anode )
    //~ {
      //~ printf( "  node : %p\n", anode );
      //~ anode = anode->get_next();
    //~ }
    //~ printf( "_rna_list[LAGGING]  (first : %p last : %p)\n", _rna_list[LAGGING]->get_first(), _rna_list[LAGGING]->get_last() );
    //~ anode = _rna_list[LAGGING]->get_first();
    //~ while ( anode )
    //~ {
      //~ printf( "  node : %p\n", anode );
      //~ anode = anode->get_next();
    //~ }
    //~ // </DEBUG>

    extracted_promoters->merge( tmp_list );

    delete tmp_list;
  }
}

void ae_genetic_unit::extract_leading_promoters_starting_after( int32_t pos, ae_list<ae_rna*>* extracted_promoters )
{
  assert( pos >= 0 && pos < _dna->get_length() );

  ae_list_node<ae_rna*>* rna_node = NULL;

  // Find the first promoters in the interval
  ae_list_node<ae_rna*>* node_first = NULL;
  rna_node = _rna_list[LEADING]->get_first();
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos )
  {
    rna_node = rna_node->get_next();
  }
  node_first = rna_node;

  if ( node_first != NULL )
  {
    // Extract the promoters (remove them from the individual's list and put the in extracted_promoters)
    ae_list<ae_rna*>* tmp_list = _rna_list[LEADING]->extract_ending_sublist( node_first );
    extracted_promoters->merge( tmp_list );
    delete tmp_list;
  }
}

void ae_genetic_unit::extract_leading_promoters_starting_before( int32_t pos, ae_list<ae_rna*>* extracted_promoters )
{
  assert( pos >= 0 && pos < _dna->get_length() );

  ae_list_node<ae_rna*>* rna_node = NULL;

  // Find the last promoters in the interval
  rna_node = _rna_list[LEADING]->get_first();

  // Find the last promoters in the interval
  ae_list_node<ae_rna*>* node_last = NULL;
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos )
  {
    node_last = rna_node;
    rna_node  = rna_node->get_next();
  }

  if ( node_last != NULL )
  {
    // Extract the promoters (remove them from the individual's list and put the in extracted_promoters)
    ae_list<ae_rna*>* tmp_list = _rna_list[LEADING]->extract_starting_sublist( node_last );
    extracted_promoters->merge( tmp_list );
    delete tmp_list;
  }
}

void ae_genetic_unit::extract_lagging_promoters_starting_before( int32_t pos, ae_list<ae_rna*>* extracted_promoters )
{
  assert( pos >= 0 && pos < _dna->get_length() );

  ae_list_node<ae_rna*>* rna_node = NULL;

  // Find the first promoters in the interval
  ae_list_node<ae_rna*>* node_first = NULL;
  rna_node = _rna_list[LAGGING]->get_first();
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos )
  {
    rna_node  = rna_node->get_next();
  }
  node_first = rna_node;

  if ( node_first != NULL )
  {
    // Extract the promoters (remove them from the individual's list and put the in extracted_promoters)
    ae_list<ae_rna*>* tmp_list = _rna_list[LAGGING]->extract_ending_sublist( node_first );
    extracted_promoters->merge( tmp_list );
    delete tmp_list;
  }
}

void ae_genetic_unit::extract_lagging_promoters_starting_after( int32_t pos, ae_list<ae_rna*>* extracted_promoters )
{
  assert( pos >= 0 && pos < _dna->get_length() );

  ae_list_node<ae_rna*>* rna_node = NULL;

  // Find the last promoters in the interval
  rna_node = _rna_list[LAGGING]->get_first();

  // Find the last promoters in the interval
  ae_list_node<ae_rna*>* node_last = NULL;
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos )
  {
    node_last = rna_node;
    rna_node  = rna_node->get_next();
  }

  if ( node_last != NULL )
  {
    // Extract the promoters (remove them from the individual's list and put the in extracted_promoters)
    ae_list<ae_rna*>* tmp_list = _rna_list[LAGGING]->extract_starting_sublist( node_last );
    extracted_promoters->merge( tmp_list );
    delete tmp_list;
  }
}

/*!
  \brief Shift all the promoters in <promoters_to_shift> by <delta_pos> in a sequence of length <seq_length>.

  Every promoter in double stranded list <promoters_to_shift> will be shifted by <delta_pos>,
  then a modulo <seq_length> will be applied
*/
/*static*/ void ae_genetic_unit::shift_promoters( ae_list<ae_rna*>** promoters_to_shift, int32_t delta_pos, int32_t seq_length )
{
  ae_list_node<ae_rna*>* rna_node  = NULL;

  // -- LEADING --
  rna_node = promoters_to_shift[LEADING]->get_first();
  while ( rna_node != NULL )
  {
    rna_node->get_obj()->shift_position( delta_pos, seq_length );
    rna_node = rna_node->get_next();
  }

  // -- LAGGING --
  rna_node = promoters_to_shift[LAGGING]->get_first();
  while ( rna_node != NULL )
  {
    rna_node->get_obj()->shift_position( delta_pos, seq_length );
    rna_node = rna_node->get_next();
  }
}

/*!
  \brief Insert promoters in double stranded list <promoters_to_insert> into <this->_rna_list>.

  The promoters in <promoters_to_insert> must already be at their rightful position according to <this>
  and the positions of the promoters from <promoters_to_insert> and <this->_rna_list> must not be interlaced
  i.e. no promoter in <this->_rna_list> must have a position in [first_prom_to_insert->pos ; last_prom_to_insert->pos]
*/
void ae_genetic_unit::insert_promoters( ae_list<ae_rna*>** promoters_to_insert )
{
  ae_list_node<ae_rna*>* rna_node            = NULL;
  ae_list_node<ae_rna*>* rna_node_to_insert  = NULL;

  // -- LEADING --
  if ( promoters_to_insert[LEADING]->get_nb_elts() > 0 )
  {
    // Get to the right position in individual's list (first promoter after the inserted segment)
    int32_t pos_last_leading_to_insert = promoters_to_insert[LEADING]->get_last()->get_obj()->get_promoter_pos();
    rna_node = _rna_list[LEADING]->get_first();
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_last_leading_to_insert )
    {
      rna_node = rna_node->get_next();
    }

    // Insert the promoters in the individual's RNA list
    rna_node_to_insert = promoters_to_insert[LEADING]->get_first();
    while ( rna_node_to_insert != NULL )
    {
      if ( rna_node != NULL )
      {
        _rna_list[LEADING]->add_before( rna_node_to_insert->get_obj(), rna_node );
      }
      else
      {
        _rna_list[LEADING]->add( rna_node_to_insert->get_obj() );
      }

      rna_node_to_insert = rna_node_to_insert->get_next();
    }
  }

  // -- LAGGING --
  if ( promoters_to_insert[LAGGING]->get_nb_elts() > 0 )
  {
    // Get to the right position in individual's list (first promoter after the inserted segment)
    int32_t pos_last_lagging_to_insert = promoters_to_insert[LAGGING]->get_last()->get_obj()->get_promoter_pos();
    rna_node = _rna_list[LAGGING]->get_first();
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos_last_lagging_to_insert )
    {
      rna_node = rna_node->get_next();
    }

    // Insert the promoters in the individual's RNA list
    rna_node_to_insert = promoters_to_insert[LAGGING]->get_first();
    while ( rna_node_to_insert != NULL )
    {
      if ( rna_node != NULL )
      {
        _rna_list[LAGGING]->add_before( rna_node_to_insert->get_obj(), rna_node );
      }
      else
      {
        _rna_list[LAGGING]->add( rna_node_to_insert->get_obj() );
      }

      rna_node_to_insert = rna_node_to_insert->get_next();
    }
  }
}

/*!
  Insert promoters in double stranded list <promoters_to_insert> into <this->_rna_list> at position <pos>

  The promoters in <promoters_to_insert> must be at their rightful position according to a stand-alone sequence
  (i.e. at a RELATIVE position). Their position will be updated automatically.
*/
void ae_genetic_unit::insert_promoters_at( ae_list<ae_rna*>** promoters_to_insert, int32_t pos )
{
  ae_list_node<ae_rna*>* rna_node = NULL;
  ae_list_node<ae_rna*>* inserted_rna_node = NULL;

  // -- LEADING --
  if ( promoters_to_insert[LEADING]->get_nb_elts() > 0 )
  {
    // Get to the right position in individual's list (first promoter after the inserted segment)
    rna_node = _rna_list[LEADING]->get_first();
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos )
    {
      rna_node = rna_node->get_next();
    }

    // Insert the promoters in the individual's RNA list
    inserted_rna_node = promoters_to_insert[LEADING]->get_first();
    while ( inserted_rna_node != NULL )
    {
      // Update promoter position
      inserted_rna_node->get_obj()->shift_position( pos, _dna->get_length() );

      // Insert
      if ( rna_node != NULL )
      {
        _rna_list[LEADING]->add_before( inserted_rna_node->get_obj(), rna_node );
      }
      else
      {
        _rna_list[LEADING]->add( inserted_rna_node->get_obj() );
      }

      inserted_rna_node = inserted_rna_node->get_next();
    }
  }

  // -- LAGGING --
  if ( promoters_to_insert[LAGGING]->get_nb_elts() > 0 )
  {
    // Get to the right position in individual's list (first promoter after the inserted segment)
    rna_node = _rna_list[LAGGING]->get_first();
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos )
    {
      rna_node = rna_node->get_next();
    }

    // Insert the duplicated promoters in the individual's RNA list
    inserted_rna_node = promoters_to_insert[LAGGING]->get_first();
    while ( inserted_rna_node != NULL )
    {
      // Update promoter position
      inserted_rna_node->get_obj()->shift_position( pos, _dna->get_length() );

      // Insert
      if ( rna_node != NULL )
      {
        _rna_list[LAGGING]->add_before( inserted_rna_node->get_obj(), rna_node );
      }
      else
      {
        _rna_list[LAGGING]->add( inserted_rna_node->get_obj() );
      }

      inserted_rna_node = inserted_rna_node->get_next();
    }
  }
}


/*!
  \brief Remove the RNAs of the LEADING strand whose starting positions lie in [pos_1 ; pos_2[
*/
void ae_genetic_unit::remove_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2 )
{
  assert( pos_1 >= 0 && pos_1 < _dna->get_length() && pos_2 >= 0 && pos_2 <= _dna->get_length() );

  if ( pos_1 > pos_2 )
  {
    remove_leading_promoters_starting_after( pos_1 );
    remove_leading_promoters_starting_before( pos_2 );
  }
  else
  {
    ae_list_node<ae_rna*>* rna_node  = _rna_list[LEADING]->get_first();

    // Go to first RNA after pos_1
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_1 )
    {
      rna_node  = rna_node->get_next();
    }

    // Delete RNAs until we pass pos_2 (or we reach the end of the list)
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_2 )
    {
      ae_list_node<ae_rna*>* next_node = rna_node->get_next();
      //~ printf( "remove LEADING promoter at [%"PRId32", %"PRId32"]\n", rna_node->get_obj()->get_promoter_pos(),
              //~ ae_utils::mod( rna_node->get_obj()->get_promoter_pos() + PROM_SIZE, _dna->get_length() ) );

      _rna_list[LEADING]->remove( rna_node, true, true );

      rna_node = next_node;
    }
  }
}


/*!
  \brief Remove the RNAs of the LAGGING strand whose starting positions lie in [pos_1 ; pos_2[

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
void ae_genetic_unit::remove_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2 )
{
  assert( pos_1 >= 0 && pos_1 <= _dna->get_length() && pos_2 >= 0 && pos_2 <= _dna->get_length() );

  if ( pos_1 == _dna->get_length() ) pos_1 = 0;
  if ( pos_2 == 0 )                  pos_2 = _dna->get_length();

  if ( pos_1 > pos_2 )
  {
    remove_lagging_promoters_starting_after( pos_1 );
    remove_lagging_promoters_starting_before( pos_2 );
  }
  else
  {
    ae_list_node<ae_rna*>* rna_node  = _rna_list[LAGGING]->get_first();

    // Go to first RNA before pos_2
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos_2 )
    {
      rna_node  = rna_node->get_next();
    }

    // Delete RNAs until we pass pos_1 (or we reach the end of the list)
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos_1 )
    {
      ae_list_node<ae_rna*>* next_node = rna_node->get_next();

      //~ printf( "remove LAGGING promoter at [%"PRId32", %"PRId32"]\n", rna_node->get_obj()->get_promoter_pos(),
              //~ ae_utils::mod( rna_node->get_obj()->get_promoter_pos() - PROM_SIZE, _dna->get_length() ) );
      _rna_list[LAGGING]->remove( rna_node, true, true );

      rna_node = next_node;
    }
  }
}


/*!
  \brief Remove the promoters from the LEADING strand whose starting positions are < pos
*/
void ae_genetic_unit::remove_leading_promoters_starting_before( int32_t pos )
{
  assert( pos >= 0 && pos < _dna->get_length() );

  ae_list_node<ae_rna*>* rna_node  = _rna_list[LEADING]->get_first();

  // Delete RNAs until we reach pos (or we reach the end of the list)
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos )
  {
    ae_list_node<ae_rna*>* next_node = rna_node->get_next();

    _rna_list[LEADING]->remove( rna_node, true, true );

    rna_node = next_node;
  }
}


/*!
  \brief Remove the promoters from the LAGGING strand whose starting positions are < pos

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
void ae_genetic_unit::remove_lagging_promoters_starting_before( int32_t pos )
{
  assert( pos >= 0 && pos < _dna->get_length() );

  ae_list_node<ae_rna*>* rna_node  = _rna_list[LAGGING]->get_last();

  // Delete RNAs until we reach pos (or we reach the beginning of the list )
  ae_list_node<ae_rna*>* prev_node = NULL;
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos )
  {
    prev_node = rna_node->get_prev();

    _rna_list[LAGGING]->remove( rna_node, true, true );

    rna_node = prev_node;
  }
}


/*!
  \brief Remove the promoters from the LEADING strand whose starting positions are >= pos
*/
void ae_genetic_unit::remove_leading_promoters_starting_after( int32_t pos )
{
  assert( pos >= 0 && pos < _dna->get_length() );

  ae_list_node<ae_rna*>* rna_node = _rna_list[LEADING]->get_last();

  // Delete RNAs until we pass pos ( or we reach the beginning of the list)
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos )
  {
    ae_list_node<ae_rna*>* prev_node = rna_node->get_prev();
    //~ printf( "remove LEADING promoter at [%"PRId32", %"PRId32"]\n", rna_node->get_obj()->get_promoter_pos(),
            //~ ae_utils::mod( rna_node->get_obj()->get_promoter_pos() + PROM_SIZE, _dna->get_length() ) );

    _rna_list[LEADING]->remove( rna_node, true, true );

    rna_node = prev_node;
  }
}


/*!
  \brief Remove the promoters from the LAGGING strand whose starting positions are >= pos

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
void ae_genetic_unit::remove_lagging_promoters_starting_after( int32_t pos )
{
  assert( pos < _dna->get_length() && pos >= 0 );

  ae_list_node<ae_rna*>* rna_node  = _rna_list[LAGGING]->get_first();

  // Delete RNAs until we pass pos (or we reach the end of the list)
  ae_list_node<ae_rna*>* next_node = NULL;
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos )
  {
    next_node = rna_node->get_next();

    //~ printf( "remove LAGGING promoter at [%"PRId32", %"PRId32"]\n", rna_node->get_obj()->get_promoter_pos(),
            //~ ae_utils::mod( rna_node->get_obj()->get_promoter_pos() - PROM_SIZE, _dna->get_length() ) );
    _rna_list[LAGGING]->remove( rna_node, true, true );

    rna_node  = next_node;
  }
}


/*!
  \brief Look for new promoters on the LEADING strand whose starting positions would lie in [pos_1 ; pos_2[
*/
void ae_genetic_unit::look_for_new_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2 )
{
  assert( pos_1 >= 0 && pos_1 < _dna->get_length() && pos_2 >= 0 && pos_2 < _dna->get_length() );

  // When pos_1 > pos_2, we will perform the search in 2 steps.
  // As positions  0 and _dna->get_length() are equivalent, it's preferable to
  // keep 0 for pos_1 and _dna->get_length() for pos_2.
  //~ if ( pos_2 == 0 ) pos_2 = _dna->get_length();


  if ( pos_1 >= pos_2 )
  {
    look_for_new_leading_promoters_starting_after( pos_1 );
    look_for_new_leading_promoters_starting_before( pos_2 );
  }
  else
  {
    int8_t dist; // Hamming distance of the sequence from the promoter consensus

    for ( int32_t i = pos_1 ; i < pos_2 ; i++ )
    {
      if ( is_promoter( LEADING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
      {
        //~ char tmp[255];
        //~ memcpy( tmp, &_dna->get_data()[i], PROM_SIZE * sizeof(char) );
        //~ printf( "new promoter found on the LEADING strand at position %"PRId32" : %s\n", i, tmp );

        // Look for the right place to insert the new promoter in the list
        ae_list_node<ae_rna*>* rna_node  = _rna_list[LEADING]->get_first();

        while( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < i )
        {
          rna_node = rna_node->get_next();
        }

        if ( rna_node == NULL )
        {
          // Add at the end of the list
          #ifndef __REGUL
            _rna_list[LEADING]->add( new ae_rna( this, LEADING, i, dist ) );
          #else
            _rna_list[LEADING]->add( new ae_rna_R( this, LEADING, i, dist ) );
          #endif
        }
        else if ( rna_node->get_obj()->get_promoter_pos() != i ) // If not already in list
        {
          // Add before rna_node
          #ifndef __REGUL
            _rna_list[LEADING]->add_before( new ae_rna( this, LEADING, i, dist ), rna_node );
          #else
            _rna_list[LEADING]->add_before( new ae_rna_R( this, LEADING, i, dist ), rna_node );
          #endif
        }
      }
    }
  }
}


/*!
  \brief Look for new promoters on the LAGGIN strand whose starting positions would lie in [pos_1 ; pos_2[

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
void ae_genetic_unit::look_for_new_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2 )
{
  assert( pos_1 >= 0 && pos_1 < _dna->get_length() && pos_2 >= 0 && pos_2 < _dna->get_length() );

  // When pos_1 > pos_2, we will perform the search in 2 steps.
  // As positions  0 and _dna->get_length() are equivalent, it's preferable to
  // keep 0 for pos_1 and _dna->get_length() for pos_2.
  //~ if ( pos_1 == _dna->get_length() ) pos_1 = 0;
  //~ if ( pos_2 == 0 )                  pos_2 = _dna->get_length();

  if ( pos_1 >= pos_2 )
  {
    look_for_new_lagging_promoters_starting_after( pos_1 );
    look_for_new_lagging_promoters_starting_before( pos_2 );
  }
  else
  {
    int8_t dist; // Hamming distance of the sequence from the promoter consensus
    ae_list_node<ae_rna*>* rna_node  = _rna_list[LAGGING]->get_first();

    for ( int32_t i = pos_2 - 1 ; i >= pos_1 ; i-- )
    {
      if ( is_promoter( LAGGING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
      {
        assert ( i >= 0 && i < _dna->get_length() );

        // Look for the right place to insert the new promoter in the list
        while( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() > i )
        {
          rna_node = rna_node->get_next();
        }

        if ( rna_node == NULL )
        {
          // Add at the end of the list
          #ifndef __REGUL
            _rna_list[LAGGING]->add( new ae_rna( this, LAGGING, i, dist ) );
          #else
            _rna_list[LAGGING]->add( new ae_rna_R( this, LAGGING, i, dist ) );
          #endif
        }
        else if ( rna_node->get_obj()->get_promoter_pos() != i ) // If not already in list
        {
          // Add before rna_node
          #ifndef __REGUL
            _rna_list[LAGGING]->add_before( new ae_rna( this, LAGGING, i, dist ), rna_node );
          #else
            _rna_list[LAGGING]->add_before( new ae_rna_R( this, LAGGING, i, dist ), rna_node );
          #endif
        }
      }
    }
  }
}


/*!
  \brief Look for new promoters on the LEADING strand whose starting positions would be >= pos
*/
void ae_genetic_unit::look_for_new_leading_promoters_starting_after( int32_t pos )
{
  assert( pos >= 0 && pos < _dna->get_length() );


  // Hamming distance of the sequence from the promoter consensus
  int8_t dist;

  // rna list node used to find the new promoter's place in the list
  ae_list_node<ae_rna*>* rna_node  = _rna_list[LEADING]->get_first();


  for ( int32_t i = pos ; i < _dna->get_length() ; i++ )
  {
    if ( is_promoter( LEADING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
    {
      //~ char tmp[255];
      //~ memcpy( tmp, &_dna->get_data()[i], PROM_SIZE * sizeof(char) );
      //~ printf( "new promoter found on the LEADING strand at position %"PRId32" : %s\n", i, tmp );

      // Look for the right place to insert the new promoter in the list
      while( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < i )
      {
        rna_node = rna_node->get_next();
      }

      if ( rna_node == NULL )
      {
        // Add at the end of the list
        #ifndef __REGUL
          _rna_list[LEADING]->add( new ae_rna( this, LEADING, i, dist ) );
        #else
          _rna_list[LEADING]->add( new ae_rna_R( this, LEADING, i, dist ) );
        #endif
      }
      else if ( rna_node->get_obj()->get_promoter_pos() != i ) // If not already in list
      {
        // Add before rna_node
        #ifndef __REGUL
          _rna_list[LEADING]->add_before( new ae_rna( this, LEADING, i, dist ), rna_node );
        #else
          _rna_list[LEADING]->add_before( new ae_rna_R( this, LEADING, i, dist ), rna_node );
        #endif
      }
    }
  }
}


/*!
  \brief Look for new promoters on the LAGGING strand whose starting positions would be >= pos

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
void ae_genetic_unit::look_for_new_lagging_promoters_starting_after( int32_t pos )
{
  assert( pos >= 0 && pos < _dna->get_length() );


  // Hamming distance of the sequence from the promoter consensus
  int8_t dist;

  // rna list node used to find the new promoter's place in the list
  ae_list_node<ae_rna*>* rna_node  = _rna_list[LAGGING]->get_first();


  for ( int32_t i = _dna->get_length() - 1 ; i >= pos ; i-- )
  {
    if ( is_promoter( LAGGING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
    {
      assert ( i >= 0 && i < _dna->get_length() );

      // Look for the right place to insert the new promoter in the list
      while( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() > i )
      {
        rna_node = rna_node->get_next();
      }

      if ( rna_node == NULL )
      {
        // Add at the end of the list
        #ifndef __REGUL
          _rna_list[LAGGING]->add( new ae_rna( this, LAGGING, i, dist ) );
        #else
          _rna_list[LAGGING]->add( new ae_rna_R( this, LAGGING, i, dist ) );
        #endif
      }
      else if ( rna_node->get_obj()->get_promoter_pos() != i ) // If not already in list
      {
        // Add before rna_node
        #ifndef __REGUL
          _rna_list[LAGGING]->add_before( new ae_rna( this, LAGGING, i, dist ), rna_node );
        #else
          _rna_list[LAGGING]->add_before( new ae_rna_R( this, LAGGING, i, dist ), rna_node );
        #endif
      }
    }
  }
}


/*!
  \brief Look for new promoters on the LEADING strand whose starting positions would be < pos
*/
void ae_genetic_unit::look_for_new_leading_promoters_starting_before( int32_t pos )
{
  assert( pos >= 0 && pos < _dna->get_length() );


  // Hamming distance of the sequence from the promoter consensus
  int8_t dist;

  // rna list node used to find the new promoter's place in the list
  ae_list_node<ae_rna*>* rna_node  = _rna_list[LEADING]->get_first();


  for ( int32_t i = 0 ; i < pos ; i++ )
  {
    if ( is_promoter( LEADING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
    {
      //~ char tmp[255];
      //~ memcpy( tmp, &_dna->get_data()[i], PROM_SIZE * sizeof(char) );
      //~ printf( "new promoter found on the LEADING strand at position %"PRId32" : %s\n", i, tmp );

      // Look for the right place to insert the new promoter in the list
      while( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < i )
      {
        rna_node = rna_node->get_next();
      }

      if ( rna_node == NULL )
      {
        // Add at the end of the list
        #ifndef __REGUL
          _rna_list[LEADING]->add( new ae_rna( this, LEADING, i, dist ) );
        #else
          _rna_list[LEADING]->add( new ae_rna_R( this, LEADING, i, dist ) );
        #endif
      }
      else if ( rna_node->get_obj()->get_promoter_pos() != i ) // If not already in list
      {
        // Add before rna_node
        #ifndef __REGUL
          _rna_list[LEADING]->add_before( new ae_rna( this, LEADING, i, dist ), rna_node );
        #else
          _rna_list[LEADING]->add_before( new ae_rna_R( this, LEADING, i, dist ), rna_node );
        #endif
      }
    }
  }
}


/*!
  \brief Look for new promoters on the LAGGING strand whose starting positions would be < pos

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
void ae_genetic_unit::look_for_new_lagging_promoters_starting_before( int32_t pos )
{
  assert( pos >= 0 && pos < _dna->get_length() );

  // Hamming distance of the sequence from the promoter consensus
  int8_t dist;

  // rna list node used to find the new promoter's place in the list
  ae_list_node<ae_rna*>* rna_node  = _rna_list[LAGGING]->get_first();


  for ( int32_t i = pos - 1 ; i >= 0 ; i-- )
  {
    if ( is_promoter( LAGGING, i, dist ) ) // dist takes the hamming distance of the sequence from the consensus
    {
      assert ( i >= 0 && i < _dna->get_length() );

      // Look for the right place to insert the new promoter in the list
      while( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() > i )
      {
        rna_node = rna_node->get_next();
      }

      if ( rna_node == NULL )
      {
        // Add at the end of the list
        #ifndef __REGUL
          _rna_list[LAGGING]->add( new ae_rna( this, LAGGING, i, dist ) );
        #else
          _rna_list[LAGGING]->add( new ae_rna_R( this, LAGGING, i, dist ) );
        #endif
      }
      else if ( rna_node->get_obj()->get_promoter_pos() != i ) // If not already in list
      {
        // Add before rna_node
        #ifndef __REGUL
          _rna_list[LAGGING]->add_before( new ae_rna( this, LAGGING, i, dist ), rna_node );
        #else
          _rna_list[LAGGING]->add_before( new ae_rna_R( this, LAGGING, i, dist ), rna_node );
        #endif
      }
    }
  }
}


/*!
  \brief Shift (by delta_post) the positions of the promoters from the LEADING strand whose starting positions are >= pos.
*/
void ae_genetic_unit::move_all_leading_promoters_after( int32_t pos, int32_t delta_pos )
{
  ae_list_node<ae_rna*>* rna_node  = _rna_list[LEADING]->get_first();

  // Go to first RNA after pos
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos )
  {
    //~ printf( "don't move leading promoter at position %"PRId32"\n", rna_node->get_obj()->get_promoter_pos() );
    rna_node = rna_node->get_next();
  }

  // Update all the remaining RNAs
  while ( rna_node != NULL )
  {
    //~ printf( "move leading promoter at position %"PRId32"\n", rna_node->get_obj()->get_promoter_pos() );
    rna_node->get_obj()->shift_position( delta_pos, _dna->get_length() );
    //~ printf( "new position : %"PRId32"\n", rna_node->get_obj()->get_promoter_pos() );

    rna_node = rna_node->get_next();
  }
}


/*!
  \brief Shift (by delta_post) the positions of the promoters from the LAGGING strand whose starting positions are >= pos.

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
void ae_genetic_unit::move_all_lagging_promoters_after( int32_t pos, int32_t delta_pos )
{
  ae_list_node<ae_rna*>* rna_node  = _rna_list[LAGGING]->get_first();
  ae_rna*       rna;

  // Update RNAs until we pass pos (or we reach the end of the list)
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() >= pos )
  {
    rna = rna_node->get_obj();
    if ( rna->get_promoter_pos() < pos ) break;

    rna->shift_position( delta_pos, _dna->get_length() );

    rna_node = rna_node->get_next();
  }
}


/*!
  \brief Copy (into new_promoter_list) the promoters from the LEADING strand whose starting positions lie in [pos_1 ; pos_2[
*/
void ae_genetic_unit::copy_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>* new_promoter_list )
{
  // 1) Go to first RNA to copy
  ae_list_node<ae_rna*>* rna_node  = _rna_list[LEADING]->get_first();
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_1 )
  {
    //~ printf( "don't move leading promoter at position %"PRId32"\n", rna_node->get_obj()->get_promoter_pos() );
    rna_node = rna_node->get_next();
  }

  // 2) Copy RNAs
  if ( pos_1 < pos_2 )
  {
    // Copy from pos_1 to pos_2
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_2 )
    {
      new_promoter_list->add( new ae_rna( this, *(rna_node->get_obj()) ) );

      rna_node = rna_node->get_next();
    }
  }
  else
  {
    // Copy from pos_1 to the end of the list
    while ( rna_node != NULL )
    {
      new_promoter_list->add( new ae_rna( this, *(rna_node->get_obj()) ) );

      rna_node = rna_node->get_next();
    }

    // Copy from the beginning of the list to pos_2
    rna_node  = _rna_list[LEADING]->get_first();
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_2 )
    {
      new_promoter_list->add( new ae_rna( this, *(rna_node->get_obj()) ) );

      rna_node = rna_node->get_next();
    }
  }
}


/*!
  \brief Copy (into new_promoter_list) the promoters from the LAGGING strand whose starting positions lie in [pos_1 ; pos_2[

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
void ae_genetic_unit::copy_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>* new_promoter_list )
{
  // Go to first RNA to copy
  ae_list_node<ae_rna*>* rna_node  = _rna_list[LAGGING]->get_last();
  while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_1 )
  {
    rna_node  = rna_node->get_prev();
  }

  // Copy RNAs
  if ( pos_1 < pos_2 )
  {
    // Copy from pos_1 to pos_2
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_2 )
    {
      new_promoter_list->add_front( new ae_rna( this, *(rna_node->get_obj()) ) );

      rna_node = rna_node->get_prev();
    }
  }
  else
  {
    // Copy from pos_1 to the beginning of the list (we are going backwards)
    while ( rna_node != NULL )
    {
      new_promoter_list->add_front( new ae_rna( this, *(rna_node->get_obj()) ) );

      rna_node = rna_node->get_prev();
    }

    // Copy from the end of the list to pos_2 (we are going backwards)
    rna_node  = _rna_list[LAGGING]->get_last();
    while ( rna_node != NULL && rna_node->get_obj()->get_promoter_pos() < pos_2 )
    {
      new_promoter_list->add_front( new ae_rna( this, *(rna_node->get_obj()) ) );

      rna_node = rna_node->get_prev();
    }
  }
}

void ae_genetic_unit::save( gzFile backup_file ) const
{
  _dna->save( backup_file );
  gzwrite( backup_file, &_min_gu_length, sizeof(_min_gu_length) );
  gzwrite( backup_file, &_max_gu_length, sizeof(_max_gu_length) );
}

int32_t ae_genetic_unit::get_nb_terminators( void )
{
  int32_t nb_term = 0;

  if ( _dna->get_length() >= TERM_SIZE )
  {
    for ( int32_t i = 0 ; i < _dna->get_length() ; i++ )
    {
      if ( is_terminator( LEADING, i ) )  // No need to count on both the LEADING and the LAGGING strand
                                          // as terminators are "shared"
      {
        nb_term++;
      }
    }
  }

  return nb_term;
}

#ifdef DEBUG
  void ae_genetic_unit::assert_promoters( void )
  {
    // Check that the lists are ordered correctly
    assert_promoters_order();

    // Make a backup of the genetic unit's lists of RNAs
    ae_list<ae_rna*>** old_rna_list = _rna_list;

    _rna_list           = new ae_list<ae_rna*>*[2];
    _rna_list[LEADING]  = new ae_list<ae_rna*>();
    _rna_list[LAGGING]  = new ae_list<ae_rna*>();

    locate_promoters();

    // Compare lists
    ae_list_node<ae_rna*>*  node_old  = NULL;
    ae_list_node<ae_rna*>*  node_new  = NULL;
    ae_rna* rna_old   = NULL;
    ae_rna* rna_new   = NULL;

    for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ )
    {
      node_old  = old_rna_list[strand]->get_first();
      node_new  = _rna_list[strand]->get_first();

      while ( node_old != NULL || node_new != NULL )
      {
        if( node_old == NULL || node_new == NULL )
        {
          printf( "****************************** NB_ELTS problem ******************************\n" );
          printf( "should be : \n" );
          print_rnas( _rna_list );
          printf( "is : \n" );
          print_rnas( old_rna_list );
          printf( "****************************************************************************\n" );
          printf( "  genome length : %" PRId32 "\n", _dna->get_length() );
          assert( node_old != NULL && node_new != NULL );
        }

        rna_old = node_old->get_obj();
        rna_new = node_new->get_obj();

        if ( rna_old->get_strand() != rna_new->get_strand() )
        {
          printf( "****************************** STRAND problem ******************************\n" );
          printf( "should be : \n" );
          print_rnas( _rna_list );
          printf( "is : \n" );
          print_rnas( old_rna_list );
          printf( "****************************************************************************\n" );
          printf( "  %" PRId32 " (%s) : %f    vs    %" PRId32 " (%s) : %f\n",
                  rna_old->get_promoter_pos(), rna_old->get_strand() == LEADING ? "LEADING" : "LAGGING", rna_old->get_basal_level(),
                  rna_new->get_promoter_pos(), rna_new->get_strand() == LEADING ? "LEADING" : "LAGGING", rna_new->get_basal_level() );
          printf( "  genome length : %" PRId32 "\n", _dna->get_length() );
          assert( rna_old->get_strand() == rna_new->get_strand() );
        }

        if ( rna_old->get_promoter_pos() != rna_new->get_promoter_pos() )
        {
          printf( "***************************** POSITION problem *****************************\n" );
          printf( "should be : \n" );
          print_rnas( _rna_list );
          printf( "is : \n" );
          print_rnas( old_rna_list );
          printf( "****************************************************************************\n" );
          printf( "  %" PRId32 " (%s) : %f    vs    %" PRId32 " (%s) : %f\n",
                  rna_old->get_promoter_pos(), rna_old->get_strand() == LEADING ? "LEADING" : "LAGGING", rna_old->get_basal_level(),
                  rna_new->get_promoter_pos(), rna_new->get_strand() == LEADING ? "LEADING" : "LAGGING", rna_new->get_basal_level() );
          printf( "  genome length : %" PRId32 "\n", _dna->get_length() );
          assert( rna_old->get_promoter_pos() == rna_new->get_promoter_pos()  );
        }

        if ( rna_old->get_basal_level() != rna_new->get_basal_level() )
        {
          printf( "*************************** BASAL LEVEL problem ****************************\n" );
          printf( "should be : \n" );
          print_rnas( _rna_list );
          printf( "is : \n" );
          print_rnas( old_rna_list );
          printf( "****************************************************************************\n" );
          printf( "  %" PRId32 " (%s) : %f    vs    %" PRId32 " (%s) : %f\n",
                  rna_old->get_promoter_pos(), rna_old->get_strand() == LEADING ? "LEADING" : "LAGGING", rna_old->get_basal_level(),
                  rna_new->get_promoter_pos(), rna_new->get_strand() == LEADING ? "LEADING" : "LAGGING", rna_new->get_basal_level() );
          printf( "  genome length : %" PRId32 "\n", _dna->get_length() );
          assert( rna_old->get_basal_level() == rna_new->get_basal_level() );
        }

        node_old = node_old->get_next();
        node_new = node_new->get_next();
      }
    }

    _rna_list[LEADING]->erase( true );
    _rna_list[LAGGING]->erase( true );
    delete _rna_list[LEADING];
    delete _rna_list[LAGGING];
    delete [] _rna_list;

    _rna_list = old_rna_list;
  }

  void ae_genetic_unit::assert_promoters_order( void )
  {
    ae_list_node<ae_rna*>* node1 = NULL;
    ae_list_node<ae_rna*>* node2 = NULL;
    ae_rna* rna1 = NULL;
    ae_rna* rna2 = NULL;

    for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ )
    {
      if ( _rna_list[strand]->get_nb_elts() >= 2 )
      {
        node1 = _rna_list[strand]->get_first();
        node2 = node1->get_next();

        while ( node2 != NULL )
        {
          rna1 = node1->get_obj();
          rna2 = node2->get_obj();

          if ( strand == LEADING )
          {
            if( rna1->get_promoter_pos() >= rna2->get_promoter_pos() )
            {
              printf( "********************** ORDER problem (LEADING) ***********************\n" );
              print_rnas();
              printf( "****************************************************************************\n" );
              assert( rna1->get_promoter_pos() < rna2->get_promoter_pos() );
            }
          }
          else
          {
            if( rna1->get_promoter_pos() <= rna2->get_promoter_pos() )
            {
              printf( "*********************** ORDER problem (LAGGING) ***********************\n" );
              print_rnas();
              printf( "****************************************************************************\n" );
              assert( rna1->get_promoter_pos() > rna2->get_promoter_pos() );
            }
          }

          node1 = node2;
          node2 = node2->get_next();
        }
      }
    }
  }
#endif

/*!
  \brief Retrieve for each base if it belongs or not to coding RNA

  \return Boolean table of sequence length size with for each base if it belongs or not to coding RNA
*/
bool* ae_genetic_unit::is_belonging_to_coding_RNA( void )
{
  int32_t genome_length = _dna->get_length();
  bool* belongs_to_coding_RNA = new bool[genome_length];
  memset( belongs_to_coding_RNA,0, genome_length );

  // Parse RNA lists and mark the corresponding bases as coding (only for the coding RNAs)
  for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ )
  {
    ae_list_node<ae_rna*>* rna_node  = _rna_list[strand]->get_first();
    ae_rna* rna = NULL;
    while ( rna_node != NULL )
    {
      rna = rna_node->get_obj();

      int32_t first;
      int32_t last;

      if ( strand == LEADING )
      {
        first = rna->get_promoter_pos();
        last  = rna->get_last_transcribed_pos();
      }
      else // ( strand == LAGGING )
      {
        last  = rna->get_promoter_pos();
        first = rna->get_last_transcribed_pos();
      }

      if (not rna->get_transcribed_proteins().empty()) // coding RNA
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
      rna_node = rna_node->get_next();
    }
  }
  return(belongs_to_coding_RNA);
}

/*!
  \brief Remove the bases that are not in coding RNA

  Remove the bases that are not in coding RNA and test at each loss that fitness is not changed
*/
void ae_genetic_unit::remove_non_coding_bases( void)
{
  Environment* env = _exp_m->get_env() ;

  reset_expression();
  locate_promoters();
  _distance_to_target_computed        = false;
  _fitness_computed                   = false;
  compute_phenotypic_contribution();
  if(_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
  {
    compute_distance_to_target(env);
  }
  compute_fitness(env);
  double initial_fitness = get_fitness();

  int32_t genome_length = _dna->get_length();
  bool* belongs_to_coding_RNA = is_belonging_to_coding_RNA();

  int32_t non_coding_bases_nb = get_nb_bases_in_0_coding_RNA();
  int32_t start = 0;
  int32_t end = 0;


  for ( int32_t i = genome_length-1 ; i >= 0 ; i-- )
  {
    if ( belongs_to_coding_RNA[i] == false )
    {
      end = i+1;
      while((i-1) >=0 && belongs_to_coding_RNA[(i-1)] == false )
      {
        i--;
      }
      start = i;
      _dna->remove(start,end);

      locate_promoters();
      reset_expression();
      _distance_to_target_computed        = false;
      _fitness_computed                   = false;
      compute_phenotypic_contribution();
      if(_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
      {
        compute_distance_to_target(env);
      }
      compute_fitness(env);
      assert(get_fitness()==initial_fitness);

      _non_coding_computed = false;
      non_coding_bases_nb = get_nb_bases_in_0_coding_RNA();
    }
  }

  locate_promoters();
  reset_expression();
  _distance_to_target_computed        = false;
  _fitness_computed                   = false;
  compute_phenotypic_contribution();
  if(_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
  {
    compute_distance_to_target(env);
  }
  compute_fitness(env);
  assert(get_fitness()==initial_fitness);

  _non_coding_computed = false;
  non_coding_bases_nb = get_nb_bases_in_0_coding_RNA();
  assert(non_coding_bases_nb==0);

  delete [] belongs_to_coding_RNA;
}

/*!
  \brief Double the bases that are not in coding RNA

  Double the bases that are not in coding RNA by addition of random bases and test at each addition that fitness is not changed
*/
void ae_genetic_unit::double_non_coding_bases(void)
{
  Environment* env = _exp_m->get_env() ;

  reset_expression();
  locate_promoters();
  _distance_to_target_computed        = false;
  _fitness_computed                   = false;
  compute_phenotypic_contribution();
  if(_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
  {
    compute_distance_to_target(env);
  }
  compute_fitness(env);
  double initial_fitness = get_fitness();

  int32_t genome_length = _dna->get_length();
  bool* belongs_to_coding_RNA = is_belonging_to_coding_RNA();

  _non_coding_computed = false;
  int32_t inital_non_coding_bases_nb = get_nb_bases_in_0_coding_RNA();
  int32_t start = 0;
  int32_t end = 0;
  int32_t length = 0;
  int32_t pos = 0;
  char * random_portion = NULL;
  bool insertion_ok = false;
  int32_t non_coding_bases_nb_before_fitness = get_nb_bases_in_0_coding_RNA();
  int32_t non_coding_bases_nb_after_fitness = get_nb_bases_in_0_coding_RNA();

  int32_t non_coding_bases_nb = get_nb_bases_in_0_coding_RNA();

  for ( int32_t i = genome_length-1 ; i >= 0 ; i-- )
  {
    if ( belongs_to_coding_RNA[i] == false )
    {
      end = i+1;
      while((i-1) >=0 && belongs_to_coding_RNA[(i-1)] == false )
      {
        i--;
      }
      start = i;
      length = end-start;

      insertion_ok = false;
      while(not insertion_ok)
      {
        random_portion = new char [length+1];
        for ( int32_t j = 0 ; j < length ; j++ )
        {
          random_portion[j] = '0' + _indiv->get_mut_prng()->random( NB_BASE );
        }
        random_portion[length] = 0;
        pos = _indiv->get_mut_prng()->random( length )+start;
        _dna->insert(pos, random_portion);

        _non_coding_computed = false;
        non_coding_bases_nb_before_fitness = get_nb_bases_in_0_coding_RNA();

        locate_promoters();
        reset_expression();
        _distance_to_target_computed        = false;
        _fitness_computed                   = false;
        compute_phenotypic_contribution();
        if(_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
        {
          compute_distance_to_target(env);
        }
        compute_fitness(env);
        assert(get_fitness()==initial_fitness);

        _non_coding_computed = false;
        non_coding_bases_nb_after_fitness = get_nb_bases_in_0_coding_RNA();

        if (non_coding_bases_nb_before_fitness != non_coding_bases_nb_after_fitness)
        {
          _dna->remove(pos, pos + length);
        }
        else
        {
          insertion_ok = true;
        }

      }
      _non_coding_computed = false;

      delete [] random_portion;
      random_portion = NULL;
    }
  }

  locate_promoters();
  reset_expression();
  _distance_to_target_computed        = false;
  _fitness_computed                   = false;
  compute_phenotypic_contribution();
  if(_exp_m->get_output_m()->get_compute_phen_contrib_by_GU())
  {
    compute_distance_to_target(env);
  }
  compute_fitness(env);
  assert(get_fitness()==initial_fitness);

  _non_coding_computed = false;
  non_coding_bases_nb = get_nb_bases_in_0_coding_RNA();
  assert(non_coding_bases_nb == 2*inital_non_coding_bases_nb);

  delete [] belongs_to_coding_RNA;
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
void ae_genetic_unit::compute_nb_of_affected_genes(const ae_mutation * mut, int & nb_genes_at_breakpoints, int & nb_genes_in_segment, int & nb_genes_in_replaced_segment)
{
  nb_genes_at_breakpoints = 0;
  nb_genes_in_segment = 0;
  nb_genes_in_replaced_segment = 0;
  int32_t genlen = get_seq_length(); // length of the genetic unit, in bp

  int32_t pos0 = -1, pos1 = -1, pos2 = -1, pos3 = -1, mutlength = -1;
  int32_t pos1donor = -1, pos2donor = -1, pos3donor = -1;
  char * seq = NULL;
  int32_t first, last;
  bool invert = false;
  ae_sense sense = DIRECT;
  MutationType type = mut->get_mut_type();
  switch(type)
    {
    case SWITCH:
      {
        mut->get_infos_point_mutation(&pos0);
        mutlength = 1;
        break;
      }
    case S_INS:
      {
        mut->get_infos_small_insertion(&pos0, &mutlength);
        break;
      }
    case S_DEL:
      {
        mut->get_infos_small_deletion(&pos0, &mutlength);
        break;
      }
    case DUPL:
      {
        mut->get_infos_duplication(&pos1, &pos2, &pos0);
        // pos2 is actually not included in the segment, the real end of the segment is pos2 - 1
        pos2 = ae_utils::mod(pos2 - 1, genlen);
        mutlength = mut->get_length();
        break;
      }
    case DEL:
      {
        mut->get_infos_deletion(&pos1, &pos2);
        pos2 = ae_utils::mod(pos2 - 1, genlen);
        mutlength = mut->get_length();
        break;
      }
    case TRANS:
      {
        mut->get_infos_translocation(&pos1, &pos2, &pos3, &pos0, &invert);
        pos2 = ae_utils::mod(pos2 - 1, genlen);
        mutlength = mut->get_length();
        break;
      }
    case INV:
      {
        mut->get_infos_inversion(&pos1, &pos2);
        pos2 = ae_utils::mod(pos2 - 1, genlen);
        mutlength = mut->get_length();
        break;
      }
    case INSERT:
      {
        mut->get_infos_insertion(&pos0, &mutlength);
        seq = new char[mutlength+1];
        mut->get_sequence_insertion(seq);
        // Make a temporary genetic unit and translate it to count how many genes are on the inserted segment
        ae_genetic_unit * tmpunit = new ae_genetic_unit( _indiv, seq, mutlength);
        tmpunit->do_transcription();
        tmpunit->do_translation();
        nb_genes_in_segment = tmpunit->get_nb_coding_RNAs();
        delete tmpunit;
        seq = NULL;
        break;
      }
    case INS_HT:
      {
        mut->get_infos_ins_HT(&pos1donor, &pos2donor, &pos0, &pos3donor, &sense, &mutlength );
        seq = new char[mutlength+1];
        mut->get_sequence_ins_HT(seq);

        // Make a temporary genetic unit and translate it to count how many genes were on the exogenote
        ae_genetic_unit * tmpunit = new ae_genetic_unit( _indiv, seq, mutlength);
        tmpunit->do_transcription();
        tmpunit->do_translation();
        nb_genes_in_segment = tmpunit->get_nb_coding_RNAs();

        // Check whether the pos3donor breakpoint killed one or several of those genes, in that case decrement nb_genes_in_segment
        ae_list_node<ae_rna*>* a_node = tmpunit->_rna_list[LEADING]->get_first();
        ae_rna * rna = NULL;
        while (a_node != NULL)
          {
            rna = (ae_rna *) a_node->get_obj();
            a_node = a_node->get_next();
            if (rna->is_coding() == false) continue;

            first = rna->get_promoter_pos();
            last = rna->get_last_transcribed_pos();

            if (breakpoint_inside_gene(pos3donor, first, last)) nb_genes_in_segment --;
          }
        a_node = tmpunit->_rna_list[LAGGING]->get_first();
        rna = NULL;
        while (a_node != NULL)
          {
            rna = (ae_rna *) a_node->get_obj();
            a_node = a_node->get_next();
            if (rna->is_coding() == false) continue;

            last = rna->get_promoter_pos();
            first = rna->get_last_transcribed_pos();

            if (breakpoint_inside_gene(pos3donor, first, last)) nb_genes_in_segment --;
          }
        delete tmpunit;
        seq = NULL;
        break;
      }
    case REPL_HT:
      {
        mut->get_infos_repl_HT(&pos1, &pos1donor, &pos2, &pos2donor, &sense, &mutlength );
        pos2 = ae_utils::mod(pos2 - 1, genlen);
        seq = new char[mutlength+1];  // seq in the sequence in the donor
        mut->get_sequence_repl_HT(seq);

        // Make a temporary genetic unit and translate it to count how many genes were on the donor segment
        ae_genetic_unit * tmpunit = new ae_genetic_unit( _indiv, seq, mutlength);
        tmpunit->do_transcription();
        tmpunit->do_translation();
        nb_genes_in_segment = tmpunit->get_nb_coding_RNAs();

        // Remove the genes that overlap ori in this temporary unit, as the transferred segment was actually linear
        ae_list_node<ae_rna*>* a_node = tmpunit->_rna_list[LEADING]->get_first();
        ae_rna * rna = NULL;
        while (a_node != NULL)
          {
            rna = (ae_rna *) a_node->get_obj();
            a_node = a_node->get_next();
            if (rna->is_coding() == false) continue;

            first = rna->get_promoter_pos();
            last = rna->get_last_transcribed_pos();
            if (first > last) nb_genes_in_segment --;
          }
        a_node = tmpunit->_rna_list[LAGGING]->get_first();
        rna = NULL;
        while (a_node != NULL)
          {
            rna = (ae_rna *) a_node->get_obj();
            a_node = a_node->get_next();
            if (rna->is_coding() == false) continue;

            first = rna->get_last_transcribed_pos();
            last = rna->get_promoter_pos();
            if (first > last) nb_genes_in_segment --;
          }

        delete tmpunit;
        seq = NULL;
        break;
      }
    default:
      {
        fprintf(stderr, "Error: unknown mutation type in ae_genetic_unit::compute_nb_of_affected_genes.\n");
      }
    }

  ae_list_node<ae_rna*>* a_node = _rna_list[LEADING]->get_first();
  ae_rna * rna = NULL;
  while (a_node != NULL)
    {
      rna = (ae_rna *) a_node->get_obj();
      a_node = a_node->get_next();
      if (rna->is_coding() == false) continue;

      first = rna->get_promoter_pos();
      last = rna->get_last_transcribed_pos();

      switch(type)
        {
        case SWITCH:
          {
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case S_INS:
          {
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case S_DEL:
          {
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case DUPL:
          {
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_segment ++;
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case DEL:
          {
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_segment ++;
            if (breakpoint_inside_gene(pos1, first, last)) nb_genes_at_breakpoints ++;
            else if (breakpoint_inside_gene(pos2, first, last)) nb_genes_at_breakpoints ++;     // else because the gene must not be counted twice if both p1 and p2 are in the same gene
            break;
          }
        case TRANS:
          {
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_segment ++;
            if (breakpoint_inside_gene(pos1, first, last)) nb_genes_at_breakpoints ++;   // beginning of the excised segment
            else if (breakpoint_inside_gene(pos2, first, last)) nb_genes_at_breakpoints ++;   // end of the excised segment
            else if (breakpoint_inside_gene(pos3, first, last)) nb_genes_at_breakpoints ++;   // breakpoint inside the segment for the reinsertion
            else if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;   // reinsertion point in the genetic unit
            break;
          }
        case INV:
          {
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_segment ++;
            if (breakpoint_inside_gene(pos1, first, last)) nb_genes_at_breakpoints ++;
            else if (breakpoint_inside_gene(pos2, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case INSERT:
          {
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case INS_HT:
          {
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case REPL_HT:
          {
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_replaced_segment ++; // the whole gene sequence was replaced by the donor DNA
            if (breakpoint_inside_gene(pos1, first, last)) nb_genes_at_breakpoints ++;  // the gene was disrupted by the breakpoint p1
            else if (breakpoint_inside_gene(pos2, first, last)) nb_genes_at_breakpoints ++;  // the gene was disrupted by the breakpoint p2
            break;
          }
        default:
          // Only simple mutation types are considered.
          break;
        }

    }


  a_node = _rna_list[LAGGING]->get_first();
  while (a_node != NULL)
    {
      rna = (ae_rna *) a_node->get_obj();
      a_node = a_node->get_next();
      if (rna->is_coding() == false) continue;

      first = rna->get_last_transcribed_pos();
      last = rna->get_promoter_pos();

     switch(type)
        {
        case SWITCH:
          {
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case S_INS:
          {
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case S_DEL:
          {
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case DUPL:
          {
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_segment ++;
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case DEL:
          {
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_segment ++;
            if (breakpoint_inside_gene(pos1, first, last)) nb_genes_at_breakpoints ++;
            else if (breakpoint_inside_gene(pos2, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case TRANS:
          {
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_segment ++;
            if (breakpoint_inside_gene(pos1, first, last)) nb_genes_at_breakpoints ++;   // beginning of the excised segment
            else if (breakpoint_inside_gene(pos2, first, last)) nb_genes_at_breakpoints ++;   // end of the excised segment
            else if (breakpoint_inside_gene(pos3, first, last)) nb_genes_at_breakpoints ++;   // breakpoint inside the segment for the reinsertion
            else if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;   // reinsertion point in the genetic unit
            break;
          }
        case INV:
          {
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_segment ++;
            if (breakpoint_inside_gene(pos1, first, last)) nb_genes_at_breakpoints ++;
            else if (breakpoint_inside_gene(pos2, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case INSERT:
          {
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case INS_HT:
          {
            if (breakpoint_inside_gene(pos0, first, last)) nb_genes_at_breakpoints ++;
            break;
          }
        case REPL_HT:
          {
            if (gene_totally_in_segment(pos1, pos2, first, last)) nb_genes_in_replaced_segment ++; // the whole gene sequence was replaced by the donor DNA
            if (breakpoint_inside_gene(pos1, first, last)) nb_genes_at_breakpoints ++;  // the gene was disrupted by the breakpoint p1
            else if (breakpoint_inside_gene(pos2, first, last)) nb_genes_at_breakpoints ++;  // the gene was disrupted by the breakpoint p2
            break;
          }
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
void ae_genetic_unit::init_statistical_data( void ) // TODO : integrate into compute_statistical_data
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
