/**
 * @file    ae_individual.cpp
 * @section AUTHORS
 *          Guillaume Beslon, Carole Knibbe, Virginie Lefort, David Parsons
 *          <david.parsons@liris.cnrs.fr>
 *
 *          LIRIS lab : http://liris.cnrs.fr/
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * @section DESCRIPTION
 *
 * ...
 */
 
 
/** @class ae_individual
 *
 * Brief description
 *
 * Long description
 * to be written
 */

// =================================================================
//                              Libraries
// =================================================================
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <algorithm>


// =================================================================
//                            Project Files
// =================================================================
#include <ae_codon.h>
#include <ae_exp_setup.h>
#include <ae_exp_manager.h>
#include <ae_grid_cell.h>
#include <ae_genetic_unit.h>
#include <ae_population.h>
#include <ae_vis_a_vis.h>
#include <ae_utils.h>

#ifdef __NO_X
  #ifdef __REGUL
    #include <ae_individual_R.h>
  #else
    #include <ae_individual.h>
  #endif
#elif defined __X11
  #ifdef __REGUL
    #include <ae_individual_R_X11.h>
  #else
    #include <ae_individual_X11.h>
  #endif
#endif

//##############################################################################
//                                                                             #
//                             Class ae_individual                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
/*!

*/
ae_individual::ae_individual( ae_exp_manager* exp_m,
                              ae_jumping_mt* mut_prng,
                              ae_jumping_mt* stoch_prng,
                              ae_params_mut* param_mut,
                              double w_max,
                              int32_t min_genome_length,
                              int32_t max_genome_length,
                              bool allow_plasmids,
                              int32_t plasmid_minimal_length,
                              int32_t id,
                              int32_t age )
{
  // Experiment manager
  _exp_m = exp_m;
  
  // PRNGs
  _mut_prng   = mut_prng;
  _stoch_prng = stoch_prng;
  
  // Replication Report
  _replic_report = NULL;
  
  // ID and rank of the indiv ; "age" of the strain
  set_id( id );
  _rank = -1; // TODO: UNRANKED
  _age = age;
  
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
  _genetic_unit_list = new ae_list<ae_genetic_unit*>();
  
  _protein_list = new ae_list<ae_protein*>();
  _rna_list     = new ae_list<ae_rna*>();
  
  // Generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for ( int8_t i = 0 ; i < 5 ; i++ )
  {
    _int_probes[i]    = 0;
    _double_probes[i] = 0;
  }
  
  // Mutation rates etc...
  _mut_params = new ae_params_mut( *param_mut );
  
  // Artificial chemistry
  _w_max = w_max;
      
  // Genome size constraints
  _min_genome_length = min_genome_length;
  _max_genome_length = max_genome_length;
  
  // Plasmids settings
  _allow_plasmids         = allow_plasmids;
  _plasmid_minimal_length = plasmid_minimal_length;
  
  
  
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
ae_individual::ae_individual( ae_exp_manager* exp_m, gzFile backup_file )
{
  _exp_m = exp_m;
  
  // Retrieve the "age" of the strain
  gzread( backup_file, &_age, sizeof(_age) );
  
  // Retrieve the PRNGs
  #ifdef DISTRIBUTED_PRNG
    _mut_prng   = new ae_jumping_mt( backup_file );
    _stoch_prng = new ae_jumping_mt( backup_file );
  #else
    _mut_prng   = exp_m->get_pop()->get_mut_prng();
    _stoch_prng = exp_m->get_pop()->get_stoch_prng();
  #endif
  
  // Retreive id and rank
  gzread( backup_file, &_id,    sizeof(_id) );
  gzread( backup_file, &_rank,  sizeof(_rank) );
  
  // Retrieve spatial coordinates
  int16_t x;
  int16_t y;
  gzread( backup_file, &x, sizeof(x) );
  gzread( backup_file, &y, sizeof(y) );
  if ( _exp_m->is_spatially_structured() )
  {
    set_grid_cell( _exp_m->get_grid_cell( x, y ) );
  }
  
  // Retreive generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  gzread( backup_file, _int_probes,     5 * sizeof(*_int_probes) );
  gzread( backup_file, _double_probes,  5 * sizeof(*_double_probes) );
  
  // Retrieve mutational parameters
  _mut_params = new ae_params_mut( backup_file );
  
  // ------------------------------------------------- Phenotypic stochasticity
  gzread( backup_file, &_with_stochasticity, sizeof(_with_stochasticity) );
  
  // Retrieve artificial chemistry parameters
  gzread( backup_file, &_w_max, sizeof(_w_max) );
      
  // Retrieve genome size constraints
  gzread( backup_file, &_min_genome_length, sizeof(_min_genome_length) );
  gzread( backup_file, &_max_genome_length, sizeof(_max_genome_length) );
  
  // Retrieve plasmids settings
  int8_t tmp_allow_plasmids;
  gzread( backup_file, &tmp_allow_plasmids, sizeof(tmp_allow_plasmids) );
  _allow_plasmids = tmp_allow_plasmids ? 1 : 0;
  if ( _allow_plasmids )
  {
    gzread( backup_file, &_plasmid_minimal_length, sizeof(_plasmid_minimal_length) );
  }
  
  // Retreive genetic units
  _genetic_unit_list = new ae_list<ae_genetic_unit*>();
  int16_t nb_gen_units;
  gzread( backup_file, &nb_gen_units,  sizeof(nb_gen_units) );
  
  for ( int16_t i = 0 ; i < nb_gen_units ; i++ )
  {
    _genetic_unit_list->add( new ae_genetic_unit( this, backup_file ) );
  }
  
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
  
  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }


  // Replication report, protein list and rna list initialization
  _replic_report = NULL;  // We are reloading from a backup, the replication report for the loaded generation 
                          // is already in the tree file corresponding to this generation and needs not be re-written.
                          // NB : If the replication report is needed in future development, it will have to be
                          // loaded from the tree file.
  
  _protein_list = new ae_list<ae_protein*>();
  _rna_list     = new ae_list<ae_rna*>();
  
  
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
  
  _placed_in_population = false;
  
  
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
ae_individual::ae_individual( const ae_individual &model )
{
  _exp_m = model._exp_m;
  
  // PRNGs
  #ifdef DISTRIBUTED_PRNG
    _mut_prng   = new ae_jumping_mt( *(model._mut_prng) );
    _stoch_prng = new ae_jumping_mt( *(model._stoch_prng) );
  #else
    _mut_prng   = model._mut_prng;
    _stoch_prng = model._stoch_prng;
  #endif
  
  _age  = model._age;
  
  _id   = model._id;
  _rank = model._rank;
  
  _evaluated                    = false;//model._evaluated;
  _transcribed                  = false;//model._transcribed;
  _translated                   = false;//model._translated;
  _folded                       = false;//model._folded;
  _phenotype_computed           = model._phenotype_computed;
  
  // Artificial chemistry parameters
  _w_max = model._w_max;
  
  // The distance to target and what results from it depend on the environment
  // and must hence be recomputed with the (possibly different) environment.
  _distance_to_target_computed  = false;
  _fitness_computed             = false;
  _statistical_data_computed    = model._statistical_data_computed;
  
  _non_coding_computed          = model._non_coding_computed;
  _modularity_computed          = model._modularity_computed;
  
  _placed_in_population         = false;
  _grid_cell = NULL;
  
  // Create an empty list of genetic units
  _genetic_unit_list = new ae_list<ae_genetic_unit*>();
  
  
  // Copy model's genetic units
  ae_list_node<ae_genetic_unit*>* gen_unit_node = model._genetic_unit_list->get_first();
  while ( gen_unit_node != NULL )
  {
    _genetic_unit_list->add( new ae_genetic_unit( this, *(gen_unit_node->get_obj()) ) );
    
    gen_unit_node = gen_unit_node->get_next();
  }
  
  
  // Copy phenotype
  if ( _phenotype_computed )
  {
    _phenotype_activ  = new ae_fuzzy_set( *(model._phenotype_activ) );
    _phenotype_inhib  = new ae_fuzzy_set( *(model._phenotype_inhib) );
    _phenotype        = new ae_phenotype( this, *(model._phenotype) );
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
  
  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
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
  
  
  // We don't copy the generation report since we are creating a clone
  // We could create a new (empty) replic report but for now, it is not needed
  _replic_report = NULL;
  
  _protein_list = new ae_list<ae_protein*>();
  _rna_list     = new ae_list<ae_rna*>();
  
  // Generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for ( int8_t i = 0 ; i < 5 ; i++ )
  {
    _int_probes[i]    = model._int_probes[i];
    _double_probes[i] = model._double_probes[i];
  }
  
  // Mutation rates etc...
  _mut_params = new ae_params_mut( *(model._mut_params) );
  
      
  // Genome size constraints
  _min_genome_length = model._min_genome_length;
  _max_genome_length = model._max_genome_length;
  
  // Plasmids settings
  _allow_plasmids         = model._allow_plasmids;
  _plasmid_minimal_length = model._plasmid_minimal_length;
  
  evaluate();
}

/*!
  This constructor creates a new individual with the same genome as it's parent.
  The location of promoters will be copied but no further process will be performed.

  The phenotype and the fitness are not set, neither is the statistical data.
*/
ae_individual::ae_individual( ae_individual* const parent, int32_t id, ae_jumping_mt* mut_prng, ae_jumping_mt* stoch_prng )
{
  _exp_m = parent->_exp_m;
  
  // PRNGs
  _mut_prng   = mut_prng;
  _stoch_prng = stoch_prng;
  
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
  
  // Artificial chemistry parameters
  _w_max = parent->_w_max;
  
  // Create new genetic units with their DNA copied from here
  // NOTE : The RNA lists (one per genetic unit) will also be copied so that we don't
  // need to look for promoters on the whole genome
  _genetic_unit_list = new ae_list<ae_genetic_unit*>();
  ae_list_node<ae_genetic_unit*>* gen_unit_node = parent->_genetic_unit_list->get_first();
  while ( gen_unit_node != NULL )
  {
    _genetic_unit_list->add( new ae_genetic_unit( this, gen_unit_node->get_obj() ) );
    
    gen_unit_node = gen_unit_node->get_next();
  }
  
  _phenotype_activ  = NULL;
  _phenotype_inhib  = NULL;
  _phenotype        = NULL;
  
  // Initialize all the fitness-related stuff
  _dist_to_target_by_segment  = NULL;
  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];
  
  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }

  // Create a new replication report to store mutational events
  if ( _exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL )
  {
    _replic_report = new ae_replication_report( this, parent );
    
    // TODO: remove this after checking it is the old way
    //_exp_m->get_output_m()->get_tree()->set_replic_report( _id, _replic_report );
  }
  else
  {
    _replic_report = NULL;
  }
  
  // Create protein and RNA access lists
  _protein_list = new ae_list<ae_protein*>();
  _rna_list     = new ae_list<ae_rna*>();
  
  // Generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for ( int8_t i = 0 ; i < 5 ; i++ )
  {
    _int_probes[i]    = parent->_int_probes[i];
    _double_probes[i] = parent->_double_probes[i];
  }
  
  // Mutation rates etc...
  _mut_params = new ae_params_mut( *(parent->_mut_params) );
      
  // Genome size constraints
  _min_genome_length = parent->_min_genome_length;
  _max_genome_length = parent->_max_genome_length;
  
  // Plasmids settings
  _allow_plasmids         = parent->_allow_plasmids;
  _plasmid_minimal_length = parent->_plasmid_minimal_length;
  
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
ae_individual::~ae_individual( void )
{
  #ifdef DISTRIBUTED_PRNG
    delete _mut_prng;
    delete _stoch_prng;
  #endif
  
  // The _replic_report pointer is destroyed, but not the report itself,
  // it will be deleted later, when the tree is written on disk and emptied.
  
  assert( _protein_list != NULL );
  _protein_list->erase( false );
  delete _protein_list;

  assert( _rna_list != NULL );
  _rna_list->erase( false );
  delete _rna_list;

  assert( _genetic_unit_list != NULL );
  _genetic_unit_list->erase( true );
  // When the unit is destoyed, its dna is destroyed too, thus the pointer 
  // to the ae_dna_replication_report is destroyed. But the 
  // dna_replic_report object itself is not deleted, its address is 
  // still contained in the global replic_report object in the tree.
  delete _genetic_unit_list;
  
  delete _phenotype_activ;
  delete _phenotype_inhib;
  delete _phenotype;
  
  if ( _dist_to_target_by_segment != NULL ) delete [] _dist_to_target_by_segment;
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
/*!
  TODO
*/
void ae_individual::set_grid_cell( ae_grid_cell* grid_cell )
{
  _placed_in_population = true;
  _grid_cell = grid_cell;
  if ( _grid_cell->get_individual() != this )
  {
    _grid_cell->set_individual( this );
  }
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_individual::compute_phenotype( void )
{
  if ( _phenotype_computed ) return; // Phenotype has already been computed, nothing to do.
  _phenotype_computed = true;
  
  // Make sure the transcription, translation and folding stages have taken place
  do_transcription_translation_folding();
  
  
  // We will use two fuzzy sets :
  //   * _phenotype_activ for the proteins realising a set of functions
  //   * _phenotype_inhib for the proteins inhibitting a set of functions
  // The phenotype will then be given by the sum of these 2 fuzzy sets
  _phenotype_activ = new ae_fuzzy_set();
  _phenotype_inhib = new ae_fuzzy_set();
  
  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit;
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->get_obj();
    
    _phenotype_activ->add( gen_unit->get_activ_contribution() );
    _phenotype_inhib->add( gen_unit->get_inhib_contribution() );
    
    gen_unit_node = gen_unit_node->get_next();
  }
  
  _phenotype_activ->add_upper_bound(   Y_MAX );
  _phenotype_inhib->add_lower_bound( - Y_MAX );
  
  _phenotype = new ae_phenotype();
  _phenotype->add( _phenotype_activ );
  _phenotype->add( _phenotype_inhib );
  _phenotype->add_lower_bound( Y_MIN );
  _phenotype->simplify();
}

void ae_individual::compute_distance_to_target( ae_environment* envir )
// Compute the areas between the phenotype and the environment for each environmental segment
// If the environment is not segmented, the total area is computed
{
  if ( _distance_to_target_computed ) return; // _distance_to_target has already been computed, nothing to do.
  _distance_to_target_computed = true;
  
  if ( ! _phenotype_computed ) compute_phenotype();
  
  // Compute the difference between the (whole) phenotype and the environment
  ae_fuzzy_set* delta = new ae_fuzzy_set( *_phenotype );
  delta->sub( envir );
  
  if ( envir->get_nb_segments() == 1 )
  {
    _dist_to_target_by_feature[METABOLISM] = delta->get_geometric_area();
  }
  else // Environment is segmented
  {
    ae_env_segment** segments = envir->get_segments();
    _dist_to_target_by_segment = new double [envir->get_nb_segments()];
    for ( int8_t i = 0 ; i < envir->get_nb_segments() ; i++ )
    {
      _dist_to_target_by_segment[i] = 0;
    }
    
    // TODO : We should take into account that we compute the areas in order (from the leftmost segment, rightwards)
    //   => We shouldn't parse the whole list of points on the left of the segment we are considering (we have 
    //      already been through them!)
    
    for ( int16_t i = 0 ; i < envir->get_nb_segments() ; i++ )
    {
      _dist_to_target_by_segment[i] = delta->get_geometric_area( segments[i]->start, segments[i]->stop );
      _dist_to_target_by_feature[segments[i]->feature] += _dist_to_target_by_segment[i];
    }
  }
  
  delete delta;
}

/*!
  Computation of a "proper" fitness value (one that increases when the individual is fitter)

  Computation of a "proper" fitness value (one that increases when the individual is fitter)
  The behaviour of this function depends on many parameters and most notably on whether it is
  a "composite" fitness or not, and on the selection scheme.
*/
void ae_individual::compute_fitness( ae_environment* envir )
{
  if ( _fitness_computed ) return; // Fitness has already been computed, nothing to do.
  _fitness_computed = true;

  #ifdef NORMALIZED_FITNESS
    for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
    {
      if (envir->get_area_by_feature(i)==0.)
      {
        _fitness_by_feature[i] = 0.;
      }
      else
      {
        _fitness_by_feature[i] =  ( envir->get_area_by_feature(i) - _dist_to_target_by_feature[i] ) / envir->get_area_by_feature(i);
        if ( (_fitness_by_feature[i] < 0.) && (i != METABOLISM) ) // non-metabolic fitness can NOT be lower than zero (we do not want individual to secrete a negative quantity of public good)
        {
          _fitness_by_feature[i] = 0.;
        }
      }
    }
    
    if (! _placed_in_population)
    {
      _fitness = _fitness_by_feature[METABOLISM];
    }
    else
    {
      _fitness =  _fitness_by_feature[METABOLISM] * ( 1 + _exp_m->get_secretion_contrib_to_fitness() * ( _grid_cell->get_compound_amount() - _exp_m->get_secretion_cost() * _fitness_by_feature[SECRETION] ) );
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
        _fitness_by_feature[SECRETION] =  exp( - _exp_m->get_selection_pressure() * _dist_to_target_by_feature[SECRETION] )
                                        - exp( - _exp_m->get_selection_pressure() * envir->get_area_by_feature(SECRETION) );
        
        if ( _fitness_by_feature[i] < 0 )
        {
          _fitness_by_feature[i] = 0;
        }         
      }
      else
      {
        _fitness_by_feature[i] = exp( - _exp_m->get_selection_pressure() * _dist_to_target_by_feature[i] );
      }  
    }
    
    // Calculate combined, total fitness here!
    // Multiply the contribution of metabolism and the amount of compound in the environment 
    if ( ! _placed_in_population )
    { 
      _fitness =  _fitness_by_feature[METABOLISM] ; 
    }
    else
    {   
      _fitness =  _fitness_by_feature[METABOLISM]
                  *  ( 1 + _exp_m->get_secretion_contrib_to_fitness() * _grid_cell->get_compound_amount()
                         - _exp_m->get_secretion_cost() * _fitness_by_feature[SECRETION] ); 
    }
  #endif
}

void ae_individual::reevaluate( ae_environment* envir /*= NULL*/ )
{
  // useful for post-treatment programs that replay mutations
  // on a single individual playing the role of the successive
  // ancestors

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
  
  if ( envir == NULL ) envir = _exp_m->get_env();

  
  ae_list_node<ae_genetic_unit*>* unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit * unit = NULL;
  while (unit_node != NULL)
  {
    unit = unit_node->get_obj();
    unit->reset_expression();
    unit_node = unit_node->get_next();
  }

  if ( _phenotype_activ != NULL )
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
  delete [] _dist_to_target_by_segment;

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }
  
  //For each RNA / individual / genetic_unit delete proteins it knows
  //Deleting the protein it self is made only once
  
  _rna_list->erase( false );
  _protein_list->erase( false ); 
  
  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_list_node<ae_rna*>* rna_node = NULL;
  ae_rna*           rna = NULL;
  ae_genetic_unit*  gen_unit = NULL;
  
  while ( gen_unit_node != NULL )
    {
      gen_unit = gen_unit_node->get_obj();

      rna_node = (gen_unit->get_rna_list()[LEADING])->get_first();
      while(rna_node !=NULL)
      {
        rna = rna_node->get_obj();
        rna->get_transcribed_proteins()->erase( false );
        rna_node = rna_node->get_next();
      }
      rna_node = (gen_unit->get_rna_list()[LAGGING])->get_first();
      while(rna_node !=NULL)
      {
        rna = rna_node->get_obj();
        rna->get_transcribed_proteins()->erase( false );
        rna_node = rna_node->get_next();
      }
      
      (gen_unit->get_protein_list()[LEADING])->erase( true );
      (gen_unit->get_protein_list()[LAGGING])->erase( true );
      
      gen_unit_node = gen_unit_node->get_next();
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
  
  evaluate( envir );
}

void ae_individual::do_transcription( void )
{
  if ( _transcribed == true ) return; // Transcription has already been performed, nothing to do.
  _transcribed = true;
  
  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->get_obj();
    
    gen_unit->do_transcription();
    _rna_list->add_list( gen_unit->get_rna_list()[LEADING] );
    _rna_list->add_list( gen_unit->get_rna_list()[LAGGING] );
     
    gen_unit_node = gen_unit_node->get_next();
  }
}

void ae_individual::do_translation( void )
{
  if ( _translated == true ) return; // ARNs have already been translated, nothing to do.
  _translated = true;
  if ( _transcribed == false ) do_transcription();
  
  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->get_obj();
    
    gen_unit->do_translation();
    _protein_list->add_list( gen_unit->get_protein_list()[LEADING] );
    _protein_list->add_list( gen_unit->get_protein_list()[LAGGING] );
     
    gen_unit_node = gen_unit_node->get_next();
  }
}

void ae_individual::do_folding( void )
{
  if ( _folded == true ) return; // Proteins have already been folded, nothing to do.
  _folded = true;
  if ( _translated == false ) do_translation();
  
  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->get_obj();
    
    gen_unit->compute_phenotypic_contribution();
     
    gen_unit_node = gen_unit_node->get_next();
  }
}

void ae_individual::evaluate( ae_environment* envir /*= NULL*/ )
{
  if ( _evaluated == true ) return; // Individual has already been evaluated, nothing to do.
  _evaluated = true;
  
  if ( envir == NULL ) envir = _exp_m->get_env();
  
  // ----------------------------------------------------------------------
  // 1) Transcription - Translation - Folding
  // ----------------------------------------------------------------------
  do_transcription_translation_folding();
  //~ print_rna_list();
  //~ print_protein_list();
  
  // ----------------------------------------------------------------------
  // 2) Compute phenotype and compare it to the environment => fitness
  // ----------------------------------------------------------------------
  compute_phenotype();
  compute_distance_to_target( envir );
  compute_fitness( envir );
  
  if ( _exp_m->get_output_m()->get_compute_phen_contrib_by_GU() )
  { 
    ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
    ae_genetic_unit*  gen_unit = NULL;
    
    while ( gen_unit_node != NULL )
    {
      gen_unit = gen_unit_node->get_obj();
      gen_unit->compute_distance_to_target( envir );
      gen_unit->compute_fitness( envir );
      gen_unit_node = gen_unit_node->get_next();
    }
  }
}


void ae_individual::inject_GU( ae_individual* donor )
{  
  // Add the GU at the end of the list
  ae_list_node<ae_genetic_unit*>* gen_unit_node = donor->get_genetic_unit_list()->get_last();
  ae_genetic_unit* temp_GU = new ae_genetic_unit( this, *(gen_unit_node->get_obj()) );  
  _genetic_unit_list->add( temp_GU );
}

void ae_individual::inject_2GUs( ae_individual* partner )
{  
  //We swap GUs from the end of the list
  ae_list_node<ae_genetic_unit*>* gen_unit_node_1 = partner->get_genetic_unit_list()->get_last();
  ae_list_node<ae_genetic_unit*>* gen_unit_node_2 = _genetic_unit_list->get_last();
  
  ae_genetic_unit* temp_GU_1 = new ae_genetic_unit( this, *(gen_unit_node_1->get_obj()) );  
  ae_genetic_unit* temp_GU_2 = new ae_genetic_unit( this, *(gen_unit_node_2->get_obj()) );  
  
  _genetic_unit_list->add( temp_GU_1 );
  partner->get_genetic_unit_list()->add( temp_GU_2 );
}

void ae_individual::compute_statistical_data( void )
{
  if ( _statistical_data_computed ) return; // Statistical data has already been computed, nothing to do.
  _statistical_data_computed = true;
  
  if ( not _phenotype_computed )
  {
    compute_phenotype();
  }
  
  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit      = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->get_obj();

    _total_genome_size                 += gen_unit->get_dna()->get_length();
    _nb_coding_RNAs                    += gen_unit->get_nb_coding_RNAs();
    _nb_non_coding_RNAs                += gen_unit->get_nb_non_coding_RNAs();
    _overall_size_coding_RNAs          += gen_unit->get_overall_size_coding_RNAs();
    _overall_size_non_coding_RNAs      += gen_unit->get_overall_size_non_coding_RNAs();
    _nb_genes_activ                    += gen_unit->get_nb_genes_activ();
    _nb_genes_inhib                    += gen_unit->get_nb_genes_inhib();
    _nb_functional_genes               += gen_unit->get_nb_functional_genes();
    _nb_non_functional_genes           += gen_unit->get_nb_non_functional_genes();
    _overall_size_functional_genes     += gen_unit->get_overall_size_functional_genes();
    _overall_size_non_functional_genes += gen_unit->get_overall_size_non_functional_genes();
    
    gen_unit_node = gen_unit_node->get_next();
  }
  
  if ( _replic_report != NULL )
  {
    // Finalize statistical data in the replication report
    _replic_report->signal_end_of_replication();
  }
}

void ae_individual::compute_non_coding( void )
{
  if ( _non_coding_computed ) return; // _non_coding has already been computed, nothing to do.
  _non_coding_computed = true;
  
  _nb_bases_in_0_CDS                  = 0;
  _nb_bases_in_0_functional_CDS       = 0;
  _nb_bases_in_0_non_functional_CDS   = 0;
  _nb_bases_in_0_RNA                  = 0;
  _nb_bases_in_0_coding_RNA           = 0;
  _nb_bases_in_0_non_coding_RNA       = 0;
  _nb_bases_in_neutral_regions        = 0;
  _nb_neutral_regions                 = 0;

  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->get_obj();
    
    _nb_bases_in_0_CDS                  += gen_unit->get_nb_bases_in_0_CDS();
    _nb_bases_in_0_functional_CDS       += gen_unit->get_nb_bases_in_0_functional_CDS();
    _nb_bases_in_0_non_functional_CDS   += gen_unit->get_nb_bases_in_0_non_functional_CDS();
    _nb_bases_in_0_RNA                  += gen_unit->get_nb_bases_in_0_RNA();
    _nb_bases_in_0_coding_RNA           += gen_unit->get_nb_bases_in_0_coding_RNA();
    _nb_bases_in_0_non_coding_RNA       += gen_unit->get_nb_bases_in_0_non_coding_RNA();
    _nb_bases_in_neutral_regions        += gen_unit->get_nb_bases_in_neutral_regions();
    _nb_neutral_regions                 += gen_unit->get_nb_neutral_regions();
    gen_unit_node = gen_unit_node->get_next();
  }
}

void ae_individual::save( gzFile backup_file ) const
{
  // Write the "age" of the strain
  gzwrite( backup_file, &_age, sizeof(_age) );
  
  #ifdef DISTRIBUTED_PRNG
    // Write the PRNG's state
    _mut_prng->save( backup_file );
    _stoch_prng->save( backup_file );
  #endif
  
  // Write id and rank
  gzwrite( backup_file, &_id,   sizeof(_id) );
  gzwrite( backup_file, &_rank, sizeof(_rank) );
  
  // Write the position of the individual (or -1, -1 if there is no spatial structure)
  int16_t x = -1;
  int16_t y = -1;
  if ( _exp_m->is_spatially_structured() )
  {
    x = _grid_cell->get_x();
    y = _grid_cell->get_y();
  }
  gzwrite( backup_file, &x, sizeof(x) );
  gzwrite( backup_file, &y, sizeof(y) );
  
  // Write generic probes
  gzwrite( backup_file, _int_probes,    5 * sizeof(*_int_probes) );
  gzwrite( backup_file, _double_probes, 5 * sizeof(*_double_probes) );
  
  // Write mutational parameters
  _mut_params->save( backup_file );
  
  // ------------------------------------------------- Phenotypic stochasticity
  gzwrite( backup_file, &_with_stochasticity, sizeof(_with_stochasticity) );
  
  // Write artificial chemistry parameters
  gzwrite( backup_file, &_w_max, sizeof(_w_max) );
      
  // Write genome size constraints
  gzwrite( backup_file, &_min_genome_length, sizeof(_min_genome_length) );
  gzwrite( backup_file, &_max_genome_length, sizeof(_max_genome_length) );
  
  // Write plasmids settings
  int8_t tmp_allow_plasmids = _allow_plasmids;
  gzwrite( backup_file, &tmp_allow_plasmids, sizeof(tmp_allow_plasmids) );
  if ( _allow_plasmids )
  {
    gzwrite( backup_file, &_plasmid_minimal_length, sizeof(_plasmid_minimal_length) );
  }
  
  // Write genetic units
  int16_t nb_gen_units = _genetic_unit_list->get_nb_elts();
  gzwrite( backup_file, &nb_gen_units,  sizeof(nb_gen_units) );
  
  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit;
  for ( int16_t i = 0 ; i < nb_gen_units ; i++ )
  {
    gen_unit = gen_unit_node->get_obj();
    
    gen_unit->save( backup_file );
    
    gen_unit_node = gen_unit_node->get_next();
  }
}

int32_t ae_individual::get_nb_terminators( void )
{
  int32_t nb_term = 0;
  
  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit      = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->get_obj();
    
    nb_term += gen_unit->get_nb_terminators();
    
    gen_unit_node = gen_unit_node->get_next();
  }
  
  return nb_term;
}


/*!
  \brief Compute reproduction statistics and statistics about the offsprings of the current individual

  * Make nb_children replications of the current individual.
  * For each replication, determine if the offsprings is neutral, beneficial or deleterious by comparison of fitness with the current individual (the parent)
  * If statistics about offsprings are required (offsprings_statistics != NULL), fitness mean, fitness variance, size mean, size variance, functional gene number mean, 
  functional gene number variance fo the nb_children offsprings are computed
  * If informations about each children are required ( replication_file != NULL), fitness, genome_size, nb of functional genes, number of coding bases, number of transcribed 
  but not translated bases, number of non transcribed bases of each offsprings are written in replication_file
  
  \param nb_children              number of replications made to have the statistics
  \param reproduction_statistics  statistics about the replications (proportion of neutral offsprings, proportion of beneficial offsprings, proportion of deleterious offsprings)
  \param offsprings_statistics    statistics about the nb_children offsprings (fitness mean, fitness variance, size mean, size variance, functional gene number mean,
                                    functional gene number variance) compute if not null
  \param replication_file         file with informations about each children of the current individual (fitness, genome_size, nb of functional genes, number of coding bases, 
                                    number of transcribed but not translated bases, number of non transcribed bases) if not null
*/
void ae_individual::compute_experimental_f_nu( int32_t nb_children, double* reproduction_statistics, double* offsprings_statistics, FILE* replication_file)
{
  double initial_fitness = get_fitness();
 
  if(reproduction_statistics != NULL) 
  {
    reproduction_statistics[0] = 0; // proportion of neutral offsprings
    reproduction_statistics[1] = 0; // proportion of beneficial offsprings
    reproduction_statistics[2] = 0; // proportion of deleterious offsprings
  }
  else
  {
    printf( "%s:%d: error: reproduction_statistics was not initialized\n", __FILE__, __LINE__);
    exit( EXIT_FAILURE );
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
	      
  double fitness_child = 0;
		
  // replicate this individual to create 'nb_children' children 
  ae_individual * child = NULL;
  
  int32_t genome_size = 0;
  int32_t nb_functional_genes = 0;
  int32_t nb_bases_in_0_functional_CDS = 0;
  int32_t nb_bases_in_0_coding_RNA = 0;

  for (int i = 0; i < nb_children; i++)
  {
    child = _exp_m->get_exp_s()->get_sel()->do_replication( this, _id );
    fitness_child = child->get_fitness(); 

    if ( fabs(initial_fitness - fitness_child) < 1e-10*std::max(initial_fitness, fitness_child) )
    { 
      reproduction_statistics[0] += 1; 
    }
    else if ( fitness_child > initial_fitness )
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
      fprintf(replication_file, "%le %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32"\n", fitness_child, genome_size, nb_functional_genes, genome_size-nb_bases_in_0_functional_CDS, 
                                nb_bases_in_0_functional_CDS-nb_bases_in_0_coding_RNA, nb_bases_in_0_coding_RNA);
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


/*!
  \brief Compute reproduction theoretical proportion of neutral offsprings

  Compute the theoretical proportion of neutral offsprings given the Carole's formula, based on the mutations and rearrangement rates and not on multiple replications
  
  \return theoretical proportion of neutral offsprings
*/
double ae_individual::compute_theoritical_f_nu( void )
{
  // We first have to collect information about genome structure.
  // Abbreviations are chosen according to Carole's formula.
  // Please notice that compared to the formula we have the beginning
  // and ends of neutral regions instead of 'functional regions'
  ae_genetic_unit* chromosome = get_genetic_unit_list()->get_first()->get_obj();
  int32_t L       = chromosome->get_dna()->get_length();
  int32_t N_G     = chromosome->get_nb_neutral_regions(); // which is not exactly Carole's original definition
  int32_t* b_i    = chromosome->get_beginning_neutral_regions();
  int32_t* e_i    = chromosome->get_end_neutral_regions();
  int32_t lambda  = chromosome->get_nb_bases_in_neutral_regions();
  int32_t l       = L - lambda; // nb bases in 'functional regions'
  
  int32_t* lambda_i = NULL;  // nb bases in ith neutral region
  if ( N_G > 0 ) // all the chromosome may be functional
  {
    lambda_i = new int32_t [ N_G ];
  
    for ( int32_t i = 0; i < N_G - 1; i++)
    {
      lambda_i[i] = e_i[i] - b_i[i] + 1;
    }
    if ( b_i[N_G-1] > e_i[N_G-1] ) // last neutral region is overlapping on the beginning of chromosome
    {
      lambda_i[N_G-1] = ( e_i[N_G-1] + L ) - b_i[N_G-1] + 1;
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
  Fv  = pow( 1 - get_point_mutation_rate()  * ( 1 - nu_local_mutation ), L);
  Fv *= pow( 1 - get_small_insertion_rate() * ( 1 - nu_local_mutation ), L);
  Fv *= pow( 1 - get_small_deletion_rate()  * ( 1 - nu_local_mutation ), L);

  // inversion ~ two local mutations
  double nu_inversion = nu_local_mutation * nu_local_mutation;
  Fv *= pow( 1 - get_inversion_rate()       * ( 1 - nu_inversion )     , L);

  // translocation ~ inversion + insertion (mathematically)
  Fv *= pow( 1 - get_translocation_rate()   * ( 1 - nu_inversion * nu_local_mutation ), L);

  // long deletion
  double nu_deletion = 0; // if N_G == 0, a deletion is always not neutral
  for ( int32_t i = 0; i < N_G; i++) { nu_deletion += lambda_i[i] * (lambda_i[i] + 1); }
  nu_deletion /= ((double) 2*L*L);
  Fv *= pow( 1 - get_deletion_rate()        * ( 1 - nu_deletion )      , L);

  // duplication ~ big deletion + insertion
  Fv *= pow( 1 - get_duplication_rate()     * ( 1 - nu_deletion * nu_local_mutation ), L);
  
  if ( lambda_i != NULL ) delete [] lambda_i;

  return Fv;
}

/*!
  \brief Remove the bases that are not in coding RNA

  Remove the bases that are not in coding RNA and test at each loss that fitness is not changed
*/
void ae_individual::remove_non_coding_bases( void)
{
  //reevaluate(_exp_m->get_env());
  //double initial_fitness = get_fitness();
  
  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->get_obj();
    gen_unit->remove_non_coding_bases(); 
    
    //reevaluate(_exp_m->get_env()); 
    //assert(get_fitness()==initial_fitness);
    gen_unit_node = gen_unit_node->get_next();
  }
  _non_coding_computed = false;
  assert(get_nb_bases_in_0_coding_RNA()==0);
}

/*!
  \brief Double the bases that are not in coding RNA

  Double the bases that are not in coding RNA by addition of random bases and test at each addition that fitness is not changed
*/
void ae_individual::double_non_coding_bases(void)
{
  //reevaluate(_exp_m->get_env());
  //double initial_fitness = get_fitness();
  int32_t initial_non_coding_base_nb = get_nb_bases_in_0_coding_RNA();
  
  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->get_obj();
    gen_unit->double_non_coding_bases(); 
    
    //reevaluate(_exp_m->get_env()); 
    //assert(get_fitness()==initial_fitness);
    gen_unit_node = gen_unit_node->get_next();
  }
  _non_coding_computed = false;
  assert(get_nb_bases_in_0_coding_RNA()==2*initial_non_coding_base_nb);
}
// =================================================================
//                           Protected Methods
// =================================================================
void ae_individual::make_protein_list( void )
{
  assert( _protein_list != NULL );
  
  // Clean list
  _protein_list->erase( false );
  
  // Make a copy of each genetic unit's protein list
  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit      = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->get_obj();
    
    _protein_list->add_list( gen_unit->get_protein_list()[LEADING] );
    _protein_list->add_list( gen_unit->get_protein_list()[LAGGING] );
    
    gen_unit_node = gen_unit_node->get_next();
  }
}

void ae_individual::make_rna_list( void )
{
  assert( _rna_list != NULL );
  
  // Clean list
  _rna_list->erase( false );
  
  // Make a copy of each genetic unit's rna list
  ae_list_node<ae_genetic_unit*>* gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = gen_unit_node->get_obj();
    
    _rna_list->add_list( gen_unit->get_rna_list()[LEADING] );
    _rna_list->add_list( gen_unit->get_rna_list()[LAGGING] );
    
    gen_unit_node = gen_unit_node->get_next();
  }
}
