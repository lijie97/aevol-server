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
  Creates a new individual with a random genome.

  Promoters will be looked for on the whole genome but no further process
  will be performed.
  The phenotype and the fitness are not set, neither is the statistical data.
*/
/*ae_individual::ae_individual( void )
{
  _prng = NULL;
  _age  = 0;
  
  _index_in_population  = -1;
  _rank_in_population   = -1;
  
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
  
  // Create a list of genetic units with only the main chromosome
  _genetic_unit_list = new ae_list<ae_genetic_unit*>();
  
  if ( _exp_m->get_allow_plasmids() )
  {
    // Create both a chromosome and a plasmid, separately 
    if ( _exp_m->get_plasmid_initial_gene() == 1 ) 
    {
      _genetic_unit_list->add( new ae_genetic_unit( this, _exp_m->get_initial_genome_length() ) );
      _genetic_unit_list->add( new ae_genetic_unit( this, _exp_m->get_plasmid_initial_length() ) );
    }
    // create a plasmid that is identical to the chromosome
    else // _exp_m->get_plasmid_initial_gene() == 2
    {
      _genetic_unit_list->add( new ae_genetic_unit( this, _exp_m->get_initial_genome_length() ) );
      _genetic_unit_list->add( new ae_genetic_unit( this, *(_genetic_unit_list->get_first()->get_obj()) ) );
    }
  
  }
  else
  {
    _genetic_unit_list->add( new ae_genetic_unit( this, _exp_m->get_initial_genome_length() ) );
  }
  
  // Create empty fuzzy sets for activation and inhibition
  _phenotype_activ  = NULL;
  _phenotype_inhib  = NULL;
  _phenotype        = NULL;
  
  // Initialize all the fitness-related stuff
  if ( ae_common::sim->get_env()->is_segmented() )
  {
    int16_t nb_segments = ae_common::sim->get_env()->get_nb_segments();
    _dist_to_target_by_segment = new double [nb_segments];
    
    for ( int16_t i = 0 ; i < nb_segments ; i++ )
    {
      _dist_to_target_by_segment[i] = 0;
    }
  }
  else
  {
    _dist_to_target_by_segment = NULL;
  }
  
  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];
  
  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }
  
  // We are at generation 0, individual has no replication to report
  _replic_report = NULL;
  
  _protein_list = new ae_list<ae_protein*>();
  _rna_list     = new ae_list<ae_rna*>();
  
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for ( int8_t i = 0 ; i < 5 ; i++ )
  {
    _int_probes[i]    = 0;
    _double_probes[i] = 0;
  }
  
  // Initialize statistical data // TODO : => function
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
}*/

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
}

/*!
  Creates a new individual with a genome read in from a file.
 
  Promoters will be looked for on the whole genome but no further process
  will be performed.
  The phenotype and the fitness are not set, neither is the statistical data.
*/
/*ae_individual::ae_individual( char* organism_file_name )
{
  _prng = NULL;
  
  _index_in_population  = -1;
  _rank_in_population   = -1;
  
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
  
  // Create a list of genetic units with only the main chromosome
  _genetic_unit_list = new ae_list<ae_genetic_unit*>();
  
  if ( ae_common::params->get_allow_plasmids() )
  {
    printf( "ERROR, not implemented for plasmids yet %s:%d\n", __FILE__, __LINE__ );
  }
  else
  {
    _genetic_unit_list->add( new ae_genetic_unit( this, organism_file_name ) );
  }
  
  // Create empty fuzzy sets for activation and inhibition
  _phenotype_activ  = NULL;
  _phenotype_inhib  = NULL;
  _phenotype        = NULL;
  
  // Initialize all the fitness-related stuff
  if ( ae_common::sim->get_env()->is_segmented() )
  {
    int16_t nb_segments = ae_common::sim->get_env()->get_nb_segments();
    _dist_to_target_by_segment = new double [nb_segments];
    
    for ( int16_t i = 0 ; i < nb_segments ; i++ )
    {
      _dist_to_target_by_segment[i] = 0;
    }
  }
  else
  {
    _dist_to_target_by_segment = NULL;
  }
  
  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];
  
  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }
  
  // We are at generation 0, individual has no replication to report
  _replic_report = NULL;
  
  _protein_list = new ae_list<ae_protein*>();
  _rna_list     = new ae_list<ae_rna*>();
  
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for ( int8_t i = 0 ; i < 5 ; i++ )
  {
    _int_probes[i]    = 0;
    _double_probes[i] = 0;
  }
  
  // Initialize statistical data // TODO : => function
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
  
  _modularity = -1;
}*/


/*!
  Creates a new individual with the given genome.
 
  Promoters will be looked for on the whole genome but no further process
  will be performed.
  The phenotype and the fitness are not set, neither is the statistical data.
 
  WARNING : 
    genome will be used directly which means the caller must not delete it
*/
//ae_individual::ae_individual( char* genome, int32_t genome_size, int32_t age /*= 0*/ )
/*{
  _age  = age;
  _prng = NULL;
  
  _index_in_population  = -1;
  _rank_in_population   = -1;
  
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
  
  // Create a list of genetic units with only the main chromosome
  _genetic_unit_list = new ae_list<ae_genetic_unit*>();
  
  if ( ae_common::params->get_allow_plasmids() )
  {
    printf( "ERROR, functionality not implemented %s:%d\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
  else
  {
    _genetic_unit_list->add( new ae_genetic_unit( this, genome, genome_size ) );
  }
  
  // Create empty fuzzy sets for activation and inhibition
  _phenotype_activ  = NULL;
  _phenotype_inhib  = NULL;
  _phenotype        = NULL;
  
  // Initialize all the fitness-related stuff
  if ( ae_common::sim->get_env()->is_segmented() )
  {
    int16_t nb_segments = ae_common::sim->get_env()->get_nb_segments();
    _dist_to_target_by_segment = new double [nb_segments];
    
    for ( int16_t i = 0 ; i < nb_segments ; i++ )
    {
      _dist_to_target_by_segment[i] = 0;
    }
  }
  else
  {
    _dist_to_target_by_segment = NULL;
  }
  
  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];
  
  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }
  
  // Individual has no replication to report
  _replic_report = NULL;
  
  _protein_list = new ae_list<ae_protein*>();
  _rna_list     = new ae_list<ae_rna*>();
  
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for ( int8_t i = 0 ; i < 5 ; i++ )
  {
    _int_probes[i]    = 0;
    _double_probes[i] = 0;
  }
  
  // Initialize statistical data // TODO : => function
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
}*/



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
  ae_genetic_unit*  gen_unit;
  
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

/*void ae_individual::compute_modularity( void )
{
  if ( _modularity_computed ) return; // _modularity has already been computed, nothing to do.
  _modularity_computed = true;
  
  double av_dist_between_genes              = 0.0;
  double av_dist_between_interacting_genes  = 0.0; // Genes which corresponding phenotypic triangles overlap
  int32_t nb_gene_pairs                     = 0;
  int32_t nb_interacting_gene_pairs         = 0;
  
  ae_list_node<ae_protein*>* gene_node_1 = NULL;
  ae_list_node<ae_protein*>* gene_node_2 = NULL;
  ae_protein*   gene_1 = NULL;
  ae_protein*   gene_2 = NULL;
  
  bool overlap;
  int32_t dist;
  
  // ------------------------------------------------
  // Compute the distances between each pair of genes
  // ------------------------------------------------
  
  // 1) LEADING genes with X (LEADING or LAGGING)
  gene_node_1 = _protein_list[LEADING]->get_first();
  while ( gene_node_1 != NULL )
  {
    gene_1 = gene_node_1->get_obj();
    
    if ( gene_1->is_degenerated() )
    {
      gene_node_1 = gene_node_1->get_next();
      continue;
    }
    
    // A) Distance with the other LEADING genes
    gene_node_2 = gene_node_1->get_next();
    
    while ( gene_node_2 != NULL )
    {
      gene_2 = gene_node_2->get_obj();
      
      if ( gene_2->is_degenerated() )
      {
        gene_node_2 = gene_node_2->get_next();
        continue;
      }
      
      int32_t f1 = gene_1->get_first_translated_pos();
      int32_t l1 = gene_1->get_last_translated_pos();
      int32_t f2 = gene_2->get_first_translated_pos();
      int32_t l2 = gene_2->get_last_translated_pos();
      
      // Check overlapping of the 2 genes
      overlap = ( f1 > l1 || l1 >= f2 || ( f2 > l2 && l2 >= f1 ) );
      
      // Compute distance between genes
      if ( overlap )
      {
        dist = 0;
      }
      else
      {
        dist = utils::min( f2 - l1 - 1, utils::mod( f1 - l2 - 1, _genome->get_length() ) );
      }
      
      // Check whether the phenotypic triangles of the genes overlap
      if ( fabs( gene_1->get_mean() - gene_2->get_mean() ) < gene_1->get_width() + gene_2->get_width() )
      {
        nb_interacting_gene_pairs++;
        av_dist_between_interacting_genes += dist;
      }
      nb_gene_pairs++;
      av_dist_between_genes += dist;
      
      
      gene_node_2 = gene_node_2->get_next();
    }
    
    // B) Distance with the LAGGING genes
    gene_node_2 = _protein_list[LAGGING]->get_first();
    
    while ( gene_node_2 != NULL )
    {
      gene_2 = gene_node_2->get_obj();
      
      if ( gene_2->is_degenerated() )
      {
        gene_node_2 = gene_node_2->get_next();
        continue;
      }
      
      int32_t f1 = gene_1->get_first_translated_pos();
      int32_t l1 = gene_1->get_last_translated_pos();
      int32_t f2 = gene_1->get_first_translated_pos();
      int32_t l2 = gene_1->get_last_translated_pos();
      
      // Check overlapping of the 2 genes
      // As gene_1 is on the leading strand and gene_2 on the lagging strand,
      // they will overlap unless the 4 points f1, l1, f2, l2 are ordered as follows :
      // f1, l1, l2, f2. Because the genome is circular, there are 4 possible cases : 
      // f1, l1, l2, f2  ;  f2, f1, l1, l2  ;  l2, f2, f1, l1  ;  l1, l2, f2, f1
      // We can express this as "genes will overlap unless they admit 3 of the following conditions :
      // (f1 < l1), (l1 < l2), (l2 < f2), (f2 < f1)"
      // We will hence test all 4 conditions and count how many are true
      int8_t trues = 0;
      if ( f1 < l1 ) trues++;
      if ( l1 < l2 ) trues++;
      if ( l2 < f2 ) trues++;
      if ( f2 < f1 ) trues++;
      overlap = ( trues != 3 );
      
      // Compute distance between genes
      if ( overlap )
      {
        dist = 0;
      }
      else
      {
        dist = utils::min( utils::mod( l2 - l1 - 1, _genome->get_length() ),
                           utils::mod( f1 - f2 - 1, _genome->get_length() ) );
      }
      
      // Check whether the phenotypic triangles of the genes overlap
      if ( fabs( gene_1->get_mean() - gene_2->get_mean() ) < gene_1->get_width() + gene_2->get_width() )
      {
        nb_interacting_gene_pairs++;
        av_dist_between_interacting_genes += dist;
      }
      nb_gene_pairs++;
      av_dist_between_genes += dist;
      
      
      gene_node_2 = gene_node_2->get_next();
    }
    
    gene_node_1 = gene_node_1->get_next();
  }
  
  // 2) LAGGING genes with X (LAGGING)
  gene_node_1 = _protein_list[LAGGING]->get_first();
  while ( gene_node_1 != NULL )
  {
    gene_1 = gene_node_1->get_obj();
    
    if ( gene_1->is_degenerated() )
    {
      gene_node_1 = gene_node_1->get_next();
      continue;
    }
    
    gene_node_2 = gene_node_1->get_next();
    while ( gene_node_2 != NULL )
    {
      gene_2 = gene_node_2->get_obj();
      
      if ( gene_2->is_degenerated() )
      {
        gene_node_2 = gene_node_2->get_next();
        continue;
      }
      
      int32_t f1 = gene_1->get_first_translated_pos();
      int32_t l1 = gene_1->get_last_translated_pos();
      int32_t f2 = gene_2->get_first_translated_pos();
      int32_t l2 = gene_2->get_last_translated_pos();
      
      // Check overlapping of the 2 genes
      overlap = ( f2 >= l1 || l1 > f1 || ( l2 <= f1 && l2 > f2 ) );
      
      // Compute distance between genes
      if ( overlap )
      {
        dist = 0;
      }
      else
      {
        dist = utils::min( utils::mod( l1 - f2 - 1, _genome->get_length() ),
                           utils::mod( l2 - f1 - 1, _genome->get_length() ) );
      }
      
      // Check whether the phenotypic triangles of the genes overlap
      if ( fabs( gene_1->get_mean() - gene_2->get_mean() ) < gene_1->get_width() + gene_2->get_width() )
      {
        nb_interacting_gene_pairs++;
        av_dist_between_interacting_genes += dist;
      }
      nb_gene_pairs++;
      av_dist_between_genes += dist;
      
      
      gene_node_2 = gene_node_2->get_next();
    }
    
    gene_node_1 = gene_node_1->get_next();
  }
  
  
  if ( nb_interacting_gene_pairs > 0 && av_dist_between_interacting_genes > 0 )
  {
    av_dist_between_genes /= (double)nb_gene_pairs;
    av_dist_between_interacting_genes /= (double)nb_interacting_gene_pairs;
  }
  else
  {
    // If the value is irrelevant, replace it with the neutral value (1)
    av_dist_between_interacting_genes = av_dist_between_genes = 1;
  }
  
  _modularity = av_dist_between_interacting_genes / av_dist_between_genes;
}*/

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

double ae_individual::compute_experimental_f_nu( int32_t nb_children, double* neutral_or_better /*=NULL*/ )
{
  double initial_fitness = get_fitness();
  double Fv = 0;
  if ( neutral_or_better != NULL ) { *neutral_or_better = 0; }

  // ------------------------------------------
  //      Simulate fitness degradation
  // ------------------------------------------
	      
  double fitness_child = 0;
		
  // replicate this individual to create 'nb_children' children 
  ae_individual * child = NULL;
  for (int i = 0; i < nb_children; i++)
    {
      child = _exp_m->get_exp_s()->get_sel()->do_replication( this, _id );
      fitness_child = child->get_fitness(); // child is automatically evaluated
      //count neutral offspring
      if ( fabs(initial_fitness - fitness_child) < 1e-15 )
      { 
        Fv += 1; 
        if ( neutral_or_better != NULL ) { *neutral_or_better += 1; }
      }
      else if ( (neutral_or_better != NULL) && (fitness_child > initial_fitness) )
      {
        *neutral_or_better += 1;
      }
      delete child;
    }
  //compute Fv
  Fv /= nb_children;
  if ( neutral_or_better != NULL ) { *neutral_or_better /= nb_children; }
  return Fv;
}

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
