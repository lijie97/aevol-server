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
#include <ae_simulation.h>
#include <ae_population.h>
#include <ae_codon.h>
#include <ae_common.h>
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
/**
 * Creates a new individual with a random genome.
 *
 * Promoters will be looked for on the whole genome but no further process
 * will be performed.
 * The phenotype and the fitness are not set, neither is the statistical data.
 */
ae_individual::ae_individual( void )
{
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
  _genetic_unit_list = new ae_list();
  
  if ( ae_common::allow_plasmids )
  {
    // Create both a chromosome and a plasmid, separately 
    if ( ae_common::plasmid_initial_gene == 1 ) 
    {
      _genetic_unit_list->add( new ae_genetic_unit( this, ae_common::initial_genome_length ) );
      _genetic_unit_list->add( new ae_genetic_unit( this, ae_common::plasmid_initial_length ) );
    }
    // create a plasmid that is identical to the chromosome
    else // ae_common::plasmid_initial_gene == 2
    {
      _genetic_unit_list->add( new ae_genetic_unit( this, ae_common::initial_genome_length ) );
      _genetic_unit_list->add( new ae_genetic_unit( this, *((ae_genetic_unit*)_genetic_unit_list->get_first()->get_obj()) ) );
    }
  
  }
  else
  {
    _genetic_unit_list->add( new ae_genetic_unit( this, ae_common::initial_genome_length ) );
  }
  
  // Create empty fuzzy sets for activation and inhibition
  _phenotype_activ = new ae_fuzzy_set();
  _phenotype_inhib = new ae_fuzzy_set();
  _phenotype = NULL;
  
  // Initialize all the fitness-related stuff
  if ( ae_common::env_axis_is_segmented )
  {
    _dist_to_target_segment = new double [ae_common::env_axis_nb_segments];
    
    for ( int16_t i = 0 ; i < ae_common::env_axis_nb_segments ; i++ )
    {
      _dist_to_target_segment[i] = 0;
    }
  }
  else
  {
    _dist_to_target_segment = NULL;
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
  
  _protein_list = new ae_list();
  _rna_list     = new ae_list();
  
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
}

// Copy constructor
ae_individual::ae_individual( const ae_individual &model )
{
  _index_in_population  = model._index_in_population;
  _rank_in_population   = model._rank_in_population;
  
  _evaluated                    = model._evaluated;
  _transcribed                  = model._transcribed;
  _translated                   = model._translated;
  _folded                       = model._folded;
  _phenotype_computed           = model._phenotype_computed;
  _distance_to_target_computed  = model._distance_to_target_computed;
  _fitness_computed             = model._fitness_computed;
  _statistical_data_computed    = model._statistical_data_computed;
  _non_coding_computed          = model._non_coding_computed;
  _modularity_computed          = model._modularity_computed;
  
  _placed_in_population         = model._placed_in_population;
  
  // Create an empty list of genetic units
  _genetic_unit_list = new ae_list();
  
  // Copy model's genetic units
  ae_list_node* gen_unit_node = model._genetic_unit_list->get_first();
  while ( gen_unit_node != NULL )
  {
    _genetic_unit_list->add( new ae_genetic_unit( this, *((ae_genetic_unit*)gen_unit_node->get_obj()) ) );
    
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
    _phenotype_activ  = new ae_fuzzy_set();;
    _phenotype_inhib  = new ae_fuzzy_set();;
    _phenotype        = NULL;
  }
  
  // Copy fitness-related stuff
  if ( ae_common::env_axis_is_segmented )
  {
    _dist_to_target_segment = new double [ae_common::env_axis_nb_segments];
    
    for ( int16_t i = 0 ; i < ae_common::env_axis_nb_segments ; i++ )
    {
      _dist_to_target_segment[i] = model._dist_to_target_segment[i];
    }
  }
  else
  {
    _dist_to_target_segment = model._dist_to_target_segment;
  }
  
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
  
  _modularity = model._modularity;
  
  // We don't copy the generation report
  // One of the reasons is that we use this at generation 0 and _replic_report of the model is NULL
  _replic_report = NULL;
  
  _protein_list = new ae_list();
  _rna_list     = new ae_list();
  
  // Generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for ( int8_t i = 0 ; i < 5 ; i++ )
  {
    _int_probes[i]    = model._int_probes[i];
    _double_probes[i] = model._double_probes[i];
  }
}

/**
 * This constructor creates a new individual with the same genome as it's parent.
 * The location of promoters will be copied but no further process will be performed.
 *
 * The phenotype and the fitness are not set, neither is the statistical data.
 */
ae_individual::ae_individual( ae_individual* const parent, int32_t index )
{
  _index_in_population  = index;
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

  // Create new genetic units with their DNA copied from here
  // NOTE : The RNA lists (one per genetic unit) will also be copied so that we don't
  // need to look for promoters on the whole genome
  _genetic_unit_list = new ae_list();
  ae_list_node* gen_unit_node = parent->_genetic_unit_list->get_first();
  while ( gen_unit_node != NULL )
  {
    _genetic_unit_list->add( new ae_genetic_unit( this, (ae_genetic_unit*)gen_unit_node->get_obj() ) );
    
    gen_unit_node = gen_unit_node->get_next();
  }
  
  // Create empty fuzzy sets for activation and inhibition
  _phenotype_activ = new ae_fuzzy_set();
  _phenotype_inhib = new ae_fuzzy_set();
  _phenotype = NULL;
  
  // Initialize all the fitness-related stuff
  if ( ae_common::env_axis_is_segmented )
  {
    _dist_to_target_segment = new double [ae_common::env_axis_nb_segments];
    
    for ( int16_t i = 0 ; i < ae_common::env_axis_nb_segments ; i++ )
    {
      _dist_to_target_segment[i] = 0;
    }
  }
  else
  {
    _dist_to_target_segment = NULL;
  }
  
  _dist_to_target_by_feature  = new double [NB_FEATURES];
  _fitness_by_feature         = new double [NB_FEATURES];
  
  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }

  
  if ( ae_common::record_tree && ae_common::tree_mode == NORMAL )
  {
    _replic_report = new ae_replication_report( this, parent );
  }
  else
  {
    _replic_report = NULL;
  }
  
  _protein_list = new ae_list();
  _rna_list     = new ae_list();
  
  // Generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for ( int8_t i = 0 ; i < 5 ; i++ )
  {
    _int_probes[i]    = 0;
    _double_probes[i] = 0;
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
  
  _modularity = -1;
}

/**
 * This constructor retreives an individual from a backup file.
 * Since this generation has already been processed, no unnecessary calculation will be done.
 * However, promoters will be looked for on the whole genome.
 *
 * No transcription, translation or other process of that kind is performed.
 */
ae_individual::ae_individual( gzFile* backup_file )
{
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
  
  // Retreive index in population
  gzread( backup_file, &_index_in_population, sizeof(_index_in_population) );
  
  // The rank is not in the backup, it can be interpolated because individuals are written
  // ordered by their increasing rank (<=> increasing fitness)
  
  // Retreive generic probes
  _int_probes     = new int32_t[5];
  _double_probes  = new double[5];
  for ( int8_t i = 0 ; i < 5 ; i++ )
  {
    gzread( backup_file, &_int_probes[i],     sizeof(_int_probes[i]) );
    gzread( backup_file, &_double_probes[i],  sizeof(_double_probes[i]) );
  }
  
  // Retreive genetic units
  _genetic_unit_list = new ae_list();
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
  _phenotype_activ = new ae_fuzzy_set();
  _phenotype_inhib = new ae_fuzzy_set();
  _phenotype = NULL;
  
  
  // Initialize all the fitness-related stuff
  if ( ae_common::env_axis_is_segmented )
  {
    _dist_to_target_segment = new double [ae_common::env_axis_nb_segments];
    
    for ( int16_t i = 0 ; i < ae_common::env_axis_nb_segments ; i++ )
    {
      _dist_to_target_segment[i] = 0;
    }
  }
  else
  {
    _dist_to_target_segment = NULL;
  }
  
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
  
  _protein_list = new ae_list();
  _rna_list     = new ae_list();
  
  
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
  
  _modularity = -1;
}

// =================================================================
//                             Destructors
// =================================================================
ae_individual::~ae_individual( void )
{
  // The _replic_report pointer is destroyed, but not the report itself,
  // it will be deleted later, when the tree is written on disk and emptied.
  
  assert( _protein_list != NULL );
  _protein_list->erase( NO_DELETE );
  delete _protein_list;

  assert( _rna_list != NULL );
  _rna_list->erase( NO_DELETE );
  delete _rna_list;

  _genetic_unit_list->erase( DELETE_OBJ );
  // When the unit is destoyed, its dna is destroyed too, thus the pointer 
  // to the ae_dna_replication_report is destroyed. But the 
  // dna_replic_report object itself is not deleted, its address is 
  // still contained in the global replic_report object in the tree.
  delete _genetic_unit_list;
  
  if ( _phenotype_activ != NULL ) delete _phenotype_activ;
  if ( _phenotype_inhib != NULL ) delete _phenotype_inhib;
  if ( _phenotype != NULL ) delete _phenotype;
  
  if ( _dist_to_target_segment != NULL )     delete [] _dist_to_target_segment;
  if ( _dist_to_target_by_feature != NULL ) delete [] _dist_to_target_by_feature;
  
  if ( _fitness_by_feature != NULL ) delete [] _fitness_by_feature;
  
  // Generic probes
  delete [] _int_probes;
  delete [] _double_probes;
}

// =================================================================
//                            Public Methods
// =================================================================
/**
 * This method creates a new individual, replicating "this".
 *
 * The new individual is fully evaluated
 * Statistical data is computed
 * The x and y parameters (default 0) are used only for spacialized populations
 */
ae_individual* ae_individual::do_replication( int32_t index, int16_t x, int16_t y )
{
  ae_individual* new_indiv;

  // Copy parent
  #ifdef __NO_X
    #ifndef __REGUL
      new_indiv = new ae_individual( this, index );
    #else
      new_indiv = new ae_individual_R( dynamic_cast<ae_individual_R*>(this), index );
    #endif
  #elif defined __X11
    #ifndef __REGUL
      new_indiv = new ae_individual_X11( dynamic_cast<ae_individual_X11*>(this), index );
    #else
      new_indiv = new ae_individual_R_X11( dynamic_cast<ae_individual_R_X11*>(this), index );
    #endif
  #endif
  
  // Update the new individual's location on the grid
  // this is needed if the population is structured
  if ( ae_common::pop_structure == true )
  {
    new_indiv->set_grid_cell( ae_common::sim->get_pop()->get_pop_grid()[x][y] );
    new_indiv->set_placed_in_population( true );
    ae_common::sim->get_pop()->get_pop_grid()[x][y]->set_individual( new_indiv );
  }

  // Perform mutations for each genetic unit
  int32_t number_GU=ae_common::sim->alea->random((int32_t) 2);// It draws a number between 0 and 1.
  
  ae_list_node*     gen_unit_node = new_indiv->_genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit      = NULL;
  if(number_GU<1) // if the number is 0 : we begin with the chromosome
	{
		while ( gen_unit_node != NULL )
		{
			gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();
			
			gen_unit->get_dna()->perform_mutations();
			
			if ( new_indiv->_replic_report != NULL )
			{
				new_indiv->_replic_report->get_dna_replic_reports()->add( gen_unit->get_dna()->get_replic_report() );
			}

			gen_unit_node = gen_unit_node->get_next();
		}

	}  
  else  // if the number is 1 : we begin with the plasmid
	{
		gen_unit_node = gen_unit_node->get_next();  // we take the second GU
		while ( gen_unit_node != NULL )
		{
			gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();
			
			gen_unit->get_dna()->perform_mutations();
			
			if ( new_indiv->_replic_report != NULL )
			{
				new_indiv->_replic_report->get_dna_replic_reports()->add( gen_unit->get_dna()->get_replic_report() );
			}

			gen_unit_node = gen_unit_node->get_next();
		}
		
		gen_unit_node = new_indiv->_genetic_unit_list->get_first(); // We make the operation on the first GU
		gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();
		
		gen_unit->get_dna()->perform_mutations();
		
		if ( new_indiv->_replic_report != NULL )
		{
			new_indiv->_replic_report->get_dna_replic_reports()->add( gen_unit->get_dna()->get_replic_report() );
		}

		gen_unit_node = gen_unit_node->get_next();
		
	}

  // TODO : check this
  //~ new_indiv->set_secretion_present( this->_secretion_present ); 
  
  // Evaluate new individual
  new_indiv->evaluate( ae_common::sim->get_env() );
  new_indiv->compute_statistical_data();
  
  if ( new_indiv->_replic_report != NULL )
  {
    // Finalize statistical data in the replication report
    new_indiv->_replic_report->signal_end_of_replication();
  }
  
  #ifdef DEBUG
    ae_common::sim->get_logs()->flush();
    new_indiv->assert_promoters();
  #endif
  
  return new_indiv;
}

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
  
  ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit;
  while ( gen_unit_node != NULL )
  {
    gen_unit = (ae_genetic_unit*)gen_unit_node->get_obj();
    
    _phenotype_activ->add( gen_unit->get_activ_contribution() );
    _phenotype_inhib->add( gen_unit->get_inhib_contribution() );
    
    gen_unit_node = gen_unit_node->get_next();
  }
  
  _phenotype_activ->add_upper_bound( MAX_Y );
  _phenotype_inhib->add_lower_bound( -MAX_Y );
  
  _phenotype = new ae_phenotype();
  _phenotype->add( _phenotype_activ );
  _phenotype->add( _phenotype_inhib );
  _phenotype->add_lower_bound( MIN_Y );
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
  
  if ( ! ae_common::env_axis_is_segmented )
  {
    _dist_to_target_by_feature[METABOLISM] = delta->get_geometric_area();
  }
  else // Environment is segmented
  {
    ae_env_segment** segments  = envir->get_segments();
    
    // TODO : We should take into account that we compute the areas in order (from the leftmost segment, rightwards)
    //   => We shouldn't parse the whole list of points on the left of the segment we are considering (we have 
    //      already been through them!)
    
    for ( int16_t i = 0 ; i < ae_common::env_axis_nb_segments ; i++ )
    {
      _dist_to_target_segment[i] = delta->get_geometric_area( segments[i]->start, segments[i]->stop );
      _dist_to_target_by_feature[segments[i]->feature] += _dist_to_target_segment[i];
    }
  }
  
  delete delta;
}
    
void ae_individual::compute_fitness( ae_environment* envir )
// Computation of a "proper" fitness value (one that increases when the individual is fitter)
// The behaviour of this function depends on many parameters and most notably on whether it is
// a "composite" fitness or not, and on the selection scheme.
{
  if ( _fitness_computed ) return; // Fitness has already been computed, nothing to do.
  _fitness_computed = true;
  
  if ( ! ae_common::composite_fitness )
  {
    if ( ae_common::selection_scheme == FITNESS_PROPORTIONATE )
    {
      _fitness_by_feature[METABOLISM] = exp( -ae_common::selection_pressure * _dist_to_target_by_feature[METABOLISM] );
    }
    else // if ( ae_common::selection_scheme == XXX_RANKING )
    {
      // For ranking schemes, we use a faster measure, without exponential computation 
      _fitness_by_feature[METABOLISM] = ((MAX_X - MIN_X) * (MAX_Y - MIN_Y)) - _dist_to_target_by_feature[METABOLISM];
    }
    _fitness = _fitness_by_feature[METABOLISM];
  }
  else // Fitness is composite (e.g. metabolism + secretion)
  {
    for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
    {
      if ( i == SECRETION )
      {
        _fitness_by_feature[SECRETION] = exp( -ae_common::selection_pressure * _dist_to_target_by_feature[SECRETION] )
                                        - exp( -ae_common::selection_pressure
                                              * envir->get_area_by_feature(SECRETION) );
        
        if ( _fitness_by_feature[i] < 0 )
        {
          _fitness_by_feature[i] = 0;
        }         
      }
      else
      {
        _fitness_by_feature[i] = exp( -ae_common::selection_pressure * _dist_to_target_by_feature[i] );
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
      _fitness =  _fitness_by_feature[METABOLISM] * 
                  ( 1 + ae_common::secretion_fitness_contrib * _grid_cell->get_compound_amount()
                      - ae_common::secretion_cost * _fitness_by_feature[SECRETION] ); 
    }
  }
}


void ae_individual::reevaluate( ae_environment* envir )
{
  // useful for post-treatment programs, who replay mutations
  // on a single individual playing the role of the successive
  // ancestors
  _transcribed                  = false;
  _translated                   = false;
  _folded                       = false;
  _evaluated                    = false;
  _phenotype_computed           = false;
  _distance_to_target_computed  = false;
  _fitness_computed             = false;
  _statistical_data_computed    = false;
  _non_coding_computed          = false;
  _modularity_computed          = false;

  
  ae_list_node * unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit * unit = NULL;
  while (unit_node != NULL)
  {
    unit = (ae_genetic_unit *) unit_node->get_obj();
    
    unit->reset_expression();
    
    unit_node = unit_node->get_next();
  }

  
  
  if (_phenotype_activ != NULL)
  {
    delete _phenotype_activ;
    _phenotype_activ = new ae_fuzzy_set();
  }

  if (_phenotype_inhib != NULL)
  {
    delete _phenotype_inhib;
    _phenotype_inhib = new ae_fuzzy_set();
  }

  if (_phenotype != NULL)
  {
    delete _phenotype;
    _phenotype = NULL;
  }

  // Initialize all the fitness-related stuff
  if ( ae_common::env_axis_is_segmented )
  {
    for ( int16_t i = 0 ; i < ae_common::env_axis_nb_segments ; i++ )
    {
      _dist_to_target_segment[i] = 0;
    }
  }
  else
  {
    _dist_to_target_segment = NULL;
  }
  
  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _dist_to_target_by_feature[i] = 0.0;
    _fitness_by_feature[i]        = 0.0;
  }
  
  _rna_list->erase( NO_DELETE );
  _protein_list->erase( NO_DELETE );
  


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
  
  _modularity = -1;

  
  evaluate( envir );
}


void ae_individual::evaluate( ae_environment* envir )
{
  if ( _evaluated == true ) return; // Individual has already been evaluated, nothing to do.
  _evaluated = true;
  
  // ----------------------------------------------------------------------
  // 1) Transcription - Translation - Folding
  // ----------------------------------------------------------------------
  do_transcription_translation_folding();
  //~ print_protein_list();
  
  // ----------------------------------------------------------------------
  // 2) Compute phenotype and compare it to the environment => fitness
  // ----------------------------------------------------------------------
  compute_phenotype();
  compute_distance_to_target( envir );
  compute_fitness( envir );
  
  if ( ae_common::compute_phen_contrib_by_GU ) 
  { 
    ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
    ae_genetic_unit*  gen_unit      = NULL;
    
    while ( gen_unit_node != NULL )
    {
      gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();
      gen_unit->compute_distance_to_target( envir );
      gen_unit->compute_fitness( envir );
      gen_unit_node = gen_unit_node->get_next();
    }
  }
}

void ae_individual::inject_GU( ae_individual* donor )
{  
  ae_list_node* gen_unit_node = donor->get_genetic_unit_list()->get_first()->get_next();
  ae_genetic_unit* temp_GU = new ae_genetic_unit( this, *((ae_genetic_unit*)gen_unit_node->get_obj()) );  
  _genetic_unit_list->add( temp_GU) ;
}


void ae_individual::compute_statistical_data( void )
{
  if ( _statistical_data_computed ) return; // Statistical data has already been computed, nothing to do.
  _statistical_data_computed = true;
  
  //~ assert(_phenotype_computed );
  
  if ( ! _phenotype_computed )
  {
    compute_phenotype();
  }
  
  ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit      = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();

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
  
  ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = (ae_genetic_unit*)gen_unit_node->get_obj();
    
    _nb_bases_in_0_CDS                  += gen_unit->get_nb_bases_in_0_CDS();
    _nb_bases_in_0_functional_CDS       += gen_unit->get_nb_bases_in_0_functional_CDS();
    _nb_bases_in_0_non_functional_CDS   += gen_unit->get_nb_bases_in_0_non_functional_CDS();
    _nb_bases_in_0_RNA                  += gen_unit->get_nb_bases_in_0_RNA();
    _nb_bases_in_0_coding_RNA           += gen_unit->get_nb_bases_in_0_coding_RNA();
    _nb_bases_in_0_non_coding_RNA       += gen_unit->get_nb_bases_in_0_non_coding_RNA();
    
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
  
  ae_list_node* gene_node_1 = NULL;
  ae_list_node* gene_node_2 = NULL;
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
    gene_1 = (ae_protein*)gene_node_1->get_obj();
    
    if ( gene_1->is_degenerated() )
    {
      gene_node_1 = gene_node_1->get_next();
      continue;
    }
    
    // A) Distance with the other LEADING genes
    gene_node_2 = gene_node_1->get_next();
    
    while ( gene_node_2 != NULL )
    {
      gene_2 = (ae_protein*)gene_node_2->get_obj();
      
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
      gene_2 = (ae_protein*)gene_node_2->get_obj();
      
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
    gene_1 = (ae_protein*)gene_node_1->get_obj();
    
    if ( gene_1->is_degenerated() )
    {
      gene_node_1 = gene_node_1->get_next();
      continue;
    }
    
    gene_node_2 = gene_node_1->get_next();
    while ( gene_node_2 != NULL )
    {
      gene_2 = (ae_protein*)gene_node_2->get_obj();
      
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

void ae_individual::write_to_backup( gzFile* backup_file )
{
  // Write index in population
  gzwrite( backup_file, &_index_in_population, sizeof(_index_in_population) );
  
  // The rank is not written, it can be interpolated because individuals are written
  // ordered by their increasing rank (<=> increasing fitness)
  
  // Write generic probes
  for ( int8_t i = 0 ; i < 5 ; i++ )
  {
    gzwrite( backup_file, &_int_probes[i],     sizeof(_int_probes[i]) );
    gzwrite( backup_file, &_double_probes[i],  sizeof(_double_probes[i]) );
  }
  
  // Write genetic units
  int16_t nb_gen_units = _genetic_unit_list->get_nb_elts();
  gzwrite( backup_file, &nb_gen_units,  sizeof(nb_gen_units) );
  
  ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit;
  for ( int16_t i = 0 ; i < nb_gen_units ; i++ )
  {
    gen_unit = (ae_genetic_unit*)gen_unit_node->get_obj();
    
    gen_unit->write_to_backup( backup_file );
    
    gen_unit_node = gen_unit_node->get_next();
  }
}

int32_t ae_individual::get_nb_terminators( void )
{
  int32_t nb_term = 0;
  
  ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit      = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();
    
    nb_term += gen_unit->get_nb_terminators();
    
    gen_unit_node = gen_unit_node->get_next();
  }
  
  return nb_term;
}


// =================================================================
//                           Protected Methods
// =================================================================
void ae_individual::make_protein_list( void )
{
  assert( _protein_list != NULL );
  
  // Clean list
  _protein_list->erase( NO_DELETE );
  
  // Make a copy of each genetic unit's protein list
  ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit      = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();
    
    _protein_list->add_list( gen_unit->get_protein_list()[LEADING] );
    _protein_list->add_list( gen_unit->get_protein_list()[LAGGING] );
    
    gen_unit_node = gen_unit_node->get_next();
  }
}

void ae_individual::make_rna_list( void )
{
  assert( _rna_list != NULL );
  
  // Clean list
  _rna_list->erase( NO_DELETE );
  
  // Make a copy of each genetic unit's rna list
  ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit      = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();
    
    _rna_list->add_list( gen_unit->get_rna_list()[LEADING] );
    _rna_list->add_list( gen_unit->get_rna_list()[LAGGING] );
    
    gen_unit_node = gen_unit_node->get_next();
  }
}
