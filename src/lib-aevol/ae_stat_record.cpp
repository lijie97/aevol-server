//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004  LIRIS.
// Web: https://liris.cnrs.fr/
// E-mail: carole.knibbe@liris.cnrs.fr
// Original Authors : Guillaume Beslon, Carole Knibbe, Virginie Lefort
//                    David Parsons
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330.0, Boston, MA  0.02111-130.07, USA.
//*****************************************************************************


/** \class
 *  \brief
 */
 
 
// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include <ae_stat_record.h>
#include <ae_list.h>
#include <ae_simulation.h>
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_genetic_unit.h>
#include <ae_replication_report.h>
#include <ae_dna_replic_report.h>





//##############################################################################
//                                                                             #
//                             Class ae_stat_record                            #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_stat_record::ae_stat_record( void )
{
  initialize_data();
}

/* If used for post-treatments, num_gener is mandatory */
ae_stat_record::ae_stat_record( ae_individual* indiv, chrom_or_gen_unit chrom_or_gu, bool compute_non_coding, int32_t num_gener )
{
  initialize_data();
  _record_type = INDIV;
    
  // ---------------
  // Simulation data
  // ---------------
  _num_gener = ( num_gener == -1 ) ? ae_common::sim->get_num_gener() : num_gener;
  _pop_size = 0; // The pop_size value is irrelevent when dealing with a single individual. It is present for column alignment.
  
  #ifdef __REGUL
    // TODO
  #endif  
    
  ae_list_node*     gen_unit_node;
  ae_genetic_unit*  gen_unit;
  
  
  // TODO : These conditions are not well managed!!!
  if ( indiv->get_genetic_unit_list()->get_nb_elts() == 1 ) // One single Genetic Unit
  {
    // -------------------------------------------------
    // Compute statistical data for the given individual
    // -------------------------------------------------
    ae_replication_report* replic_report = NULL;
    if ( _num_gener > 0 && ae_common::record_tree )
    {
      replic_report = indiv->get_replic_report();
    }
    gen_unit = (ae_genetic_unit*) indiv->get_genetic_unit_list()->get_first()->get_obj();
    
    // Metabolic error stats
    _metabolic_error = (double) indiv->get_dist_to_target_by_feature( METABOLISM );
    _parent_metabolic_error = ( replic_report != NULL ) ? replic_report->get_parent_metabolic_error() : 0.0;
    
    // Fitness
    _fitness = indiv->get_fitness();
  
    // Secretion stats
    if ( ae_common::use_secretion )
    {
       _secretion_error   = (double) indiv->get_dist_to_target_by_feature( SECRETION );
       _compound_secreted = (double) indiv->get_fitness_by_feature(SECRETION);
       _compound_amount   = (double) indiv->get_grid_cell()->get_compound_amount();
       _parent_secretion_error = 0.0;
  
      if ( replic_report != NULL )
      {
        _parent_secretion_error = replic_report->get_parent_secretion_error();
      }
    }
    else
    {
      _secretion_error   = 0.0;
      _compound_secreted = 0.0;
      _compound_amount   = 0.0;
      _parent_secretion_error = 0.0;
    }
    

    
    // Genes and RNA stats
    _amount_of_dna               = gen_unit->get_dna()->get_length();
    _nb_coding_rnas              = gen_unit->get_nb_coding_RNAs();
    _nb_non_coding_rnas          = gen_unit->get_nb_non_coding_RNAs();
    _av_size_coding_rnas         = gen_unit->get_av_size_coding_RNAs();
    _av_size_non_coding_rnas     = gen_unit->get_av_size_non_coding_RNAs();
    _nb_functional_genes         = gen_unit->get_nb_functional_genes();
    _nb_non_functional_genes     = gen_unit->get_nb_non_functional_genes();
    _av_size_functional_gene     = gen_unit->get_av_size_functional_genes();
    _av_size_non_functional_gene = gen_unit->get_av_size_non_functional_genes();


    // Non coding stats
    if ( compute_non_coding )
    {
      _nb_bases_in_0_CDS                = gen_unit->get_nb_bases_in_0_CDS();
      _nb_bases_in_0_functional_CDS     = gen_unit->get_nb_bases_in_0_functional_CDS();
      _nb_bases_in_0_non_functional_CDS = gen_unit->get_nb_bases_in_0_non_functional_CDS();
      _nb_bases_in_0_RNA                = gen_unit->get_nb_bases_in_0_RNA();
      _nb_bases_in_0_coding_RNA         = gen_unit->get_nb_bases_in_0_coding_RNA();
      _nb_bases_in_0_non_coding_RNA     = gen_unit->get_nb_bases_in_0_non_coding_RNA();
      
      _nb_bases_non_essential                     = gen_unit->get_nb_bases_non_essential();
      _nb_bases_non_essential_including_nf_genes  = gen_unit->get_nb_bases_non_essential_including_nf_genes();
    }
    
    // Mutation stats
    if ( replic_report != NULL )
    {
      _nb_mut    = gen_unit->get_dna()->get_replic_report()->get_nb_small_mutations();
      _nb_rear   = gen_unit->get_dna()->get_replic_report()->get_nb_rearrangements();
      _nb_switch = gen_unit->get_dna()->get_replic_report()->get_nb_switch();
      _nb_indels = gen_unit->get_dna()->get_replic_report()->get_nb_indels();
      _nb_dupl   = gen_unit->get_dna()->get_replic_report()->get_nb_duplications();
      _nb_del    = gen_unit->get_dna()->get_replic_report()->get_nb_deletions();
      _nb_trans  = gen_unit->get_dna()->get_replic_report()->get_nb_translocations();
      _nb_inv    = gen_unit->get_dna()->get_replic_report()->get_nb_inversions();
      
      // Rearrangement rate stats
      int32_t parent_genome_size = replic_report->get_parent_genome_size();
      _dupl_rate  = _nb_dupl  / parent_genome_size;
      _del_rate   = _nb_del   / parent_genome_size;
      _trans_rate = _nb_trans / parent_genome_size;
      _inv_rate   = _nb_inv   / parent_genome_size;
      
      //~ printf( "_nb_dupl : %"PRId32"\n_nb_del : %"PRId32"\n_nb_trans : %"PRId32"\n_nb_inv : %"PRId32"\n",
              //~ (int32_t) _nb_dupl, (int32_t) _nb_del, (int32_t) _nb_trans, (int32_t) _nb_inv );
      //~ printf( "parent genome size : %"PRId32"\n", parent_genome_size );
      //~ printf( "_dupl_rate : %lf\n_del_rate : %lf\n_trans_rate : %lf\n_inv_rate : %lf\n",
              //~ _dupl_rate, _del_rate, _trans_rate, _inv_rate );
      //~ getchar();
      
      _mean_align_score = replic_report->get_mean_align_score();
    }
  }
  else if ( chrom_or_gu == ALL_GU )
  {
    // -------------------------------------------------
    // Compute statistical data for the given individual
    // -------------------------------------------------
    ae_replication_report* replic_report = NULL;
    if ( _num_gener > 0 && ae_common::record_tree )
    {
      replic_report = indiv->get_replic_report();
    }
    
    // Metabolic error stats
    _metabolic_error = (double) indiv->get_dist_to_target_by_feature( METABOLISM );
    _parent_metabolic_error = ( replic_report != NULL ) ? replic_report->get_parent_metabolic_error() : 0.0;
    
    // Fitness
    _fitness = indiv->get_fitness();
  
    // Secretion stats
    if ( ae_common::use_secretion )
    {
       _secretion_error = (double) indiv->get_dist_to_target_by_feature( SECRETION );
       _compound_secreted = (double) indiv->get_fitness_by_feature(SECRETION);
       _compound_amount   = (double) indiv->get_grid_cell()->get_compound_amount();
       _parent_secretion_error = 0.0;
  
      if ( _num_gener > 0 && ae_common::record_tree )
      {
        _parent_secretion_error = replic_report->get_parent_secretion_error();
      }
    }
    else
    {
      _secretion_error   = 0.0;
      _compound_secreted = 0.0;
      _compound_amount   = 0.0;
      _parent_secretion_error = 0.0;
    }

    
    
    gen_unit_node = indiv->get_genetic_unit_list()->get_first();
    while ( gen_unit_node != NULL )
    {
      gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();

      // Genes and RNA stats
      _amount_of_dna               += gen_unit->get_dna()->get_length();
      _nb_coding_rnas              += gen_unit->get_nb_coding_RNAs();
      _nb_non_coding_rnas          += gen_unit->get_nb_non_coding_RNAs();
      _av_size_coding_rnas         += gen_unit->get_av_size_coding_RNAs();
      _av_size_non_coding_rnas     += gen_unit->get_av_size_non_coding_RNAs();
      _nb_functional_genes         += gen_unit->get_nb_functional_genes();
      _nb_non_functional_genes     += gen_unit->get_nb_non_functional_genes();
      _av_size_functional_gene     += gen_unit->get_av_size_functional_genes();
      _av_size_non_functional_gene += gen_unit->get_av_size_non_functional_genes();

      // Non coding stats
      if ( compute_non_coding )
      {
        _nb_bases_in_0_CDS                += gen_unit->get_nb_bases_in_0_CDS();
        _nb_bases_in_0_functional_CDS     += gen_unit->get_nb_bases_in_0_functional_CDS();
        _nb_bases_in_0_non_functional_CDS += gen_unit->get_nb_bases_in_0_non_functional_CDS();
        _nb_bases_in_0_RNA                += gen_unit->get_nb_bases_in_0_RNA();
        _nb_bases_in_0_coding_RNA         += gen_unit->get_nb_bases_in_0_coding_RNA();
        _nb_bases_in_0_non_coding_RNA     += gen_unit->get_nb_bases_in_0_non_coding_RNA();
        
        _nb_bases_non_essential                     += gen_unit->get_nb_bases_non_essential();
        _nb_bases_non_essential_including_nf_genes  += gen_unit->get_nb_bases_non_essential_including_nf_genes();
      }
      
      // Mutation stats
      if ( replic_report != NULL )
      {
        _nb_mut    += gen_unit->get_dna()->get_replic_report()->get_nb_small_mutations();
        _nb_rear   += gen_unit->get_dna()->get_replic_report()->get_nb_rearrangements();
        _nb_switch += gen_unit->get_dna()->get_replic_report()->get_nb_switch();
        _nb_indels += gen_unit->get_dna()->get_replic_report()->get_nb_indels();
        _nb_dupl   += gen_unit->get_dna()->get_replic_report()->get_nb_duplications();
        _nb_del    += gen_unit->get_dna()->get_replic_report()->get_nb_deletions();
        _nb_trans  += gen_unit->get_dna()->get_replic_report()->get_nb_translocations();
        _nb_inv    += gen_unit->get_dna()->get_replic_report()->get_nb_inversions();
      }

     gen_unit_node = gen_unit_node->get_next();
    }
    
    // Rearrangement rate stats
    if ( replic_report != NULL )
    {
      int32_t parent_genome_size = replic_report->get_parent_genome_size();
      _dupl_rate  = _nb_dupl  / parent_genome_size;
      _del_rate   = _nb_del   / parent_genome_size;
      _trans_rate = _nb_trans / parent_genome_size;
      _inv_rate   = _nb_inv   / parent_genome_size;
      _mean_align_score = replic_report->get_mean_align_score();
    }
  }
  else // => We have a multi-GU individual and we want only the main chromosome or only the plasmids
  // WARNING (TODO) As it is coded, this will work only if there is ONE SINGLE PLASMID!
  {
    if ( chrom_or_gu == PLASMIDS )
    {
      gen_unit_node = indiv->get_genetic_unit_list()->get_first()->get_next();
    }
    if ( chrom_or_gu == CHROM) 
    {
      gen_unit_node = indiv->get_genetic_unit_list()->get_first();
    }

    gen_unit = (ae_genetic_unit*)gen_unit_node->get_obj();
    
    // -------------------------------------------------
    // Compute statistical data for the given individual
    // -------------------------------------------------
    ae_replication_report*  replic_report = NULL;
    if ( _num_gener > 0 && ae_common::record_tree )
    {
      replic_report = indiv->get_replic_report();
    }
    
    // Metabolic error stats
    _metabolic_error = (double) indiv->get_dist_to_target_by_feature( METABOLISM );
    _parent_metabolic_error = ( replic_report != NULL ) ? replic_report->get_parent_metabolic_error() : 0.0;
    
    // Fitness
    _fitness = indiv->get_fitness();
  
    // Secretion stats
    if ( ae_common::use_secretion )
    {
       _secretion_error = (double) gen_unit->get_dist_to_target_by_feature( SECRETION );
       _compound_secreted = (double) gen_unit->get_fitness_by_feature(SECRETION);
       _compound_amount   = (double) indiv->get_grid_cell()->get_compound_amount();
       _parent_secretion_error = 0.0;
  
      if ( _num_gener > 0 && ae_common::record_tree )
      {
        _parent_secretion_error = replic_report->get_parent_secretion_error();
      }
    }
    else
    {
      _secretion_error   = 0.0;
      _compound_secreted = 0.0;
      _compound_amount   = 0.0;
      _parent_secretion_error = 0.0;
    }
    
      // Genes and RNA stats
    _amount_of_dna               = gen_unit->get_dna()->get_length();
    _nb_coding_rnas              = gen_unit->get_nb_coding_RNAs();
    _nb_non_coding_rnas          = gen_unit->get_nb_non_coding_RNAs();
    _av_size_coding_rnas         = gen_unit->get_av_size_coding_RNAs();
    _av_size_non_coding_rnas     = gen_unit->get_av_size_non_coding_RNAs();
    _nb_functional_genes         = gen_unit->get_nb_functional_genes();
    _nb_non_functional_genes     = gen_unit->get_nb_non_functional_genes();
    _av_size_functional_gene     = gen_unit->get_av_size_functional_genes();
    _av_size_non_functional_gene = gen_unit->get_av_size_non_functional_genes();
    
      // Non coding stats
    if ( compute_non_coding )
    {
      _nb_bases_in_0_CDS                  = gen_unit->get_nb_bases_in_0_CDS();
      _nb_bases_in_0_functional_CDS       = gen_unit->get_nb_bases_in_0_functional_CDS();
      _nb_bases_in_0_non_functional_CDS   = gen_unit->get_nb_bases_in_0_non_functional_CDS();
      _nb_bases_in_0_RNA                  = gen_unit->get_nb_bases_in_0_RNA();
      _nb_bases_in_0_coding_RNA           = gen_unit->get_nb_bases_in_0_coding_RNA();
      _nb_bases_in_0_non_coding_RNA       = gen_unit->get_nb_bases_in_0_non_coding_RNA();
      
      _nb_bases_non_essential                     = gen_unit->get_nb_bases_non_essential();
      _nb_bases_non_essential_including_nf_genes  = gen_unit->get_nb_bases_non_essential_including_nf_genes();
    }
    
    // Mutation stats
    if ( _num_gener > 0 && ae_common::record_tree )
    {
      _nb_mut    = gen_unit->get_dna()->get_replic_report()->get_nb_small_mutations();
      _nb_rear   = gen_unit->get_dna()->get_replic_report()->get_nb_rearrangements();
      _nb_switch = gen_unit->get_dna()->get_replic_report()->get_nb_switch();
      _nb_indels = gen_unit->get_dna()->get_replic_report()->get_nb_indels();
      _nb_dupl   = gen_unit->get_dna()->get_replic_report()->get_nb_duplications();
      _nb_del    = gen_unit->get_dna()->get_replic_report()->get_nb_deletions();
      _nb_trans  = gen_unit->get_dna()->get_replic_report()->get_nb_translocations();
      _nb_inv    = gen_unit->get_dna()->get_replic_report()->get_nb_inversions();
    }
    
    // Rearrangement rate stats
    if ( replic_report != NULL )
    {
      int32_t parent_genome_size = replic_report->get_parent_genome_size();
      _dupl_rate  = _nb_dupl  / parent_genome_size;
      _del_rate   = _nb_del   / parent_genome_size;
      _trans_rate = _nb_trans / parent_genome_size;
      _inv_rate   = _nb_inv   / parent_genome_size;
      _mean_align_score = replic_report->get_mean_align_score();
    }
  }
}

ae_stat_record::ae_stat_record( ae_population* pop, chrom_or_gen_unit chrom_or_gu )
{
  initialize_data();
  
  _record_type = POP;
  
  // ---------------
  // Simulation data
  // ---------------
  _num_gener = (double) ae_common::sim->get_num_gener();
  _pop_size  = (double) pop->get_nb_indivs();

  // ------------------------------------------------------------------
  // Compute statistical data for the each individual in the population
  // ------------------------------------------------------------------
  ae_list_node*   indiv_node  = pop->get_indivs()->get_first();
  ae_individual*  indiv       = NULL;

  while ( indiv_node != NULL )
  {
    indiv = (ae_individual*) indiv_node->get_obj();
    ae_stat_record* indiv_stat_record = new ae_stat_record( indiv, chrom_or_gu, false );
    this->add( indiv_stat_record );
    delete indiv_stat_record;
    
    indiv_node = indiv_node->get_next();
  }
  
  
  // ------------------------------------------------------------------
  // Divide every accumulator by the number of indivs in the population
  // ------------------------------------------------------------------
  this->divide( _pop_size );
}

ae_stat_record::ae_stat_record( const ae_stat_record &model )
{
  _num_gener = model._num_gener;
  _pop_size  = model._pop_size;
  
  _metabolic_error         = model._metabolic_error;
  _parent_metabolic_error  = model._parent_metabolic_error;
  
  _secretion_error         = model._secretion_error;
  _parent_secretion_error  = model._parent_secretion_error;
  
  _compound_secreted = model._compound_secreted;
  _compound_amount   = model._compound_amount;
  
  _fitness = model._fitness;
  
  _amount_of_dna                = model._amount_of_dna;
  _nb_coding_rnas               = model._nb_coding_rnas;
  _nb_non_coding_rnas           = model._nb_non_coding_rnas;
  _av_size_coding_rnas          = model._av_size_coding_rnas;
  _av_size_non_coding_rnas      = model._av_size_non_coding_rnas;
  _nb_functional_genes          = model._nb_functional_genes;
  _nb_non_functional_genes      = model._nb_non_functional_genes;
  _av_size_functional_gene      = model._av_size_functional_gene;
  _av_size_non_functional_gene  = model._av_size_non_functional_gene;

  _nb_mut    = model._nb_mut;
  _nb_rear   = model._nb_rear;
  _nb_switch = model._nb_switch;
  _nb_indels = model._nb_indels;
  _nb_dupl   = model._nb_dupl;
  _nb_del    = model._nb_del;
  _nb_trans  = model._nb_trans;
  _nb_inv    = model._nb_inv;
  
  _nb_bases_in_0_CDS                = model._nb_bases_in_0_CDS;
  _nb_bases_in_0_functional_CDS     = model._nb_bases_in_0_functional_CDS;
  _nb_bases_in_0_non_functional_CDS = model._nb_bases_in_0_non_functional_CDS;
  _nb_bases_in_0_RNA                = model._nb_bases_in_0_RNA;
  _nb_bases_in_0_coding_RNA         = model._nb_bases_in_0_coding_RNA;
  _nb_bases_in_0_non_coding_RNA     = model._nb_bases_in_0_non_coding_RNA;
      
  _nb_bases_non_essential                     = model._nb_bases_non_essential;
  _nb_bases_non_essential_including_nf_genes  = model._nb_bases_non_essential_including_nf_genes;
    
  #ifdef __REGUL
    _nb_influences                 = model._nb_influences;
    _nb_enhancing_influences       = model._nb_enhancing_influences;
    _nb_operating_influences       = model._nb_operating_influences;
    _av_value_influences           = model._av_value_influences;
    _av_value_enhancing_influences = model._av_value_enhancing_influences;
    _av_value_operating_influences = model._av_value_operating_influences;
  #endif
}

// =================================================================
//                             Destructors
// =================================================================
ae_stat_record::~ae_stat_record( void )
{
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_stat_record::initialize_data( void )
{
  _num_gener = 0.0;
  _pop_size  = 0.0;
  
  _metabolic_error         = 0.0;
  _parent_metabolic_error  = 0.0;

  _secretion_error         = 0.0;
  _parent_secretion_error  = 0.0;
  
  _compound_secreted = 0.0;
  _compound_amount   = 0.0;
  
  _fitness = 0.0;
  
  _amount_of_dna               = 0.0;
  _nb_coding_rnas              = 0.0;
  _nb_non_coding_rnas          = 0.0;
  _av_size_coding_rnas         = 0.0;
  _av_size_non_coding_rnas     = 0.0;
  _nb_functional_genes         = 0.0;
  _nb_non_functional_genes     = 0.0;
  _av_size_functional_gene     = 0.0;
  _av_size_non_functional_gene = 0.0;

  _nb_mut    = 0.0;
  _nb_rear   = 0.0;
  _nb_switch = 0.0;
  _nb_indels = 0.0;
  _nb_dupl   = 0.0;
  _nb_del    = 0.0;
  _nb_trans  = 0.0;
  _nb_inv    = 0.0;
  
  _dupl_rate  = 0.0;
  _del_rate   = 0.0;
  _trans_rate = 0.0;
  _inv_rate   = 0.0;
  _mean_align_score = 0.0;
  
  _nb_bases_in_0_CDS                = 0.0;
  _nb_bases_in_0_functional_CDS     = 0.0;
  _nb_bases_in_0_non_functional_CDS = 0.0;
  _nb_bases_in_0_RNA                = 0.0;
  _nb_bases_in_0_coding_RNA         = 0.0;
  _nb_bases_in_0_non_coding_RNA     = 0.0;
    
  _nb_bases_non_essential                     = 0.0;
  _nb_bases_non_essential_including_nf_genes  = 0.0;
    
  #ifdef __REGUL
    _nb_influences                 = 0.0;
    _nb_enhancing_influences       = 0.0;
    _nb_operating_influences       = 0.0;
    _av_value_influences           = 0.0;
    _av_value_enhancing_influences = 0.0;
    _av_value_operating_influences = 0.0;
  #endif
}

void ae_stat_record::write_to_file( FILE* stat_file, stats_type stat_type_to_print) const
{
  if ( _record_type == INDIV )
  {
    if ( stat_type_to_print == FITNESS_STATS )
    {
      fprintf( stat_file, "%"PRId32" %"PRId32" %lf %lf %lf %lf %lf %lf %lf %lf", 
              (int32_t) _num_gener,
              (int32_t) _pop_size,
              _fitness,              
              _amount_of_dna, 
              _metabolic_error,
              _parent_metabolic_error,
              _secretion_error,
              _parent_secretion_error,
              _compound_secreted,
              _compound_amount );

      #ifdef __REGUL
        fprintf(  stat_file, " %"PRId32" %"PRId32" %"PRId32" %lf %lf %lf",
              (int32_t) _nb_influences,
              (int32_t) _nb_enhancing_influences,
              (int32_t) _nb_operating_influences,
              _av_value_influences,
              _av_value_enhancing_influences,
              _av_value_operating_influences );
      #endif
    }
    if ( stat_type_to_print == MUTATION_STATS )
    {
      fprintf(  stat_file, "%"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32"", 
              (int32_t) _num_gener,
              (int32_t) _nb_mut,
              (int32_t) _nb_rear,
              (int32_t) _nb_switch,
              (int32_t) _nb_indels,
              (int32_t) _nb_dupl,
              (int32_t) _nb_del,
              (int32_t) _nb_trans,
              (int32_t) _nb_inv );

    }
    if ( stat_type_to_print == GENES_STATS )
    {
      fprintf(  stat_file, "%"PRId32" %"PRId32" %"PRId32" %lf %lf %"PRId32" %"PRId32" %lf %lf ",
              (int32_t) _num_gener,
              (int32_t) _nb_coding_rnas,
              (int32_t) _nb_non_coding_rnas,
              _av_size_coding_rnas,
              _av_size_non_coding_rnas,
              (int32_t) _nb_functional_genes,
              (int32_t) _nb_non_functional_genes,
              _av_size_functional_gene,
              _av_size_non_functional_gene );
    }
    if ( stat_type_to_print == BP_STATS )
    {
      fprintf(  stat_file, "%"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32"",
              (int32_t) _num_gener,
              (int32_t) _nb_bases_in_0_CDS,
              (int32_t) _nb_bases_in_0_functional_CDS,
              (int32_t) _nb_bases_in_0_non_functional_CDS,
              (int32_t) _nb_bases_in_0_RNA,
              (int32_t) _nb_bases_in_0_coding_RNA,
              (int32_t) _nb_bases_in_0_non_coding_RNA,
              (int32_t) _nb_bases_non_essential,
              (int32_t) _nb_bases_non_essential_including_nf_genes );
    }
    if ( stat_type_to_print == REAR_STATS )
    {
      fprintf(  stat_file, "%"PRId32" %lf %lf %lf %lf %lf",
              (int32_t) _num_gener,
              _dupl_rate,
              _del_rate,
              _trans_rate,
              _inv_rate,
              _mean_align_score );
    }
  }
  else // if _record_type == POP
  {
   if ( stat_type_to_print == FITNESS_STATS )
    {
      fprintf(  stat_file, "%"PRId32" %"PRId32" %lf %lf %lf %lf %lf %lf %lf %lf", 
              (int32_t) _num_gener,
              (int32_t) _pop_size,
              _fitness,              
              _amount_of_dna, 
              _metabolic_error,
              _parent_metabolic_error,
              _secretion_error,
              _parent_secretion_error,
              _compound_secreted,
              _compound_amount);

      #ifdef __REGUL
        fprintf(  stat_file, " %"PRId32" %"PRId32" %"PRId32" %lf %lf %lf",
              (int32_t) _nb_influences,
              (int32_t) _nb_enhancing_influences,
              (int32_t) _nb_operating_influences,
              _av_value_influences,
              _av_value_enhancing_influences,
              _av_value_operating_influences );
      #endif
    }
    if ( stat_type_to_print == MUTATION_STATS )
    {        
      fprintf(  stat_file, "%"PRId32" %lf %lf %lf %lf %lf %lf %lf %lf",
              (int32_t) _num_gener,
              _nb_mut,
              _nb_rear,
              _nb_switch,
              _nb_indels,
              _nb_dupl,
              _nb_del,
              _nb_trans,
              _nb_inv );

    }
    if ( stat_type_to_print == GENES_STATS )
    {
      fprintf(  stat_file, "%"PRId32" %lf %lf %lf %lf %lf %lf %lf %lf",
              (int32_t) _num_gener,
              _nb_coding_rnas,
              _nb_non_coding_rnas,
              _av_size_coding_rnas,
              _av_size_non_coding_rnas,
              _nb_functional_genes,
              _nb_non_functional_genes,
              _av_size_functional_gene,
              _av_size_non_functional_gene );
    }
    if ( stat_type_to_print == BP_STATS )
    {
     // TO DO (if needed) : base-pair stats for all individuals, not just for the best one. 
     //
     // fprintf(  stat_file, "%"PRId32" %lf %lf %lf %lf %lf %lf %lf %lf",
     //         (int32_t)_num_gener,
     //         _nb_bases_in_0_CDS,
     //         _nb_bases_in_0_functional_CDS,
     //         _nb_bases_in_0_non_functional_CDS,
     //         _nb_bases_in_0_RNA,
     //         _nb_bases_in_0_coding_RNA,
     //         _nb_bases_in_0_non_coding_RNA,
     //         _nb_bases_non_essential,
     //         _nb_bases_non_essential_including_nf_genes );
    }
    if ( stat_type_to_print == REAR_STATS )
    {
      fprintf(  stat_file, "%"PRId32" %lf %lf %lf %lf %lf",
              (int32_t) _num_gener,
              _dupl_rate,
              _del_rate,
              _trans_rate,
              _inv_rate,
              _mean_align_score );
    }
  }
  
  fprintf( stat_file, "\n" );
}

void ae_stat_record::divide( double divisor )
{
  // NB : _num_gener and pop_size are global values and are not to be divided.
  
  _fitness                 /= divisor;
  
  _metabolic_error         /= divisor;
  _parent_metabolic_error  /= divisor;

  _secretion_error         /= divisor;
  _parent_secretion_error  /= divisor;
  
  _compound_secreted       /= divisor;
  _compound_amount         /= divisor;
  
  _amount_of_dna               /= divisor;
  _nb_coding_rnas              /= divisor;
  _nb_non_coding_rnas          /= divisor;
  _av_size_coding_rnas         /= divisor;
  _av_size_non_coding_rnas     /= divisor;
  _nb_functional_genes          /= divisor;
  _nb_non_functional_genes      /= divisor;
  _av_size_functional_gene      /= divisor;
  _av_size_non_functional_gene  /= divisor;

  _nb_mut    /= divisor;
  _nb_rear   /= divisor;
  _nb_switch /= divisor;
  _nb_indels /= divisor;
  _nb_dupl   /= divisor;
  _nb_del    /= divisor;
  _nb_trans  /= divisor;
  _nb_inv    /= divisor;
  
  //~ printf( "PREFINAL %lf %lf %lf %lf\n", _dupl_rate, _del_rate, _trans_rate, _inv_rate );
  _dupl_rate  /= divisor;
  _del_rate   /= divisor;
  _trans_rate /= divisor;
  _inv_rate   /= divisor;
  //~ printf( "FINAL %lf %lf %lf %lf\n", _dupl_rate, _del_rate, _trans_rate, _inv_rate );
  //~ getchar();
  _mean_align_score /= divisor;
  
  _nb_bases_in_0_CDS                /= divisor;
  _nb_bases_in_0_functional_CDS      /= divisor;
  _nb_bases_in_0_non_functional_CDS  /= divisor;
  _nb_bases_in_0_RNA                /= divisor;
  _nb_bases_in_0_coding_RNA         /= divisor;
  _nb_bases_in_0_non_coding_RNA     /= divisor;
    
  _nb_bases_non_essential                     /= divisor;
  _nb_bases_non_essential_including_nf_genes  /= divisor;
    
  #ifdef __REGUL
    _nb_influences                 /= divisor;
    _nb_enhancing_influences       /= divisor;
    _nb_operating_influences       /= divisor;
    _av_value_influences           /= divisor;
    _av_value_enhancing_influences /= divisor;
    _av_value_operating_influences /= divisor;
  #endif
}

void ae_stat_record::add( ae_stat_record* to_add )
{
  // NB : _num_gener and pop_size are global values and are not to be summed.
  
  _fitness                 += to_add->_fitness;
  
  _metabolic_error         += to_add->_metabolic_error;
  _parent_metabolic_error  += to_add->_parent_metabolic_error;

  _secretion_error         += to_add->_secretion_error;
  _parent_secretion_error  += to_add->_parent_secretion_error;
  
  _compound_secreted       += to_add->_compound_secreted;
  _compound_amount         += to_add->_compound_amount;
  
  _amount_of_dna               += to_add->_amount_of_dna;
  _nb_coding_rnas              += to_add->_nb_coding_rnas;
  _nb_non_coding_rnas          += to_add->_nb_non_coding_rnas;
  _av_size_coding_rnas         += to_add->_av_size_coding_rnas;
  _av_size_non_coding_rnas     += to_add->_av_size_non_coding_rnas;
  _nb_functional_genes         += to_add->_nb_functional_genes;
  _nb_non_functional_genes     += to_add->_nb_non_functional_genes;
  _av_size_functional_gene     += to_add->_av_size_functional_gene;
  _av_size_non_functional_gene += to_add->_av_size_non_functional_gene;

  _nb_mut    += to_add->_nb_mut;
  _nb_rear   += to_add->_nb_rear;
  _nb_switch += to_add->_nb_switch;
  _nb_indels += to_add->_nb_indels;
  _nb_dupl   += to_add->_nb_dupl;
  _nb_del    += to_add->_nb_del;
  _nb_trans  += to_add->_nb_trans;
  _nb_inv    += to_add->_nb_inv;
    
  _dupl_rate  += to_add->_dupl_rate;
  _del_rate   += to_add->_del_rate;
  _trans_rate += to_add->_trans_rate;
  _inv_rate   += to_add->_inv_rate;
  //~ printf( "%lf %lf %lf %lf\n", to_add->_dupl_rate, to_add->_del_rate, to_add->_trans_rate, to_add->_inv_rate );
  _mean_align_score += to_add->_mean_align_score;
  
  _nb_bases_in_0_CDS                += to_add->_nb_bases_in_0_CDS;
  _nb_bases_in_0_functional_CDS     += to_add->_nb_bases_in_0_functional_CDS;
  _nb_bases_in_0_non_functional_CDS += to_add->_nb_bases_in_0_non_functional_CDS;
  _nb_bases_in_0_RNA                += to_add->_nb_bases_in_0_RNA;
  _nb_bases_in_0_coding_RNA         += to_add->_nb_bases_in_0_coding_RNA;
  _nb_bases_in_0_non_coding_RNA     += to_add->_nb_bases_in_0_non_coding_RNA;
    
  _nb_bases_non_essential                     += to_add->_nb_bases_non_essential;
  _nb_bases_non_essential_including_nf_genes  += to_add->_nb_bases_non_essential_including_nf_genes;
    
  #ifdef __REGUL
    _nb_influences                 += to_add->_nb_influences;
    _nb_enhancing_influences       += to_add->_nb_enhancing_influences;
    _nb_operating_influences       += to_add->_nb_operating_influences;
    _av_value_influences           += to_add->_av_value_influences;
    _av_value_enhancing_influences += to_add->_av_value_enhancing_influences;
    _av_value_operating_influences += to_add->_av_value_operating_influences;
  #endif
}

// =================================================================
//                           Protected Methods
// =================================================================
