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
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//*****************************************************************************

/** \class
 *  \brief
 */

#ifndef __AE_COMMON_H__
#define  __AE_COMMON_H__

// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>

// =================================================================
//                            Project Files
// =================================================================
#include <ae_enums.h>
#include <ae_macros.h>
#include <ae_rand_mt.h>
#ifdef __REGUL
  #include <ae_array_short.h>
  #include <ae_matrix_double.h>
#endif

// =================================================================
//                          Class declarations
// =================================================================
class ae_simulation;
#ifdef __X11
  class ae_simulation_X11;
#endif
class ae_list;






class ae_common
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================

    // =================================================================
    //                             Destructors
    // =================================================================

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    // Backup functions
    static void write_to_backup( gzFile* backup_file );
    static void read_from_backup( gzFile* backup_file, bool verbose = false );
    static void print_to_file( void );

    // Binding functions
    static          void    init_binding_matrix( void );
    static  inline  double  get_binding( int16_t i, int16_t j );
    
    static void clean( void ); // Release memory

    // =================================================================
    //                           Public Attributes
    // =================================================================
    static ae_simulation* sim;

#ifdef __X11
    static ae_simulation_X11* sim_display;
#endif

    // PseudoRandom Number Generator
    static int32_t seed;
    static int32_t env_seed; // Seed for the environment's random generator (environmental variation)

    // Initial conditions
    static int32_t  initial_genome_length;
    static int32_t  min_genome_length;
    static int32_t  max_genome_length;
    static int8_t   init_method;
    static int32_t  nb_generations;
    static double   min_w;
    static double   max_w;
    
    // Statistics collection
    static int32_t      backup_step;
    static int32_t      big_backup_step;
    static int32_t      tree_step;
    static bool         record_tree;
    static ae_tree_mode tree_mode;    
    static bool         more_stats;   // TODO : explain
    static int32_t      dump_period;  // TODO : explain

    // Population  size, structure, and other properties
    static int32_t  init_pop_size;
    static bool     pop_structure;
    static int16_t  grid_x; 
    static int16_t  grid_y;
    static int32_t  migration_number; // TODO : explain

    // Mutation rates
    static double   point_mutation_rate;
    static double   small_insertion_rate;
    static double   small_deletion_rate;
    static int16_t  max_indel_size;




    // Rearrangements and Transfer
    static bool   with_4pts_trans;
    static bool   with_alignments;
    static bool   with_transfer;
    static double transfer_ins_rate;
    static double transfer_repl_rate;

    // Rearrangement rates (without alignements)
    static double duplication_rate;
    static double deletion_rate;
    static double translocation_rate;
    static double inversion_rate;

    // Rearrangement rates (with alignements)
    static double neighbourhood_rate;
    static double duplication_proportion;
    static double deletion_proportion;
    static double translocation_proportion;
    static double inversion_proportion;

    // Alignements
    static ae_align_fun_shape align_fun_shape;
    static double  align_sigm_lambda;
    static int16_t align_sigm_mean;
    static int16_t align_lin_min;
    static int16_t align_lin_max;
    
    static int16_t align_max_shift;     // Maximum shift of one seq on the other
    static int16_t align_w_zone_h_len;  // Work zone half length
    static int16_t align_match_bonus;   // Corresponding residues match bonus
    static int16_t align_mismatch_cost; // Corresponding residues mismatch cost




    // Selection
    static ae_selection_scheme  selection_scheme;
    static double               selection_pressure;
    
    // Environment
    static ae_list    env_gaussians;
    static ae_list    env_custom_points;
    static int16_t    env_sampling;
    
    static ae_env_var env_var_method;
    static double     env_sigma;
    static int32_t    env_tau;
    
    static bool                 env_axis_is_segmented;        // true if environment is segmented    
    static int16_t              env_axis_nb_segments;
    static double*              env_axis_segment_boundaries;  // (sorted) including MIN_X and MAX_X
    static ae_env_axis_feature* env_axis_features;            // Corresponding features
    static bool                 env_separate_segments;        // Automatically separate segments
    
    static bool                 composite_fitness;            // true if the are multiple components of fitness,
                                                              // i.e. 2 or more non-neutral features
    
    // Secretion
    static bool   use_secretion;
    static double secretion_fitness_contrib;      // (0,1)
    // proportion that difuses into each cell, every generation, 0 for no difusion
    static double secretion_difusion_prop;
    // proportion of secreted substance that degrades every generation generation
    static double secretion_degradation_prop;
    // starting configuration of secretion grid; 0, all are 0; 1, point source of secreted compund
    static double secretion_init;
    // cost of secreting the compound, as a proportion of the amount secreted
    static double secretion_cost;
    
    
    // Plasmids
    static bool     allow_plasmids;
    static int32_t  plasmid_initial_length;
    static int32_t  plasmid_initial_gene;
    static int32_t  plasmid_minimal_length;
    static double   prob_horiz_trans;
    static bool     compute_phen_contrib_by_GU;
    
    
    // Translation cost
    static double translation_cost;
    
    
    // Log files
    static int8_t logs;

#ifdef __REGUL
    // Regulation factors
    static  double  hill_shape_n;
    static  double  hill_shape_theta;
    static  double  hill_shape;

    // Degradation equation
    static  double  degradation_rate;
    static  double  degradation_step;

    // Individual life
    static  int16_t         individual_evaluation_nbr;
    static  ae_array_short* individual_evaluation_dates;

    // Binding matrix
    static  double  binding_zeros_percentage;

    // Heredity
    static  bool    with_heredity;
    static  double  protein_presence_limit;
#endif

  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_common( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_common( const ae_common &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
#ifdef __REGUL
    // Binding matrix
    static  ae_matrix_double* _binding_matrix;
#endif

};

// =====================================================================
//                          Accessors definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
#ifdef __REGUL
double ae_common::get_binding( int16_t i, int16_t j )
{
 return ae_common::_binding_matrix->get_value( i, j );
}
#endif

#endif // __AE_COMMON_H__
