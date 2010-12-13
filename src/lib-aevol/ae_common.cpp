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


// =================================================================
//                              Libraries
// =================================================================
#include <stdlib.h>
#include <math.h>

// =================================================================
//                            Project Files
// =================================================================
#include <ae_common.h>
#include <ae_simulation.h>
#include <ae_list.h>
#include <ae_gaussian.h>
#include <ae_point_2d.h>

//##############################################################################
//                                                                             #
//                               Class ae_common                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

ae_simulation*  ae_common::sim  = NULL;

#ifdef __X11
ae_simulation_X11* ae_common::sim_display = NULL;
#endif

// PseudoRandom Number Generator
int32_t ae_common::seed      = 0;
int32_t ae_common::env_seed  = 0;


// Initial conditions
int32_t ae_common::initial_genome_length  = 5000;
int32_t ae_common::min_genome_length      = 10;
int32_t ae_common::max_genome_length      = 10000000;
int8_t  ae_common::init_method            = ONE_GOOD_GENE | CLONE;
int32_t ae_common::nb_generations         = 1000000;
double  ae_common::min_w                  = 0.0;
double  ae_common::max_w                  = 0.033333333;

// Statistics collection
int32_t       ae_common::backup_step      = 500;
int32_t       ae_common::big_backup_step  = 10000;
int32_t       ae_common::tree_step        = 100;
bool          ae_common::record_tree      = true;
ae_tree_mode  ae_common::tree_mode        = NORMAL;
bool          ae_common::more_stats       = false;
int32_t       ae_common::dump_period      = 0;

// Population size, structure, and other properties
int32_t ae_common::init_pop_size     = 1000;
bool    ae_common::pop_structure     = false;
int16_t ae_common::grid_x            = 25; 
int16_t ae_common::grid_y            = 40; 
int32_t ae_common::migration_number  = 0; 

// Mutation rates
double  ae_common::point_mutation_rate    = 1e-5;
double  ae_common::small_insertion_rate   = 1e-5;
double  ae_common::small_deletion_rate    = 1e-5;
int16_t ae_common::max_indel_size         = 6;

// Rearrangements and Transfer
bool ae_common::with_4pts_rears   = true;
bool ae_common::with_alignments   = false;
bool ae_common::with_transfer     = false;

// Rearrangement rates (without alignements)
double ae_common::duplication_rate      = 5e-5;
double ae_common::deletion_rate         = 5e-5;
double ae_common::translocation_rate    = 5e-5;
double ae_common::inversion_rate        = 5e-5;

// Rearrangement rates (with alignements)
double ae_common::neighbourhood_rate        = 5e-5;
double ae_common::duplication_proportion    = 0.3;
double ae_common::deletion_proportion       = 0.3;
double ae_common::translocation_proportion  = 0.3;
double ae_common::inversion_proportion      = 0.3;

// Alignements
int16_t ae_common::align_min_score      = 0;
int16_t ae_common::align_max_score      = 100;
int16_t ae_common::align_max_shift      = 20; // Maximum shift of one seq on the other
int16_t ae_common::align_w_zone_h_len   = 50; // Work zone half length
int16_t ae_common::align_match_bonus    = 1;  // Corresponding residues match bonus
int16_t ae_common::align_mismatch_cost  = 2;  // Corresponding residues mismatch cost

// Selection
ae_selection_scheme ae_common::selection_scheme    = RANK_EXPONENTIAL;
double              ae_common::selection_pressure  = 0.998;
    
// Environment
ae_list     ae_common::env_gaussians;
ae_list     ae_common::env_custom_points;
int16_t     ae_common::env_sampling = 300;

// Environment variation
ae_env_var  ae_common::env_var_method = NONE;
double      ae_common::env_sigma      = 0.01;
int32_t     ae_common::env_tau        = 1000;

// Environment x-axis segmentation
bool                  ae_common::env_axis_is_segmented        = false;
int16_t               ae_common::env_axis_nb_segments         = 1;
double*               ae_common::env_axis_segment_boundaries  = NULL;
ae_env_axis_feature*  ae_common::env_axis_features            = NULL;
bool                  ae_common::env_separate_segments        = false;

bool                  ae_common::composite_fitness            = false;


// Secretion
bool    ae_common::use_secretion                = false;
double  ae_common::secretion_fitness_contrib    = 0;
double  ae_common::secretion_difusion_prop      = 0;
double  ae_common::secretion_degradation_prop   = 0;
double  ae_common::secretion_init               = 0;
double  ae_common::secretion_cost               = 0;

// Plasmids
bool    ae_common::allow_plasmids             = false;
int32_t ae_common::plasmid_initial_length     = 1000;
int32_t ae_common::plasmid_initial_gene       = 0;
int32_t ae_common::plasmid_minimal_length     = 40;
double  ae_common::prob_horiz_trans           = 0;
bool    ae_common::compute_phen_contrib_by_GU = false;

// Translation cost
double ae_common::translation_cost  = 0;

// Log files
int8_t ae_common::logs = 0;


#ifdef __REGUL
// Regulation factors
double ae_common::hill_shape_n      = 4;
double ae_common::hill_shape_theta  = 0.5;
double ae_common::hill_shape        = pow( hill_shape_theta, hill_shape_n );

// Degradation equation
double ae_common::degradation_rate  = 1;
double ae_common::degradation_step  = 0.1;

// Individual life
int16_t         ae_common::individual_evaluation_nbr   = 1;
ae_array_short* ae_common::individual_evaluation_dates = NULL;

// Binding matrix
double            ae_common::binding_zeros_percentage = 0.75;
ae_matrix_double* ae_common::_binding_matrix          = NULL;

// Heredity
bool    ae_common::with_heredity          = false;
double  ae_common::protein_presence_limit = 1e-5;
#endif

// =================================================================
//                             Constructors
// =================================================================

// =================================================================
//                             Destructors
// =================================================================

// =================================================================
//                            Public Methods
// =================================================================
void ae_common::write_to_backup( gzFile* backup_file )
{
  // PseudoRandom Number Generator
  gzwrite( backup_file, &seed,                        sizeof(seed)                      );
  gzwrite( backup_file, &env_seed,                    sizeof(env_seed)                  );

  // Initial conditions
  gzwrite( backup_file, &initial_genome_length,       sizeof(initial_genome_length)     );
  gzwrite( backup_file, &min_genome_length,           sizeof(min_genome_length)         );
  gzwrite( backup_file, &max_genome_length,           sizeof(max_genome_length)         );
  gzwrite( backup_file, &init_method,                 sizeof(init_method)               );
  gzwrite( backup_file, &min_w,                       sizeof(min_w)                     );
  gzwrite( backup_file, &max_w,                       sizeof(max_w)                     );

  // Statistics collection 
  gzwrite( backup_file, &backup_step,                 sizeof(backup_step)               );
  gzwrite( backup_file, &big_backup_step,             sizeof(big_backup_step)           );
  gzwrite( backup_file, &tree_step,                   sizeof(tree_step)                 );
  int8_t tmp_record_tree = record_tree? 1 : 0;
  gzwrite( backup_file, &tmp_record_tree,             sizeof(tmp_record_tree)           );
  int8_t tmp_tree_mode = tree_mode;
  gzwrite( backup_file, &tmp_tree_mode,               sizeof(tmp_tree_mode)             );
  int8_t tmp_more_stats = more_stats? 1 : 0;
  gzwrite( backup_file, &tmp_more_stats,              sizeof(tmp_more_stats)            );
  gzwrite( backup_file, &dump_period,                 sizeof(dump_period)               );

  
  // Population size, structure, and other properties 
  gzwrite( backup_file, &init_pop_size,               sizeof(init_pop_size)             );
  int8_t tmp_pop_structure = pop_structure? 1 : 0;
  gzwrite( backup_file, &tmp_pop_structure,           sizeof(tmp_pop_structure)         );
  gzwrite( backup_file, &grid_x,                      sizeof(grid_x)                    );
  gzwrite( backup_file, &grid_y,                      sizeof(grid_y)                    );
  gzwrite( backup_file, &migration_number,            sizeof(migration_number)          );
  
  // Mutation rates
  gzwrite( backup_file, &point_mutation_rate,         sizeof(point_mutation_rate)       );
  gzwrite( backup_file, &small_insertion_rate,        sizeof(small_insertion_rate)      );
  gzwrite( backup_file, &small_deletion_rate,         sizeof(small_deletion_rate)       );
  gzwrite( backup_file, &max_indel_size,              sizeof(max_indel_size)            );

  // Rearrangements and Transfer
  int8_t tmp_with_4pts_rears = with_4pts_rears? 1 : 0;
  gzwrite( backup_file, &tmp_with_4pts_rears,         sizeof(tmp_with_4pts_rears)       );
  int8_t tmp_with_alignments = with_alignments? 1 : 0;
  gzwrite( backup_file, &tmp_with_alignments,         sizeof(tmp_with_alignments)       );
  int8_t tmp_with_transfer = with_transfer? 1 : 0;
  gzwrite( backup_file, &tmp_with_transfer,           sizeof(tmp_with_transfer)         );

  // Rearrangement rates (without alignements)
  gzwrite( backup_file, &duplication_rate,            sizeof(duplication_rate)          );
  gzwrite( backup_file, &deletion_rate,               sizeof(deletion_rate)             );
  gzwrite( backup_file, &translocation_rate,          sizeof(translocation_rate)        );
  gzwrite( backup_file, &inversion_rate,              sizeof(inversion_rate)            );

  // Rearrangement rates (with alignements)
  gzwrite( backup_file, &neighbourhood_rate,          sizeof(neighbourhood_rate)        );
  gzwrite( backup_file, &duplication_proportion,      sizeof(duplication_proportion)    );
  gzwrite( backup_file, &deletion_proportion,         sizeof(deletion_proportion)       );
  gzwrite( backup_file, &translocation_proportion,    sizeof(translocation_proportion)  );
  gzwrite( backup_file, &inversion_proportion,        sizeof(inversion_proportion)      );

  // Alignements
  gzwrite( backup_file, &align_min_score,             sizeof(align_min_score)           );
  gzwrite( backup_file, &align_max_score,             sizeof(align_max_score)           );
  gzwrite( backup_file, &align_max_shift,             sizeof(align_max_shift)           );
  gzwrite( backup_file, &align_w_zone_h_len,          sizeof(align_w_zone_h_len)        );
  gzwrite( backup_file, &align_match_bonus,           sizeof(align_match_bonus)         );
  gzwrite( backup_file, &align_mismatch_cost,         sizeof(align_mismatch_cost)       );

  // Selection
  int8_t tmp_selection_scheme = selection_scheme;
  gzwrite( backup_file, &tmp_selection_scheme,        sizeof(tmp_selection_scheme)      );
  gzwrite( backup_file, &selection_pressure,          sizeof(selection_pressure)        );
  
  // Environment gaussians
  int16_t nb_gaussians = env_gaussians.get_nb_elts();
  gzwrite( backup_file, &nb_gaussians, sizeof(nb_gaussians) );
  //~ printf( "write %d bytes (nb_gaussians : %"PRId16")\n", sizeof(nb_gaussians), nb_gaussians );
  ae_list_node* gaussian_node = env_gaussians.get_first();
  ae_gaussian*  gaussian;
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    gaussian = ( ae_gaussian* ) gaussian_node->get_obj();

    gaussian->write_to_backup( backup_file );

    gaussian_node = gaussian_node->get_next();
  }
  
  // Environment custom points
  int16_t nb_points = env_custom_points.get_nb_elts();
  gzwrite( backup_file, &nb_points, sizeof(nb_points) );
  //~ printf( "write %d bytes (nb_points : %"PRId16")\n", sizeof(nb_points), nb_points );
  ae_list_node* point_node = env_custom_points.get_first();
  ae_point_2d*  point;
  for ( int16_t i = 0 ; i < nb_points ; i++ )
  {
    point = ( ae_point_2d* ) point_node->get_obj();

    point->write_to_backup( backup_file );

    point_node = point_node->get_next();
  }
  
  // Environment sampling
  gzwrite( backup_file, &env_sampling, sizeof(env_sampling) );
  //~ printf( "write %d bytes (env_sampling : %"PRId16")\n", sizeof(env_sampling), env_sampling );
  
  // Environment variation
  int8_t tmp_env_var_method = env_var_method;
  gzwrite( backup_file, &tmp_env_var_method,          sizeof(tmp_env_var_method)  );
  gzwrite( backup_file, &env_sigma,                   sizeof(env_sigma)           );
  gzwrite( backup_file, &env_tau,                     sizeof(env_tau)             );
  
  
  // Environment x-axis segmentation
  int8_t tmp_env_axis_is_segmented = env_axis_is_segmented? 1 : 0;
  gzwrite( backup_file, &tmp_env_axis_is_segmented, sizeof(tmp_env_axis_is_segmented) );
  
  if ( env_axis_is_segmented )
  {
    gzwrite( backup_file, &env_axis_nb_segments, sizeof(env_axis_nb_segments) );
    
    // Write segment boundaries
    for ( int16_t i = 0 ; i < env_axis_nb_segments + 1 ; i++ )
    {
      gzwrite( backup_file, &env_axis_segment_boundaries[i], sizeof(env_axis_segment_boundaries[i]) );
    }
    
    // Write segment features
    for ( int16_t i = 0 ; i < env_axis_nb_segments ; i++ )
    {
      int8_t tmp_env_axis_features = env_axis_features[i];
      gzwrite( backup_file, &tmp_env_axis_features, sizeof(tmp_env_axis_features) );
    }
    
    int8_t tmp_env_separate_segments = env_separate_segments? 1 : 0;
    gzwrite( backup_file, &tmp_env_separate_segments,   sizeof(tmp_env_separate_segments) );
  }
  
  // Is the fitness composed of various parts? (e.g. metabolism and secretion)
  int8_t tmp_composite_fitness = composite_fitness? 1 : 0;
  gzwrite( backup_file, &tmp_composite_fitness, sizeof(tmp_composite_fitness) );

  
  // Secretion
  int8_t tmp_use_secretion = use_secretion? 1 : 0;
  gzwrite( backup_file, &tmp_use_secretion,           sizeof(tmp_use_secretion)           );
  gzwrite( backup_file, &secretion_fitness_contrib,   sizeof(secretion_fitness_contrib)   );
  gzwrite( backup_file, &secretion_difusion_prop,     sizeof(secretion_difusion_prop)     );
  gzwrite( backup_file, &secretion_degradation_prop,  sizeof(secretion_degradation_prop)  );
  gzwrite( backup_file, &secretion_init,              sizeof(secretion_init)              );
  gzwrite( backup_file, &secretion_cost,              sizeof(secretion_cost)              );

  // Plasmids
  int8_t tmp_allow_plasmids = allow_plasmids? 1 : 0;
  gzwrite( backup_file, &tmp_allow_plasmids,              sizeof(tmp_allow_plasmids)            );
  gzwrite( backup_file, &plasmid_initial_length,          sizeof(plasmid_initial_length)        );
  gzwrite( backup_file, &plasmid_initial_gene,            sizeof(plasmid_initial_gene)          );
  gzwrite( backup_file, &plasmid_minimal_length,          sizeof(plasmid_minimal_length)        );
  gzwrite( backup_file, &prob_horiz_trans,                sizeof(prob_horiz_trans)              );
  int8_t tmp_compute_phen_contrib_by_GU = compute_phen_contrib_by_GU? 1 : 0;
  gzwrite( backup_file, &tmp_compute_phen_contrib_by_GU, sizeof(tmp_compute_phen_contrib_by_GU) );
  
  // Translation cost
  gzwrite( backup_file, &translation_cost, sizeof(translation_cost) );
  
  // Log files
  gzwrite( backup_file, &logs, sizeof(logs) );
  

  // R-AEVOL specific
  #ifdef __REGUL
    // Regulation factors
    gzwrite( backup_file, &hill_shape_n,                sizeof(hill_shape_n)      );
    gzwrite( backup_file, &hill_shape_theta,            sizeof(hill_shape_theta)  );
    gzwrite( backup_file, &hill_shape,                  sizeof(hill_shape)        );

    // Degradation rate
    gzwrite( backup_file, &degradation_rate,            sizeof(degradation_rate)  );
    gzwrite( backup_file, &degradation_step,            sizeof(degradation_step)  );

    // Individual life
    gzwrite( backup_file, &individual_evaluation_nbr,   sizeof(individual_evaluation_nbr)  );
    individual_evaluation_dates->write_to_backup( backup_file );

    // Binding matrix
    gzwrite( backup_file, &binding_zeros_percentage,    sizeof(binding_zeros_percentage)    );
    _binding_matrix->write_to_backup( backup_file );
    
    // Heredity
    int8_t tmp_with_heredity = with_heredity? 1 : 0;
    gzwrite( backup_file, &tmp_with_heredity,           sizeof(tmp_with_heredity)           );
    gzwrite( backup_file, &protein_presence_limit,      sizeof(protein_presence_limit)      );
  #endif
}

void ae_common::read_from_backup( gzFile* backup_file, bool verbose )
{
  // PseudoRandom Number Generator
  gzread( backup_file, &seed,                         sizeof(seed)                      );
  gzread( backup_file, &env_seed,                     sizeof(env_seed)                  );

  // Initial conditions
  gzread( backup_file, &initial_genome_length,        sizeof(initial_genome_length)     );
  gzread( backup_file, &min_genome_length,            sizeof(min_genome_length)         );
  gzread( backup_file, &max_genome_length,            sizeof(max_genome_length)         );
  gzread( backup_file, &init_method,                  sizeof(init_method)               );
  gzread( backup_file, &min_w,                        sizeof(min_w)                     );
  gzread( backup_file, &max_w,                        sizeof(max_w)                     );

  // Statistics collection
  gzread( backup_file, &backup_step,                  sizeof(backup_step)               );
  gzread( backup_file, &big_backup_step,              sizeof(big_backup_step)           );
  gzread( backup_file, &tree_step,                    sizeof(tree_step)                 );
  int8_t tmp_record_tree;
  gzread( backup_file, &tmp_record_tree,              sizeof(tmp_record_tree)           );
  record_tree = (tmp_record_tree != 0);
  int8_t tmp_tree_mode;
  gzread( backup_file, &tmp_tree_mode,                sizeof(tmp_tree_mode)             );
  tree_mode = (ae_tree_mode) tmp_tree_mode;
  int8_t tmp_more_stats;
  gzread( backup_file, &tmp_more_stats,               sizeof(tmp_more_stats)            );
  more_stats = (tmp_more_stats != 0);
  gzread( backup_file, &dump_period,                  sizeof(dump_period)               );
  
  // Population size, structure, and other properties 
  gzread( backup_file, &init_pop_size,                sizeof(init_pop_size)             );
  int8_t tmp_pop_structure;
  gzread( backup_file, &tmp_pop_structure,            sizeof(tmp_pop_structure)         );
  pop_structure = (tmp_pop_structure!=0);
  gzread( backup_file, &grid_x,                       sizeof(grid_x)                    );
  gzread( backup_file, &grid_y,                       sizeof(grid_y)                    );
  gzread( backup_file, &migration_number,             sizeof(migration_number)          );
  
  // Mutation rates
  gzread( backup_file, &point_mutation_rate,          sizeof(point_mutation_rate)       );
  gzread( backup_file, &small_insertion_rate,         sizeof(small_insertion_rate)      );
  gzread( backup_file, &small_deletion_rate,          sizeof(small_deletion_rate)       );
  gzread( backup_file, &max_indel_size,               sizeof(max_indel_size)            );

  // Rearrangements and Transfer
  int8_t tmp_with_4pts_rears;
  gzread( backup_file, &tmp_with_4pts_rears,          sizeof(tmp_with_4pts_rears)       );
  with_4pts_rears = (tmp_with_4pts_rears!=0);
  int8_t tmp_with_alignments;
  gzread( backup_file, &tmp_with_alignments,          sizeof(tmp_with_alignments)       );
  with_alignments = (tmp_with_alignments!=0);
  int8_t tmp_with_transfer;
  gzread( backup_file, &tmp_with_transfer,            sizeof(tmp_with_transfer)         );
  with_transfer = (tmp_with_transfer!=0);

  // Rearrangement rates (without alignements)
  gzread( backup_file, &duplication_rate,             sizeof(duplication_rate)          );
  gzread( backup_file, &deletion_rate,                sizeof(deletion_rate)             );
  gzread( backup_file, &translocation_rate,           sizeof(translocation_rate)        );
  gzread( backup_file, &inversion_rate,               sizeof(inversion_rate)            );

  // Rearrangement rates (with alignements)
  gzread( backup_file, &neighbourhood_rate,           sizeof(neighbourhood_rate)        );
  gzread( backup_file, &duplication_proportion,       sizeof(duplication_proportion)    );
  gzread( backup_file, &deletion_proportion,          sizeof(deletion_proportion)       );
  gzread( backup_file, &translocation_proportion,     sizeof(translocation_proportion)  );
  gzread( backup_file, &inversion_proportion,         sizeof(inversion_proportion)      );

  // Alignements
  gzread( backup_file, &align_min_score,             sizeof(align_min_score)            );
  gzread( backup_file, &align_max_score,             sizeof(align_max_score)            );
  gzread( backup_file, &align_max_shift,             sizeof(align_max_shift)            );
  gzread( backup_file, &align_w_zone_h_len,          sizeof(align_w_zone_h_len)         );
  gzread( backup_file, &align_match_bonus,           sizeof(align_match_bonus)          );
  gzread( backup_file, &align_mismatch_cost,         sizeof(align_mismatch_cost)        );

  // Selection
  if ( verbose )
  {
    printf( "    Loading selection data\n" );
  }
  int8_t tmp_selection_scheme;
  gzread( backup_file, &tmp_selection_scheme,         sizeof(tmp_selection_scheme)      );
  selection_scheme = (ae_selection_scheme) tmp_selection_scheme;
  gzread( backup_file, &selection_pressure,           sizeof(selection_pressure)        );
  
  // Environment gaussians
  if ( verbose )
  {
    printf( "    Loading environmental data\n" );
  }
  int16_t nb_gaussians;
  gzread( backup_file, &nb_gaussians, sizeof(nb_gaussians) );
  //~ printf( "read %d bytes (nb_gaussians : %"PRId16")\n", sizeof(nb_gaussians), nb_gaussians );
  
  if(env_gaussians.is_empty() == false)  
  {
    env_gaussians.erase( DELETE_OBJ );
  }
  
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    env_gaussians.add( new ae_gaussian( backup_file ) );
  }
    
  
  // Environment custom points
  int16_t nb_points;
  gzread( backup_file, &nb_points, sizeof(nb_points) );
  //~ printf( "read %d bytes (nb_points : %"PRId16")\n", sizeof(nb_points), nb_points );

  if(env_custom_points.is_empty() == false)  
  {
    env_custom_points.erase( DELETE_OBJ );
  }

  for ( int16_t i = 0 ; i < nb_points ; i++ )
  {
    env_custom_points.add( new ae_point_2d( backup_file ) );
  }
  
  
  // Environment sampling
  gzread( backup_file, &env_sampling,                sizeof(env_sampling)      );
  //~ printf( "read %d bytes (env_sampling : %"PRId16")\n", sizeof(env_sampling), env_sampling );
  
  
  // Environment variation
  int8_t tmp_env_var_method;
  gzread( backup_file, &tmp_env_var_method,          sizeof(tmp_env_var_method) );
  env_var_method = (ae_env_var) tmp_env_var_method;
  gzread( backup_file, &env_sigma,                   sizeof(env_sigma)          );
  gzread( backup_file, &env_tau,                     sizeof(env_tau)            );
  
  
  // Environment x-axis segmentation
  if ( verbose )
  {
    printf( "    Loading x-axis segmentation data\n" );
  }
  int8_t tmp_env_axis_is_segmented;
  gzread( backup_file, &tmp_env_axis_is_segmented, sizeof(tmp_env_axis_is_segmented) );
  env_axis_is_segmented = (tmp_env_axis_is_segmented != 0);
  
  if ( env_axis_is_segmented )
  {
    gzread( backup_file, &env_axis_nb_segments, sizeof(env_axis_nb_segments) );
    env_axis_segment_boundaries = new double [env_axis_nb_segments + 1];
    for ( int16_t i = 0 ; i < env_axis_nb_segments + 1 ; i++ )
    {
      gzread( backup_file, &env_axis_segment_boundaries[i], sizeof(env_axis_segment_boundaries[i]) );
    }
    env_axis_features = new ae_env_axis_feature [env_axis_nb_segments];
    for ( int16_t i = 0 ; i < env_axis_nb_segments ; i++ )
    {
      int8_t tmp_env_axis_features;
      gzread( backup_file, &tmp_env_axis_features, sizeof(tmp_env_axis_features) );
      env_axis_features[i] = (ae_env_axis_feature) tmp_env_axis_features;
    }
    int8_t tmp_env_separate_segments;
    gzread( backup_file, &tmp_env_separate_segments,   sizeof(tmp_env_separate_segments) );
    env_separate_segments = (tmp_env_separate_segments!=0);
  }
  
  int8_t tmp_composite_fitness;
  gzread( backup_file, &tmp_composite_fitness, sizeof(tmp_composite_fitness) );
  composite_fitness = (tmp_composite_fitness != 0);

  // Secretion
  if ( verbose )
  {
    printf( "    Loading secretion data\n" );
  }
  int8_t tmp_use_secretion;
  gzread( backup_file, &tmp_use_secretion,            sizeof(tmp_use_secretion)           );
  use_secretion = (tmp_use_secretion!=0);
  gzread( backup_file, &secretion_fitness_contrib,    sizeof(secretion_fitness_contrib)   );
  gzread( backup_file, &secretion_difusion_prop,      sizeof(secretion_difusion_prop)     );
  gzread( backup_file, &secretion_degradation_prop,   sizeof(secretion_degradation_prop)  );
  gzread( backup_file, &secretion_init,               sizeof(secretion_init)              );
  gzread( backup_file, &secretion_cost,               sizeof(secretion_cost)              );

  // Plasmids
  if ( verbose )
  {
    printf( "    Loading plasmid management data\n" );
  }
  int8_t tmp_allow_plasmids;
  gzread( backup_file, &tmp_allow_plasmids,              sizeof(tmp_allow_plasmids)              );
  allow_plasmids = (tmp_allow_plasmids!=0);
  gzread( backup_file, &plasmid_initial_length,          sizeof(plasmid_initial_length)         );
  gzread( backup_file, &plasmid_initial_gene,            sizeof(plasmid_initial_gene)           );
  gzread( backup_file, &plasmid_minimal_length,          sizeof(plasmid_minimal_length)         );
  gzread( backup_file, &prob_horiz_trans,                sizeof(prob_horiz_trans)               );
  int8_t tmp_compute_phen_contrib_by_GU;
  gzread( backup_file, &tmp_compute_phen_contrib_by_GU,  sizeof(tmp_compute_phen_contrib_by_GU) );
  compute_phen_contrib_by_GU = (tmp_compute_phen_contrib_by_GU!=0);
  
  
  // Translation cost
  gzread( backup_file, &translation_cost, sizeof(translation_cost) );
  
  
  // Log files
  gzread( backup_file, &logs, sizeof(logs) );

#ifdef __REGUL
  // Regulation factors
  gzread( backup_file, &hill_shape_n,                 sizeof(hill_shape_n)      );
  gzread( backup_file, &hill_shape_theta,             sizeof(hill_shape_theta)  );
  gzread( backup_file, &hill_shape,                   sizeof(hill_shape)        );
  
  // Degradation rate
  gzread( backup_file, &degradation_rate,             sizeof(degradation_rate)  );
  gzread( backup_file, &degradation_step,             sizeof(degradation_step)  );
  
  // Individual life
  gzread( backup_file, &individual_evaluation_nbr,    sizeof(individual_evaluation_nbr)    );
  individual_evaluation_dates  = new ae_array_short( individual_evaluation_nbr );
  individual_evaluation_dates->read_from_backup( backup_file );


  // Binding matrix
  gzread( backup_file, &binding_zeros_percentage,    sizeof(binding_zeros_percentage)    );
  _binding_matrix = new ae_matrix_double( MAX_QUADON, MAX_CODON );
  _binding_matrix->read_from_backup( backup_file );

      
  // Heredity
  int8_t tmp_with_heredity;
  gzread( backup_file, &tmp_with_heredity,           sizeof(tmp_with_heredity)           );
  with_heredity = (tmp_with_heredity!=0);
  gzread( backup_file, &protein_presence_limit,      sizeof(protein_presence_limit)      );
#endif
}

void ae_common::print_to_file( void )
{
  FILE* param_out = fopen( "param.out", "w" );
  assert( param_out );
  
  // PseudoRandom Number Generator
  fprintf( param_out, "seed :                      %"PRId32"\n", seed                    );
  fprintf( param_out, "env_seed :                  %"PRId32"\n", env_seed                );

  // Initial conditions
  fprintf( param_out, "initial_genome_length :     %"PRId32"\n", initial_genome_length   );
  fprintf( param_out, "min_genome_length :         %"PRId32"\n", min_genome_length       );
  fprintf( param_out, "max_genome_length :         %"PRId32"\n", max_genome_length       );
  
  fprintf( param_out, "init_method :              " );
  if ( init_method & ONE_GOOD_GENE )
  {
    fprintf( param_out, " ONE_GOOD_GENE" );
  }
  if ( init_method & CLONE )
  {
    fprintf( param_out, " CLONE" );
  }
  if ( init_method & WITH_INS_SEQ )
  {
    fprintf( param_out, " WITH_INS_SEQ" );
  }
  fprintf( param_out, "\n" );
  
  fprintf( param_out, "MIN_W :               %lf\n", min_w );
  fprintf( param_out, "MAX_W :               %lf\n", max_w );

  // Statistics collection
  fprintf( param_out, "backup_step :               %"PRId32"\n",  backup_step );
  fprintf( param_out, "tree_step :                 %"PRId32"\n",  tree_step );
  fprintf( param_out, "record_tree :               %s\n",  record_tree? "true" : "false" );
  switch ( tree_mode )
  {
    case LIGHT :
    {
      fprintf( param_out, "tree_mode   :               LIGHT\n" );
      break;
    }
    case NORMAL :
    {
      fprintf( param_out, "tree_mode   :               NORMAL\n" );
      break;
    }
    default :
    {
      fprintf( param_out, "tree_mode   :               UNKNOWN\n" );
      break;
    }
  }
  fprintf( param_out, "more_stats :                %s\n",  more_stats? "true" : "false"  );
  fprintf( param_out, "dump_period :               %"PRId32"\n",  dump_period            );

  // Population size, structure, and other properties 
  fprintf( param_out, "init_pop_size :             %"PRId32"\n", init_pop_size             );
  fprintf( param_out, "pop_structure :             %s\n",  pop_structure? "true" : "false" );
  fprintf( param_out, "grid_x :                    %"PRId16"\n", grid_x                    ); 
  fprintf( param_out, "grid_y :                    %"PRId16"\n", grid_y                    ); 
  fprintf( param_out, "migration_number :          %"PRId32"\n", migration_number          ); 
  
  // Mutation rates
  fprintf( param_out, "point_mutation_rate :       %e\n",  point_mutation_rate        );
  fprintf( param_out, "small_insertion_rate :      %e\n",  small_insertion_rate       );
  fprintf( param_out, "small_deletion_rate :       %e\n",  small_deletion_rate        );
  fprintf( param_out, "max_indel_size :            %"PRId16"\n", max_indel_size       );

  // Rearrangements and Transfer
  fprintf( param_out, "with_4pts_rears :           %s\n",  with_4pts_rears? "true" : "false" );
  fprintf( param_out, "with_alignments :           %s\n",  with_alignments? "true" : "false" );
  fprintf( param_out, "with_transfer :             %s\n",  with_transfer? "true" : "false"   );

  // Rearrangement rates (without alignements)
  fprintf( param_out, "duplication_rate :          %e\n",  duplication_rate           );
  fprintf( param_out, "deletion_rate :             %e\n",  deletion_rate              );
  fprintf( param_out, "translocation_rate :        %e\n",  translocation_rate         );
  fprintf( param_out, "inversion_rate :            %e\n",  inversion_rate             );

  // Rearrangement rates (with alignements)
  fprintf( param_out, "neighbourhood_rate :        %e\n",  neighbourhood_rate         );
  fprintf( param_out, "duplication_proportion :    %e\n",  duplication_proportion     );
  fprintf( param_out, "deletion_proportion :       %e\n",  deletion_proportion        );
  fprintf( param_out, "translocation_proportion :  %e\n",  translocation_proportion   );
  fprintf( param_out, "inversion_proportion :      %e\n",  inversion_proportion       );

  // Alignements
  fprintf( param_out, "align_min_score :           %"PRId16"\n", align_min_score            );
  fprintf( param_out, "align_max_score :           %"PRId16"\n", align_max_score            );
  fprintf( param_out, "align_max_shift :           %"PRId16"\n", align_max_shift            );
  fprintf( param_out, "align_w_zone_h_len :        %"PRId16"\n", align_w_zone_h_len         );
  fprintf( param_out, "align_match_bonus :         %"PRId16"\n", align_match_bonus          );
  fprintf( param_out, "align_mismatch_cost :       %"PRId16"\n", align_mismatch_cost        );

  // Selection
  switch ( selection_scheme )
  {
    case RANK_LINEAR :
    {
      fprintf( param_out, "selection_scheme   :        RANK_LINEAR\n" );
      break;
    }
    case RANK_EXPONENTIAL :
    {
      fprintf( param_out, "selection_scheme   :        RANK_EXPONENTIAL\n" );
      break;
    }
    case FITNESS_PROPORTIONATE :
    {
      fprintf( param_out, "selection_scheme   :        FITNESS_PROPORTIONATE\n" );
      break;
    }
    case FITTEST :
    {
      fprintf( param_out, "selection_scheme   :        FITTEST\n" );
      break;
    }
    default :
    {
      fprintf( param_out, "selection_scheme   :        UNKNOWN\n" );
      break;
    }
  }
  fprintf( param_out, "selection_pressure :        %e\n",  selection_pressure         );
  
  // Environment
  // TODO
  //~ static ae_list    env_gaussians;
  //~ static ae_list    env_custom_points;
  //~ static int16_t    env_sampling;
  
  //~ static ae_env_var env_var_method;
  //~ static double     env_sigma;
  //~ static int32_t    env_tau;

  //~ static bool                 env_axis_is_segmented;
  //~ static int16_t              env_axis_nb_segments;
  //~ static double*              env_axis_segment_boundaries;
  //~ static ae_env_axis_feature* env_axis_features;
  //~ static bool                 env_separate_segments;
  
  //~ static bool                 composite_fitness;
  
  
  // Secretion
  fprintf( param_out, "use_secretion :              %s\n", use_secretion? "true" : "false" );
  fprintf( param_out, "secretion_fitness_contrib :  %e\n", secretion_fitness_contrib       );
  fprintf( param_out, "secretion_difusion_prop :    %e\n", secretion_difusion_prop         );
  fprintf( param_out, "secretion_degradation_prop : %e\n", secretion_degradation_prop      );
  fprintf( param_out, "secretion_init :             %e\n", secretion_init                  );
  fprintf( param_out, "secretion_cost :             %e\n", secretion_cost                  );

  // Plasmids
  fprintf( param_out, "allow_plasmids :             %s\n", allow_plasmids? "true" : "false"              );
  fprintf( param_out, "plasmid_initial_length :     %"PRId32"\n", plasmid_initial_length                 );
  fprintf( param_out, "plasmid_initial_gene :       %"PRId32"\n", plasmid_initial_gene                   );
  fprintf( param_out, "plasmid_minimal_length :     %"PRId32"\n", plasmid_minimal_length                 );
  fprintf( param_out, "prob_horiz_trans :           %e\n", prob_horiz_trans                              );
  fprintf( param_out, "compute_phen_contrib_by_GU : %s\n", compute_phen_contrib_by_GU? "true" : "false"  );

  // Translation cost
  fprintf( param_out, "translation_cost :          %e\n",  translation_cost           );

  #ifdef __REGUL
    fprintf( param_out, "\n********** RAEVOL SPECIFIC PARAMETERS **********\n" );
  
    // Regulation factors
    fprintf( param_out, "hill_shape_n :                %e\n",                 hill_shape_n  );
    fprintf( param_out, "hill_shape_theta :            %e\n",             hill_shape_theta  );
    fprintf( param_out, "hill_shape :                  %e\n",                   hill_shape  );

    // Degradation rate
    fprintf( param_out, "degradation_rate :            %e\n",             degradation_rate  );
    fprintf( param_out, "degradation_step :            %e\n",             degradation_step  );

    // Individual life
    fprintf( param_out, "individual_evaluation_nbr :   %"PRId16"\n",    individual_evaluation_nbr  );
    fprintf( param_out, "individual_evaluation_dates : ");
    individual_evaluation_dates->print_to_file( param_out );

    // Binding matrix
    fprintf( param_out, "binding_zeros_percentage :    %e\n",     binding_zeros_percentage    );
    fprintf( param_out, "binding_matrix :\n");
    _binding_matrix->print_to_file( param_out );

    // Heredity
    fprintf( param_out, "with_heredity :               %s\n", with_heredity? "true" : "false"  );
    fprintf( param_out, "protein_presence_limit :      %e\n", protein_presence_limit           );
  #endif
  
  fclose( param_out );
}

void ae_common::clean( void )
{
  if (sim != NULL)
  {
    delete sim;
    sim = NULL;
  }
  
  env_gaussians.erase( DELETE_OBJ );
  env_custom_points.erase( DELETE_OBJ );
  
  if ( env_axis_features != NULL )
  {
    delete [] env_axis_features;
    env_axis_features = NULL;
  }

  if ( env_axis_segment_boundaries != NULL )
  {
    delete [] env_axis_segment_boundaries;
    env_axis_segment_boundaries = NULL;
  }
  
  #ifdef __REGUL
    if ( individual_evaluation_dates != NULL )
    {
      delete individual_evaluation_dates;
    }
    if ( _binding_matrix != NULL )
    {
      delete _binding_matrix;
      _binding_matrix = NULL;
    }
  #endif
}

#ifdef __REGUL
void ae_common::init_binding_matrix( void )
{
  _binding_matrix = new ae_matrix_double( MAX_QUADON, MAX_CODON );

  for( int8_t i = 0; i < MAX_QUADON; i++ )  // i for the quadons
  {
    for( int8_t j = 0; j < MAX_CODON; j++ )  // j for the codons
    {
      if( sim->alea->random() > binding_zeros_percentage)
      {
        _binding_matrix->set_value( i, j, sim->alea->random() );
      }
      else
      {
        _binding_matrix->set_value( i, j, 0 );
      }
    }
  }
}
#endif

// =================================================================
//                           Protected Methods
// =================================================================
