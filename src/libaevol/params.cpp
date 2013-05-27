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


/*! \class
    \brief
*/
 
 
// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include <params.h>




//##############################################################################
//                                                                             #
//                                Class params                                 #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
params::params( void )
{
  // ---------------------------------- Number of generations to be simulated
  _nb_gener = 0;
  
  // ----------------------------------------- PseudoRandom Number Generators
  _seed           = 0;
  _mut_seed       = 0;
  _stoch_seed     = 0;
  _env_var_seed   = 0;
  _env_noise_seed = 0;
  
  // ------------------------------------------------------------ Constraints
  _min_genome_length  = 10;
  _max_genome_length  = 10000000;
  _w_max              = 0.033333333;
    
  // ----------------------------------------------------- Initial conditions
  _initial_genome_length  = 5000;
  _init_method            = ONE_GOOD_GENE & CLONE;
  _init_pop_size          = 1000;
    
  // ------------------------------------------------------------ Environment
  _env_gaussians      = NULL;
  _env_custom_points  = NULL;
  _env_sampling       = 300;
    
  // ---------------------------------------- Environment x-axis segmentation
  _env_axis_is_segmented        = false;
  _env_axis_nb_segments         = 1;
  _env_axis_segment_boundaries  = NULL;
  _env_axis_features            = NULL;
  _env_axis_separate_segments   = false;
  
  // -------------------------------------------------- Environment variation
  _env_var_method = NO_VAR;
  _env_var_sigma  = 0;
  _env_var_tau    = 0;
  
  // ------------------------------------------------------ Environment noise
  _env_noise_method       = NO_NOISE;
  _env_noise_alpha        = 0;
  _env_noise_sigma        = 0;
  _env_noise_prob         = 0;
  _env_noise_sampling_log = 0;

  // --------------------------------------------------------- Mutation rates
  _point_mutation_rate  = 1e-5;
  _small_insertion_rate = 1e-5;
  _small_deletion_rate  = 1e-5;
  _max_indel_size       = 6;

  // -------------------------------------------- Rearrangements and Transfer
  _with_4pts_trans  = true;
  _with_alignments  = false;
  _with_HT          = false;
  _HT_ins_rate      = 0.0;
  _HT_repl_rate     = 0.0;

  // ------------------------------ Rearrangement rates (without alignements)
  _duplication_rate   = 5e-5;
  _deletion_rate      = 5e-5;
  _translocation_rate = 5e-5;
  _inversion_rate     = 5e-5;

  // --------------------------------- Rearrangement rates (with alignements)
  _neighbourhood_rate       = 5e-5;
  _duplication_proportion   = 0.3;
  _deletion_proportion      = 0.3;
  _translocation_proportion = 0.3;
  _inversion_proportion     = 0.3;

  // ------------------------------------------------------------ Alignements
  _align_fun_shape    = SIGMOID;
  _align_sigm_lambda  = 4;
  _align_sigm_mean    = 50;
  _align_lin_min      = 0;
  _align_lin_max      = 100;
  
  _align_max_shift      = 20;
  _align_w_zone_h_len   = 50;
  _align_match_bonus    = 1;
  _align_mismatch_cost  = 2;

  // ----------------------------------------------- Phenotypic Stochasticity
  _with_stochasticity = false;

  // -------------------------------------------------------------- Selection
  _selection_scheme   = RANK_EXPONENTIAL;
  _selection_pressure = 0.998;
    
  // ------------------------------------------------------ Spatial structure
  _spatially_structured       = false;
  int16_t  _grid_width        = 0;
  int16_t  _grid_height       = 0;
  int32_t  _migration_number  = 0;
  
  // -------------------------------------------------------------- Secretion
  _with_secretion               = false;
  _secretion_contrib_to_fitness = 0;
  _secretion_diffusion_prop     = 0;
  _secretion_degradation_prop   = 0;
  _secretion_cost               = 0;
  _secretion_init               = 0;
  
  // --------------------------------------------------------------- Plasmids
  _allow_plasmids             = false;
  _plasmid_initial_length     = 1000;
  _plasmid_initial_gene       = 0;
  _plasmid_minimal_length     = 40;
  _prob_plasmid_HT            = 0;
  _tune_donor_ability         = 0;
  _tune_recipient_ability     = 0;
  _donor_cost                 = 0;
  _recipient_cost             = 0;
  _compute_phen_contrib_by_GU = false;
  _swap_GUs         = false;
  
  // ------------------------------------------------------- Translation cost
  _translation_cost = 0;
    
  // ---------------------------------------------------------------- Outputs
  _stats            = 0;
  _delete_old_stats = false;
  
  // Backups
  _backup_step      = 500;
  _big_backup_step  = 10000;
  
  // Tree
  _record_tree  = false;
  _tree_step    = 100;
  _tree_mode    = NORMAL;
  
  // Dumps
  _make_dumps = false;
  _dump_step  = 1000;
  
  // Logs
  _logs = 0;
  
  // Other
  _more_stats = false;

  #ifdef __REGUL
    // ------------------------------------------------------- Binding matrix
    _binding_zeros_percentage = 75;
  #endif
}

// =================================================================
//                             Destructors
// =================================================================
params::~params( void )
{
  // these are allocated during reading of the parameter file
  if ( _env_axis_segment_boundaries != NULL ) delete [] _env_axis_segment_boundaries;
  if ( _env_axis_features != NULL ) delete [] _env_axis_features;
}

// =================================================================
//                            Public Methods
// =================================================================
void params::print_to_file( FILE* file )
{
  // ------------------------------------------------------------ Constraints
  fprintf( file, "\nConstraints ---------------------------------------------\n" );
  fprintf( file, "min_genome_length :          %"PRId32"\n", _min_genome_length       );
  fprintf( file, "max_genome_length :          %"PRId32"\n", _max_genome_length       );
  fprintf( file, "W_MAX :                      %f\n",        _w_max                   );
  
  // --------------------------------------------------------- Mutation rates
  fprintf( file, "\nMutation rates ------------------------------------------\n" );
  fprintf( file, "point_mutation_rate :        %e\n",  _point_mutation_rate        );
  fprintf( file, "small_insertion_rate :       %e\n",  _small_insertion_rate       );
  fprintf( file, "small_deletion_rate :        %e\n",  _small_deletion_rate        );
  fprintf( file, "max_indel_size :             %"PRId16"\n", _max_indel_size       );

  // -------------------------------------------- Rearrangements and Transfer
  fprintf( file, "\nRearrangements and Transfer -----------------------------\n" );
  fprintf( file, "with_4pts_trans :            %s\n",  _with_4pts_trans? "true" : "false" );
  fprintf( file, "with_alignments :            %s\n",  _with_alignments? "true" : "false" );
  fprintf( file, "with_HT :                    %s\n",  _with_HT? "true" : "false"   );
  fprintf( file, "HT_ins_rate :                %e\n",  _HT_ins_rate );
  fprintf( file, "HT_repl_rate :               %e\n",  _HT_repl_rate );

  // ---------------------------------------------------- Rearrangement rates
  if ( _with_alignments )
  {
    fprintf( file, "\nRearrangement rates (with alignements) ------------------\n" );
    fprintf( file, "neighbourhood_rate :         %e\n",  _neighbourhood_rate         );
    fprintf( file, "duplication_proportion :     %e\n",  _duplication_proportion     );
    fprintf( file, "deletion_proportion :        %e\n",  _deletion_proportion        );
    fprintf( file, "translocation_proportion :   %e\n",  _translocation_proportion   );
    fprintf( file, "inversion_proportion :       %e\n",  _inversion_proportion       );
  }
  else
  {
    fprintf( file, "\nRearrangement rates (without alignements) ----------------\n" );
    fprintf( file, "duplication_rate :           %e\n",  _duplication_rate           );
    fprintf( file, "deletion_rate :              %e\n",  _deletion_rate              );
    fprintf( file, "translocation_rate :         %e\n",  _translocation_rate         );
    fprintf( file, "inversion_rate :             %e\n",  _inversion_rate             );
  }

  // ------------------------------------------------------------ Alignements
  fprintf( file, "\nAlignements ---------------------------------------------\n" );
  fprintf( file, "align_fun_shape :            %"PRId16"\n", _align_fun_shape       );
  fprintf( file, "align_sigm_lambda :          %f\n",        _align_sigm_lambda     );
  fprintf( file, "align_sigm_mean :            %"PRId16"\n", _align_sigm_mean       );
  fprintf( file, "align_lin_min :              %"PRId16"\n", _align_lin_min         );
  fprintf( file, "align_lin_max :              %"PRId16"\n", _align_lin_max         );
  fprintf( file, "align_max_shift :            %"PRId16"\n", _align_max_shift       );
  fprintf( file, "align_w_zone_h_len :         %"PRId16"\n", _align_w_zone_h_len    );
  fprintf( file, "align_match_bonus :          %"PRId16"\n", _align_match_bonus     );
  fprintf( file, "align_mismatch_cost :        %"PRId16"\n", _align_mismatch_cost   );

  // -------------------------------------------------------------- Selection
  fprintf( file, "\nSelection -----------------------------------------------\n" );
  switch ( _selection_scheme )
  {
    case RANK_LINEAR :
    {
      fprintf( file, "selection_scheme :           RANK_LINEAR\n" );
      break;
    }
    case RANK_EXPONENTIAL :
    {
      fprintf( file, "selection_scheme :           RANK_EXPONENTIAL\n" );
      break;
    }
    case FITNESS_PROPORTIONATE :
    {
      fprintf( file, "selection_scheme :           FITNESS_PROPORTIONATE\n" );
      break;
    }
    case FITTEST :
    {
      fprintf( file, "selection_scheme :           FITTEST\n" );
      break;
    }
    default :
    {
      fprintf( file, "selection_scheme :           UNKNOWN\n" );
      break;
    }
  }
  fprintf( file, "selection_pressure :         %e\n",  _selection_pressure );
  
  
  // -------------------------------------------------------------- Secretion
  fprintf( file, "\nSecretion -----------------------------------------------\n" );
  fprintf( file, "with_secretion :                %s\n", _with_secretion? "true" : "false" );
  fprintf( file, "secretion_contrib_to_fitness :  %e\n", _secretion_contrib_to_fitness    );
  fprintf( file, "secretion_diffusion_prop :      %e\n", _secretion_diffusion_prop        );
  fprintf( file, "secretion_degradation_prop :    %e\n", _secretion_degradation_prop      );
  fprintf( file, "secretion_cost :                %e\n", _secretion_cost                  );

  // --------------------------------------------------------------- Plasmids
  fprintf( file, "\nPlasmids ------------------------------------------------\n" );
  fprintf( file, "allow_plasmids :             %s\n", _allow_plasmids? "true" : "false"              );
  fprintf( file, "plasmid_minimal_length :     %"PRId32"\n", _plasmid_minimal_length                 );
  fprintf( file, "prob_plasmid_HT :            %e\n", _prob_plasmid_HT                               );
  fprintf( file, "tune_donor_ability :         %e\n", _tune_donor_ability                            );
  fprintf( file, "tune_recipient_ability :     %e\n", _tune_recipient_ability                        );
  fprintf( file, "donor_cost :                 %e\n", _donor_cost                                    );
  fprintf( file, "recipient_cost :             %e\n", _recipient_cost                                );
  fprintf( file, "compute_phen_contrib_by_GU : %s\n", _compute_phen_contrib_by_GU? "true" : "false"  );
  fprintf( file, "swap_GUs :                   %s\n",  _swap_GUs? "true" : "false"   );

  // ------------------------------------------------------- Translation cost
  fprintf( file, "\nTranslation cost ----------------------------------------\n" );
  fprintf( file, "translation_cost :           %e\n",  _translation_cost           );
}

// =================================================================
//                           Protected Methods
// =================================================================

// =================================================================
//                          Non inline accessors
// =================================================================
