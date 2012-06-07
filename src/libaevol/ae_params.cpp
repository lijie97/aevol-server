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
#include <ae_params.h>




//##############################################################################
//                                                                             #
//                               Class ae_params                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_params::ae_params( void )
{
  // ------------------------------------------------------------ Constraints
  _min_genome_length = 10;
  _max_genome_length = 10000000;
  _min_w = 0.0;
  _max_w = 0.033333333;

  // --------------------------------------------------------- Mutation rates
  _point_mutation_rate  = 1e-5;
  _small_insertion_rate = 1e-5;
  _small_deletion_rate  = 1e-5;
  _max_indel_size = 6;

  // -------------------------------------------- Rearrangements and Transfer
  _with_4pts_trans = true;
  _with_alignments = false;
  _with_transfer   = false;
  _swap_GUs        = false;
  _transfer_ins_rate  = 0.0;
  _transfer_repl_rate = 0.0;

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

  // -------------------------------------------------------------- Selection
  _selection_scheme   = RANK_EXPONENTIAL;
  _selection_pressure = 0.998;
  
  // -------------------------------------------------------------- Secretion
  _use_secretion              = false;
  _secretion_fitness_contrib  = 0;
  _secretion_diffusion_prop   = 0;
  _secretion_degradation_prop = 0;
  _secretion_cost             = 0;
  
  // --------------------------------------------------------------- Plasmids
  _allow_plasmids             = false;
  _plasmid_minimal_length     = 40;
  _prob_horiz_trans           = 0;
  _nb_horiz_trans             = 1;
  _compute_phen_contrib_by_GU = false;
  
  // ------------------------------------------------------- Translation cost
  _translation_cost = 0;
}

// =================================================================
//                             Destructors
// =================================================================
ae_params::~ae_params( void )
{
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_params::write_to_backup( gzFile* backup_file )
{
  // ------------------------------------------------------------ Constraints
  gzwrite( backup_file, &_min_genome_length, sizeof(_min_genome_length) );
  gzwrite( backup_file, &_max_genome_length, sizeof(_max_genome_length) );
  gzwrite( backup_file, &_min_w,             sizeof(_min_w)             );
  gzwrite( backup_file, &_max_w,             sizeof(_max_w)             );
  
  // --------------------------------------------------------- Mutation rates
  gzwrite( backup_file, &_point_mutation_rate,  sizeof(_point_mutation_rate)  );
  gzwrite( backup_file, &_small_insertion_rate, sizeof(_small_insertion_rate) );
  gzwrite( backup_file, &_small_deletion_rate,  sizeof(_small_deletion_rate)  );
  gzwrite( backup_file, &_max_indel_size,       sizeof(_max_indel_size)       );

  // -------------------------------------------- Rearrangements and Transfer
  int8_t tmp_with_4pts_trans = _with_4pts_trans? 1 : 0;
  gzwrite( backup_file, &tmp_with_4pts_trans, sizeof(tmp_with_4pts_trans) );
  int8_t tmp_with_alignments = _with_alignments? 1 : 0;
  gzwrite( backup_file, &tmp_with_alignments, sizeof(tmp_with_alignments) );
  int8_t tmp_with_transfer = _with_transfer? 1 : 0;
  gzwrite( backup_file, &tmp_with_transfer, sizeof(tmp_with_transfer) );
  int8_t tmp_swap_GUs = _swap_GUs? 1 : 0;
  gzwrite( backup_file, &tmp_swap_GUs, sizeof(tmp_swap_GUs) );
  gzwrite( backup_file, &_transfer_ins_rate,          sizeof(_transfer_ins_rate)  );
  gzwrite( backup_file, &_transfer_repl_rate,         sizeof(_transfer_repl_rate) );

  // ------------------------------ Rearrangement rates (without alignements)
  gzwrite( backup_file, &_duplication_rate,   sizeof(_duplication_rate)          );
  gzwrite( backup_file, &_deletion_rate,      sizeof(_deletion_rate)             );
  gzwrite( backup_file, &_translocation_rate, sizeof(_translocation_rate)        );
  gzwrite( backup_file, &_inversion_rate,     sizeof(_inversion_rate)            );

  // --------------------------------- Rearrangement rates (with alignements)
  gzwrite( backup_file, &_neighbourhood_rate,       sizeof(_neighbourhood_rate)        );
  gzwrite( backup_file, &_duplication_proportion,   sizeof(_duplication_proportion)    );
  gzwrite( backup_file, &_deletion_proportion,      sizeof(_deletion_proportion)       );
  gzwrite( backup_file, &_translocation_proportion, sizeof(_translocation_proportion)  );
  gzwrite( backup_file, &_inversion_proportion,     sizeof(_inversion_proportion)      );

  // ------------------------------------------------------------ Alignements
  int8_t tmp_align_fun_shape = _align_fun_shape;
  gzwrite( backup_file, &tmp_align_fun_shape,  sizeof(tmp_align_fun_shape)           );
  gzwrite( backup_file, &_align_sigm_lambda,   sizeof(_align_sigm_lambda)         );
  gzwrite( backup_file, &_align_sigm_mean,     sizeof(_align_sigm_mean)           );
  gzwrite( backup_file, &_align_lin_min,       sizeof(_align_lin_min)             );
  gzwrite( backup_file, &_align_lin_max,       sizeof(_align_lin_max)             );
  gzwrite( backup_file, &_align_max_shift,     sizeof(_align_max_shift)           );
  gzwrite( backup_file, &_align_w_zone_h_len,  sizeof(_align_w_zone_h_len)        );
  gzwrite( backup_file, &_align_match_bonus,   sizeof(_align_match_bonus)         );
  gzwrite( backup_file, &_align_mismatch_cost, sizeof(_align_mismatch_cost)       );

  // -------------------------------------------------------------- Selection
  int8_t tmp_selection_scheme = _selection_scheme;
  gzwrite( backup_file, &tmp_selection_scheme, sizeof(tmp_selection_scheme) );
  gzwrite( backup_file, &_selection_pressure,  sizeof(_selection_pressure)  );
  
  // -------------------------------------------------------------- Secretion
  int8_t tmp_use_secretion = _use_secretion? 1 : 0;
  gzwrite( backup_file, &tmp_use_secretion,           sizeof(tmp_use_secretion)           );
  gzwrite( backup_file, &_secretion_fitness_contrib,  sizeof(_secretion_fitness_contrib)  );
  gzwrite( backup_file, &_secretion_diffusion_prop,   sizeof(_secretion_diffusion_prop)   );
  gzwrite( backup_file, &_secretion_degradation_prop, sizeof(_secretion_degradation_prop) );
  gzwrite( backup_file, &_secretion_cost,             sizeof(_secretion_cost)             );
  
  // --------------------------------------------------------------- Plasmids
  int8_t tmp_allow_plasmids = _allow_plasmids? 1 : 0;
  gzwrite( backup_file, &tmp_allow_plasmids,      sizeof(tmp_allow_plasmids)      );
  gzwrite( backup_file, &_plasmid_minimal_length, sizeof(_plasmid_minimal_length) );
  gzwrite( backup_file, &_prob_horiz_trans,       sizeof(_prob_horiz_trans)       );
  gzwrite( backup_file, &_nb_horiz_trans,         sizeof(_nb_horiz_trans)         );
  int8_t tmp_compute_phen_contrib_by_GU = _compute_phen_contrib_by_GU? 1 : 0;
  gzwrite( backup_file, &tmp_compute_phen_contrib_by_GU, sizeof(tmp_compute_phen_contrib_by_GU) );
  
  // ------------------------------------------------------- Translation cost
  gzwrite( backup_file, &_translation_cost, sizeof(_translation_cost) );
}

void ae_params::read_from_backup( gzFile* backup_file, bool verbose )
{
  // ------------------------------------------------------------ Constraints
  gzread( backup_file, &_min_genome_length, sizeof(_min_genome_length) );
  gzread( backup_file, &_max_genome_length, sizeof(_max_genome_length) );
  gzread( backup_file, &_min_w,             sizeof(_min_w)             );
  gzread( backup_file, &_max_w,             sizeof(_max_w)             );
  
  // --------------------------------------------------------- Mutation rates
  gzread( backup_file, &_point_mutation_rate,  sizeof(_point_mutation_rate)  );
  gzread( backup_file, &_small_insertion_rate, sizeof(_small_insertion_rate) );
  gzread( backup_file, &_small_deletion_rate,  sizeof(_small_deletion_rate)  );
  gzread( backup_file, &_max_indel_size,       sizeof(_max_indel_size)       );

  // -------------------------------------------- Rearrangements and Transfer
  int8_t tmp_with_4pts_trans;
  gzread( backup_file, &tmp_with_4pts_trans, sizeof(tmp_with_4pts_trans) );
  _with_4pts_trans = (tmp_with_4pts_trans!=0);
  int8_t tmp_with_alignments;
  gzread( backup_file, &tmp_with_alignments, sizeof(tmp_with_alignments) );
  _with_alignments = (tmp_with_alignments!=0);
  int8_t tmp_with_transfer;
  gzread( backup_file, &tmp_with_transfer, sizeof(tmp_with_transfer) );
  _with_transfer = (tmp_with_transfer!=0);
  int8_t tmp_swap_GUs;
  gzread( backup_file, &tmp_swap_GUs, sizeof(tmp_swap_GUs) );
  _swap_GUs = (tmp_swap_GUs!=0);
  gzread( backup_file, &_transfer_ins_rate,  sizeof(_transfer_ins_rate)  );
  gzread( backup_file, &_transfer_repl_rate, sizeof(_transfer_repl_rate) );

  // ------------------------------ Rearrangement rates (without alignements)
  gzread( backup_file, &_duplication_rate,   sizeof(_duplication_rate)   );
  gzread( backup_file, &_deletion_rate,      sizeof(_deletion_rate)      );
  gzread( backup_file, &_translocation_rate, sizeof(_translocation_rate) );
  gzread( backup_file, &_inversion_rate,     sizeof(_inversion_rate)     );

  // --------------------------------- Rearrangement rates (with alignements)
  gzread( backup_file, &_neighbourhood_rate,       sizeof(_neighbourhood_rate)       );
  gzread( backup_file, &_duplication_proportion,   sizeof(_duplication_proportion)   );
  gzread( backup_file, &_deletion_proportion,      sizeof(_deletion_proportion)      );
  gzread( backup_file, &_translocation_proportion, sizeof(_translocation_proportion) );
  gzread( backup_file, &_inversion_proportion,     sizeof(_inversion_proportion)     );

  // ------------------------------------------------------------ Alignements
  int8_t tmp_align_fun_shape;
  gzread( backup_file, &tmp_align_fun_shape, sizeof(tmp_align_fun_shape) );
  _align_fun_shape = (ae_align_fun_shape) tmp_align_fun_shape;
  gzread( backup_file, &_align_sigm_lambda,   sizeof(_align_sigm_lambda)   );
  gzread( backup_file, &_align_sigm_mean,     sizeof(_align_sigm_mean)     );
  gzread( backup_file, &_align_lin_min,       sizeof(_align_lin_min)       );
  gzread( backup_file, &_align_lin_max,       sizeof(_align_lin_max)       );
  gzread( backup_file, &_align_max_shift,     sizeof(_align_max_shift)     );
  gzread( backup_file, &_align_w_zone_h_len,  sizeof(_align_w_zone_h_len)  );
  gzread( backup_file, &_align_match_bonus,   sizeof(_align_match_bonus)   );
  gzread( backup_file, &_align_mismatch_cost, sizeof(_align_mismatch_cost) );
  
  // -------------------------------------------------------------- Selection
  if ( verbose )
  {
    printf( "    Loading selection parameters\n" );
  }
  int8_t tmp_selection_scheme;
  gzread( backup_file, &tmp_selection_scheme, sizeof(tmp_selection_scheme) );
  _selection_scheme = (ae_selection_scheme) tmp_selection_scheme;
  gzread( backup_file, &_selection_pressure, sizeof(_selection_pressure) );
  

  // -------------------------------------------------------------- Secretion
  if ( verbose )
  {
    printf( "    Loading secretion parameters\n" );
  }
  int8_t tmp_use_secretion;
  gzread( backup_file, &tmp_use_secretion, sizeof(tmp_use_secretion) );
  _use_secretion = (tmp_use_secretion!=0);
  gzread( backup_file, &_secretion_fitness_contrib,  sizeof(_secretion_fitness_contrib)  );
  gzread( backup_file, &_secretion_diffusion_prop,   sizeof(_secretion_diffusion_prop)   );
  gzread( backup_file, &_secretion_degradation_prop, sizeof(_secretion_degradation_prop) );
  gzread( backup_file, &_secretion_cost,             sizeof(_secretion_cost)             );

  // --------------------------------------------------------------- Plasmids
  if ( verbose )
  {
    printf( "    Loading plasmid management parameters\n" );
  }
  int8_t tmp_allow_plasmids;
  gzread( backup_file, &tmp_allow_plasmids, sizeof(tmp_allow_plasmids) );
  _allow_plasmids = (tmp_allow_plasmids!=0);
  gzread( backup_file, &_plasmid_minimal_length, sizeof(_plasmid_minimal_length) );
  gzread( backup_file, &_prob_horiz_trans,       sizeof(_prob_horiz_trans)       );
  gzread( backup_file, &_nb_horiz_trans,         sizeof(_nb_horiz_trans)         );
  
  int8_t tmp_compute_phen_contrib_by_GU;
  gzread( backup_file, &tmp_compute_phen_contrib_by_GU, sizeof(tmp_compute_phen_contrib_by_GU) );
  _compute_phen_contrib_by_GU = (tmp_compute_phen_contrib_by_GU != 0);
  
  // ------------------------------------------------------- Translation cost
  gzread( backup_file, &_translation_cost, sizeof(_translation_cost) );
}

void ae_params::print_to_file( FILE* file )
{
  // ------------------------------------------------------------ Constraints
  fprintf( file, "\nConstraints ---------------------------------------------\n" );
  fprintf( file, "min_genome_length :          %"PRId32"\n", _min_genome_length       );
  fprintf( file, "max_genome_length :          %"PRId32"\n", _max_genome_length       );
  fprintf( file, "MIN_W :                      %f\n",        _min_w                   );
  fprintf( file, "MAX_W :                      %f\n",        _max_w                   );
  
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
  fprintf( file, "with_transfer :              %s\n",  _with_transfer? "true" : "false"   );
  fprintf( file, "swap_GUs :                   %s\n",  _swap_GUs? "true" : "false"   );
  fprintf( file, "transfer_ins_rate :          %e\n",  _transfer_ins_rate );
  fprintf( file, "transfer_repl_rate :         %e\n",  _transfer_repl_rate );

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
  fprintf( file, "use_secretion :              %s\n", _use_secretion? "true" : "false" );
  fprintf( file, "secretion_fitness_contrib :  %e\n", _secretion_fitness_contrib       );
  fprintf( file, "secretion_diffusion_prop :   %e\n", _secretion_diffusion_prop        );
  fprintf( file, "secretion_degradation_prop : %e\n", _secretion_degradation_prop      );
  fprintf( file, "secretion_cost :             %e\n", _secretion_cost                  );

  // --------------------------------------------------------------- Plasmids
  fprintf( file, "\nPlasmids ------------------------------------------------\n" );
  fprintf( file, "allow_plasmids :             %s\n", _allow_plasmids? "true" : "false"              );
  fprintf( file, "plasmid_minimal_length :     %"PRId32"\n", _plasmid_minimal_length                 );
  fprintf( file, "prob_horiz_trans :           %e\n", _prob_horiz_trans                              );
  fprintf( file, "nb_horiz_trans :             %d\n", _nb_horiz_trans                                );
  fprintf( file, "compute_phen_contrib_by_GU : %s\n", _compute_phen_contrib_by_GU? "true" : "false"  );

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
