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


/*! \class ae_params_mut
    \brief
*/
 
 
// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include <ae_params_mut.h>




//##############################################################################
//                                                                             #
//                             Class ae_params_mut                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_params_mut::ae_params_mut( void )
{
  // --------------------------------------------------------- Mutation rates
  _point_mutation_rate  = 0.0;
  _small_insertion_rate = 0.0;
  _small_deletion_rate  = 0.0;
  _max_indel_size       = 0;

  // -------------------------------------------- Rearrangements and Transfer
  _with_4pts_trans  = false;
  _with_alignments  = false;
  _with_HT          = false;
  _swap_GUs         = false;
  _HT_ins_rate      = 0.0;
  _HT_repl_rate     = 0.0;

  // ------------------------------ Rearrangement rates (without alignements)
  _duplication_rate   = 0.0;
  _deletion_rate      = 0.0;
  _translocation_rate = 0.0;
  _inversion_rate     = 0.0;

  // --------------------------------- Rearrangement rates (with alignements)
  _neighbourhood_rate = 0.0;
  
  _duplication_proportion   = 0.0;
  _deletion_proportion      = 0.0;
  _translocation_proportion = 0.0;
  _inversion_proportion     = 0.0;
}

ae_params_mut::ae_params_mut( const ae_params_mut& model )
{
  // --------------------------------------------------------- Mutation rates
  _point_mutation_rate  = model._point_mutation_rate;
  _small_insertion_rate = model._small_insertion_rate;
  _small_deletion_rate  = model._small_deletion_rate;
  _max_indel_size       = model._max_indel_size;

  // -------------------------------------------- Rearrangements and Transfer
  _with_4pts_trans  = model._with_4pts_trans;
  _with_alignments  = model._with_alignments;
  _with_HT          = model._with_HT;
  _swap_GUs         = model._swap_GUs;
  _HT_ins_rate      = model._HT_ins_rate;
  _HT_repl_rate     = model._HT_repl_rate;

  // ------------------------------ Rearrangement rates (without alignements)
  _duplication_rate   = model._duplication_rate;
  _deletion_rate      = model._deletion_rate;
  _translocation_rate = model._translocation_rate;
  _inversion_rate     = model._inversion_rate;

  // --------------------------------- Rearrangement rates (with alignements)
  _neighbourhood_rate = model._neighbourhood_rate;
  
  _duplication_proportion   = model._duplication_proportion;
  _deletion_proportion      = model._deletion_proportion;
  _translocation_proportion = model._translocation_proportion;
  _inversion_proportion     = model._inversion_proportion;
}

ae_params_mut::ae_params_mut( gzFile backup_file )
{
  // --------------------------------------------------------- Mutation rates
  gzread( backup_file, &_point_mutation_rate,  sizeof(_point_mutation_rate) );
  gzread( backup_file, &_small_insertion_rate, sizeof(_small_insertion_rate) );
  gzread( backup_file, &_small_deletion_rate,  sizeof(_small_deletion_rate) );
  gzread( backup_file, &_max_indel_size,       sizeof(_max_indel_size) );

  // -------------------------------------------- Rearrangements and Transfer
  int8_t tmp;
  gzread( backup_file, &tmp, sizeof(tmp) );
  _with_4pts_trans = (tmp != 0);
  gzread( backup_file, &tmp, sizeof(tmp) );
  _with_alignments = (tmp != 0);
  gzread( backup_file, &tmp, sizeof(tmp) );
  _with_HT = (tmp != 0);
  gzread( backup_file, &tmp, sizeof(tmp) );
  _swap_GUs = (tmp != 0);
  gzread( backup_file, &_HT_ins_rate,  sizeof(_HT_ins_rate) );
  gzread( backup_file, &_HT_repl_rate, sizeof(_HT_repl_rate) );

  // ------------------------------ Rearrangement rates (without alignements)
  gzread( backup_file, &_duplication_rate,   sizeof(_duplication_rate) );
  gzread( backup_file, &_deletion_rate,      sizeof(_deletion_rate) );
  gzread( backup_file, &_translocation_rate, sizeof(_translocation_rate) );
  gzread( backup_file, &_inversion_rate,     sizeof(_inversion_rate) );

  // --------------------------------- Rearrangement rates (with alignements)
  gzread( backup_file, &_neighbourhood_rate,       sizeof(_neighbourhood_rate) );
  gzread( backup_file, &_duplication_proportion,   sizeof(_duplication_proportion) );
  gzread( backup_file, &_deletion_proportion,      sizeof(_deletion_proportion) );
  gzread( backup_file, &_translocation_proportion, sizeof(_translocation_proportion) );
  gzread( backup_file, &_inversion_proportion,     sizeof(_inversion_proportion) );
}

// =================================================================
//                             Destructors
// =================================================================
ae_params_mut::~ae_params_mut( void )
{
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_params_mut::save( gzFile backup_file ) const
{
  // --------------------------------------------------------- Mutation rates
  gzwrite( backup_file, &_point_mutation_rate,  sizeof(_point_mutation_rate) );
  gzwrite( backup_file, &_small_insertion_rate, sizeof(_small_insertion_rate) );
  gzwrite( backup_file, &_small_deletion_rate,  sizeof(_small_deletion_rate) );
  gzwrite( backup_file, &_max_indel_size,       sizeof(_max_indel_size) );

  // -------------------------------------------- Rearrangements and Transfer
  int8_t tmp = _with_4pts_trans ? 1 : 0;
  gzwrite( backup_file, &tmp,  sizeof(tmp) );
  tmp = _with_alignments ? 1 : 0;
  gzwrite( backup_file, &tmp,  sizeof(tmp) );
  tmp = _with_HT ? 1 : 0;
  gzwrite( backup_file, &tmp,  sizeof(tmp) );
  tmp = _swap_GUs ? 1 : 0;
  gzwrite( backup_file, &tmp,  sizeof(tmp) );
  gzwrite( backup_file, &_HT_ins_rate,  sizeof(_HT_ins_rate) );
  gzwrite( backup_file, &_HT_repl_rate, sizeof(_HT_repl_rate) );

  // ------------------------------ Rearrangement rates (without alignements)
  gzwrite( backup_file, &_duplication_rate,   sizeof(_duplication_rate) );
  gzwrite( backup_file, &_deletion_rate,      sizeof(_deletion_rate) );
  gzwrite( backup_file, &_translocation_rate, sizeof(_translocation_rate) );
  gzwrite( backup_file, &_inversion_rate,     sizeof(_inversion_rate) );

  // --------------------------------- Rearrangement rates (with alignements)
  gzwrite( backup_file, &_neighbourhood_rate,       sizeof(_neighbourhood_rate) );
  gzwrite( backup_file, &_duplication_proportion,   sizeof(_duplication_proportion) );
  gzwrite( backup_file, &_deletion_proportion,      sizeof(_deletion_proportion) );
  gzwrite( backup_file, &_translocation_proportion, sizeof(_translocation_proportion) );
  gzwrite( backup_file, &_inversion_proportion,     sizeof(_inversion_proportion) );
}

// =================================================================
//                           Protected Methods
// =================================================================

// =================================================================
//                          Non inline accessors
// =================================================================
