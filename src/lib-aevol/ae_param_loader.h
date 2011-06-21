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


#ifndef __AE_PARAM_LOADER_H__
#define  __AE_PARAM_LOADER_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>


// =================================================================
//                            Project Files
// =================================================================
#include <ae_common.h>


// =================================================================
//                          Class declarations
// =================================================================
class ae_simulation;

typedef enum
{
  MIN_TRIANGLE_WIDTH,
  MAX_TRIANGLE_WIDTH,
  ENV_AXIS_SEGMENTS,
  ENV_AXIS_FEATURES,
  ENV_SEPARATE_SEGMENTS,
  BACKUP_STEP,
  BIG_BACKUP_STEP,
  RECORD_TREE,
  TREE_MODE,
  TREE_STEP,
  MORE_STATS,
  DUMP_PERIOD,
  INITIAL_GENOME_LENGTH,
  MIN_GENOME_LENGTH,
  MAX_GENOME_LENGTH,
  INIT_POP_SIZE,
  NB_GENER,
  POP_STRUCTURE,
  MIGRATION_NUMBER, 
  GRID_X,
  GRID_Y,
  INIT_METHOD,
  POINT_MUTATION_RATE,
  SMALL_INSERTION_RATE,
  SMALL_DELETION_RATE,
  MAX_INDEL_SIZE,
  WITH_4PTS_TRANS,
  WITH_ALIGNMENTS,
  DUPLICATION_RATE,
  DELETION_RATE,
  TRANSLOCATION_RATE,
  INVERSION_RATE,
  WITH_TRANSFER,
  TRANSFER_INS_RATE,
  TRANSFER_REPL_RATE,
  NEIGHBOURHOOD_RATE,
  DUPLICATION_PROPORTION,
  DELETION_PROPORTION,
  TRANSLOCATION_PROPORTION,
  INVERSION_PROPORTION,
  ALIGN_FUNCTION,
  ALIGN_MAX_SHIFT,
  ALIGN_W_ZONE_H_LEN,
  ALIGN_MATCH_BONUS,
  ALIGN_MISMATCH_COST,
  SELECTION_SCHEME,
  SELECTION_PRESSURE,
  SEED,
  TRANSLATION_COST,
  ENV_ADD_POINT,
  ENV_ADD_GAUSSIAN,
  ENV_SAMPLING,
  ENV_VARIATION,
  ENV_SEED,
  SECRETION_FITNESS_CONTRIB,
  SECRETION_DIFUSION_PROP,
  SECRETION_DEGRADATION_PROP, 
  SECRETION_INITIAL,
  SECRETION_COST,
  ALLOW_PLASMIDS,
  PLASMID_INITIAL_LENGTH,
  PLASMID_INITIAL_GENE,
  PLASMID_MINIMAL_LENGTH,
  PROB_HORIZONTAL_TRANS,
  NB_HORIZONTAL_TRANS,
  COMPUTE_PHEN_CONTRIB_BY_GU,
  LOG,

#ifdef __REGUL
  HILL_SHAPE_N,
  HILL_SHAPE_THETA,
  DEGRADATION_RATE,
  DEGRADATION_STEP,
  INDIVIDUAL_EVALUATION_DATES,
  INDIVIDUAL_LIFE_TIME,
  BINDING_ZEROS_PERCENTAGE,
  WITH_HEREDITY,
  PROTEIN_PRESENCE_LIMIT,
#endif

  UNDEFINED
} ae_keywd;


class f_line
{
  public :
    f_line( void );

    ae_keywd get_keywd( void );

    int16_t nb_words;
    char    words[50][255];
};




class ae_param_loader
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_param_loader( void );
    ae_param_loader( const char* file_name );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_param_loader( void );

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    static void format_line( f_line*, char*, bool* );
    static void interpret_line( f_line* line, int32_t cur_line );
    void load( void );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ ae_param_loader( void )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    ae_param_loader( const ae_param_loader &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                           Protected Methods
    // =================================================================
    f_line* get_line( void );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    FILE*   param_in;
    int32_t cur_line;
};


// =====================================================================
//                          Accessors definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

#endif // __AE_PARAM_LOADER_H__
