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


#ifndef __param_loader_H__
#define  __param_loader_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>


// =================================================================
//                            Project Files
// =================================================================
#include <params.h>
#include <ae_params_mut.h>
#include <ae_rand_mt.h>


// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;
class ae_environment;
class ae_individual;

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
  GRID_WIDTH,
  GRID_HEIGHT,
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
  SWAP_GUS,
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
  SECRETION_CONTRIB_TO_FITNESS,
  SECRETION_DIFFUSION_PROP,
  SECRETION_DEGRADATION_PROP, 
  SECRETION_INITIAL,
  SECRETION_COST,
  ALLOW_PLASMIDS,
  PLASMID_INITIAL_LENGTH,
  PLASMID_INITIAL_GENE,
  PLASMID_MINIMAL_LENGTH,
  PROB_PLASMID_HT,
  NB_PLASMID_HT,
  COMPUTE_PHEN_CONTRIB_BY_GU,
  LOG,
  GENERATION_OVERLOAD,

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




class param_loader
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    param_loader( const char* file_name );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~param_loader( void );

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void read_file( void );
    void load( ae_exp_manager* exp_m, bool verbose = false );
    
    f_line* get_line( void ); 
    
    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    param_loader( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    param_loader( const param_loader &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                           Protected Methods
    // =================================================================
    static void format_line( f_line*, char*, bool* );
    void interpret_line( f_line* line, int32_t cur_line );
    ae_individual* create_random_individual( ae_exp_manager* exp_m, ae_params_mut* param_mut, int32_t id ) const;
    ae_individual* create_random_individual_with_good_gene( ae_exp_manager* exp_m, ae_params_mut* param_mut, int32_t id ) const;
    ae_individual* create_clone( ae_individual* dolly, int32_t id ) const;
    
    

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    //~ ae_exp_manager* _exp_m;
    
    ae_rand_mt* _alea;
    
    char*   _param_file_name;
    FILE*   _param_file;
    
    params* _param_values;
    
    int32_t _cur_line;
};


// =====================================================================
//                          Accessors definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

#endif // __param_loader_H__
