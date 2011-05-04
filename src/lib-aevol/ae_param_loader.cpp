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
#include <string.h>
#include <assert.h>

// =================================================================
//                            Project Files
// =================================================================
#include <ae_common.h>
#include <ae_param_loader.h>

#include <ae_rand_mt.h>
#include <ae_simulation.h>
#include <ae_gaussian.h>
#include <ae_env_segment.h>
#include <ae_point_2d.h>

#ifdef __REGUL
  #include <ae_array_short.h>
#endif

#define INPUT_FILE_NAME   "param.in"
#define OUTPUT_FILE_NAME  "param.out"




//##############################################################################
//                                                                             #
//                           Class ae_param_loader                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_param_loader::ae_param_loader( void )
{
  param_in  = fopen( "param.in",  "r" );
  if ( param_in == NULL )
  {
    printf( "ERROR : couldn't open file param.in\n" );
    exit( EXIT_FAILURE );
  }
  assert( param_in );
  cur_line = 0;
}

// =================================================================
//                             Destructors
// =================================================================
ae_param_loader::~ae_param_loader( void )
{
  fclose( param_in );
}

// =================================================================
//                            Public Methods
// =================================================================

void ae_param_loader::interpret_line( f_line* line, int32_t cur_line )
{
  // Interpret line
  switch ( line->get_keywd() )
  {
    case UNDEFINED :
    {
      printf( "ERROR in param file \"%s\" on line %"PRId32" : undefined key word \"%s\"\n",
              INPUT_FILE_NAME, cur_line, line->words[0] );
      exit( EXIT_FAILURE );
    }
    case MIN_TRIANGLE_WIDTH :
    {
      MIN_W = atof( line->words[1] );
      break;
    }
    case MAX_TRIANGLE_WIDTH :
    {
      MAX_W = atof( line->words[1] );
      break;
    }
    case ENV_AXIS_SEGMENTS :
    {
      // Check consistency
      if ( ae_common::env_axis_segment_boundaries != NULL )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : ENV_AXIS_SEGMENT keyword found twice.\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE );
      }
      if ( ae_common::env_axis_features != NULL && ae_common::env_axis_nb_segments != line->nb_words )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : number of segment boundaries and number of features don't match.\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE );
      }
      
      // OK, proceed
      ae_common::env_axis_is_segmented = true;
      ae_common::env_axis_nb_segments = line->nb_words; // in fact line->nb_words - 1 + 1
      ae_common::env_axis_segment_boundaries = new double [ae_common::env_axis_nb_segments + 1];
      
      ae_common::env_axis_segment_boundaries[0] = MIN_X;
      ae_common::env_axis_segment_boundaries[ae_common::env_axis_nb_segments] = MAX_X;
      for ( int16_t i = 1 ; i < ae_common::env_axis_nb_segments ; i++ )
      {
        ae_common::env_axis_segment_boundaries[i] = atof( line->words[i] );
      }
      
      break;
    }
    case ENV_AXIS_FEATURES :
    {
      // Check consistency
      if ( ae_common::env_axis_features != NULL )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : ENV_AXIS_FEATURES keyword found twice.\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE );
      }
      if ( ae_common::env_axis_segment_boundaries != NULL && ae_common::env_axis_nb_segments != line->nb_words - 1 )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : number of segment boundaries and number of features don't match.\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE );
      }
      
      // OK, proceed
      ae_common::env_axis_is_segmented = true;
      ae_common::env_axis_nb_segments = line->nb_words - 1;
      ae_common::env_axis_features = new ae_env_axis_feature[ae_common::env_axis_nb_segments];
      
      for ( int16_t i = 0 ; i < ae_common::env_axis_nb_segments ; i++ )
      {
        if ( strcmp( line->words[i+1], "NEUTRAL" ) == 0 )
        {
          ae_common::env_axis_features[i] = NEUTRAL;
        }
        else if ( strcmp( line->words[i+1], "METABOLISM" ) == 0 )
        {
          ae_common::env_axis_features[i] = METABOLISM;
        }
        else if ( strcmp( line->words[i+1], "SECRETION" ) == 0 )
        {
          ae_common::use_secretion      = true;
          ae_common::composite_fitness  = true;
          ae_common::env_axis_features[i] = SECRETION;
        }
        else if ( strcmp( line->words[i+1], "TRANSFER" ) == 0 )
        {
          ae_common::env_axis_features[i] = TRANSFER;
        }
        else
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown axis feature \"%s\".\n",
                  INPUT_FILE_NAME, cur_line, line->words[1] );
          exit( EXIT_FAILURE );
        }
      }
      
      break;
    }
    case ENV_SEPARATE_SEGMENTS :
    {
      ae_common::env_separate_segments = true;
      break;
    }
    case RECORD_TREE :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        ae_common::record_tree  = true;
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        ae_common::record_tree = false;
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown tree recording option (use true/false).\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE ); 
      }        
      break;
    }
    case TREE_MODE :
    {
      if ( strcmp( line->words[1], "light" ) == 0 )
      {
        ae_common::tree_mode  = LIGHT;
      }
      else if ( strcmp( line->words[1], "normal" ) == 0 )
      {
        ae_common::tree_mode  = NORMAL;
      }       
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown tree mode option (use normal/light).\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE );
      }        
      break;
    }
    case MORE_STATS :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        ae_common::more_stats  = true;
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        ae_common::more_stats = false;
      }       
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown more stats option (use true/false).\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE ); 
      }        
      break;
    }
    case DUMP_PERIOD :
    {
      ae_common::dump_period = atol( line->words[1] );
      break;
    }
    case BACKUP_STEP :
    {
      ae_common::backup_step = atol( line->words[1] );
      break;
    }
    case BIG_BACKUP_STEP :
    {
      ae_common::big_backup_step = atol( line->words[1] );
      break;
    }
    case TREE_STEP :
    {
      ae_common::tree_step = atol( line->words[1] );
      break;
    }
    case INITIAL_GENOME_LENGTH :
    {
      ae_common::initial_genome_length = atol( line->words[1] );
      break;
    }
    case MIN_GENOME_LENGTH :
    {
      if ( strncmp( line->words[1], "NONE", 4 ) == 0 )
      {
        ae_common::min_genome_length = 1; // Must not be 0
      }
      else
      {
        ae_common::min_genome_length = atol( line->words[1] );
        if ( ae_common::min_genome_length == 0 )
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : MIN_GENOME_LENGTH must be > 0.\n",
                  INPUT_FILE_NAME, cur_line );
          exit( EXIT_FAILURE ); 
        }
      }
      break;
    }
    case MAX_GENOME_LENGTH :
    {
      if ( strncmp( line->words[1], "NONE", 4 ) == 0 )
      {
        ae_common::max_genome_length = INT32_MAX;
      }
      else
      {
        ae_common::max_genome_length = atol( line->words[1] );
      }
      break;
    }
    case INIT_POP_SIZE :
    {
      ae_common::init_pop_size = atol( line->words[1] );
      break;      
    }
    case POP_STRUCTURE :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        ae_common::pop_structure  = true;
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        ae_common::pop_structure = false;
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown population structure (use true/false).\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE ); 
      }
      break;
    }
    case GRID_X :
    {
      ae_common::grid_x = atol( line->words[1] );
      break;
    }
    case GRID_Y :
    {
      ae_common::grid_y = atol( line->words[1] );
      break;      
    }
    case MIGRATION_NUMBER :
    {
      ae_common::migration_number = atoi( line->words[1] );
      break;
    }
    case INIT_METHOD :
    {
      ae_common::init_method = 0;
      
      for ( int8_t i = 1 ; i < line->nb_words ; i++ )
      {
        if ( strcmp( line->words[i], "ONE_GOOD_GENE" ) == 0 )
        {
          ae_common::init_method |= ONE_GOOD_GENE;
        }
        else if ( strcmp( line->words[i], "CLONE" ) == 0 )
        {
          ae_common::init_method |= CLONE;
        }
        else if ( strcmp( line->words[i], "WITH_INS_SEQ" ) == 0 )
        {
          ae_common::init_method |= WITH_INS_SEQ;
        }   
        else
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown initialization method %s.\n",
                  INPUT_FILE_NAME, cur_line, line->words[1] );
          exit( EXIT_FAILURE ); 
        }
      }
      break;
    }
    case POINT_MUTATION_RATE :
    {
      ae_common::point_mutation_rate = atof( line->words[1] );
      break;
    }
    case SMALL_INSERTION_RATE :
    {
      ae_common::small_insertion_rate = atof( line->words[1] );
      break;
    }
    case SMALL_DELETION_RATE :
    {
      ae_common::small_deletion_rate = atof( line->words[1] );
      break;
    }
    case MAX_INDEL_SIZE :
    {
      ae_common::max_indel_size = atol( line->words[1] );
      break;
    }
    case DUPLICATION_RATE :
    {
      ae_common::duplication_rate = atof( line->words[1] );
      break;
    }
    case DELETION_RATE :
    {
      ae_common::deletion_rate = atof( line->words[1] );
      break;
    }
    case TRANSLOCATION_RATE :
    {
      ae_common::translocation_rate = atof( line->words[1] );
      break;
    }
    case INVERSION_RATE :
    {
      ae_common::inversion_rate = atof( line->words[1] );
      break;
    }
    case NEIGHBOURHOOD_RATE :
    {
      ae_common::neighbourhood_rate = atof( line->words[1] );
      break;
    }
    case DUPLICATION_PROPORTION :
    {
      ae_common::duplication_proportion = atof( line->words[1] );
      break;
    }
    case DELETION_PROPORTION :
    {
      ae_common::deletion_proportion = atof( line->words[1] );
      break;
    }
    case TRANSLOCATION_PROPORTION :
    {
      ae_common::translocation_proportion = atof( line->words[1] );
      break;
    }
    case INVERSION_PROPORTION :
    {
      ae_common::inversion_proportion = atof( line->words[1] );
      break;
    }
    case ALIGN_FUNCTION :
    {
      if ( line->nb_words != 2 && line->nb_words != 4 )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : incorrect number of parameters for keyword \"%s\".\n",
                INPUT_FILE_NAME, cur_line, line->words[0] );
        exit( EXIT_FAILURE );
      }
      
      if ( strcmp( line->words[1], "LINEAR" ) == 0 )
      {
        ae_common::align_fun_shape = LINEAR;
        
        if ( line->nb_words == 4 )
        {
          ae_common::align_lin_min  = atol( line->words[2] );
          ae_common::align_lin_max  = atol( line->words[3] );
        }
      }
      else if ( strcmp( line->words[1], "SIGMOID" ) == 0 )
      {
        ae_common::align_fun_shape = SIGMOID;
        
        if ( line->nb_words == 4 )
        {
          ae_common::align_sigm_lambda  = atol( line->words[2] );
          ae_common::align_sigm_mean    = atol( line->words[3] );
        }
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown align function shape \"%s\".\n",
                INPUT_FILE_NAME, cur_line, line->words[1] );
        exit( EXIT_FAILURE );
      }
      break;
    }
    case ALIGN_MAX_SHIFT :
    {
      ae_common::align_max_shift = atol( line->words[1] );
      break;
    }
    case ALIGN_W_ZONE_H_LEN :
    {
      ae_common::align_w_zone_h_len = atol( line->words[1] );
      break;
    }
    case ALIGN_MATCH_BONUS :
    {
      ae_common::align_match_bonus = atol( line->words[1] );
      break;
    }
    case ALIGN_MISMATCH_COST :
    {
      ae_common::align_mismatch_cost = atol( line->words[1] );
      break;
    }
    case SELECTION_SCHEME :
    {
      if ( strncmp( line->words[1], "lin", 3 ) == 0 )
      {
        ae_common::selection_scheme = RANK_LINEAR;
      }
      else if ( strncmp( line->words[1], "exp", 3 ) == 0 )
      {
        ae_common::selection_scheme = RANK_EXPONENTIAL;
      }
      else if ( strncmp( line->words[1], "fitness", 7 ) == 0 )
      {
        ae_common::selection_scheme = FITNESS_PROPORTIONATE;
      }
      else if ( strcmp( line->words[1], "fittest" ) == 0 )
      {
        ae_common::selection_scheme = FITTEST;
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown selection scheme \"%s\".\n",
                INPUT_FILE_NAME, cur_line, line->words[1] );
        exit( EXIT_FAILURE );
      }
      break;
    }
    case SELECTION_PRESSURE :
    {
      ae_common::selection_pressure = atof( line->words[1] );
      break;
    }
    case SEED :
    {
      ae_common::seed = atol( line->words[1] );
      break;
    }
    case NB_GENER :
    {
      ae_common::nb_generations = atol( line->words[1] );
      break;
    }
    case WITH_4PTS_TRANS :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        ae_common::with_4pts_trans = true;
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        printf( "ERROR: 3 points rear hasn't been implemented yet\n" );
        exit( EXIT_FAILURE );
        ae_common::with_4pts_trans = false;
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown 4pts_trans option (use true/false).\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE ); 
      }
      break;
    }
    case WITH_ALIGNMENTS :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        ae_common::with_alignments = true;
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        ae_common::with_alignments = false;
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown alignement option (use true/false).\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE ); 
      }
      break;
    }
    case WITH_TRANSFER :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        ae_common::with_transfer = true;
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        ae_common::with_transfer = false;
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown transfer option (use true/false).\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE ); 
      }
      break;
    }
    case TRANSFER_INS_RATE :
    {
      ae_common::transfer_ins_rate  = atof( line->words[1] );
      break;
    }
    case TRANSFER_REPL_RATE :
    {
      ae_common::transfer_repl_rate  = atof( line->words[1] );
      break;
    }
    case TRANSLATION_COST :
    {
      ae_common::translation_cost = atof( line->words[1] );
      break;
    }
    case ENV_ADD_POINT :
    {
      ae_common::env_custom_points.add( new ae_point_2d(  atof( line->words[1] ), atof( line->words[2] ) ) );
      break;
    }
    case ENV_ADD_GAUSSIAN :
    {
      ae_common::env_gaussians.add( new ae_gaussian(  atof( line->words[1] ),
                                                      atof( line->words[2] ),
                                                      atof( line->words[3] ) ) );
      break;
    }
    case ENV_SAMPLING :
    {
      ae_common::env_sampling = atoi( line->words[1] );
      break;
    }
    case ENV_VARIATION :
    {
      if ( strcmp( line->words[1], "none" ) == 0 )
      {
        ae_common::env_var_method = NONE;
      }
      else if ( strcmp( line->words[1], "autoregressive_mean_variation" ) == 0 )
      {
        ae_common::env_var_method = AUTOREGRESSIVE_MEAN_VAR;
        ae_common::env_sigma      = atof( line->words[2] );
        ae_common::env_tau        = atol( line->words[3] );
      }
      else if ( strcmp( line->words[1], "add_local_gaussians" ) == 0 )
      {
        ae_common::env_var_method = LOCAL_GAUSSIANS_VAR;
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown environment variation method.\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE );
      }
      break;
    }
    case ENV_SEED :
    {
      static bool seed_already_set = false;
      if ( seed_already_set == true )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : duplicate entry for ENV_SEED.\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE );
      }
      ae_common::env_seed = atoi( line->words[1] );
      seed_already_set = true;
      break;
    }
    case SECRETION_FITNESS_CONTRIB :
    {
      ae_common::secretion_fitness_contrib = atof( line->words[1] );
      break; 
    }
    case SECRETION_DIFUSION_PROP :
    {
      ae_common::secretion_difusion_prop = atof( line->words[1] );
      break; 
    }
    case SECRETION_DEGRADATION_PROP :
    {
      ae_common::secretion_degradation_prop = atof( line->words[1] );
      if ( ae_common::secretion_degradation_prop > 1 || ae_common::secretion_degradation_prop < 0 ) 
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : degradation must be in (0,1).\n", INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE );
      }
      break; 
    }
    case SECRETION_INITIAL :
    {
      ae_common::secretion_init = atof( line->words[1] );
      break; 
    }
    case SECRETION_COST :
    {
      ae_common::secretion_cost = atof( line->words[1] );
      break; 
    }
    case ALLOW_PLASMIDS :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        ae_common::allow_plasmids = true;
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        ae_common::allow_plasmids = false;
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown allow_plasmids option (use true/false).\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE ); 
      }
      break; 
    }
    case PLASMID_INITIAL_LENGTH :
    {
      ae_common::plasmid_initial_length = atoi( line->words[1] );
      break; 
    }
    case PLASMID_INITIAL_GENE :
    {
      ae_common::plasmid_initial_gene = atoi( line->words[1] );
      break; 
    }
    case PLASMID_MINIMAL_LENGTH :
    {
      ae_common::plasmid_minimal_length = atoi( line->words[1] );
      break; 
    }
    case PROB_HORIZONTAL_TRANS :
    {
      ae_common::prob_horiz_trans = atof( line->words[1] );
      break; 
    }
    case COMPUTE_PHEN_CONTRIB_BY_GU :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        ae_common::compute_phen_contrib_by_GU = true;
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        ae_common::compute_phen_contrib_by_GU = false;
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown compute_phen_contrib_by_GU option (use true/false).\n",
                INPUT_FILE_NAME, cur_line );
        exit( EXIT_FAILURE ); 
      }
      break; 
    }
    case LOG :
    {
      ae_common::logs = 0;
      
      for ( int8_t i = 1 ; i < line->nb_words ; i++ )
      {
        if ( strcmp( line->words[i], "TRANSFER" ) == 0 )
        {
          ae_common::logs |= LOG_TRANSFER;
        }
        else if ( strcmp( line->words[i], "REAR" ) == 0 )
        {
          ae_common::logs |= LOG_REAR;
        }
        else if ( strcmp( line->words[i], "BARRIER" ) == 0 )
        {
          ae_common::logs |= LOG_BARRIER;
        }
        else if ( strcmp( line->words[i], "LOADS" ) == 0 )
        {
          ae_common::logs |= LOG_LOADS;
        }   
        else
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown log option %s.\n",
                  INPUT_FILE_NAME, cur_line, line->words[1] );
          exit( EXIT_FAILURE ); 
        }
      }
      break;
    }
    
    #ifdef __REGUL
      case HILL_SHAPE_N :
      {
        ae_common::hill_shape_n = atof( line->words[1] );
        break;
      }
      case HILL_SHAPE_THETA :
      {
        ae_common::hill_shape_theta = atof( line->words[1] );
        break;
      }
      case DEGRADATION_RATE :
      {
        ae_common::degradation_rate = atof( line->words[1] );
        break;
      }
      case DEGRADATION_STEP :
      {
        ae_common::degradation_step = atof( line->words[1] );
	// we have to test if 1/degradation_step is an integer
        int i = static_cast<int>(1/ae_common::degradation_step);
        if( 1/ae_common::degradation_step != i )
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : DEGRADATION STEP\n", INPUT_FILE_NAME, cur_line );
          printf( "This step have to divide 1.\n");
          exit( EXIT_FAILURE );
        }
        break;
      }
      case INDIVIDUAL_EVALUATION_DATES :
      {
        ae_common::individual_evaluation_nbr    = line->nb_words - 1;
        if( !ae_common::individual_evaluation_nbr )
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : number of dates\n", INPUT_FILE_NAME, cur_line );
          exit( EXIT_FAILURE );
        }
        ae_common::individual_evaluation_dates  = new ae_array_short( ae_common::individual_evaluation_nbr );
        for( int16_t i = 0; i < ae_common::individual_evaluation_nbr; i++)
        {
          ae_common::individual_evaluation_dates->set_value( i, atoi( line->words[1 + i] ) );
        }
        ae_common::individual_evaluation_dates->sort();
        break;
      }
      case BINDING_ZEROS_PERCENTAGE :
      {
        ae_common::binding_zeros_percentage = atof( line->words[1] );
        break;
      }
      case WITH_HEREDITY :
      {
        if ( strncmp( line->words[1], "true", 4 ) == 0 )
        {
          ae_common::with_heredity = true;
        }
        else if ( strncmp( line->words[1], "false", 5 ) == 0 )
        {
          ae_common::with_heredity = false;
        }
        else
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown with_heredity option (use true/false).\n",
                  INPUT_FILE_NAME, cur_line );
          exit( EXIT_FAILURE ); 
        }
        
        break;
      }
      case PROTEIN_PRESENCE_LIMIT :
      {
        ae_common::protein_presence_limit = atof( line->words[1] );
        break;
      }
    #endif
    
    
    default :
      printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown error\n", INPUT_FILE_NAME, cur_line );
      exit( EXIT_FAILURE );
      break;
  }
}

void ae_param_loader::load( void )
{
  // Thi rewind are only necessary when using multiple param files
  rewind( param_in );

  cur_line = 0;
  f_line* line;

  while ( ( line = get_line() ) != NULL ) // TODO : write line = new f_line( param_in ) => f_line::f_line( char* )
  {
    interpret_line( line, cur_line );
    delete line;
  }
  
  
  //  printf("Params loaded !\n");
// 
//   // Check mandatory parameters
//   #ifdef __REGUL
//     if( !ae_common::individual_evaluation_dates )
//     {
//       printf( "ERROR in param file \"%s\" : you must specifiy dates of evaluation\n", INPUT_FILE_NAME );
//       exit( EXIT_FAILURE );
//     }
//     if( !ae_common::individual_life_time )
//     {
//       printf( "ERROR in param file \"%s\" : you must specifiy a life time for individuals\n", INPUT_FILE_NAME );
//       exit( EXIT_FAILURE );
//     }
//   #endif
}



// =================================================================
//                           Protected Methods
// =================================================================

void ae_param_loader::format_line( f_line* formated_line, char* line, bool* line_is_interpretable )
{
  int16_t i = 0;
  int16_t j;

  // Parse line
  while ( line[i] != '\n' && line[i] != '\0' )
  {
    //      printf("i = %d\n", i);
    j = 0;

    // Flush white spaces and tabs
    while ( line[i] == ' ' || line[i] == 0x09 ) i++; // 0x09 is the ASCII code for TAB

    // Check comments
    if ( line[i] == '#' ) break;

    // If we got this far, there is content in the line
    *line_is_interpretable = true;

    // Parse word
    while ( line[i] != ' '  && line[i] != '\n' && line[i] != '\0' )
    {
      formated_line->words[formated_line->nb_words][j++] = line[i++];
    }

    // Add '\0' at end of word
    formated_line->words[formated_line->nb_words++][j] = '\0';
  }
}

f_line* ae_param_loader::get_line( void )
{
  char line[255];
  f_line* formated_line = new f_line();

  bool found_interpretable_line = false; // Found line that is neither a comment nor empty

  while ( !feof( param_in ) && !found_interpretable_line )
  {
    if ( !fgets( line, 255, param_in ) )
    {
      delete formated_line;
      return NULL;
    }
    cur_line++;
    format_line( formated_line, line, &found_interpretable_line );
  }

  if ( found_interpretable_line )
  {
    return formated_line;
  }
  else
  {
    delete formated_line;
    return NULL;
  }
}




f_line::f_line( void )
{
  nb_words = 0;
}

ae_keywd f_line::get_keywd( void )
{
  if ( !strcmp( words[0], "MIN_TRIANGLE_WIDTH" ) )          return MIN_TRIANGLE_WIDTH;
  if ( !strcmp( words[0], "MAX_TRIANGLE_WIDTH" ) )          return MAX_TRIANGLE_WIDTH;
  if ( !strcmp( words[0], "ENV_AXIS_SEGMENTS" ) )           return ENV_AXIS_SEGMENTS;
  if ( !strcmp( words[0], "ENV_AXIS_FEATURES" ) )           return ENV_AXIS_FEATURES;
  if ( !strcmp( words[0], "ENV_SEPARATE_SEGMENTS" ) )       return ENV_SEPARATE_SEGMENTS;
  if ( !strcmp( words[0], "ALLOW_PLASMIDS" ) )              return ALLOW_PLASMIDS;
  if ( !strcmp( words[0], "INITIAL_GENOME_LENGTH" ) )       return INITIAL_GENOME_LENGTH;
  if ( !strcmp( words[0], "MIN_GENOME_LENGTH" ) )           return MIN_GENOME_LENGTH;
  if ( !strcmp( words[0], "MAX_GENOME_LENGTH" ) )           return MAX_GENOME_LENGTH;
  if ( !strcmp( words[0], "BACKUP_STEP" ) )                 return BACKUP_STEP;
  if ( !strcmp( words[0], "BIG_BACKUP_STEP" ) )             return BIG_BACKUP_STEP;
  if ( !strcmp( words[0], "TREE_STEP" ) )                   return TREE_STEP;
  if ( !strcmp( words[0], "RECORD_TREE" ) )                 return RECORD_TREE;
  if ( !strcmp( words[0], "TREE_MODE") )                    return TREE_MODE;
  if ( !strcmp( words[0], "MORE_STATS") )                   return MORE_STATS;
  if ( !strcmp( words[0], "DUMP_PERIOD" ) )                 return DUMP_PERIOD;
  if ( !strcmp( words[0], "INIT_POP_SIZE" ) )               return INIT_POP_SIZE;
  if ( !strcmp( words[0], "NB_GENER" ) )                    return NB_GENER;
  if ( !strcmp( words[0], "POP_STRUCTURE" ) )               return POP_STRUCTURE;
  if ( !strcmp( words[0], "GRID_X" ) )                      return GRID_X;
  if ( !strcmp( words[0], "GRID_Y" ) )                      return GRID_Y;
  if ( !strcmp( words[0], "MIGRATION_NUMBER" ) )            return MIGRATION_NUMBER;
  if ( !strcmp( words[0], "INIT_METHOD" ) )                 return INIT_METHOD;
  if ( !strcmp( words[0], "POINT_MUTATION_RATE" ) )         return POINT_MUTATION_RATE;
  if ( !strcmp( words[0], "SMALL_INSERTION_RATE" ) )        return SMALL_INSERTION_RATE;
  if ( !strcmp( words[0], "SMALL_DELETION_RATE" ) )         return SMALL_DELETION_RATE;
  if ( !strcmp( words[0], "MAX_INDEL_SIZE" ) )              return MAX_INDEL_SIZE;
  if ( !strcmp( words[0], "WITH_4PTS_TRANS" ) )             return WITH_4PTS_TRANS;
  if ( !strcmp( words[0], "WITH_ALIGNMENTS" ) )             return WITH_ALIGNMENTS;
  if ( !strcmp( words[0], "DUPLICATION_RATE" ) )            return DUPLICATION_RATE;
  if ( !strcmp( words[0], "DELETION_RATE" ) )               return DELETION_RATE;
  if ( !strcmp( words[0], "TRANSLOCATION_RATE" ) )          return TRANSLOCATION_RATE;
  if ( !strcmp( words[0], "INVERSION_RATE" ) )              return INVERSION_RATE;
  if ( !strcmp( words[0], "WITH_TRANSFER" ) )               return WITH_TRANSFER;
  if ( !strcmp( words[0], "TRANSFER_INS_RATE" ) )           return TRANSFER_INS_RATE;
  if ( !strcmp( words[0], "TRANSFER_REPL_RATE" ) )          return TRANSFER_REPL_RATE;
  if ( !strcmp( words[0], "NEIGHBOURHOOD_RATE" ) )          return NEIGHBOURHOOD_RATE;
  if ( !strcmp( words[0], "DUPLICATION_PROPORTION" ) )      return DUPLICATION_PROPORTION;
  if ( !strcmp( words[0], "DELETION_PROPORTION" ) )         return DELETION_PROPORTION;
  if ( !strcmp( words[0], "TRANSLOCATION_PROPORTION" ) )    return TRANSLOCATION_PROPORTION;
  if ( !strcmp( words[0], "INVERSION_PROPORTION" ) )        return INVERSION_PROPORTION;
  if ( !strcmp( words[0], "ALIGN_FUNCTION" ) )              return ALIGN_FUNCTION;
  if ( !strcmp( words[0], "ALIGN_MAX_SHIFT" ) )             return ALIGN_MAX_SHIFT;
  if ( !strcmp( words[0], "ALIGN_W_ZONE_H_LEN" ) )          return ALIGN_W_ZONE_H_LEN;
  if ( !strcmp( words[0], "ALIGN_MATCH_BONUS" ) )           return ALIGN_MATCH_BONUS;
  if ( !strcmp( words[0], "ALIGN_MISMATCH_COST" ) )         return ALIGN_MISMATCH_COST;
  if ( !strcmp( words[0], "SELECTION_SCHEME" ) )            return SELECTION_SCHEME;
  if ( !strcmp( words[0], "SELECTION_PRESSURE" ) )          return SELECTION_PRESSURE;
  if ( !strcmp( words[0], "SEED" ) )                        return SEED;
  if ( !strcmp( words[0], "TRANSLATION_COST" ) )            return TRANSLATION_COST;
  if ( !strcmp( words[0], "ENV_ADD_POINT" ) )               return ENV_ADD_POINT;
  if ( !strcmp( words[0], "ENV_ADD_GAUSSIAN" ) )            return ENV_ADD_GAUSSIAN;
  if ( !strcmp( words[0], "ENV_SAMPLING" ) )                return ENV_SAMPLING;
  if ( !strcmp( words[0], "ENV_VARIATION" ) )               return ENV_VARIATION;
  if ( !strcmp( words[0], "ENV_SEED" ) )                    return ENV_SEED;
  if ( !strcmp( words[0], "SECRETION_FITNESS_CONTRIB" ) )   return SECRETION_FITNESS_CONTRIB;
  if ( !strcmp( words[0], "SECRETION_DIFUSION_PROP" ) )     return SECRETION_DIFUSION_PROP;
  if ( !strcmp( words[0], "SECRETION_DEGRADATION_PROP" ) )  return SECRETION_DEGRADATION_PROP;
  if ( !strcmp( words[0], "SECRETION_INITIAL" ) )           return SECRETION_INITIAL;
  if ( !strcmp( words[0], "SECRETION_COST" ) )              return SECRETION_COST;
  if ( !strcmp( words[0], "ALLOW_PLASMIDS" ) )              return ALLOW_PLASMIDS;
  if ( !strcmp( words[0], "PLASMID_INITIAL_LENGTH" ) )      return PLASMID_INITIAL_LENGTH;
  if ( !strcmp( words[0], "PLASMID_INITIAL_GENE" ) )        return PLASMID_INITIAL_GENE;
  if ( !strcmp( words[0], "PLASMID_MINIMAL_LENGTH" ) )      return PLASMID_MINIMAL_LENGTH;
  if ( !strcmp( words[0], "PROB_HORIZONTAL_TRANS" ) )       return PROB_HORIZONTAL_TRANS;
  if ( !strcmp( words[0], "COMPUTE_PHEN_CONTRIB_BY_GU" ) )  return COMPUTE_PHEN_CONTRIB_BY_GU;
  if ( !strcmp( words[0], "LOG" ) )                         return LOG;
 
  #ifdef __REGUL
    if ( !strcmp( words[0], "HILL_SHAPE_N" ) )                return HILL_SHAPE_N;
    if ( !strcmp( words[0], "HILL_SHAPE_THETA" ) )            return HILL_SHAPE_THETA;
    if ( !strcmp( words[0], "DEGRADATION_RATE" ) )            return DEGRADATION_RATE;
    if ( !strcmp( words[0], "DEGRADATION_STEP" ) )            return DEGRADATION_STEP;
    if ( !strcmp( words[0], "INDIVIDUAL_EVALUATION_DATES" ) ) return INDIVIDUAL_EVALUATION_DATES;
    if ( !strcmp( words[0], "INDIVIDUAL_LIFE_TIME" ) )        return INDIVIDUAL_LIFE_TIME;
    if ( !strcmp( words[0], "BINDING_ZEROS_PERCENTAGE" ) )    return BINDING_ZEROS_PERCENTAGE;
    if ( !strcmp( words[0], "WITH_HEREDITY" ) )               return WITH_HEREDITY;
    if ( !strcmp( words[0], "PROTEIN_PRESENCE_LIMIT" ) )      return PROTEIN_PRESENCE_LIMIT;
  #endif
  
  return UNDEFINED;
}
