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
#include <param_loader.h>

#include <ae_exp_manager.h>
#include <ae_exp_setup.h>
#include <ae_output_manager.h>
#include <ae_population.h>
#include <ae_individual.h>

#include <ae_jumping_mt.h>
#include <ae_gaussian.h>
#include <ae_env_segment.h>
#include <ae_point_2d.h>

//~ #ifdef __X11
  //~ #include <ae_individual_X11.h>
//~ #endif

#ifdef __REGUL
  #include <ae_array_short.h>
#endif


// =================================================================
//                          Class declarations
// =================================================================
class ae_environment;




//##############################################################################
//                                                                             #
//                             Class param_loader                              #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================


param_loader::param_loader( const char* file_name )
{
  _param_file_name = strdup( file_name );
  _param_file  = fopen( _param_file_name,  "r" );
  
  if ( _param_file == NULL )
  {
    printf( "ERROR : couldn't open file %s\n", file_name );
    exit( EXIT_FAILURE );
  }
  
  assert( _param_file );
  _cur_line = 0;
  
  read_file();
}

// =================================================================
//                             Destructors
// =================================================================
param_loader::~param_loader( void )
{
  free( _param_file_name );
  fclose( _param_file );
  
  delete _param_values;
}

// =================================================================
//                            Public Methods
// =================================================================

void param_loader::interpret_line( f_line* line, int32_t _cur_line )
{
  // Interpret line
  switch ( line->get_keywd() )
  {
    case UNDEFINED :
    {
      printf( "ERROR in param file \"%s\" on line %"PRId32" : undefined key word \"%s\"\n",
              _param_file_name, _cur_line, line->words[0] );
      exit( EXIT_FAILURE );
    }
    case MIN_TRIANGLE_WIDTH :
    {
      printf( "ERROR in param file \"%s\" on line %"PRId32": this option is no longer supported.\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE ); 
      break;
    }
    case MAX_TRIANGLE_WIDTH :
    {
      _param_values->set_w_max( atof( line->words[1] ) );
      break;
    }
    case ENV_AXIS_SEGMENTS :
    {
      // Check consistency
      if ( _param_values->get_env_axis_segment_boundaries() != NULL )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : ENV_AXIS_SEGMENT keyword found twice.\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE );
      }
      if ( _param_values->get_env_axis_features() != NULL && _param_values->get_env_axis_nb_segments() != line->nb_words )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : number of segment boundaries and number of features don't match.\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE );
      }
      
      // OK, proceed
      _param_values->set_env_axis_is_segmented( true );
      _param_values->set_env_axis_nb_segments( line->nb_words ); // in fact line->nb_words - 1 + 1
      
      double* env_axis_segment_boundaries = new double [_param_values->get_env_axis_nb_segments() + 1];
      env_axis_segment_boundaries[0] = X_MIN;
      for ( int16_t i = 1 ; i < _param_values->get_env_axis_nb_segments() ; i++ )
      {
        env_axis_segment_boundaries[i] = atof( line->words[i] );
      }
      env_axis_segment_boundaries[_param_values->get_env_axis_nb_segments()] = X_MAX;
      _param_values->set_env_axis_segment_boundaries( env_axis_segment_boundaries );
      // do not delete env_axis_segment_boundaries, no in-depth copy !
      break;
    }
    case ENV_AXIS_FEATURES :
    {
      // Check consistency
      if ( _param_values->get_env_axis_features() != NULL )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : ENV_AXIS_FEATURES keyword found twice.\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE );
      }
      if ( _param_values->get_env_axis_segment_boundaries() != NULL && _param_values->get_env_axis_nb_segments() != line->nb_words - 1 )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : number of segment boundaries and number of features don't match.\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE );
      }
      
      // OK, proceed
      _param_values->set_env_axis_is_segmented( true );
      _param_values->set_env_axis_nb_segments( line->nb_words - 1 );
      ae_env_axis_feature* env_axis_features = new ae_env_axis_feature[_param_values->get_env_axis_nb_segments()];
      
      for ( int16_t i = 0 ; i < _param_values->get_env_axis_nb_segments() ; i++ )
      {
        if ( strcmp( line->words[i+1], "NEUTRAL" ) == 0 )
        {
          env_axis_features[i] = NEUTRAL;
        }
        else if ( strcmp( line->words[i+1], "METABOLISM" ) == 0 )
        {
          env_axis_features[i] = METABOLISM;
        }
        else if ( strcmp( line->words[i+1], "SECRETION" ) == 0 )
        {
          _param_values->set_use_secretion( true );
          env_axis_features[i] = SECRETION;
        }
        else if ( strcmp( line->words[i+1], "TRANSFER" ) == 0 )
        {
          env_axis_features[i] = TRANSFER;
        }
        else
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown axis feature \"%s\".\n",
                  _param_file_name, _cur_line, line->words[1] );
          exit( EXIT_FAILURE );
        }
      }
      _param_values->set_env_axis_features( env_axis_features );
      // do not delete env_axis_features, no in-depth copy !
      
      break;
    }
    case ENV_SEPARATE_SEGMENTS :
    {
      _param_values->set_env_axis_separate_segments( true );
      break;
    }
    case RECORD_TREE :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        _param_values->set_record_tree( true );
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        _param_values->set_record_tree( false );
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown tree recording option (use true/false).\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE ); 
      }        
      break;
    }
    case TREE_MODE :
    {
      if ( strcmp( line->words[1], "light" ) == 0 )
      {
        _param_values->set_tree_mode( LIGHT );
      }
      else if ( strcmp( line->words[1], "normal" ) == 0 )
      {
        _param_values->set_tree_mode( NORMAL );
      }       
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown tree mode option (use normal/light).\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE );
      }        
      break;
    }
    case MORE_STATS :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        _param_values->set_more_stats( true );
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        _param_values->set_more_stats( false );
      }       
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown more stats option (use true/false).\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE ); 
      }        
      break;
    }
    case DUMP_STEP :
    {
      _param_values->set_dump_step( atol( line->words[1] ) );
      break;
    }
    case BACKUP_STEP :
    {
      _param_values->set_backup_step( atol( line->words[1] ) );
      break;
    }
    case BIG_BACKUP_STEP :
    {
      _param_values->set_big_backup_step( atol( line->words[1] ) );
      break;
    }
    case TREE_STEP :
    {
      _param_values->set_tree_step( atol( line->words[1] ) );
      break;
    }
    case INITIAL_GENOME_LENGTH :
    {
      _param_values->set_initial_genome_length( atol( line->words[1] ) );
      break;
    }
    case MIN_GENOME_LENGTH :
    {
      if ( strncmp( line->words[1], "NONE", 4 ) == 0 )
      {
        _param_values->set_min_genome_length( 1 ); // Must not be 0
      }
      else
      {
        _param_values->set_min_genome_length( atol( line->words[1] ) );
        if ( _param_values->get_min_genome_length() == 0 )
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : MIN_GENOME_LENGTH must be > 0.\n",
                  _param_file_name, _cur_line );
          exit( EXIT_FAILURE ); 
        }
      }
      break;
    }
    case MAX_GENOME_LENGTH :
    {
      if ( strncmp( line->words[1], "NONE", 4 ) == 0 )
      {
        _param_values->set_max_genome_length( INT32_MAX );
      }
      else
      {
        _param_values->set_max_genome_length( atol( line->words[1] ) );
      }
      break;
    }
    case INIT_POP_SIZE :
    {
      _param_values->set_init_pop_size( atol( line->words[1] ) );
      break;      
    }
    case POP_STRUCTURE :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        _param_values->set_spatially_structured( true );
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        _param_values->set_spatially_structured( false );
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown population structure (use true/false).\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE ); 
      }
      break;
    }
    case GRID_WIDTH :
    {
      _param_values->set_grid_width( atol( line->words[1] ) );
      break;
    }
    case GRID_HEIGHT :
    {
      _param_values->set_grid_height( atol( line->words[1] ) );
      break;      
    }
    case MIGRATION_NUMBER :
    {
      _param_values->set_migration_number( atol( line->words[1] ) );
      break;
    }
    case INIT_METHOD :
    {
      int8_t init_method = _param_values->get_init_method();
      
      for ( int8_t i = 1 ; i < line->nb_words ; i++ )
      {
        if ( strcmp( line->words[i], "ONE_GOOD_GENE" ) == 0 )
        {
          init_method |= ONE_GOOD_GENE;
        }
        else if ( strcmp( line->words[i], "CLONE" ) == 0 )
        {
          init_method |= CLONE;
        }
        else if ( strcmp( line->words[i], "WITH_INS_SEQ" ) == 0 )
        {
          init_method |= WITH_INS_SEQ;
        }   
        else
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown initialization method %s.\n",
                  _param_file_name, _cur_line, line->words[1] );
          exit( EXIT_FAILURE ); 
        }
      }
      _param_values->set_init_method( init_method );
      break;
    }
    case POINT_MUTATION_RATE :
    {
      _param_values->set_point_mutation_rate( atof( line->words[1] ) );
      break;
    }
    case SMALL_INSERTION_RATE :
    {
      _param_values->set_small_insertion_rate( atof( line->words[1] ) );
      break;
    }
    case SMALL_DELETION_RATE :
    {
      _param_values->set_small_deletion_rate( atof( line->words[1] ) );
      break;
    }
    case MAX_INDEL_SIZE :
    {
      _param_values->set_max_indel_size( atol( line->words[1] ) );
      break;
    }
    case DUPLICATION_RATE :
    {
      _param_values->set_duplication_rate( atof( line->words[1] ) );
      break;
    }
    case DELETION_RATE :
    {
      _param_values->set_deletion_rate( atof( line->words[1] ) );
      break;
    }
    case TRANSLOCATION_RATE :
    {
      _param_values->set_translocation_rate( atof( line->words[1] ) );
      break;
    }
    case INVERSION_RATE :
    {
      _param_values->set_inversion_rate( atof( line->words[1] ) );
      break;
    }
    case NEIGHBOURHOOD_RATE :
    {
      _param_values->set_neighbourhood_rate( atof( line->words[1] ) );
      break;
    }
    case DUPLICATION_PROPORTION :
    {
      _param_values->set_duplication_proportion( atof( line->words[1] ) );
      break;
    }
    case DELETION_PROPORTION :
    {
      _param_values->set_deletion_proportion( atof( line->words[1] ) );
      break;
    }
    case TRANSLOCATION_PROPORTION :
    {
      _param_values->set_translocation_proportion( atof( line->words[1] ) );
      break;
    }
    case INVERSION_PROPORTION :
    {
      _param_values->set_inversion_proportion( atof( line->words[1] ) );
      break;
    }
    case ALIGN_FUNCTION :
    {
      if ( line->nb_words != 2 && line->nb_words != 4 )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : incorrect number of parameters for keyword \"%s\".\n",
                _param_file_name, _cur_line, line->words[0] );
        exit( EXIT_FAILURE );
      }
      
      if ( strcmp( line->words[1], "LINEAR" ) == 0 )
      {
        _param_values->set_align_fun_shape( LINEAR );
        
        if ( line->nb_words == 4 )
        {
          _param_values->set_align_lin_min( atol( line->words[2] ) );
          _param_values->set_align_lin_max( atol( line->words[3] ) );
        }
      }
      else if ( strcmp( line->words[1], "SIGMOID" ) == 0 )
      {
        _param_values->set_align_fun_shape( SIGMOID );
        
        if ( line->nb_words == 4 )
        {
          _param_values->set_align_sigm_lambda( atol( line->words[2] ) );
          _param_values->set_align_sigm_mean( atol( line->words[3] ) );
        }
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown align function shape \"%s\".\n",
                _param_file_name, _cur_line, line->words[1] );
        exit( EXIT_FAILURE );
      }
      break;
    }
    case ALIGN_MAX_SHIFT :
    {
      _param_values->set_align_max_shift( atol( line->words[1] ) );
      break;
    }
    case ALIGN_W_ZONE_H_LEN :
    {
      _param_values->set_align_w_zone_h_len( atol( line->words[1] ) );
      break;
    }
    case ALIGN_MATCH_BONUS :
    {
      _param_values->set_align_match_bonus( atol( line->words[1] ) );
      break;
    }
    case ALIGN_MISMATCH_COST :
    {
      _param_values->set_align_mismatch_cost( atol( line->words[1] ) );
      break;
    }
    case SELECTION_SCHEME :
    {
      if ( strncmp( line->words[1], "lin", 3 ) == 0 )
      {
        _param_values->set_selection_scheme( RANK_LINEAR );
      }
      else if ( strncmp( line->words[1], "exp", 3 ) == 0 )
      {
        _param_values->set_selection_scheme( RANK_EXPONENTIAL );
      }
      else if ( strncmp( line->words[1], "fitness", 7 ) == 0 )
      {
        _param_values->set_selection_scheme( FITNESS_PROPORTIONATE );
      }
      else if ( strcmp( line->words[1], "fittest" ) == 0 )
      {
        _param_values->set_selection_scheme( FITTEST );
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown selection scheme \"%s\".\n",
                _param_file_name, _cur_line, line->words[1] );
        exit( EXIT_FAILURE );
      }
      break;
    }
    case SELECTION_PRESSURE :
    {
      _param_values->set_selection_pressure( atof( line->words[1] ) );
      break;
    }
    case SEED :
    {
      static bool seed_already_set = false;
      if ( seed_already_set )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : duplicate entry for SEED.\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE );
      }
      _param_values->set_seed( atol( line->words[1] ) );
      seed_already_set = true;
      break;
    }
    case NB_GENER :
    {
      _param_values->set_nb_gener( atol( line->words[1] ) );
      break;
    }
    case WITH_4PTS_TRANS :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        _param_values->set_with_4pts_trans( true );
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        printf( "ERROR: 3 points translocation hasn't been implemented yet\n" );
        exit( EXIT_FAILURE );
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown 4pts_trans option (use true/false).\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE ); 
      }
      break;
    }
    case WITH_ALIGNMENTS :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        _param_values->set_with_alignments( true );
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        _param_values->set_with_alignments( false );
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown alignement option (use true/false).\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE ); 
      }
      break;
    }
    case WITH_TRANSFER :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        _param_values->set_with_HT( true );
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        _param_values->set_with_HT( false );
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown transfer option (use true/false).\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE ); 
      }
      break;
    }
    case SWAP_GUS:
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        _param_values->set_swap_GUs( true );
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        _param_values->set_swap_GUs( false );
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown swap option (use true/false).\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE ); 
      }
      break;
    }
    case TRANSFER_INS_RATE :
    {
      _param_values->set_HT_ins_rate( atof( line->words[1] ) );
      break;
    }
    case TRANSFER_REPL_RATE :
    {
      _param_values->set_HT_repl_rate( atof( line->words[1] ) );
      break;
    }
    case TRANSLATION_COST :
    {
      _param_values->set_translation_cost( atof( line->words[1] ) );
      break;
    }
    case ENV_ADD_POINT :
    {
      _param_values->add_env_custom_point( new ae_point_2d(  atof( line->words[1] ), atof( line->words[2] ) ) );
      break;
    }
    case ENV_ADD_GAUSSIAN :
    {
      if ( _param_values->_env_gaussians == NULL ) _param_values->_env_gaussians = new ae_list();
      _param_values->add_env_gaussian( new ae_gaussian( atof( line->words[1] ), atof( line->words[2] ), atof( line->words[3] ) ) );
      break;
    }
    case ENV_SAMPLING :
    {
      _param_values->set_env_sampling( atoi( line->words[1] ) );
      break;
    }
    case ENV_VARIATION :
    {
      if ( strcmp( line->words[1], "none" ) == 0 )
      {
        _param_values->set_env_var_method( NONE );
      }
      else if ( strcmp( line->words[1], "autoregressive_mean_variation" ) == 0 )
      {
        _param_values->set_env_var_method( AUTOREGRESSIVE_MEAN_VAR );
        _param_values->set_env_var_sigma( atof( line->words[2] ) );
        _param_values->set_env_var_tau( atol( line->words[3] ) );
      }
      else if ( strcmp( line->words[1], "add_local_gaussians" ) == 0 )
      {
        _param_values->set_env_var_method( LOCAL_GAUSSIANS_VAR );
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown environment variation method.\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE );
      }
      break;
    }
    case ENV_SEED :
    {
      static bool env_seed_already_set = false;
      if ( env_seed_already_set )
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : duplicate entry for ENV_SEED.\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE );
      }
      _param_values->set_env_seed( atoi( line->words[1] ) );
      env_seed_already_set = true;
      break;
    }
    case SECRETION_CONTRIB_TO_FITNESS :
    {
      _param_values->set_secretion_contrib_to_fitness( atof( line->words[1] ) );
      break; 
    }
    case SECRETION_DIFFUSION_PROP :
    {
      _param_values->set_secretion_diffusion_prop( atof( line->words[1] ) );
      break; 
    }
    case SECRETION_DEGRADATION_PROP :
    {
      _param_values->set_secretion_degradation_prop( atof( line->words[1] ) );
      if ( _param_values->get_secretion_degradation_prop() > 1 || _param_values->get_secretion_degradation_prop() < 0 ) 
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : degradation must be in (0,1).\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE );
      }
      break; 
    }
    case SECRETION_INITIAL :
    {
      _param_values->set_secretion_init( atof( line->words[1] ) );
      break; 
    }
    case SECRETION_COST :
    {
      _param_values->set_secretion_cost( atof( line->words[1] ) );
      break; 
    }
    case ALLOW_PLASMIDS :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        _param_values->set_allow_plasmids( true );
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        _param_values->set_allow_plasmids( false );
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown allow_plasmids option (use true/false).\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE ); 
      }
      break; 
    }
    case PLASMID_INITIAL_LENGTH :
    {
      _param_values->set_plasmid_initial_length( atoi( line->words[1] ) );
      break; 
    }
    case PLASMID_INITIAL_GENE :
    {
      _param_values->set_plasmid_initial_gene( atoi( line->words[1] ) );
      break; 
    }
    case PLASMID_MINIMAL_LENGTH :
    {
      _param_values->set_plasmid_minimal_length( atoi( line->words[1] ) );
      break; 
    }
    case PROB_PLASMID_HT  :
    {
      _param_values->set_prob_plasmid_HT( atof( line->words[1] ) );
      break; 
    }
    case NB_PLASMID_HT :
    {
      _param_values->set_with_plasmid_HT( true );
      _param_values->set_nb_plasmid_HT( atoi( line->words[1] ) );
      break;
    }
	
    case COMPUTE_PHEN_CONTRIB_BY_GU :
    {
      if ( strncmp( line->words[1], "true", 4 ) == 0 )
      {
        _param_values->set_compute_phen_contrib_by_GU( true );
      }
      else if ( strncmp( line->words[1], "false", 5 ) == 0 )
      {
        _param_values->set_compute_phen_contrib_by_GU( false );
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown compute_phen_contrib_by_GU option (use true/false).\n",
                _param_file_name, _cur_line );
        exit( EXIT_FAILURE ); 
      }
      break; 
    }
    case LOG :
    {
      int8_t tmp_to_be_logged = 0;
      
      for ( int8_t i = 1 ; i < line->nb_words ; i++ )
      {
        if ( strcmp( line->words[i], "TRANSFER" ) == 0 )
        {
          tmp_to_be_logged |= LOG_TRANSFER;
        }
        else if ( strcmp( line->words[i], "REAR" ) == 0 )
        {
          tmp_to_be_logged |= LOG_REAR;
        }
        else if ( strcmp( line->words[i], "BARRIER" ) == 0 )
        {
          tmp_to_be_logged |= LOG_BARRIER;
        }
        else if ( strcmp( line->words[i], "LOADS" ) == 0 )
        {
          tmp_to_be_logged |= LOG_LOADS;
        }   
        else
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown log option %s.\n",
                  _param_file_name, _cur_line, line->words[1] );
          exit( EXIT_FAILURE ); 
        }
      }
      
      _param_values->set_logs( tmp_to_be_logged );
      break;
    }
    
    #ifdef __REGUL
      case HILL_SHAPE_N :
      {
        _param_values->set_hill_shape_n( atof( line->words[1] ) );
        break;
      }
      case HILL_SHAPE_THETA :
      {
        _param_values->set_hill_shape_theta( atof( line->words[1] ) );
        break;
      }
      case DEGRADATION_RATE :
      {
        _param_values->set_degradation_rate( atof( line->words[1] ) );
        break;
      }
      case DEGRADATION_STEP :
      {
        _param_values->set_degradation_step( atof( line->words[1] ) );
        // Check that 1/degradation_step is an integer
        if( 1/_param_values->get_degradation_step() != ((int) 1/_param_values->get_degradation_step()) )
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : DEGRADATION STEP\n",
                  _param_file_name, _cur_line );
          printf( "This step has to divide 1.\n");
          exit( EXIT_FAILURE );
        }
        break;
      }
      case INDIVIDUAL_EVALUATION_DATES :
      {
        _param_values->set_individual_evaluation_nbr( line->nb_words - 1 );
        if( _param_values->get_individual_evaluation_nbr() == 0 )
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : no evaluation dates provided\n",
                  _param_file_name, _cur_line );
          exit( EXIT_FAILURE );
        }
        ae_array_short* individual_evaluation_dates  = new ae_array_short( _param_values->get_individual_evaluation_nbr() );
        for( int16_t i = 0 ; i < _param_values->set_individual_evaluation_nbr() ; i++ )
        {
          individual_evaluation_dates->set_value( i, atoi( line->words[1 + i] ) );
        }
        individual_evaluation_dates->sort();
        _param_values->set_individual_evaluation_dates( individual_evaluation_dates );
        break;
      }
      case BINDING_ZEROS_PERCENTAGE :
      {
        _param_values->set_binding_zeros_percentage( atof( line->words[1] ) );
        break;
      }
      case WITH_HEREDITY :
      {
        if ( strncmp( line->words[1], "true", 4 ) == 0 )
        {
          _param_values->set_with_heredity( true );
        }
        else if ( strncmp( line->words[1], "false", 5 ) == 0 )
        {
          _param_values->set_with_heredity( false );
        }
        else
        {
          printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown with_heredity option (use true/false).\n",
                  _param_file_name, _cur_line );
          exit( EXIT_FAILURE ); 
        }
        
        break;
      }
      case PROTEIN_PRESENCE_LIMIT :
      {
        _param_values->set_protein_presence_limit( atof( line->words[1] ) );
        break;
      }
    #endif
    
    
    default :
      printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown error\n",
              _param_file_name, _cur_line );
      exit( EXIT_FAILURE );
      break;
  }
}

void param_loader::read_file( void )
{
  // The rewind is only necessary when using multiple param files
  rewind( _param_file );

  _cur_line = 0;
  f_line* line;
  
  delete _param_values;
  _param_values = new params();

  while ( ( line = get_line() ) != NULL ) // TODO : write line = new f_line( _param_file ) => f_line::f_line( char* )
  {
    interpret_line( line, _cur_line );
    delete line;
  }
}

void param_loader::load( ae_exp_manager* exp_m, bool verbose )
{
  _prng = new ae_jumping_mt( _param_values->_seed );
  
  // Create aliases (syntaxic sugars)
  ae_population*      pop       = exp_m->get_pop();
  ae_exp_setup*       exp_s     = exp_m->get_exp_s();
  ae_output_manager*  output_m  = exp_m->get_output_m();
  ae_environment*     env       = exp_s->get_env();
  ae_selection*       sel       = exp_s->get_sel();
  
  // 1) ------------------------------------- Initialize the experimental setup
  sel->set_prng( new ae_jumping_mt(*_prng) );

  // ---------------------------------------------------------------- Selection
  sel->set_selection_scheme( _param_values->_selection_scheme );
  sel->set_selection_pressure( _param_values->_selection_pressure );
  
  // ----------------------------------------------------------------- Transfer
  sel->set_with_HT( _param_values->_with_HT );
  sel->set_HT_ins_rate( _param_values->_HT_ins_rate );
  sel->set_HT_repl_rate( _param_values->_HT_repl_rate );
  sel->set_with_plasmid_HT( _param_values->_with_plasmid_HT );
  sel->set_nb_plasmid_HT( _param_values->_nb_plasmid_HT );
  sel->set_prob_plasmid_HT( _param_values->_prob_plasmid_HT );
  sel->set_swap_GUs( _param_values->_swap_GUs );
  
  // -------------------------------------------------------- Spatial structure
  sel->set_spatially_structured( _param_values->_spatially_structured );
  if ( _param_values->_spatially_structured )
  {
    sel->set_grid_size( _param_values->_grid_width, _param_values->_grid_height );
    sel->set_migration_number( _param_values->_migration_number );
  }
  
  // ---------------------------------------------------------------- Secretion
  sel->set_use_secretion( _param_values->_use_secretion );
  if ( _param_values->_use_secretion )
  {
    sel->set_secretion_contrib_to_fitness( _param_values->_secretion_contrib_to_fitness );
    sel->set_secretion_cost( _param_values->_secretion_cost );
  }
  
  
  
  // 2) ------------------------------------------------ Create the environment
  // Move the gaussian list and the list of custom points from the parameters
  // to the environment
  env->set_gaussians( _param_values->get_env_gaussians() );
  _param_values->set_env_gaussians( NULL );
  env->set_custom_points( _param_values->get_env_custom_points() );
  _param_values->set_env_custom_points( NULL );
  
  // Copy the sampling
  env->set_sampling( _param_values->get_env_sampling() );
  
  // Set the environment segmentation
  if ( _param_values->get_env_axis_is_segmented() )
  {
    env->set_segmentation( _param_values->get_env_axis_nb_segments(),
                           _param_values->get_env_axis_segment_boundaries(),
                           _param_values->get_env_axis_features(),
                           _param_values->get_env_axis_separate_segments() );
  }
  
  // Set environmental variation
  if ( _param_values->get_env_var_method() != NONE )
  {
    env->set_variation_method( _param_values->get_env_var_method() );
    env->set_prng_var( new ae_jumping_mt( _param_values->get_env_seed() ) );
    env->set_var_sigma_tau( _param_values->get_env_var_sigma(), _param_values->get_env_var_tau() );
  }
  
  // Build the environment
  env->build();
  
  if ( verbose )
  {
    printf( "Entire geometric area of the environment : %f\n", env->get_geometric_area() );
  }
  
  
  // 3) --------------------------------------------- Create the new population
  // Generate a template ae_mut_param object
  ae_params_mut* param_mut = new ae_params_mut();
  param_mut->set_point_mutation_rate( _param_values->get_point_mutation_rate() );
  param_mut->set_small_insertion_rate( _param_values->get_small_insertion_rate() );
  param_mut->set_small_deletion_rate( _param_values->get_small_deletion_rate() );
  param_mut->set_max_indel_size( _param_values->get_max_indel_size() );
  param_mut->set_with_4pts_trans( _param_values->get_with_4pts_trans() );
  param_mut->set_with_alignments( _param_values->get_with_alignments() );
  param_mut->set_with_HT( _param_values->get_with_HT() );
  param_mut->set_swap_GUs( _param_values->get_swap_GUs() );
  param_mut->set_HT_ins_rate( _param_values->get_HT_ins_rate() );
  param_mut->set_HT_repl_rate( _param_values->get_HT_repl_rate() );
  param_mut->set_duplication_rate( _param_values->get_duplication_rate() );
  param_mut->set_deletion_rate( _param_values->get_deletion_rate() );
  param_mut->set_translocation_rate( _param_values->get_translocation_rate() );
  param_mut->set_inversion_rate( _param_values->get_inversion_rate() );
  param_mut->set_neighbourhood_rate( _param_values->get_neighbourhood_rate() );
  param_mut->set_duplication_proportion( _param_values->get_duplication_proportion() );
  param_mut->set_deletion_proportion( _param_values->get_deletion_proportion() );
  param_mut->set_translocation_proportion( _param_values->get_translocation_proportion() );
  param_mut->set_inversion_proportion( _param_values->get_inversion_proportion() );
  
  ae_individual* indiv        = NULL;
  int32_t        id_new_indiv = 0;
  
  if ( _param_values->get_init_method() & ONE_GOOD_GENE )
  {
    if ( _param_values->get_init_method() & CLONE )
    {
      // Create an individual with a "good" gene (in fact, make an indiv whose
      // fitness is better than that corresponding to a flat phenotype)
      // and set its id
      indiv = create_random_individual_with_good_gene( exp_m, param_mut, id_new_indiv++ );
      
      // Add it to the list
      pop->add_indiv( indiv );      
    
      // Make the clones and add them to the list of individuals
      ae_individual* clone = NULL;
      for ( int32_t i = 1 ; i < _param_values->get_init_pop_size() ; i++ )
      {
        // Create a clone, setting its id
        clone = create_clone( indiv, id_new_indiv++ );
        
        // Add it to the list
        pop->add_indiv( clone );
      }
    }
    else // if ( ! CLONE )
    {
      for ( int32_t i = 0 ; i < _param_values->get_init_pop_size() ; i++ )
      {
        // Create an individual and set its id
        indiv = create_random_individual_with_good_gene( exp_m, param_mut, id_new_indiv++ );
        
        // Add it to the list
        pop->add_indiv( indiv );
      }
      
      pop->sort_individuals();
    }
  }
  else // if ( ! ONE_GOOD_GENE )
  {    
    if ( _param_values->get_init_method() & CLONE )
    {
      // Create a random individual and set its id
      indiv = create_random_individual( exp_m, param_mut, id_new_indiv++ );
      
      // Add it to the list
      pop->add_indiv( indiv );
      
      // Make the clones and add them to the list of individuals
      ae_individual* clone = NULL;
      for ( int32_t i = 1 ; i < _param_values->get_init_pop_size() ; i++ )
      {
        // Create a clone, setting its id
        clone = create_clone( indiv, id_new_indiv++ );
        
        // Add it to the list
        pop->add_indiv( clone );
      }
    }
    else // if ( ! CLONE )
    {
      for ( int32_t i = 0 ; i < _param_values->get_init_pop_size() ; i++ )
      {
        // Create a random individual and set its id
        indiv = create_random_individual( exp_m, param_mut, id_new_indiv++ );
        
        // Add it to the list
        pop->add_indiv( indiv );
      }
      
      pop->sort_individuals();
    }
  }
  
  // NB: Since we are only creating (initial) backup files, we need not
  // establish a link between the individuals and the spatial structure (if any)
  
  
  
  // 4) ------------------------------------------ Set the recording parameters
  output_m->set_backup_step( _param_values->get_backup_step() );
  output_m->set_big_backup_step( _param_values->get_big_backup_step() );
  
  if ( _param_values->get_record_tree() )
  {
    output_m->init_tree( _param_values->get_tree_mode(), _param_values->get_tree_step() );
  }
  
  if ( _param_values->get_make_dumps() )
  {
    output_m->set_dump_step( _param_values->get_dump_step() );
  }
}



// =================================================================
//                           Protected Methods
// =================================================================

void param_loader::format_line( f_line* formated_line, char* line, bool* line_is_interpretable )
{
  int16_t i = 0;
  int16_t j;

  // Parse line
  while ( line[i] != '\n' && line[i] != '\0' && line[i] != '\r' )
  {
    j = 0;

    // Flush white spaces and tabs
    while ( line[i] == ' ' || line[i] == 0x09 ) i++; // 0x09 is the ASCII code for TAB

    // Check comments
    if ( line[i] == '#' ) break;

    // If we got this far, there is content in the line
    *line_is_interpretable = true;

    // Parse word
    while ( line[i] != ' '  && line[i] != '\n' && line[i] != '\0' && line[i] != '\r' )
    {
      formated_line->words[formated_line->nb_words][j++] = line[i++];
    }

    // Add '\0' at end of word if it's not empty (line ending with space or tab)
    if ( j != 0 )
    {
      formated_line->words[formated_line->nb_words++][j] = '\0';
    }
  }
}

f_line* param_loader::get_line( void )
{
  char line[255];
  f_line* formated_line = new f_line();

  bool found_interpretable_line = false; // Found line that is neither a comment nor empty

  while ( !feof( _param_file ) && !found_interpretable_line )
  {
    if ( !fgets( line, 255, _param_file ) )
    {
      delete formated_line;
      return NULL;
    }
    _cur_line++;
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
  if ( !strcmp( words[0], "DUMP_STEP" ) )                   return DUMP_STEP;
  if ( !strcmp( words[0], "INIT_POP_SIZE" ) )               return INIT_POP_SIZE;
  if ( !strcmp( words[0], "NB_GENER" ) )                    return NB_GENER;
  if ( !strcmp( words[0], "POP_STRUCTURE" ) )               return POP_STRUCTURE;
  if ( !strcmp( words[0], "GRID_WIDTH" ) )                  return GRID_WIDTH;
  if ( !strcmp( words[0], "GRID_HEIGHT" ) )                 return GRID_HEIGHT;
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
  if ( !strcmp( words[0], "SWAP_GUS" ) )                    return SWAP_GUS;
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
  if ( !strcmp( words[0], "SECRETION_CONTRIB_TO_FITNESS") ) return SECRETION_CONTRIB_TO_FITNESS;
  if ( !strcmp( words[0], "SECRETION_DIFFUSION_PROP" ) )    return SECRETION_DIFFUSION_PROP;
  if ( !strcmp( words[0], "SECRETION_DEGRADATION_PROP" ) )  return SECRETION_DEGRADATION_PROP;
  if ( !strcmp( words[0], "SECRETION_INITIAL" ) )           return SECRETION_INITIAL;
  if ( !strcmp( words[0], "SECRETION_COST" ) )              return SECRETION_COST;
  if ( !strcmp( words[0], "ALLOW_PLASMIDS" ) )              return ALLOW_PLASMIDS;
  if ( !strcmp( words[0], "PLASMID_INITIAL_LENGTH" ) )      return PLASMID_INITIAL_LENGTH;
  if ( !strcmp( words[0], "PLASMID_INITIAL_GENE" ) )        return PLASMID_INITIAL_GENE;
  if ( !strcmp( words[0], "PLASMID_MINIMAL_LENGTH" ) )      return PLASMID_MINIMAL_LENGTH;
  if ( !strcmp( words[0], "PROB_PLASMID_HT" ) )             return PROB_PLASMID_HT;
  if ( !strcmp( words[0], "HT_PLASMID_HT" ) )               return NB_PLASMID_HT;
  if ( !strcmp( words[0], "COMPUTE_PHEN_CONTRIB_BY_GU" ) )  return COMPUTE_PHEN_CONTRIB_BY_GU;
  if ( !strcmp( words[0], "LOG" ) )                         return LOG;
  if ( !strcmp( words[0], "GENERATION_OVERLOAD" ))          return GENERATION_OVERLOAD;
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

ae_individual* param_loader::create_random_individual( ae_exp_manager* exp_m, ae_params_mut* param_mut, int32_t id ) const
{
  // Generate a random genome
  char * random_genome = new char [_param_values->_initial_genome_length + 1];
  for ( int32_t i = 0 ; i < _param_values->_initial_genome_length ; i++ )
  {
    random_genome[i] = '0' + _prng->random( NB_BASE );
  }
  random_genome[_param_values->_initial_genome_length] = 0;
  
  
  // ------------------------------------------------------- Global constraints
  // Create an individual with this genome and set its id
  ae_individual* indiv = new ae_individual( exp_m,
                                            new ae_jumping_mt(*_prng),
                                            param_mut,
                                            _param_values->_w_max,
                                            _param_values->_min_genome_length,
                                            _param_values->_max_genome_length,
                                            _param_values->_allow_plasmids,
                                            _param_values->_plasmid_minimal_length,
                                            id, 0 );
  // <Graphical debug>
  //~ #ifdef __X11
    //~ indiv = new ae_individual_X11( exp_m, new ae_jumping_mt(*_prng), param_mut, _param_values->_w_max, id, 0 );
  //~ #else
    //~ indiv = new ae_individual( exp_m, new ae_jumping_mt(*_prng), param_mut, _param_values->_w_max, id, 0 );
  //~ #endif
  // </Graphical debug>
  indiv->add_GU( random_genome, _param_values->_initial_genome_length );
  
  
  // Insert a few IS in the sequence
  /*if ( ae_common::init_params->get_init_method() & WITH_INS_SEQ )
  {
    // Create a random sequence
    int32_t seq_len = 50;
    char* ins_seq = new char[seq_len+1];
    int16_t nb_insert = 50;
    int16_t nb_invert = 50;
    
    for ( int32_t i = 0 ; i < seq_len ; i++ )
    {
      ins_seq[i] = '0' + ae_common::sim->prng->random( NB_BASE );
    }
    ins_seq[seq_len] = '\0';
    
    
    // Insert the sequence at random positions
    ae_mutation* mut1 = NULL;
    for ( int16_t i = 0 ; i < nb_insert ; i++ )
    {
      mut1 = indiv->get_genetic_unit(0)->get_dna()->do_insertion( ins_seq, seq_len );
      delete mut1;
    }
    
    
    // Invert the sequence and insert it at random positions
    char* inverted_seq = new char[seq_len+1];
    for ( int32_t i = 0 ; i < seq_len ; i++ )
    {
      inverted_seq[i] = (ins_seq[seq_len-1-i] == '1') ? '0' : '1';
    }
    inverted_seq[seq_len] = '\0';
    
    for ( int16_t i = 0 ; i < nb_invert ; i++ )
    {
      mut1 = indiv->get_genetic_unit(0)->get_dna()->do_insertion( inverted_seq, seq_len );
      delete mut1;
    }
    
    delete [] ins_seq;
    delete [] inverted_seq;
  }*/
  
  // Evaluate the newly created individual
  indiv->evaluate( exp_m->get_env() );
  
  return indiv;
}

ae_individual* param_loader::create_random_individual_with_good_gene( ae_exp_manager* exp_m, ae_params_mut* param_mut, int32_t id ) const
{
  // Create a random individual and evaluate it
  ae_individual* indiv = create_random_individual( exp_m, param_mut, id );
  
  // While the created individual is not better than the flat individual (indiv whith no metabolic gene),
  // we delete it and replace it by another random individual
  double env_metabolic_area;
  
  env_metabolic_area = exp_m->get_env()->get_area_by_feature( METABOLISM );

  // If there are plasmids, make sure there is at least one metabolic gene on each genetic units
  if ( _param_values->get_allow_plasmids() ) 
  {
    if ( _param_values->get_plasmid_initial_gene() == 1 )
    {
      while ( indiv->get_genetic_unit(0)->get_dist_to_target_by_feature( METABOLISM ) >= env_metabolic_area  ||
              indiv->get_genetic_unit(1)->get_dist_to_target_by_feature( METABOLISM ) >= env_metabolic_area  )
      {
        delete indiv;
        indiv = create_random_individual( exp_m, param_mut, id );
      }
    }
    else // if  ( ae_common::init_params->get_plasmid_initial_gene() == 2 )
    {
      // here things work the same as before, but in the constructor of the individual, 
      // a single genetic unit is created and then copied from the chromosome to the plasmid
      while ( indiv->get_dist_to_target_by_feature( METABOLISM ) >= env_metabolic_area )
      {
        delete indiv;
        indiv = create_random_individual( exp_m, param_mut, id );
      }
    }
  }
  else
  {
    while ( indiv->get_dist_to_target_by_feature( METABOLISM ) >= env_metabolic_area )
    {
      delete indiv;
      indiv = create_random_individual( exp_m, param_mut, id );
    }
  }
  
  // Compute the "good" individual's statistics
  indiv->compute_statistical_data();
  
  printf( "metabolic error of the generated individual : %f (%"PRId32" gene(s))\n",
          indiv->get_dist_to_target_by_feature(METABOLISM), indiv->get_protein_list()->get_nb_elts() );
  
  return indiv;
}

ae_individual* param_loader::create_clone( ae_individual* dolly, int32_t id ) const
{
  ae_individual* indiv;
  
  indiv = new ae_individual( *dolly );
  
  //~ #ifdef __X11
    //~ indiv = new ae_individual_X11( *(dynamic_cast<ae_individual_X11*>(dolly)) );
  //~ #else
    //~ indiv = new ae_individual( *dolly );
  //~ #endif
  
  indiv->set_id( id );
  
  return indiv;
}

