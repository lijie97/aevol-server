// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
// 
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// 
//*****************************************************************************




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
#include <ae_align.h>

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
  if ( strcmp( line->words[0], "MIN_TRIANGLE_WIDTH" ) == 0 )
  {
    printf( "ERROR in param file \"%s\" on line %"PRId32": this option is no longer supported.\n", _param_file_name, _cur_line );
    exit( EXIT_FAILURE );
  }
  else if ( strcmp( line->words[0], "MAX_TRIANGLE_WIDTH" ) == 0 )
  {
    _param_values->set_w_max( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "ENV_AXIS_SEGMENTS" ) == 0 )
  {
    printf( "ERROR in param file \"%s\" on line %"PRId32": this option is no longer supported.\nCheck ENV_AXIS_FEATURES new syntax", _param_file_name, _cur_line );
    exit( EXIT_FAILURE );
  }
  else if ( strcmp( line->words[0], "ENV_AXIS_FEATURES" ) == 0 )
  {
    // Set general segmentation data
    _param_values->set_env_axis_nb_segments( line->nb_words / 2 );
    
    // Set segmentation boundaries
    double* env_axis_segment_boundaries = new double [_param_values->get_env_axis_nb_segments() + 1];
    env_axis_segment_boundaries[0] = X_MIN;
    for ( int16_t i = 1 ; i < _param_values->get_env_axis_nb_segments() ; i++ )
    {
      env_axis_segment_boundaries[i] = atof( line->words[2*i] );
    }
    env_axis_segment_boundaries[_param_values->get_env_axis_nb_segments()] = X_MAX;
    _param_values->set_env_axis_segment_boundaries( env_axis_segment_boundaries );
    // do not delete env_axis_segment_boundaries, no in-depth copy, will be deleted with _param_values
    
    // Set segment features
    ae_env_axis_feature* env_axis_features = new ae_env_axis_feature[_param_values->get_env_axis_nb_segments()];
    for ( int16_t i = 0 ; i < _param_values->get_env_axis_nb_segments() ; i++ )
    {
      if ( strcmp( line->words[2*i+1], "NEUTRAL" ) == 0 )
      {
        env_axis_features[i] = NEUTRAL;
      }
      else if ( strcmp( line->words[2*i+1], "METABOLISM" ) == 0 )
      {
        env_axis_features[i] = METABOLISM;
      }
      else if ( strcmp( line->words[2*i+1], "SECRETION" ) == 0 )
      {
        _param_values->set_with_secretion( true );
        env_axis_features[i] = SECRETION;
      }
      else if ( strcmp( line->words[2*i+1], "DONOR" ) == 0 )
      {
        env_axis_features[i] = DONOR;
      }
      else if ( strcmp( line->words[2*i+1], "RECIPIENT" ) == 0 )
      {
        env_axis_features[i] = RECIPIENT;
      }
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown axis feature \"%s\".\n",
                _param_file_name, _cur_line, line->words[2*i+1] );
        exit( EXIT_FAILURE );
      }
    }
    _param_values->set_env_axis_features( env_axis_features );
    // do not delete env_axis_features, no in-depth copy, will be deleted with _param_values 
  }
  else if ( strcmp( line->words[0], "ENV_SEPARATE_SEGMENTS" ) == 0 )
  {
    _param_values->set_env_axis_separate_segments( true );
  }
  else if ( strcmp( line->words[0], "RECORD_TREE" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "TREE_MODE" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "MORE_STATS" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "DUMP_STEP" ) == 0 )
  {
    _param_values->set_dump_step( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "BACKUP_STEP" ) == 0 )
  {
    _param_values->set_backup_step( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "BIG_BACKUP_STEP" ) == 0 )
  {
    _param_values->set_big_backup_step( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "TREE_STEP" ) == 0 )
  {
    _param_values->set_tree_step( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "INITIAL_GENOME_LENGTH" ) == 0 )
  {
    _param_values->set_initial_genome_length( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "MIN_GENOME_LENGTH" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "MAX_GENOME_LENGTH" ) == 0 )
  {
    if ( strncmp( line->words[1], "NONE", 4 ) == 0 )
    {
      _param_values->set_max_genome_length( INT32_MAX );
    }
    else
    {
      _param_values->set_max_genome_length( atol( line->words[1] ) );
    }
  }
  else if ( strcmp( line->words[0], "INIT_POP_SIZE" ) == 0 )
  {
    _param_values->set_init_pop_size( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "POP_STRUCTURE" ) == 0 )
  {
    if ( strcmp( line->words[1], "grid" ) == 0 )
    {
      _param_values->set_spatially_structured( true );
      _param_values->set_grid_width( atol( line->words[2] ) );
      _param_values->set_grid_height( atol( line->words[3] ) );
    }
    else
    {
      printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown population structure.\n", _param_file_name, _cur_line );
      exit( EXIT_FAILURE ); 
    }
  }
  else if ( strcmp( line->words[0], "MIGRATION_NUMBER" ) == 0 )
  {
    _param_values->set_migration_number( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "INIT_METHOD" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "POINT_MUTATION_RATE" ) == 0 )
  {
    _param_values->set_point_mutation_rate( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "SMALL_INSERTION_RATE" ) == 0 )
  {
    _param_values->set_small_insertion_rate( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "SMALL_DELETION_RATE" ) == 0 )
  {
    _param_values->set_small_deletion_rate( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "MAX_INDEL_SIZE" ) == 0 )
  {
    _param_values->set_max_indel_size( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "DUPLICATION_RATE" ) == 0 )
  {
    _param_values->set_duplication_rate( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "DELETION_RATE" ) == 0 )
  {
    _param_values->set_deletion_rate( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "TRANSLOCATION_RATE" ) == 0 )
  {
    _param_values->set_translocation_rate( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "INVERSION_RATE" ) == 0 )
  {
    _param_values->set_inversion_rate( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "NEIGHBOURHOOD_RATE" ) == 0 )
  {
    _param_values->set_neighbourhood_rate( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "DUPLICATION_PROPORTION" ) == 0 )
  {
    _param_values->set_duplication_proportion( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "DELETION_PROPORTION" ) == 0 )
  {
    _param_values->set_deletion_proportion( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "TRANSLOCATION_PROPORTION" ) == 0 )
  {
    _param_values->set_translocation_proportion( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "INVERSION_PROPORTION" ) == 0 )
  {
    _param_values->set_inversion_proportion( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "ALIGN_FUNCTION" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "ALIGN_MAX_SHIFT" ) == 0 )
  {
    _param_values->set_align_max_shift( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "ALIGN_W_ZONE_H_LEN" ) == 0 )
  {
    _param_values->set_align_w_zone_h_len( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "ALIGN_MATCH_BONUS" ) == 0 )
  {
    _param_values->set_align_match_bonus( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "ALIGN_MISMATCH_COST" ) == 0 )
  {
    _param_values->set_align_mismatch_cost( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "STOCHASTICITY" ) == 0 )
  {
    if ( strncmp( line->words[1], "true", 4 ) == 0 )
    {
      _param_values->set_with_stochasticity( true );
    }
  }
  else if ( strcmp( line->words[0], "SELECTION_SCHEME" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "SELECTION_PRESSURE" ) == 0 )
  {
    _param_values->set_selection_pressure( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "SEED" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "MUT_SEED" ) == 0 )
  {
    static bool mut_seed_already_set = false;
    if ( mut_seed_already_set )
    {
      printf( "ERROR in param file \"%s\" on line %"PRId32" : duplicate entry for MUT_SEED.\n",
              _param_file_name, _cur_line );
      exit( EXIT_FAILURE );
    }
    _param_values->set_mut_seed( atol( line->words[1] ) );
    mut_seed_already_set = true;
  }
  else if ( strcmp( line->words[0], "STOCH_SEED" ) == 0 )
  {
    static bool stoch_seed_already_set = false;
    if ( stoch_seed_already_set )
    {
      printf( "ERROR in param file \"%s\" on line %"PRId32" : duplicate entry for STOCH_SEED.\n",
              _param_file_name, _cur_line );
      exit( EXIT_FAILURE );
    }
    _param_values->set_stoch_seed( atol( line->words[1] ) );
    stoch_seed_already_set = true;
  }
  else if ( strcmp( line->words[0], "NB_GENER" ) == 0 )
  {
    _param_values->set_nb_gener( atol( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "WITH_4PTS_TRANS" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "WITH_ALIGNMENTS" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "WITH_TRANSFER" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "REPL_TRANSFER_WITH_CLOSE_POINTS" ) == 0 )
  {
    if ( strncmp( line->words[1], "true", 4 ) == 0 )
    {
      _param_values->set_repl_HT_with_close_points( true );
    }
    else if ( strncmp( line->words[1], "false", 5 ) == 0 )
    {
      _param_values->set_repl_HT_with_close_points( false );
    }
    else
    {
      printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown transfer option (use true/false).\n",
              _param_file_name, _cur_line );
      exit( EXIT_FAILURE ); 
    }
  }
  else if ( strcmp( line->words[0], "SWAP_GUS" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "TRANSFER_INS_RATE" ) == 0 )
  {
    _param_values->set_HT_ins_rate( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "TRANSFER_REPL_RATE" ) == 0 )
  {
    _param_values->set_HT_repl_rate( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "REPL_TRANSFER_DETACH_RATE" ) == 0 )
  {
    _param_values->set_repl_HT_detach_rate( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "TRANSLATION_COST" ) == 0 )
  {
    _param_values->set_translation_cost( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "ENV_ADD_POINT" ) == 0 )
  {
    _param_values->add_env_custom_point( new ae_point_2d(  atof( line->words[1] ), atof( line->words[2] ) ) );
  }
  else if ( strcmp( line->words[0], "ENV_ADD_GAUSSIAN" ) == 0 )
  {
    if ( _param_values->_env_gaussians == NULL ) _param_values->_env_gaussians = new ae_list<ae_gaussian*>();
    
    _param_values->add_env_gaussian( new ae_gaussian( atof( line->words[1] ), atof( line->words[2] ), atof( line->words[3] ) ) );
  }
  else if ( strcmp( line->words[0], "ENV_GAUSSIAN" ) == 0 )
  {
    if ( _param_values->_env_gaussians == NULL ) _param_values->_env_gaussians = new ae_list<ae_gaussian*>();
    
    _param_values->add_env_gaussian( new ae_gaussian( atof( line->words[1] ), atof( line->words[2] ), atof( line->words[3] ) ) );
  }
  else if ( strcmp( line->words[0], "ENV_SAMPLING" ) == 0 )
  {
    _param_values->set_env_sampling( atoi( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "ENV_VARIATION" ) == 0 )
  {
    static bool env_var_already_set = false;
    if ( env_var_already_set )
    {
      printf( "ERROR in param file \"%s\" on line %"PRId32" : duplicate entry for %s.\n",
              _param_file_name, _cur_line, line->words[0] );
      exit( EXIT_FAILURE );
    }
    env_var_already_set = true;
    
    if ( strcmp( line->words[1], "none" ) == 0 )
    {
      assert( line->nb_words == 2 );
      _param_values->set_env_var_method( NO_VAR );
    }
    else if ( strcmp( line->words[1], "autoregressive_mean_variation" ) == 0 )
    {
      assert( line->nb_words == 5 );
      _param_values->set_env_var_method( AUTOREGRESSIVE_MEAN_VAR );
      _param_values->set_env_var_sigma( atof( line->words[2] ) );
      _param_values->set_env_var_tau( atol( line->words[3] ) );
      _param_values->set_env_var_seed( atoi( line->words[4] ) );
    }
    else if ( strcmp( line->words[1], "autoregressive_height_variation" ) == 0 )
    {
      assert( line->nb_words == 5 );
      _param_values->set_env_var_method( AUTOREGRESSIVE_HEIGHT_VAR );
      _param_values->set_env_var_sigma( atof( line->words[2] ) );
      _param_values->set_env_var_tau( atol( line->words[3] ) );
      _param_values->set_env_var_seed( atoi( line->words[4] ) );
    }
    else if ( strcmp( line->words[1], "add_local_gaussians" ) == 0 )
    {
      assert( line->nb_words == 3 );
      _param_values->set_env_var_method( LOCAL_GAUSSIANS_VAR );
      _param_values->set_env_var_seed( atoi( line->words[2] ) );
    }
    else
    {
      printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown environment variation method.\n",
              _param_file_name, _cur_line );
      exit( EXIT_FAILURE );
    }
  }
  else if ( strcmp( line->words[0], "ENV_NOISE" ) == 0 )
  {
    static bool env_noise_already_set = false;
    if ( env_noise_already_set )
    {
      printf( "ERROR in param file \"%s\" on line %"PRId32" : duplicate entry for %s.\n",
              _param_file_name, _cur_line, line->words[0] );
      exit( EXIT_FAILURE );
    }
    env_noise_already_set = true;
    
    if ( strcmp( line->words[1], "none" ) == 0 )
    {
      assert( line->nb_words == 2 );
      _param_values->set_env_noise_method( NO_NOISE );
    }
    else if ( strcmp( line->words[1], "FRACTAL" ) == 0 )
    {
			assert( line->nb_words == 6 );
      _param_values->set_env_noise_method( FRACTAL );
      _param_values->set_env_noise_sampling_log( atoi( line->words[2] ) );
      _param_values->set_env_noise_sigma( atof( line->words[3] ) );
      _param_values->set_env_noise_alpha( atof( line->words[4] ) );
      _param_values->set_env_noise_prob( atof( line->words[5] ) );
      _param_values->set_env_noise_seed( atoi( line->words[6] ) );
    }
    else
    {
      printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown environment noise method.\n",
              _param_file_name, _cur_line );
      exit( EXIT_FAILURE );
    }
  }
  else if ( strcmp( line->words[0], "SECRETION_CONTRIB_TO_FITNESS" ) == 0 )
  {
    _param_values->set_secretion_contrib_to_fitness( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "SECRETION_DIFFUSION_PROP" ) == 0 )
  {
    _param_values->set_secretion_diffusion_prop( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "SECRETION_DEGRADATION_PROP" ) == 0 )
  {
    _param_values->set_secretion_degradation_prop( atof( line->words[1] ) );
    if ( _param_values->get_secretion_degradation_prop() > 1 || _param_values->get_secretion_degradation_prop() < 0 ) 
    {
      printf( "ERROR in param file \"%s\" on line %"PRId32" : degradation must be in (0,1).\n",
              _param_file_name, _cur_line );
      exit( EXIT_FAILURE );
    }
  }
  else if ( strcmp( line->words[0], "SECRETION_INITIAL" ) == 0 )
  {
    _param_values->set_secretion_init( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "SECRETION_COST" ) == 0 )
  {
    _param_values->set_secretion_cost( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "ALLOW_PLASMIDS" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "PLASMID_INITIAL_LENGTH" ) == 0 )
  {
    _param_values->set_plasmid_initial_length( atoi( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "PLASMID_INITIAL_GENE" ) == 0 )
  {
    _param_values->set_plasmid_initial_gene( atoi( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "PLASMID_MINIMAL_LENGTH" ) == 0 )
  {
    _param_values->set_plasmid_minimal_length( atoi( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "PROB_PLASMID_HT" ) == 0 )
  {
    _param_values->set_prob_plasmid_HT( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "TUNE_DONOR_ABILITY" ) == 0 )
  {
    _param_values->set_tune_donor_ability( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "TUNE_RECIPIENT_ABILITY" ) == 0 )
  {
    _param_values->set_tune_recipient_ability( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "DONOR_COST" ) == 0 )
  {
    _param_values->set_donor_cost( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "RECIPIENT_COST" ) == 0 )
  {
    _param_values->set_recipient_cost( atof( line->words[1] ) );
  }
  else if ( strcmp( line->words[0], "COMPUTE_PHEN_CONTRIB_BY_GU" ) == 0 )
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
  }
  else if ( strcmp( line->words[0], "LOG" ) == 0 )
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
      /*else if ( strcmp( line->words[i], "LOADS" ) == 0 )
      {
        tmp_to_be_logged |= LOG_LOADS;
      }   */
      else
      {
        printf( "ERROR in param file \"%s\" on line %"PRId32" : unknown log option %s.\n",
                _param_file_name, _cur_line, line->words[1] );
        exit( EXIT_FAILURE ); 
      }
    }
    
    _param_values->set_logs( tmp_to_be_logged );
  }
  #ifdef __REGUL
    else if ( strcmp( line->words[0], "HILL_SHAPE_N" ) == 0 )
    {
      _param_values->set_hill_shape_n( atof( line->words[1] ) );
    }
    else if ( strcmp( line->words[0], "HILL_SHAPE_THETA" ) == 0 )
    {
      _param_values->set_hill_shape_theta( atof( line->words[1] ) );
    }
    else if ( strcmp( line->words[0], "DEGRADATION_RATE" ) == 0 )
    {
      _param_values->set_degradation_rate( atof( line->words[1] ) );
    }
    else if ( strcmp( line->words[0], "DEGRADATION_STEP" ) == 0 )
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
    }
    else if ( strcmp( line->words[0], "INDIVIDUAL_EVALUATION_DATES" ) == 0 )
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
    }
    else if ( strcmp( line->words[0], "BINDING_ZEROS_PERCENTAGE" ) == 0 )
    {
      _param_values->set_binding_zeros_percentage( atof( line->words[1] ) );
    }
    else if ( strcmp( line->words[0], "WITH_HEREDITY" ) == 0 )
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
    }
    else if ( strcmp( line->words[0], "PROTEIN_PRESENCE_LIMIT" ) == 0 )
    {
      _param_values->set_protein_presence_limit( atof( line->words[1] ) );
    }
  #endif
  else
  {
    printf( "ERROR in param file \"%s\" on line %"PRId32" : undefined key word \"%s\"\n", _param_file_name, _cur_line, line->words[0] );
    exit( EXIT_FAILURE );
  }
}

void param_loader::read_file( void )
{
  // The rewind is only necessary when using multiple param files
  rewind( _param_file );

  _cur_line = 0;
  f_line* line;
  
  _param_values = new params();

  while ( ( line = get_line() ) != NULL ) // TODO : write line = new f_line( _param_file ) => f_line::f_line( char* )
  {
    interpret_line( line, _cur_line );
    delete line;
  }
}


void param_loader::load( ae_exp_manager* exp_m, bool verbose )
{
	// Initialize _prng
  _prng = new ae_jumping_mt( _param_values->_seed );
  
  // Initialize _mut_prng and _stoch_prng :
  // if mut_seed (respectively stoch_seed) not given in param.in, choose it at random
  if ( _param_values->_mut_seed == 0 ) {
		_param_values->set_mut_seed( _prng->random( 1000000 ) );
	}
	if ( _param_values->_stoch_seed == 0 ) {
		_param_values->set_stoch_seed( _prng->random( 1000000 ) );
	}
	_mut_prng = new ae_jumping_mt( _param_values->_mut_seed );
	_stoch_prng = new ae_jumping_mt( _param_values->_stoch_seed );
  
  // Create aliases (syntaxic sugars)
  ae_exp_setup*       exp_s     = exp_m->get_exp_s();
  ae_population*      pop       = exp_m->get_pop();
  ae_environment*     env       = exp_m->get_env();
  ae_selection*       sel       = exp_m->get_sel();
  ae_output_manager*  output_m  = exp_m->get_output_m();
  
  // If the population is spatially structured,
  // check that the population fits in the spatial structure
  if ( _param_values->_spatially_structured )
  {
    if ( _param_values->get_init_pop_size() > _param_values->_grid_width * _param_values->_grid_height )
    {
      printf( "ERROR: the spatial structure is too small to fit the entire population\n" );
      exit( EXIT_FAILURE );
    }
    else if ( _param_values->get_init_pop_size() < _param_values->_grid_width * _param_values->_grid_height )
    {
      printf( "ERROR: the current implementation does not support grids that are not filled with individuals.\n" );
      exit( EXIT_FAILURE );
    }
  }
  
  // 1) ------------------------------------- Initialize the experimental setup
  sel->set_prng( new ae_jumping_mt(*_prng) );

  // ---------------------------------------------------------------- Selection
  sel->set_selection_scheme( _param_values->_selection_scheme );
  sel->set_selection_pressure( _param_values->_selection_pressure );
  
  // ----------------------------------------------------------------- Transfer
  exp_s->set_with_HT( _param_values->_with_HT );
  exp_s->set_repl_HT_with_close_points( _param_values->get_repl_HT_with_close_points() );
  exp_s->set_HT_ins_rate( _param_values->_HT_ins_rate );
  exp_s->set_HT_repl_rate( _param_values->_HT_repl_rate );
  exp_s->set_repl_HT_detach_rate( _param_values->_repl_HT_detach_rate );
  
  // ----------------------------------------------------------------- Plasmids
  exp_s->set_with_plasmids( _param_values->_allow_plasmids );
  exp_s->set_prob_plasmid_HT( _param_values->_prob_plasmid_HT );
  exp_s->set_tune_donor_ability( _param_values->_tune_donor_ability );
  exp_s->set_tune_recipient_ability( _param_values->_tune_recipient_ability );
  exp_s->set_donor_cost( _param_values->_donor_cost );
  exp_s->set_recipient_cost( _param_values->_recipient_cost );
  exp_s->set_swap_GUs( _param_values->_swap_GUs );
  output_m->set_compute_phen_contrib_by_GU( _param_values->_compute_phen_contrib_by_GU );
  
  // -------------------------------------------------------- Spatial structure
  if ( _param_values->_spatially_structured )
  {
    ae_jumping_mt* sp_struct_prng = new ae_jumping_mt(*_prng);
    exp_m->set_spatial_structure( _param_values->_grid_width,
                                  _param_values->_grid_height,
                                  sp_struct_prng );
    ae_spatial_structure* sp_struct = exp_m->get_spatial_structure();
    sp_struct->set_secretion_degradation_prop( _param_values->_secretion_degradation_prop );
    sp_struct->set_secretion_diffusion_prop( _param_values->_secretion_diffusion_prop );
    sp_struct->set_migration_number( _param_values->_migration_number );
  }
  
  // ---------------------------------------------------------------- Secretion
  exp_s->set_with_secretion( _param_values->_with_secretion );
  exp_s->set_secretion_contrib_to_fitness( _param_values->_secretion_contrib_to_fitness );
  exp_s->set_secretion_cost( _param_values->_secretion_cost );
  
  
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
  env->set_segmentation( _param_values->get_env_axis_nb_segments(),
                         _param_values->get_env_axis_segment_boundaries(),
                         _param_values->get_env_axis_features(),
                         _param_values->get_env_axis_separate_segments() );
  
  // Set environmental variation
  if ( _param_values->get_env_var_method() != NO_VAR )
  {
    env->set_var_method( _param_values->get_env_var_method() );
    env->set_var_prng( new ae_jumping_mt( _param_values->get_env_var_seed() ) );
    env->set_var_sigma_tau( _param_values->get_env_var_sigma(), _param_values->get_env_var_tau() );
  }
  
  // Set environmental noise
  if ( _param_values->get_env_noise_method() != NO_NOISE )
  {
    env->set_noise_method( _param_values->get_env_noise_method() );
    env->set_noise_sampling_log( _param_values->get_env_noise_sampling_log() );
    env->set_noise_prng( new ae_jumping_mt( _param_values->get_env_noise_seed() ) );
    env->set_noise_alpha( _param_values->get_env_noise_alpha() );
    env->set_noise_sigma( _param_values->get_env_noise_sigma() );
    env->set_noise_prob( _param_values->get_env_noise_prob()  );
  }
  
  // Build the environment
  env->build();
  
  if ( verbose )
  {
    printf( "Entire geometric area of the environment : %f\n", env->get_geometric_area() );
  }
  
  
  // 3) --------------------------------------------- Create the new population
  pop->set_mut_prng( new ae_jumping_mt(*_mut_prng) );
  pop->set_stoch_prng( new ae_jumping_mt(*_stoch_prng) );
  
  // Generate a model ae_mut_param object
  ae_params_mut* param_mut = new ae_params_mut();
  param_mut->set_point_mutation_rate( _param_values->get_point_mutation_rate() );
  param_mut->set_small_insertion_rate( _param_values->get_small_insertion_rate() );
  param_mut->set_small_deletion_rate( _param_values->get_small_deletion_rate() );
  param_mut->set_max_indel_size( _param_values->get_max_indel_size() );
  param_mut->set_with_4pts_trans( _param_values->get_with_4pts_trans() );
  param_mut->set_with_alignments( _param_values->get_with_alignments() );
  param_mut->set_with_HT( _param_values->get_with_HT() );
  param_mut->set_repl_HT_with_close_points( _param_values->get_repl_HT_with_close_points() );
  param_mut->set_HT_ins_rate( _param_values->get_HT_ins_rate() );
  param_mut->set_HT_repl_rate( _param_values->get_HT_repl_rate() );
  param_mut->set_repl_HT_detach_rate( _param_values->get_repl_HT_detach_rate() );
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
      
      indiv->set_with_stochasticity( _param_values->get_with_stochasticity() );
      
      // Add it to the list
      pop->add_indiv( indiv );
    
      // Make the clones and add them to the list of individuals
      ae_individual* clone = NULL;
      for ( int32_t i = 1 ; i < _param_values->get_init_pop_size() ; i++ )
      {
        // Create a clone, setting its id
        clone = create_clone( indiv, id_new_indiv++ );
        
        #ifdef DISTRIBUTED_PRNG
          #error Not implemented yet !
          indiv->do_prng_jump();
        #endif
        
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

        #ifdef DISTRIBUTED_PRNG
          #error Not implemented yet !
          indiv->do_prng_jump();
        #endif
        
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
  
  // If the population is spatially structured, set each individual's position
  if ( exp_m->is_spatially_structured() )
  {
    int16_t x, y;
    int16_t x_max = exp_m->get_grid_width();
    int16_t y_max = exp_m->get_grid_height();
    ae_grid_cell* grid_cell = NULL;
    
    ae_list_node<ae_individual*>* indiv_node = pop->get_indivs()->get_first();
    ae_individual*  indiv = NULL;
    
    while ( indiv_node != NULL )
    {
      indiv = indiv_node->get_obj();
      
      do
      {
        x = _prng->random( x_max );
        y = _prng->random( y_max );
        grid_cell = exp_m->get_grid_cell( x, y );
      } while ( grid_cell->get_individual() != NULL );
      
      grid_cell->set_individual( indiv );
      
      indiv_node = indiv_node->get_next();
    }
  }
  
  
  
  // 4) ------------------------------------------ Set the recording parameters
  output_m->set_backup_step( _param_values->get_backup_step() );
  output_m->set_big_backup_step( _param_values->get_big_backup_step() );
  
  if ( _param_values->get_record_tree() )
  {
    output_m->init_tree( exp_m, _param_values->get_tree_mode(), _param_values->get_tree_step() );
  }
  
  if ( _param_values->get_make_dumps() )
  {
    output_m->set_dump_step( _param_values->get_dump_step() );
  }
  output_m->set_logs( _param_values->get_logs() );
  
  delete param_mut;
  delete _prng; // Each class that needed it has now its own copy
  delete _mut_prng;
  delete _stoch_prng;
}



// =================================================================
//                           Protected Methods
// =================================================================
/*!
  \brief Format a line by parsing it and the words inside
  
  \param formated_line the resulted formated line
  \param line original line in char*
  \param line_is_interpretable boolean with about the possible intrepretation of the line
*/
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

/*!
  \brief Get a line in a file and format it
  
  \return line (pointer)
  
  \see format_line(f_line* formated_line, char* line, bool* line_is_interpretable )
*/
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

/*!
  \brief Create an individual with random sequences
  
  \param exp_m global exp_manager
  \param param_mut mutation parameter of the newly created individual
  \param id index of newly created individual in the population
  \return clone of dolly
*/
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
  #ifdef DISTRIBUTED_PRNG
    #error Not implemented yet !
  #endif
  ae_individual* indiv = new ae_individual( exp_m,
                                            _prng,
                                            _prng,
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
  
  if (_param_values->_allow_plasmids) // We create a plasmid
  {
    char * plasmid_genome = new char [_param_values->_initial_genome_length + 1]; // As ae_dna constructor do not allocate memory but directly use the provided string, we allocate the memory here.
    if ( _param_values->get_plasmid_initial_gene() == 1 ) // Then the plasmid is generated independently from the chromosome
    {
      for ( int32_t i = 0 ; i < _param_values->_initial_genome_length ; i++ )
      {
        plasmid_genome[i] = '0' + _prng->random( NB_BASE );
      }
      plasmid_genome[_param_values->_initial_genome_length] = 0;
    }
    else // The plasmid has the same genome than the chromosome
    {
      strncpy(plasmid_genome,random_genome,_param_values->_initial_genome_length+1);
    }
    indiv->add_GU( plasmid_genome, _param_values->_initial_genome_length );
    plasmid_genome = NULL; // should not be deleted since it is now the plasmid dna
  }
  
  random_genome = NULL; // should not be deleted since it is now the chromosomal dna
  
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

/*!
  \brief Create an individual with random sequences. The individual have to have at least one good functional gene
  
  \param exp_m global exp_manager
  \param param_mut mutation parameter of the newly created individual
  \param id index of newly created individual in the population
  \return clone of dolly
*/
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
  
  //~ printf( "metabolic error of the generated individual : %f (%"PRId32" gene(s))\n",
          //~ indiv->get_dist_to_target_by_feature(METABOLISM), indiv->get_protein_list()->get_nb_elts() );
  
  return indiv;
}

/*!
  \brief Create of clone of an ae_individual 
  
  \param dolly original individual that would be cloned
  \param id index of the clone in the population
  \return clone of dolly
*/
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
  //~ printf( "metabolic error of the clonal individual : %f (%"PRId32" gene(s))\n",
          //~ indiv->get_dist_to_target_by_feature(METABOLISM), indiv->get_protein_list()->get_nb_elts());
  return indiv;
}

