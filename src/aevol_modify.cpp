//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004 LIRIS.
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
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>



// =================================================================
//                            Project Files
// =================================================================
#include <f_line.h>
#include <ae_population.h>
#ifdef __X11
  #include <ae_exp_manager_X11.h>
#else
  #include <ae_exp_manager.h>
#endif

// =================================================================
//                         Function declarations
// =================================================================
void print_help( char* prog_name );
f_line* get_line( FILE* param_file );
void format_line( f_line* formated_line, char* line, bool* line_is_interpretable );

int main( int argc, char* argv[] )
{
  // 1) Initialize command-line option variables with default values
  char* param_file_name = NULL;
  char* pop_file_name   = NULL;
  bool verbose          = false;
  int32_t num_gener = -1;  
  
  // 2) Define allowed options
  const char * options_list = "hf:g:";
  static struct option long_options_list[] = {
    { "help",     no_argument,        NULL, 'h' },
    { "file",     required_argument,  NULL, 'f' }, // Provide file with parameters to change
    { "gener",    required_argument,  NULL, 'g' },
    { 0, 0, 0, 0 }
  };
      
  // 3) Get actual values of the command-line options
  int option;
  while ( ( option = getopt_long(argc, argv, options_list, long_options_list, NULL) ) != -1 ) 
  {
    switch ( option ) 
    {
      case 'h' :
      {
        print_help( argv[0] );
        exit( EXIT_SUCCESS );
      }
      case 'f' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -f or --file : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }
        
        param_file_name = optarg;
        break;
      }
      case 'g' :
      {
        num_gener = atoi( optarg );
        break;
      }
      default :
      {
        // An error message is printed in getopt_long, we just need to exit
        exit( EXIT_FAILURE );
      }
    }
  }
  
  // 4) Check the consistancy of the command-line options
  if ( num_gener == -1 )
  {
    printf( "%s: error: You must provide a generation number.\n", argv[0] );
    exit( EXIT_FAILURE );
  }
  
  // 5) Initialize the experiment manager
  #ifndef __NO_X
    ae_exp_manager* exp_manager = new ae_exp_manager_X11();
  #else
    ae_exp_manager* exp_manager = new ae_exp_manager();
  #endif
  exp_manager->load( num_gener, false, verbose );

  // 6) Retrieve the population, the environment, the selection,...
  ae_population* pop = exp_manager->get_pop();
  ae_environment* env = exp_manager->get_env();
  ae_selection* sel = exp_manager->get_sel();
    
  // 7) Interpret and apply changes
  printf("Interpret and apply changes\n");
  FILE* param_file  = fopen( param_file_name,  "r" );
  if ( param_file == NULL )
  {
    printf( "%s:%d: error: could not open parameter file %s\n", __FILE__, __LINE__, param_file_name );
    exit( EXIT_FAILURE );
  }
  
  f_line* line;
  while ( ( line = get_line(param_file) ) != NULL ) 
  {
    if ( strcmp( line->words[0], "POINT_MUTATION_RATE" ) == 0 )
    {
      pop->set_overall_point_mutation_rate( atof( line->words[1] ) );
      printf("\tChange of overall point mutation rate to %f\n",atof( line->words[1] ));
    }
    else if ( strcmp( line->words[0], "SMALL_INSERTION_RATE" ) == 0 )
    {
      pop->set_overall_small_insertion_rate( atof( line->words[1] ) );
      printf("\tChange of overall small insertion rate to %f\n",atof( line->words[1] ));
    }
    else if ( strcmp( line->words[0], "SMALL_DELETION_RATE" ) == 0 )
    {
      pop->set_overall_small_deletion_rate( atof( line->words[1] ) );
      printf("\tChange of overall small deletion rate to %f\n",atof( line->words[1] ));
    }
    else if ( strcmp( line->words[0], "MAX_INDEL_SIZE" ) == 0 )
    {
      pop->set_overall_max_indel_size( atol( line->words[1] ) );
      printf("\tChange of overall maximum indel size to %f\n",atof( line->words[1] ));
    }
    else if ( strcmp( line->words[0], "DUPLICATION_RATE" ) == 0 )
    {
      pop->set_overall_duplication_rate( atof( line->words[1] ) );
      printf("\tChange of overall duplication rate to %f\n",atof( line->words[1] ));
    }
    else if ( strcmp( line->words[0], "DELETION_RATE" ) == 0 )
    {
      pop->set_overall_deletion_rate( atof( line->words[1] ) );
      printf("\tChange of overall deletion rate to %f\n",atof( line->words[1] ));
    }
    else if ( strcmp( line->words[0], "TRANSLOCATION_RATE" ) == 0 )
    {
      pop->set_overall_translocation_rate( atof( line->words[1] ) );
      printf("\tChange of overall translocation rate to %f\n",atof( line->words[1] ));
    }
    else if ( strcmp( line->words[0], "INVERSION_RATE" ) == 0 )
    {
      pop->set_overall_inversion_rate( atof( line->words[1] ) );
      printf("\tChange of overall inversion to %f\n",atof( line->words[1] ));
    }
    else if ( strcmp( line->words[0], "ENV_ADD_GAUSSIAN" ) == 0 )
    {
      env->add_gaussian( atof(line->words[1]), atof(line->words[2]), atof(line->words[3]));
      env->add_initial_gaussian( atof(line->words[1]), atof(line->words[2]), atof(line->words[3])); //usefull in case of autoregressive mean variation to compute delta_m
      printf("\tAddition of a gaussian with %f, %f, %f \n",atof(line->words[1]), atof(line->words[2]), atof(line->words[3]));
    }
    else if ( strcmp( line->words[0], "ENV_VARIATION" ) == 0 )
    {
      static bool env_var_already_set = false;
      if ( env_var_already_set )
      {
        printf( "%s:%d: ERROR in param file : duplicate entry for %s.\n", __FILE__, __LINE__, line->words[0] );
        exit( EXIT_FAILURE );
      }
      env_var_already_set = true;
      
      if ( strcmp( line->words[1], "none" ) == 0 )
      {
        assert( line->nb_words == 2 );
        env->set_var_method( NO_VAR );
        printf("\tNo more environmental variation\n");
      }
      else if ( strcmp( line->words[1], "autoregressive_mean_variation" ) == 0 )
      {
        assert( line->nb_words == 5 );
        env->set_var_method( AUTOREGRESSIVE_MEAN_VAR );
        env->set_var_sigma( atof( line->words[2] ) );
        env->set_var_tau( atol( line->words[3] ) );
        env->set_var_prng( new ae_jumping_mt(atoi( line->words[4])));
        printf("\tChange of environmental variation to a autoregressive mean variation with sigma=%f, tau=%ld and seed=%d\n", atof( line->words[2] ),atol( line->words[3] ),atoi( line->words[4]));
      }
      else if ( strcmp( line->words[1], "autoregressive_height_variation" ) == 0 )
      {
        assert( line->nb_words == 5 );
        env->set_var_method( AUTOREGRESSIVE_HEIGHT_VAR );
        env->set_var_sigma( atof( line->words[2] ) );
        env->set_var_tau( atol( line->words[3] ) );
        env->set_var_prng( new ae_jumping_mt(atoi( line->words[4])));
        printf("\tChange of environmental variation to a autoregressive height variation with sigma=%f, tau=%ld and seed=%d\n", atof( line->words[2] ),atol( line->words[3] ),atoi( line->words[4]));
      }
      else if ( strcmp( line->words[1], "add_local_gaussians" ) == 0 )
      {
        assert( line->nb_words == 3 );
        env->set_var_method( LOCAL_GAUSSIANS_VAR );
        env->set_var_prng( new ae_jumping_mt(atoi(line->words[2])));
        printf("\tChange of environmental variation to a local gaussians variation with seed=%d\n", atoi( line->words[2]));
      }
      else
      {
        printf( "%s:%d: ERROR in param file : unknown environment variation method.\n", __FILE__, __LINE__ );
        exit( EXIT_FAILURE );
      }
    }
    else if ( strcmp( line->words[0], "SEED" ) == 0 )
    {
      int32_t seed = atoi( line->words[1] ) ;
      
      ae_jumping_mt* prng = new ae_jumping_mt( seed );
      
      // Change prng in ae_selection 
      sel->set_prng( new ae_jumping_mt(*prng) );
      
      if( exp_manager->is_spatially_structured())
      {
        ae_spatial_structure* sp_struct = exp_manager->get_spatial_structure();
        sp_struct->set_prng(new ae_jumping_mt(*prng) );
      }
      
      // Change prng of the population
      pop->set_mut_prng( new ae_jumping_mt(*prng) );
      pop->set_stoch_prng( new ae_jumping_mt(*prng) );
      printf("\tChange of the seed to %d in selection, mutations and individuals' stochasticity \n",atoi( line->words[1] ));
    }
    
  
    delete line;
  }
  fclose( param_file );
  printf("Ok\n");
  
  // 8) Save the changements
  printf("Save the changements into backup\t");
  exp_manager->write_setup_files();
  exp_manager->save();
  printf("Ok\n");
}



f_line* get_line( FILE* param_file )
{
  char line[255];
  f_line* formated_line = new f_line();

  bool found_interpretable_line = false; 

  while ( !feof( param_file ) && !found_interpretable_line )
  {
    if ( !fgets( line, 255, param_file ) )
    {
      delete formated_line;
      return NULL;
    }
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

void format_line( f_line* formated_line, char* line, bool* line_is_interpretable )
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


void print_help( char* prog_name ) 
{
	printf( "******************************************************************************\n" );
	printf( "*                        aevol - Artificial Evolution                        *\n" );
	printf( "******************************************************************************\n" );
	printf( "Usage : change_pop -h\n" );
	printf( "   or : change_pop [-f param_file -p pop_file]\n" );
	printf( "  -h, --help       Display this screen\n" );
  printf( "  -g, --gener      Specify generation number (input)\n" );
  printf( "  -f, --file       File with parameter to change\n");
	printf( "Change a simulation as specified in the parameter file.\n" );
    printf( "(default: param_to_change.in)\n\n" );
}

