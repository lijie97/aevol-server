//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2009  LIRIS.
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
#include <errno.h>
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <sys/stat.h>



// =================================================================
//                            Project Files
// =================================================================
#ifdef __X11
  #include <ae_exp_manager_X11.h>
#else
  #include <ae_exp_manager.h>
#endif

#include <ae_macros.h>





#ifndef __NO_X
void catch_usr1( int sig_num );
#endif

#ifdef __IN2P3
void catch_usr1( int sig_num );
void catch_xcpu( int sig_num );
#endif

void print_help( char* prog_path );




int main( int argc, char* argv[] )
{
  // Set handlers for signals to be caught
  #ifndef __NO_X
    signal( SIGUSR1, catch_usr1 );
  #endif

  #ifdef __IN2P3
    signal( SIGUSR1, catch_usr1 );
    signal( SIGXCPU, catch_xcpu );
  #endif  
  
  
  // Print warning for debug mode
  #ifdef DEBUG
    printf( "aevol is being run in DEBUG mode\n" );
  #endif
  
 
  
  // =========================================================================
  //                           Get command-line options
  // =========================================================================
  // 1) Initialize command-line option variables with default values
  // 2) Define allowed options
  // 3) Get actual values of the command-line options
  // 4) Check the consistancy of the command-line options
  // 5) Set file names according to options
  // =========================================================================
  
  // -------------------------------------------------------------------------
  // 1) Initialize command-line option variables with default values
  // -------------------------------------------------------------------------
  bool use_text_files = false;
  bool pause_on_startup = false;
  bool verbose = false;
  
  int32_t num_gener = 0;
  int32_t nb_gener  = 1000;
  
  #ifndef __NO_X
    bool show_display_on_startup = true;
  #endif
  
  
  // -------------------------------------------------------------------------
  // 2) Define allowed options
  // -------------------------------------------------------------------------
  const char * options_list = "hn:r:tvwx";
  static struct option long_options_list[] = {
    { "help",     no_argument,        NULL, 'h' },
    { "nbgener",  required_argument,  NULL, 'n' }, // Number of generations to be run
    { "resume",   required_argument,  NULL, 'r' }, // Resume from generation X
    { "text",     no_argument,        NULL, 't' }, // Use text files instead of gzipped binary files
    { "verbose",  no_argument,        NULL, 'v' }, // Be verbose
    { "wait",     no_argument,        NULL, 'w' }, // Pause after loading
    { "noX",      no_argument,        NULL, 'x' }, // Don't display X outputs on start
    { 0, 0, 0, 0 }
  };
  
  
  // -------------------------------------------------------------------------
  // 3) Get actual values of the command-line options
  // -------------------------------------------------------------------------
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
      case 'n' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -n or --nbgener : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }
        
        nb_gener = atoi( optarg );
        
        break;
      }
      case 'r' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -r or --resume : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }
        
        num_gener = atol( optarg );
		  
        break;      
      }
      case 'v' :
      {
        verbose = true;
        
        break;
      }
      case 'w' :
      {
        pause_on_startup = true;
        
        break;
      }
      case 'x' :
      {
        #ifdef __NO_X
          printf( "%s: error: Program was compiled with __NO_X option, no visualisation available.\n", argv[0] );
          exit( EXIT_FAILURE );
        #else
          show_display_on_startup = false;
        #endif
        
        break;
      }
    }
  }
  
  
  // -------------------------------------------------------------------------
  // 4) Check the consistancy of the command-line options
  // -------------------------------------------------------------------------
  // NOTA: nothing to be checked for (for now)
  //
  //~ if ( env_file_name == NULL || pop_file_name == NULL )
  //~ {
    //~ printf( "%s: error: You must either resume a simulation of provide both an environment backup and a population backup.\n", argv[0] );
    //~ // NB : spatial structure is tested later, when we know if the population is structured or not
    //~ exit( EXIT_FAILURE );
  //~ }
  
  
  // -------------------------------------------------------------------------
  // 5) Set file names according to options
  // -------------------------------------------------------------------------  
  char* env_file_name       = new char[255];
  char* pop_file_name       = new char[255];
  char* exp_setup_file_name = new char[255];
  char* out_prof_file_name  = new char[255];
  char* sp_struct_file_name = new char[255];
  
  sprintf( env_file_name,       ENV_FNAME_FORMAT,       num_gener );
  sprintf( pop_file_name,       POP_FNAME_FORMAT,       num_gener );
  sprintf( exp_setup_file_name, EXP_S_FNAME_FORMAT,     num_gener );
  sprintf( out_prof_file_name,  OUT_P_FNAME_FORMAT,     num_gener );
  sprintf( sp_struct_file_name, SP_STRUCT_FNAME_FORMAT, num_gener );

  // Check existence of optional files in file system.
  // Missing files will cause the corresponding file_name variable to be nullified
  struct stat stat_buf;
  if ( stat( sp_struct_file_name, &stat_buf ) == -1 )
  {
    if ( errno == ENOENT )
    {
      delete [] sp_struct_file_name;
      sp_struct_file_name = NULL;
    }
    else
    {
      printf( "%s:%d: error: unknown error.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    }
  }
  
  
  // =================================================================
  //                          Load the simulation
  // =================================================================
  #ifndef __NO_X
    ae_exp_manager* exp_manager = new ae_exp_manager_X11();
  #else
    ae_exp_manager* exp_manager = new ae_exp_manager();
  #endif
  
  exp_manager->load( num_gener, exp_setup_file_name, out_prof_file_name, env_file_name, pop_file_name, sp_struct_file_name, verbose );
  exp_manager->set_nb_gener( nb_gener );
  
  // Make a numbered copy of each static input file (dynamic files are saved elsewhere)
  // TODO (?)
  
  delete [] exp_setup_file_name;
  delete [] out_prof_file_name;
  delete [] env_file_name;
  delete [] pop_file_name;
  delete [] sp_struct_file_name;

#ifndef __NO_X
  if ( show_display_on_startup )
  {
    ((ae_exp_manager_X11*) exp_manager)->toggle_display_on_off();
  }
#endif
  
  
  // =================================================================
  //                         Run the simulation
  // =================================================================
  exp_manager->run_evolution();
  
  
  delete exp_manager;
  //~ return EXIT_SUCCESS;
}


#ifndef __NO_X
void catch_usr1( int sig_num )
{
  signal( SIGUSR1, catch_usr1 );
  
  printf( "display on/off\n" );
  //~ exp_manager->toggle_display_on_off();
}
#endif

#ifdef __IN2P3
void catch_usr1( int sig_num )
{
  //system( "./make_big_backup.py" );
  exit( EXIT_FAILURE );
}

void catch_xcpu( int sig_num )
{
  //system( "./make_big_backup.py" );
  exit( EXIT_FAILURE );
}
#endif

void print_help( char* prog_path ) 
{
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ( prog_name = strrchr( prog_path, '/' ) ) prog_name++;
  else prog_name = prog_path;
  
	printf( "******************************************************************************\n" );
	printf( "*                        aevol - Artificial Evolution                        *\n" );
	printf( "******************************************************************************\n" );
	printf( "Usage : %s -h\n", prog_name );
  printf( "   or : %s [options]\n", prog_name );
  printf( "  -h, --help               Display this screen\n" );
  printf( "  -n, --nbgener NB_GENER   Number of generations to be run\n" );
  printf( "  -r, --resume GENERATION  Resume simulation at given generation\n" );
  printf( "  -t, --text               Use text files instead of binary files when possible\n" );
  printf( "  -v, --verbose            Be verbose\n" );
  printf( "  -w, --wait               Pause after loading\n" );
  printf( "  -x, --noX                Don't display X outputs upon start\n" );
  printf( "                           send SIGUSR1 to switch X output on/off\n" );
  printf( "Run an aevol simulation\n\n" );
}



    //~ { "setup",    required_argument,  NULL, 's' }, // Provide exp setup file
    //~ { "out_prof", required_argument,  NULL, 'o' }, // Provide output profile file
    //~ { "envir",    required_argument,  NULL, 'e' }, // Provide environment file
    //~ { "pop",      required_argument,  NULL, 'p' }, // Provide population file
    //~ { "grid",     required_argument,  NULL, 'g' }, // Provide spatial structure file

	//~ printf( "******************************************************************************\n" );
	//~ printf( "*                        aevol - Artificial Evolution                        *\n" );
	//~ printf( "******************************************************************************\n" );
	//~ printf( "Usage : aevol -h\n" );
  //~ printf( "   or : aevol [options]\n" );
  //~ printf( "  -h, --help               Display this screen\n" );
  //~ printf( "  -r, --resume GENERATION  Resume simulation at given generation\n" );
  //~ printf( "  -s, --setup FILE         Provide experimental setup file\n" );
  //~ printf( "  -o, --out_prof FILE      Provide output profile file\n" );
  //~ printf( "  -e, --envir FILE         Provide environment file\n" );
  //~ printf( "  -p, --pop FILE           Provide population file\n" );
  //~ printf( "  -g, --grid FILE          Provide spatial structure file\n" );
  //~ printf( "  -n, --nbgener NB_GENER   Number of generations to be run\n" );
  //~ printf( "  -w, --wait               Pause after loading\n" );
  //~ printf( "  -x, --noX                Don't display X outputs upon start\n" );
  //~ printf( "                           send SIGUSR1 to switch X output on/off\n" );
  //~ printf( "Run an aevol simulation\n\n" );
