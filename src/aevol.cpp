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
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>



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

void print_help( char* prog_name );




int main( int argc, char* argv[] )
{
  // Catch SIGUSR1 in handler catch_usr1. SIGUSR1 will toggle
  #ifndef __NO_X
    signal( SIGUSR1, catch_usr1 );
  #endif

  #ifdef __IN2P3
    signal( SIGUSR1, catch_usr1 );
    signal( SIGXCPU, catch_xcpu );
  #endif  
  
  #ifdef DEBUG
    printf( "aevol is being run in DEBUG mode\n" );
  #endif
  
 
  // =================================================================
  //                      Get command-line options
  // =================================================================
  //
  // 1) Initialize command-line option variables with default values
  char* exp_setup_file_name = new char[63];
  char* out_prof_file_name  = new char[63];
  strcpy( exp_setup_file_name,  "exp_setup.ae" );
  strcpy( out_prof_file_name,   "output_profile.ae" );
  char* env_file_name       = NULL;
  char* pop_file_name       = NULL;
  char* sp_struct_file_name = NULL;
  
  bool  pause_on_startup = false;
  
  int32_t num_gener = 0;
  int32_t nb_gener  = 0;
  
  
  #ifndef __NO_X
    bool show_display_on_startup = true;
  #endif
  
  // 2) Define allowed options
  const char * options_list = "hr:s:o:e:p:g:n:wx";
  static struct option long_options_list[] = {
    { "help",     no_argument,        NULL, 'h' },
    { "resume",   required_argument,  NULL, 'r' }, // Resume from generation X
    { "setup",    required_argument,  NULL, 's' }, // Provide exp setup file
    { "out_prof", required_argument,  NULL, 'o' }, // Provide output profile file
    { "envir",    required_argument,  NULL, 'e' }, // Provide environment file
    { "pop",      required_argument,  NULL, 'p' }, // Provide population file
    { "grid",     required_argument,  NULL, 'g' }, // Provide spatial structure file
    { "nbgener",  required_argument,  NULL, 'n' }, // Number of generations to be run
    { "wait",     no_argument,        NULL, 'w' }, // Pause after loading
    { "noX",      no_argument,        NULL, 'x' }, // Don't display X outputs on start
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
      case 'r' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -r or --resume : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }
        
        num_gener = atol( optarg );
        
        env_file_name       = new char[255];
        pop_file_name       = new char[255];
        sp_struct_file_name = new char[255];
        
        sprintf( env_file_name,       ENV_FNAME_FORMAT,       num_gener );
        sprintf( pop_file_name,       POP_FNAME_FORMAT,       num_gener );
        sprintf( sp_struct_file_name, SP_STRUCT_FNAME_FORMAT, num_gener );
        
        break;      
      }
      case 's' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -s or --setup : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }
        
        delete [] exp_setup_file_name;
        exp_setup_file_name = new char[strlen(optarg)+1];
        memcpy( exp_setup_file_name, optarg, strlen(optarg)+1 );
        
        break;
      }
      case 'o' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -o or --out_prof : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }
        
        delete [] out_prof_file_name;
        out_prof_file_name = new char[strlen(optarg)+1];
        memcpy( out_prof_file_name, optarg, strlen(optarg)+1 );
        
        break;
      }
      case 'e' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -e or --envir : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }
        
        env_file_name = new char[strlen(optarg)+1];
        memcpy( env_file_name, optarg, strlen(optarg)+1 );
        
        break;
      }
      case 'p' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -p or --pop : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }
        
        pop_file_name = new char[strlen(optarg)+1];
        memcpy( pop_file_name, optarg, strlen(optarg)+1 );
        
        break;
      }
      case 'g' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -g or --grid : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }
        
        sp_struct_file_name = new char[strlen(optarg)+1];
        memcpy( sp_struct_file_name, optarg, strlen(optarg)+1 );
        
        break;
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
  
  // 3) Check the consistancy of the command-line options
  if ( env_file_name == NULL || pop_file_name == NULL )
  {
    printf( "%s: error: You must either resume a simulation of provide both an environment backup and a population backup.\n", argv[0] );
    // NB : spatial structure is tested later, when we know if the population is structured or not
    exit( EXIT_FAILURE );
  }
  
  
  // =================================================================
  //                          Load the simulation
  // =================================================================
  #ifndef __NO_X
    ae_exp_manager* exp_manager = new ae_exp_manager_X11();
  #else
    ae_exp_manager* exp_manager = new ae_exp_manager();
  #endif
  
  exp_manager->load_experiment( exp_setup_file_name, out_prof_file_name, env_file_name, pop_file_name, sp_struct_file_name, true );
  
  if ( num_gener > 0 )
  {
    exp_manager->set_first_gener( num_gener );
  }
  if ( nb_gener > 0 )
  {
    exp_manager->set_nb_gener( nb_gener );
  }
  
  delete [] exp_setup_file_name;
  delete [] out_prof_file_name;
  delete [] env_file_name;
  delete [] pop_file_name;
  delete [] sp_struct_file_name;

    
#ifndef __NO_X
  if ( show_display_on_startup )
  {
    ((ae_exp_manager_X11*)exp_manager)->toggle_display_on_off();
  }
#endif
  
  // =================================================================
  //                         Run the simulation
  // =================================================================
  exp_manager->run_evolution();
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

void print_help( char* prog_name ) 
{
  printf( "\n************* aevol - Artificial Evolution ************* \n\n" );
  printf( "This program is Free Software. No Warranty.\n\
Copyright (C) 2009  LIRIS.\n\n\
Usage : %s -h\n\
   or : %s [options]\n\n\
\t-h or --help : Display this screen\n\
\t-x or --noX  : Disable X display\n\n\
Options (i : integer, d : double, s : string) :\n\n\
\t-f or --file  s   : Begin with the backup file s\n\
\t-n or --nbgen i   : Set the number of generation to do to i \n\
\t-s or --seed  i   : Set the random number seed to i \n\
\t-p or --param \"s\" : Override parameters as specified in s.\n\
\t                    s must follow the syntax of a parameter file \n\n\
If an option is not set, the programm use the default value for this parameter.\n\n\
<press enter to display examples>",
   prog_name+2, prog_name+2 );
  
  getchar();
  
  printf( "\n\nExamples :\n\
\t%s -h : \n\
\t  Print this screen\n\n\
\t%s -f file.ae -n 1000 : \n\
\t  Read the file file.ae and do 1000 new generations\n\n\
\t%s --noX -f file.ae -p \"BACKUP_STEP 50\" : \n\
\t  Read file.ae and continue, writing a backup every 50 generations.\n\
\t  Graphical display is disabled (send SIGUSR1 to process to re-enable).\n\n",
   prog_name+2, prog_name+2, prog_name+2 );
}
