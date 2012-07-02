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
//~ #ifdef __X11
  //~ #include <ae_exp_setup_X11.h>
  //~ #include <ae_population_X11.h>
//~ #elif defined __NO_X
  //~ #include <ae_exp_setup.h>
  //~ #include <ae_population.h>
//~ #else
  //~ #error You must specify a graphic option
//~ #endif

#include <ae_exp_manager.h>
//~ #include <ae_param_overloader.h>
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
  char* exp_setup_file_name = NULL;
  char* pop_file_name       = NULL;
  char* out_man_file_name   = NULL;
  int32_t num_gener = 0;
  
  //~ ae_param_overloader* param_overloader = new ae_param_overloader();
  
  
  #ifndef __NO_X
    bool show_display_on_startup   = true;
  #endif
  
  // 2) Define allowed options
  const char * options_list = "hr:e:p:o:xn:";
  static struct option long_options_list[] = {
    { "help",     no_argument,        NULL, 'h' },
    { "restart",  required_argument,  NULL, 'r' },
    { "exp",      required_argument,  NULL, 'e' },
    { "pop",      required_argument,  NULL, 'p' },
    { "out",      required_argument,  NULL, 'o' },
    { "noX",      no_argument,        NULL, 'x' },
    { "nbgener",  required_argument,  NULL, 'n' },
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
      case 'x' :
      {
        #ifdef __NO_X
          printf( "Program was compiled with __NO_X option, no visualisation available.\n" );
        #else
          show_display_on_startup = false;
        #endif
        
        break;
      }
      case 'r' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "ERROR : Option -r or --restart : missing argument.\n" );
          exit( EXIT_FAILURE );
        }
        
        num_gener = atol( optarg );
        
        exp_setup_file_name = new char[255];
        pop_file_name       = new char[255];
        out_man_file_name   = new char[255];
        
        sprintf( exp_setup_file_name, EXP_SETUP_BACKUP_FNAME_FORMAT,  num_gener );
        sprintf( pop_file_name,       POP_BACKUP_FNAME_FORMAT,        num_gener );
        sprintf( out_man_file_name,   OUT_PROF_BACKUP_FNAME_FORMAT,   num_gener );
        
        break;      
      }
      case 'n' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "ERROR : Option -n or --nbgener : missing argument.\n" );
          exit( EXIT_FAILURE );
        }
        
        char* tmp = new char[50];
        sprintf( tmp, "NB_GENER %s", optarg );
        //~ param_overloader->store_overload( tmp );
        // NOTE that the actual overloading of the parameters is done after loading the parameter file
        delete [] tmp;
        
        break;
      }
    }
  }
  
  
  // =================================================================
  //                          Load the simulation
  // =================================================================
  ae_exp_manager* exp_manager = new ae_exp_manager();
  
  exp_manager->load_experiment( exp_setup_file_name, pop_file_name, out_man_file_name, true );
  
  delete [] pop_file_name;
  delete [] exp_setup_file_name;
  //~ delete param_overloader;

    
#ifndef __NO_X
  //~ if ( show_display_on_startup )
  //~ {
    //~ ae_common::sim_display->toggle_display_on_off();
  //~ }
#endif
  
  // =================================================================
  //                         Run the simulation
  // =================================================================
  //~ if ( ae_common::nb_generations > 0 )
  //~ {
    //~ exp_manager->run_evolution();
  //~ }
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
