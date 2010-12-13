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
#include <ae_common.h>
#include <ae_param_overloader.h>
#include <ae_simulation.h>

#ifdef __X11
#include <ae_simulation_X11.h>
#endif





#ifndef __NO_X
void catch_usr1( int sig_num );
#endif

#ifdef __IN2P3
void catch_usr1( int sig_num );
#endif

void print_help( char* prog_name );




int main( int argc, char* argv[] )
{
  // Catch SIGUSR1 in handler catch_usr1. SIGUSR1 will toggle
#ifndef __NO_X
  signal( SIGUSR1, catch_usr1 );
#endif
  
  
  #ifdef DEBUG
    printf( "aevol is being run in DEBUG mode\n" );
  #endif
  
 
  // =================================================================
  //                      Get command-line options
  // =================================================================
  //
  // 1) Initialize command-line option variables with default values
  char*   initial_backup_file_name  = NULL;
  
  ae_param_overloader* param_overloader = new ae_param_overloader();
  
  
  #ifndef __NO_X
    bool show_display_on_startup   = true;
  #endif
  
  // 2) Define allowed options
  const char * options_list = "hxf:n:s:p:";
  static struct option long_options_list[] = {
    { "help",     no_argument, NULL, 'h' },
    { "noX",      no_argument, NULL, 'x' },
    { "file",     required_argument, NULL, 'f' },
    { "nbgener",  required_argument, NULL, 'n' },
    { "seed",     required_argument, NULL, 's' },
    { "param",    required_argument, NULL, 'p' },
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
      case 'f' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "ERROR : Option -f or --file : missing argument.\n" );
          exit( EXIT_FAILURE );
        }
        
        initial_backup_file_name = new char[strlen(optarg) + 1];
        sprintf( initial_backup_file_name, "%s", optarg );
        
        break;      
      }
      case 's' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "ERROR : Option -s or --seed : missing argument.\n" );
          exit( EXIT_FAILURE );
        }
        
        char* tmp = new char[50];
        sprintf( tmp, "SEED %s", optarg );
        param_overloader->store_overload( tmp );
        delete tmp;
        
        break;
      }
      case 'p' :
      {
        if( strcmp( optarg, "") == 0)
        {
          printf( "ERROR : Option -p or --param : missing agument.\n" );
          exit( EXIT_FAILURE );
        }
        
        param_overloader->store_overload( optarg );
        // NOTE that the actual overloading of the parameters is done after loading the parameter file
        
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
        param_overloader->store_overload( tmp );
        // NOTE that the actual overloading of the parameters is done after loading the parameter file
        delete [] tmp;
        
        break;
      }
    }
  }
  
  
  // =================================================================
  //                       Create the simulation
  // =================================================================
#ifdef __NO_X
  if ( initial_backup_file_name == NULL )
  {    
    // Create a new simulation
    ae_common::sim = new ae_simulation( param_overloader );
  }
  else
  {
    printf( "Loading simulation from backup file <%s>...\n", initial_backup_file_name );
    
    // Load simulation from backup
    ae_common::sim = new ae_simulation( initial_backup_file_name, true, param_overloader );
  }
#elif defined __X11
  if ( initial_backup_file_name == NULL )
  {
    
    // Create a new simulation
    ae_common::sim = ae_common::sim_display = new ae_simulation_X11( param_overloader );
  }
  else
  {
    printf( "Loading simulation from backup file <%s>...\n", initial_backup_file_name );
    
    // Load simulation from backup
    ae_common::sim = ae_common::sim_display = new ae_simulation_X11( initial_backup_file_name, true, param_overloader );
  }
#else
#error You must specify a graphic option
#endif
  
  delete [] initial_backup_file_name;
  delete param_overloader;
  ae_common::print_to_file();

    
#ifndef __NO_X
  if ( show_display_on_startup )
  {
    ae_common::sim_display->toggle_display_on_off();
  }
#endif
  
  // =================================================================
  //                         Run the simulation
  // =================================================================
  if ( ae_common::nb_generations > 0 )
  {
    ae_common::sim->run();
  }
  

  ae_common::clean(); // deletes the param overloader (among other things) 
}


#ifndef __NO_X
void catch_usr1( int sig_num )
{
  signal( SIGUSR1, catch_usr1 );
  
  printf( "display on/off\n" );
  ae_common::sim_display->toggle_display_on_off();
}
#endif

#ifdef __IN2P3
void catch_usr1( int sig_num )
{
  system( "./make_big_backup.py" );
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
\t%s -f file.ae --n 1000 : \n\
\t  Read the file file.ae and do 1000 new generations\n\n\
\t%s --noX -f file.ae -p \"BACKUP_STEP 50\" : \n\
\t  Read file.ae and continue, writing a backup every 50 generations.\n\
\t  Graphical display is disabled (send SIGUSR1 to process to re-enable).\n\n",
   prog_name+2, prog_name+2, prog_name+2 );
}
