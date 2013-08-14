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
void print_version( void );




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
  bool use_text_files   = false;
  bool pause_on_startup = false;
  bool verbose          = false;
  
  int32_t num_gener = 0;
  int32_t nb_gener  = 1000;
  
  #ifndef __NO_X
    bool show_display_on_startup = true;
  #endif
  
  
  // -------------------------------------------------------------------------
  // 2) Define allowed options
  // -------------------------------------------------------------------------
  const char * options_list = "hn:r:tvVwx";
  static struct option long_options_list[] = {
    { "help",     no_argument,        NULL, 'h' },
    { "nbgener",  required_argument,  NULL, 'n' }, // Number of generations to be run
    { "resume",   required_argument,  NULL, 'r' }, // Resume from generation X
    { "text",     no_argument,        NULL, 't' }, // Use text files instead of gzipped binary files
    { "verbose",  no_argument,        NULL, 'v' }, // Be verbose
    { "version",  no_argument,        NULL, 'V' },
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
      case 'V' :
      {
        print_version();
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
      case 't' :
      {
        use_text_files = true;
        printf( "%s: error: Option -t or --text not yet implemented.\n", argv[0] );
        exit( EXIT_FAILURE );
        
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
      default :
      {
        // An error message is printed in getopt_long, we just need to exit
        exit( EXIT_FAILURE );
      }
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
  
  exp_manager->load( num_gener, false, verbose );
  exp_manager->set_nb_gener( nb_gener );
  
  /*ae_list_node<ae_individual*>*   indiv_node_tmp = exp_manager->get_pop()->get_indivs()->get_first();
  ae_individual*  indiv_tmp;
  /*for ( int32_t i = 0 ; i < exp_manager->get_pop()->get_nb_indivs() ; i++ )
  {
    indiv_tmp = indiv_node_tmp->get_obj();
    printf("%le %"PRId32" %"PRId32"\n", indiv_tmp->get_fitness(), indiv_tmp->get_nb_bases_in_0_coding_RNA(), indiv_tmp->get_int_probes()[0]);
    indiv_node_tmp = indiv_node_tmp->get_next();
  }*/
  
  // Make a numbered copy of each static input file (dynamic files are saved elsewhere)
  // TODO (?)

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







/*!
  \brief 
  
*/
void print_help( char* prog_path ) 
{
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ( prog_name = strrchr( prog_path, '/' ) ) prog_name++;
  else prog_name = prog_path;
  
	printf( "******************************************************************************\n" );
	printf( "*                                                                            *\n" );
	printf( "*                        aevol - Artificial Evolution                        *\n" );
	printf( "*                                                                            *\n" );
	printf( "* Aevol is a simulation platform that allows one to let populations of       *\n" );
  printf( "* digital organisms evolve in different conditions and study experimentally  *\n" );
  printf( "* the mechanisms responsible for the structuration of the genome and the     *\n" );
  printf( "* transcriptome.                                                             *\n" );
	printf( "*                                                                            *\n" );
	printf( "******************************************************************************\n" );
  printf( "\n" );
	printf( "%s: run an aevol simulation.\n", prog_name );
  printf( "\n" );
	printf( "Usage : %s -h or --help\n", prog_name );
	printf( "   or : %s -V or --version\n", prog_name );
	printf( "   or : aevol_run_X11 [-r GENER] [-n NB_GENER] [-tvwx]\n", prog_name );
	printf( "\nOptions\n" );
	printf( "  -h, --help\n\tprint this help, then exit\n\n" );
	printf( "  -V, --version\n\tprint version number, then exit\n\n" );
  printf( "  -r, --resume GENER\n\tspecify generation to resume simulation at (default 0)\n\n" );
  printf( "  -n, --nbgener NB_GENER\n\tspecify number of generations to be run (default 1000)\n\n" );
  printf( "  -t, --text\n\tuse text files instead of binary files when possible\n\n" );
	printf( "  -v, --verbose\n\tbe verbose\n\n" );
  printf( "  -w, --wait\n\tpause after loading\n\n" );
  printf( "  -x, --noX\n\tdon't display X outputs upon start\n\tsend SIGUSR1 to switch X output on/off\n" );
}

/*!
  \brief 
  
*/
void print_version( void ) 
{
	printf( "aevol %s\n", VERSION );
}