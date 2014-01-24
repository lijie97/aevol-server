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
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>

#ifdef __NO_X
#error This program requires graphics libraries
#else
#include <X11/Xlib.h>
#endif

// =================================================================
//                            Project Files
// =================================================================
#include <ae_exp_manager_X11.h>
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_list.h>





void print_help( char* prog_name );


int main( int argc, char* argv[] )
{

  

  // =================================================================
  //                      Get command-line options
  // =================================================================
  //
  // 1) Initialize command-line option variables with default values
  int32_t     gener         = 0;
  
  // 2) Define allowed options
  const char * options_list = "hg:";
  static struct option long_options_list[] = {
    { "gener", 1, NULL, 'g' },
    { 0, 0, 0, 0 }
  };

  // 3) Get actual values of the command-line options
  int option;
  while ( ( option = getopt_long(argc, argv, options_list, long_options_list, NULL) ) != -1 ) 
  {
          switch ( option ) 
            {
            case 'h' :
              print_help( argv[0] );
              exit( EXIT_SUCCESS );
              break;
            case 'g' :
	      if ( strcmp( optarg, "" ) == 0 )
		{
		  printf( "%s: error: Option -g or --gener : missing argument.\n", argv[0] );
		  exit( EXIT_FAILURE );
		}
        
	      gener = atol( optarg );
              break;
            }
  }
  
  
  // =================================================================
  //                       Read the backup file
  // =================================================================

  printf( "Displaying generation %"PRId32"... ", gener );
  fflush(stdout);

  // Load simulation from backup
  ae_exp_manager_X11* exp_manager = new ae_exp_manager_X11();
  exp_manager->load( gener, false, true, false );
  printf("done\n");
  


  // =================================================================
  //                       Draw the windows
  // =================================================================


  exp_manager->toggle_display_on_off();
  exp_manager->display();
  while ( exp_manager->quit_signal_received() == false )
  {
    exp_manager->handle_events();
  }
  


  delete exp_manager;

  return EXIT_SUCCESS;
}


void print_help( char* prog_name ) 
{
  printf( "\n************* aevol - Artificial Evolution ************* \n\n" );
  printf( "This program is Free Software. No Warranty.\n" );
  printf( "Copyright (C) 2009  LIRIS.\n" );
  printf( "Usage : %s -h\n", prog_name );
  printf( "   or : %s -f file.ae\n", prog_name );
  printf( "\t-h : Display this screen\n" );
  printf( "\t-g or --gener n    : Display the population at generation n\n" );
}
