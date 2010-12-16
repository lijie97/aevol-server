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
#include <ae_common.h>
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_list.h>
#include <ae_simulation_X11.h>





void print_help( char* prog_name );





int main( int argc, char* argv[] )
{

  

  // =================================================================
  //                      Get command-line options
  // =================================================================
  //
  // 1) Initialize command-line option variables with default values

  char* backup_file_name  = NULL;
  
  // 2) Define allowed options
  const char * options_list = "hf:";
  static struct option long_options_list[] = {
    { "file", 1, NULL, 'f' },
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
            case 'f' :
              backup_file_name = new char[strlen(optarg) + 1];
              sprintf( backup_file_name, "%s", optarg );
              break;
            }
  }
  
  
  // =================================================================
  //                       Read the backup file
  // =================================================================

  if ( backup_file_name == NULL )
    {
      printf("You must specify a backup file. Please use the option -f or --file.\n");
      exit(EXIT_FAILURE);
    }
  else
    {
      printf( "Reading backup file <%s>... ", backup_file_name );
      fflush(stdout);

      // Load simulation from backup
      ae_common::sim = ae_common::sim_display = new ae_simulation_X11( backup_file_name, false );
      printf("done\n");
    }



  // =================================================================
  //                      Evaluate the individuals
  // =================================================================


  // The constructor of the ae_simulation has read the genomes of the individuals
  // and located their promoters, but has not performed the translation nor the
  // phenotype computation. We must do it now.
  
  printf( "Computing phenotypes... " );
  fflush( stdout );

  ae_list_node  * indiv_node = ((ae_common::sim->get_pop())->get_indivs())->get_first();
  ae_individual * indiv      = NULL;


  while( indiv_node != NULL )
  {
    indiv = (ae_individual *) indiv_node->get_obj();
    
    // ===============  Raevol-related remark ========================
    // evaluate is a virtual method. Thus the method actually executed
    // is not determined at the compilation but during the execution.
    // The REAL type of the object (ae_individual or ae_individual_R)
    // will be determined. Depending on that, the method actually 
    // executed will be either ae_individual::evaluate() or
    // ae_individual_R::evaluate().
    indiv->evaluate( ae_common::sim->get_env() );

    indiv_node = indiv_node->get_next();
  }

  printf( "done\n" );


  // =================================================================
  //                       Draw the windows
  // =================================================================


  ae_common::sim_display->toggle_display_on_off();
  ae_common::sim_display->display();
  while ( ae_common::sim_display->quit_signal_received() == false )
  {
    ae_common::sim_display->handle_events();
  }
  


  delete [] backup_file_name;

  delete ae_common::sim_display;

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
  printf( "\t-f or --file s    : Begin with the backup file s\n" );
}
