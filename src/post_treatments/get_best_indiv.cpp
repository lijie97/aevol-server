// This program load a backup file, extract the best individual and
// save it into a binary file. It is intended to work with
// set_simulation_from_indiv to create a population of clones of one
// selected individual with a choosed phenotype and to start a
// new simulation from this population.
//
// Example : get_best_indiv -f backup/gen_020000.ae -t best_020000.i
//


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
#include <zlib.h>


// =================================================================
//                            Project Files
// =================================================================
#include <ae_common.h>
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_simulation.h>





void print_help( char* prog_name );





int main( int argc, char* argv[] )
{
  // Initialize command-line option variables with default values
  char* backup_file_name  = NULL;
  char* indiv_file_name  = NULL;

  // Define allowed options
  const char * options_list = "hf:t:";
  static struct option long_options_list[] = {
    { "from", 1, NULL, 'f' },
    { "to", 1, NULL, 't' },
    { 0, 0, 0, 0 }
  };

  // Get actual values of the command-line options
  int option;
  while ( ( option = getopt_long(argc, argv, options_list, long_options_list, NULL) ) != -1 ) 
  {
          switch ( option ) 
            {
            case 'h' :
              print_help( argv[0] );
              exit( EXIT_SUCCESS );
              break;
            case 't' :
              indiv_file_name = new char[strlen(optarg) + 1];
              sprintf( indiv_file_name, "%s", optarg );
              break;
            case 'f' :
              backup_file_name = new char[strlen(optarg) + 1];
              sprintf( backup_file_name, "%s", optarg );
              break;			  
            }
  }
  

  if ( indiv_file_name == NULL )
  {
    printf("You must specify an indiv file. Please use the option -t or --to.\n");
    exit(EXIT_FAILURE);
  }
  
  if ( backup_file_name == NULL )
  {
    printf("You must specify a backup file. Please use the option -f or --from.\n");
    exit(EXIT_FAILURE);
  }
  

  // Load simulation from backup
  gzFile indiv_file = gzopen( indiv_file_name, "w" );
  ae_common::sim = new ae_simulation( backup_file_name, false );
  
  // Evaluate the individuals
  (ae_common::sim->get_pop())->evaluate_individuals(ae_common::sim->get_env());
  
  // Saving the best one
  ae_common::sim->get_pop()->get_best()->write_to_backup(indiv_file);
  
  // Close the file
  gzclose(indiv_file);
  
  return EXIT_SUCCESS;

}

void print_help( char* prog_name ) 
{
  printf( "\n************* aevol - Artificial Evolution ************* \n\n" );
  printf( "This program is Free Software. No Warranty.\n\n\
Usage : get_best_indiv -h\n\
   or : get_best_indiv -f source -t destination \n\n\
\t-h : Display this screen\n\
\t--from source    : Load the individuals from the backup file 'source'\n\
\t--to destination   : Save the best individual to the indiv file 'destination'\n");
}
