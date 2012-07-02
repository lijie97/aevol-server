// This programs loads an individual previously saved in a bin file
// by get_best_indiv, create a new simulation in the current
// directory (which has to contain a param.in file), and initialize
// the population of this simulation with clones of this individual.
//
// This allow us to start a new simulation from a population formed
// by clones a choosen individual instead of clones of a random and 
// almost-flat individual.
//
// Example :
//
// mkdir next_run
// cd previous_run
// get_best_indiv -f backup/gen_020000.ae -t ../next_run/best.i
// cd ../next_run
// (create a param.in file)
// set_genome_from_indiv -f best.i
// aevol -f backup/gen_000000.ae -n 20000
//
// This example with extract the best individual from previous_run
// at generation 20000, create a new simulation new_run from clones
// of this individual and run this new simulation for 20000 gen.
//
// Nothing else is kept from previous_run, all the parameters are
// read from next_run/param.in. The number of clone is read from
// INIT_POP_SIZE in this file.
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
#include <ae_experiment.h>
#include <ae_param_loader.h>





void print_help( char* prog_name );





int main( int argc, char* argv[] )
{
  // Initialize command-line option variables with default values
  char* init_file_name  = NULL;

  // Define allowed options
  const char * options_list = "hf:";
  static struct option long_options_list[] = {
    { "from", 1, NULL, 'f' },
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
            case 'f' :
              init_file_name = new char[strlen(optarg) + 1];
              sprintf( init_file_name, "%s", optarg );
              break;			  
            }
  }
  
  
  if ( init_file_name == NULL )
  {
    printf("You must specify an indiv file. Please use the option -f or --from.\n");
    exit(EXIT_FAILURE);
  }
  
  int32_t _num_gener=0;
  
  //~ // We create a new simulation
  //~ printf("Creating the new simulation.\n");
  ae_experiment* sim1 = new ae_experiment();
  //~ ae_param_loader* param_loader = new ae_param_loader( "param.in" );
  //~ sim1->load_params( param_loader, NULL );
  //~ delete param_loader;

  
  // We modify its individuals
  printf("Loading the individual.\n");  
  gzFile* init_file = (gzFile*) gzopen( init_file_name, "r" );
  ae_individual* indiv = new ae_individual( init_file );
  printf("Closing the file.\n");
  gzclose(init_file);
  
  ae_list* list_indivs = NULL;//sim1->get_pop()->get_indivs();
  int32_t nb_indivs = NULL;//sim1->get_pop()->get_nb_indivs();
  
  ae_list_node*  indiv_node = list_indivs->get_first();
  int i;
  for (i=0;i<nb_indivs;i++)
  {
    ae_individual* new_indiv = new ae_individual( *indiv );
    ae_individual* old_indiv = (ae_individual*) indiv_node->get_obj();
    delete old_indiv;
    indiv_node->set_obj(new_indiv);
    indiv_node = indiv_node->get_next();
  }
	
  if ( ae_common::pop_structure == true )
  {
    ae_list_node*  indiv_node = list_indivs->get_first();
    ae_grid_cell*** _pop_grid = NULL;//sim1->get_pop()->get_pop_grid();
    for ( int16_t x = 0 ; x < ae_common::grid_x ; x++ )
    {
      for ( int16_t y = 0 ; y < ae_common::grid_y ; y++ )
      {
        _pop_grid[x][y]->set_individual((ae_individual*) indiv_node->get_obj());
        _pop_grid[x][y]->get_individual()->set_grid_cell( _pop_grid[x][y] );
        indiv_node = indiv_node->get_next();
      }
    }
  }
  
  // Evaluate the new individuals
  ae_population* pop = NULL;
  /*sim1->get_pop()*/pop->evaluate_individuals( sim1->get_env() );
  
  // Delete the backup files created by the ae_experiment constructor
  char backup_file_name[50];
  char best_indiv_file_name[50];
  
#ifdef __REGUL
  sprintf( backup_file_name, "backup/gen_%06"PRId32".rae", _num_gener );
  sprintf( best_indiv_file_name, "backup/best_%06"PRId32".rae", _num_gener );
#else
  sprintf( backup_file_name, "backup/gen_%06"PRId32".ae", _num_gener );
  sprintf( best_indiv_file_name, "backup/best_%06"PRId32".ae", _num_gener );
#endif
  
  remove(backup_file_name);
  remove(best_indiv_file_name);
  
  // And create new ones
  printf("Writting the backup file\n");
  sim1->write_backup();
  return EXIT_SUCCESS;

}

void print_help( char* prog_name ) 
{
  printf( "\n************* aevol - Artificial Evolution ************* \n\n" );
  printf( "This program is Free Software. No Warranty.\n\n\
Usage : set_genome_from_indiv -h\n\
   or : set_genome_from_indiv -f source\n\n\
\t-h : Display this screen\n\
\t--from source    : Create a new simulation from param.in then load the individual from the file 'source'\n");
}
