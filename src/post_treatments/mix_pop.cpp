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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_exp_manager.h>
#include <ae_population.h>



// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);
void print_version( void );





int main( int argc, char* argv[] )
{
  // 1) Initialize command-line option variables with default values
  char*   input_dir = NULL;
  int32_t num_gener = -1;
  char*   pop_file_name = NULL;
  char*   output_dir = NULL;


  // 2) Define allowed options
  const char * options_list = "hVi:g:p:o:";
  static struct option long_options_list[] = {
    { "help",     no_argument,        NULL, 'h' },
    { "version",  no_argument,        NULL, 'V' },
    { "input",    required_argument,  NULL, 'i' },
    { "gener",    required_argument,  NULL, 'g' },
    { "pop",      required_argument,  NULL, 'p' },
    { "out",      required_argument,  NULL, 'o' },
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
        print_help(argv[0]);
        exit( EXIT_SUCCESS );
      }
      case 'V' :
      {
        print_version();
        exit( EXIT_SUCCESS );
      }
      case 'i' :
      {
        input_dir = new char[strlen(optarg)+1];
        strcpy( input_dir, optarg );
        break;
      }
      case 'g' :
      {
        num_gener = atol(optarg);
        break;
      }
      case 'p' :
      {
        pop_file_name = new char[strlen(optarg)+1];
        strcpy( pop_file_name, optarg );
        break;
      }
      case 'o' :
      {
        output_dir = new char[strlen(optarg)+1];
        strcpy( output_dir, optarg );
        break;
      }
      default :
      {
        // An error message is printed in getopt_long, we just need to exit
        exit( EXIT_FAILURE );
      }
    }
  }

  // 4) Check for missing mandatory arguments
  if (input_dir == NULL)
  {
    printf( "%s: error: You must provide an input directory.\n", argv[0] );
    exit( EXIT_FAILURE );
  }
  if ( num_gener == -1 )
  {
    printf( "%s: error: You must provide a generation number.\n", argv[0] );
    exit( EXIT_FAILURE );
  }
  if (pop_file_name == NULL)
  {
    printf( "%s: error: You must provide a population file.\n", argv[0] );
    exit( EXIT_FAILURE );
  }
  if (output_dir == NULL)
  {
    printf( "%s: error: You must provide an output directory.\n", argv[0] );
    exit( EXIT_FAILURE );
  }


  // Load the experiment
  ae_exp_manager* exp = new ae_exp_manager();
  exp->load(input_dir, num_gener, false, false, false);

  // Load the population from which to get replacement individuals
  ae_population* pop = new ae_population(NULL);
  pop->load(pop_file_name, false);

  int32_t pop_size = exp->get_pop()->get_nb_indivs();

  // Create a new list of individuals
  ae_list<ae_individual*>* new_indivs = NULL;

  // Extract the first half of the Initial population into the new list
  new_indivs = exp->get_pop()->get_indivs()->extract_starting_sublist(pop_size/2);

  // Add the first half of the second population to the new list
  new_indivs->merge(pop->get_indivs()->extract_starting_sublist(pop_size/2));

  // Replace the list of individuals in our experiment with our new list
  exp->get_pop()->replace_population(new_indivs);

  // Evaluate the individuals vs the environment defined in the experiment
  // and sort the population accordingly
  exp->get_pop()->evaluate_individuals(exp->get_env());
  exp->get_pop()->sort_individuals();

  // Save the new experiment
  exp->save_copy(output_dir);


  delete exp;
  delete pop;
}


/*!
  \brief

*/
void print_help(char* prog_path)
{
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) prog_name++;
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
  printf( "%s: does nothing (template file).\n", prog_name );
  printf( "\n" );
  printf( "Usage : %s -h or --help\n", prog_name );
  printf( "   or : %s -V or --version\n", prog_name );
  printf( "   or : %s [-o] [-a arg]\n", prog_name );
  printf( "\nOptions\n" );
  printf( "  -h, --help\n\tprint this help, then exit\n\n" );
  printf( "  -V, --version\n\tprint version number, then exit\n\n" );
  printf( "  -o, --opt\n\toption with no argument)\n" );
  printf( "  -a, --arg argument\n\toption with an argument\n" );
}


/*!
  \brief Print aevol version number

*/
void print_version( void )
{
  printf( "aevol %s\n", VERSION );
}








