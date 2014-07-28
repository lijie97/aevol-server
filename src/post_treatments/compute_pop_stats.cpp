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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <inttypes.h>
#include <getopt.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>




// =================================================================
//                            Project Files
// =================================================================
#include <ae_macros.h>
#include <ae_utils.h>
#include <population_statistics.h>
#include <ae_exp_manager.h>






// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);
void print_version( void );



// =====================================================================
//                         Main Function
// =====================================================================


//#define FV_FILE "fv.out"
//#define REP_FILE "replications.out"

int main( int argc, char* argv[] )
{
  // ----------------------------------------
  //     command-line option parsing
  // ----------------------------------------
  int32_t     nb_children       = 1000;
  int32_t     backup_step       = 0;
  int32_t     generation_number = -1;
  int32_t     wanted_rank       = -1;
  int32_t     wanted_index      = -1;
  int32_t     begin_generation  = 0;

  const char * options_list = "hVe:b:r:i:n:";
  static struct option long_options_list[] = {
  	{"help",        no_argument,        NULL, 'h'},
    {"version",     no_argument,        NULL, 'V' },
    {"end",         required_argument,  NULL, 'e' },
    {"begin",       required_argument,  NULL, 'b' },
    {"rank",        required_argument,  NULL, 'r'},
    {"index",       required_argument,  NULL, 'i'},
    {"nb-children", required_argument,  NULL, 'n'},
    {0, 0, 0, 0}
  };

  int option = -1;
  while((option=getopt_long(argc,argv,options_list,long_options_list,NULL))!=-1)
  {
    switch(option)
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
      case 'e' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -e or --end : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }
        generation_number = atol( optarg );
        break;
      }
      case 'n' :
        nb_children = atol(optarg);
        break;
      case 'b':
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -b or --begin : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }
        begin_generation = atol( optarg );
        break;
      }
      case 'r' :
        wanted_rank = atol(optarg);
        wanted_index = -1;
        break;
      case 'i' :
        wanted_index = atol(optarg);
        wanted_rank = -1;
        break;
    }
  }

  if(wanted_rank == -1 && wanted_index ==-1){
    wanted_rank = 1;
  }

  if ( generation_number == -1 )
  {
    printf( "%s: error: You must provide a generation number.\n", argv[0] );
    exit( EXIT_FAILURE );
  }

  analysis_type type = MULTIPLE_GENERATIONS;

  population_statistics* population_statistics_compute = new population_statistics(type, nb_children, wanted_rank, wanted_index);

  // Load simulation
  ae_exp_manager* exp_manager = new ae_exp_manager();
  exp_manager->load( begin_generation, false, true, false );
  backup_step = exp_manager->get_backup_step();

  printf("\n\n Generation : %d\n\n", begin_generation);
  population_statistics_compute->compute_reproduction_stats(exp_manager,begin_generation);
  population_statistics_compute->compute_population_stats(begin_generation);
  delete exp_manager;

  for ( int32_t i = begin_generation + backup_step  ; i <= generation_number ; i += backup_step )
    {
      printf("\n\n Generation : %d\n\n", i);
      exp_manager = new ae_exp_manager();
      exp_manager->load( i, false, true, false );
      
      population_statistics_compute->compute_reproduction_stats(exp_manager,i);
      population_statistics_compute->compute_population_stats(i);
      delete exp_manager;
    }
  
  delete population_statistics_compute;

  return EXIT_SUCCESS;
}




/*!
  \brief

*/
void print_help(char* prog_path)
{
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ( ( prog_name = strrchr( prog_path, '/' ) ) ) prog_name++;
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
  printf( "%s: computes some population and replication statistics at each available backup from begin_gener until end_gener.\n", prog_name );
  printf( "These statistics are saved in global_pop_stats.out and in files inside stats/pop_stats.\n");
  printf( "The replication statistics (information about the children_nb offsprings) of the individual of given rank or index are written at each backup.\n" );
  printf( "\n" );
  printf( "WARNING: This program should not be used for simulations run with lateral\n" );
  printf( "transfer. When an individual has more than one parent, the notion of lineage\n" );
  printf( "used here is not relevant.\n" );
  printf( "\n" );
  printf( "Usage : %s -h\n", prog_name);
  printf( "   or : %s -V or --version\n", prog_name );
  printf( "   or : %s -e GENER2 [-b GENER1] [-r RANK | -i INDEX]\n", prog_name );
  printf( "\nOptions\n" );
  printf( "  -h, --help\n\tprint this help, then exit\n\n" );
  printf( "  -V, --version\n\tprint version number, then exit\n\n" );
  printf( "  -b GENER1 or --begin GENER1 :\n" );
  printf( "\tFirst backup used to compute the statistics\n" );
  printf( "\t-e GENER2 or --end GENER2 :\n" );
  printf( "\tLast backup used to compute the statistics\n" );
  printf( "\t-i INDEX or --index INDEX :\n" );
  printf( "\tIndex of individual of whom we want information about the offsprings at each backup\n" );
  printf( "\n" );
  printf( "\t-r RANK or --rank RANK :\n" );
  printf( "\tRank of individual of whom we want information about the offsprings at each backup\n" );
}


/*!
  \brief Print aevol version number

*/
void print_version( void )
{
  printf( "aevol %s\n", VERSION );
}

