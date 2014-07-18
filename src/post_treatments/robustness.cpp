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
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <inttypes.h>
#include <getopt.h>
#include <math.h>
#include <sys/stat.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_macros.h>
#include <population_statistics.h>

#ifndef __NO_X
  #include <ae_exp_manager_X11.h>
#else
  #include <ae_exp_manager.h>
#endif



// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);
void print_version( void );





// =====================================================================
//                         Main Function
// =====================================================================


int main( int argc, char* argv[] )
{
  // ----------------------------------------
  //     command-line option parsing
  // ----------------------------------------
  int32_t nb_children     = 1000;
  int32_t wanted_rank     = -1;
  int32_t wanted_index    = -1;
  int32_t num_gener       = 0;

  const char * options_list = "hVg:n:r:i:";
  static struct option long_options_list[] = {
    {"help",          no_argument,        NULL, 'h'},
    {"version",       no_argument,        NULL, 'V'},
    {"gener",           required_argument,  NULL, 'g'},
    {"nb-children",   required_argument,  NULL, 'n'},
    {"rank",          required_argument,  NULL, 'r'},
    {"index",         required_argument,  NULL, 'i'},
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
      case 'g' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -g or --gener : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }

        num_gener = atol( optarg );
        break;
      }
      case 'n' :
        nb_children = atol(optarg);
        break;
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

  analysis_type type = ONE_GENERATION;

  population_statistics* population_statistics_compute = new population_statistics(type, nb_children, wanted_rank, wanted_index);

  // Load simulation
  #ifndef __NO_X
    ae_exp_manager* exp_manager = new ae_exp_manager_X11();
  #else
    ae_exp_manager* exp_manager = new ae_exp_manager();
  #endif
  exp_manager->load( num_gener, false, true, false );

  population_statistics_compute->compute_reproduction_stats(exp_manager, num_gener);

  delete exp_manager;
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
  if ( ( prog_name = strrchr( prog_path, '/' )) ) prog_name++;
  else prog_name = prog_path;

  printf( "\n" );
  printf( "*********************** aevol - Artificial Evolution ******************* \n" );
  printf( "*                                                                      * \n" );
  printf( "*                     Robustness post-treatment program                * \n" );
  printf( "*                                                                      * \n" );
  printf( "************************************************************************ \n" );
  printf( "\n\n" );
  printf( "This program is Free Software. No Warranty.\n" );
  printf( "\n" );
  printf( "Usage : %s -h\n", prog_name);
  printf( "   or : %s -V or --version\n", prog_name );
  printf( "   or : %s -g numgener [-n nbchildren] [-r rank | -i index]\n", prog_name);
  printf( "\n" );
  printf( "This program computes the replication statistics of all the individuals of a given generation,\n");
  printf( "like the proportion of neutral, beneficial, deleterious offsprings. This is done by simulating\n");
  printf( "\'nbchildren\' replications for each individual, with its mutation, rearrangement and transfer rates.\n" );
  printf( "Depending on those rates and genome size, there can be several events per replication.\n" );
  printf( "Those global statistics are written in %s/robustness_numgener.out, one line per individual\n", STATS_DIR );
  printf( "in the specified generation.\n\n" );
  printf( "The program also outputs detailed statistics for one of the individuals (the best one by default). \n");
  printf( "The detailed statistics for this individual are written in %s/replication_numgener.out,\n", STATS_DIR);
  printf( "with one line per simulated child of this particular individual.\n");
  printf( "\n" );
  printf( "\n" );
  printf( "\t-h or --help    : Display this help, then exit\n" );
  printf( "\n" );
  printf( "\t-V or --version : Print version number, then exit\n" );
  printf( "\n" );
  printf( "\t-g numgener or --gener numgener : \n" );
  printf( "\t                  Generation at which the statistics are computed\n" );
  printf( "\n" );
  printf( "\t-n nbchildren or --nb-children nbchildren : \n" );
  printf( "\t                  Use nbchildren replications per individual to compute its statistics. Default = 1000.\n" );
  printf( "\n" );
  printf( "\t-i index or --index index : \n" );
  printf( "\t                  Index of individual for whom we want detailed information about the simulated offsprings\n" );
  printf( "\n" );
  printf( "\t-r rank or --rank rank : \n" );
  printf( "\t                  Rank of individual for whom we want detailed information about the simulated offsprings\n" );


  printf( "\n" );
}


/*!
  \brief Print aevol version number

*/
void print_version( void )
{
  printf( "aevol %s\n", VERSION );
}
