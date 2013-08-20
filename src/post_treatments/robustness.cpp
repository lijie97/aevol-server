/** \class 
 *  \brief 
 *  \author  Carole Knibbe
 *  Input files :
 *               
 *  Output file :
 *               
 */


// =======================================================================
//                        Standard Libraries
// =======================================================================
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <inttypes.h>
#include <getopt.h>
#include <math.h>
#include <sys/stat.h>


// =======================================================================
//                        Project Libraries
// =======================================================================

#include <population_statistics.h>
#ifndef __NO_X
  #include <ae_exp_manager_X11.h>
#else
  #include <ae_exp_manager.h>
#endif



// =======================================================================
//                       Secondary Functions
// =======================================================================

void print_help( void );

// function copied from ae_individual's computation of experimental fv.
// In addition, it provides and prints information about replications
/*double compute_experimental_fv( ae_individual* indiv, int nb_children, double* neutral_or_better, FILE* replication_file );

// count how many proteins were modified after replication
int count_affected_genes( ae_individual* parent, ae_individual* child );*/


// TODO: update this function...
// reconstruct final individual from backup and lineage
// ae_individual * get_final_individual_using_dstory();     



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
  int32_t num_gener       = 100;

  const char * options_list = "he:n:r:i:"; 
  static struct option long_options_list[] = {
    {"help",          no_argument,        NULL, 'h'},
    {"end",           required_argument,  NULL, 'e' }, 
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
      case 'h' : print_help(); exit(EXIT_SUCCESS);  break;
      case 'e' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -e or --end : missing argument.\n", argv[0] );
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

void print_help( void )
{
  printf( "\n" ); 
  printf( "*********************** aevol - Artificial Evolution ******************* \n" );
  printf( "*                                                                      * \n" );
  printf( "*                     Robustness post-treatment program                * \n" );
  printf( "*                                                                      * \n" );
  printf( "************************************************************************ \n" );
  printf( "\n\n" ); 
  printf( "This program is Free Software. No Warranty.\n" );
  printf( "Copyright (C) 2009  LIRIS.\n" );
  printf( "\n" ); 
  printf( "Usage : robustness -h\n");
  printf( "or :    robustness -e end_gener [-n children_nb] [-r rank | -i index]\n");
  printf( "\n" ); 
  printf( "This program computes replication statistics at backup at end_gener like proportion of neutral, beneficial, deleterious offsprings\n" );
  printf( "and statistics about the offsprings, written in robustness_end_gener.out.\n");
  printf( "The replication statistics (information about the children_nb offsprings) of the individual of rank or index are written in replication_end_gener.out.\n");
  printf( "\n" ); 
  printf( "WARNING: This program should not be used for simulations run with lateral\n" ); 
  printf( "transfer. When an individual has more than one parent, the notion of lineage\n" ); 
  printf( "used here is not relevant.\n" );
  printf( "\n" );  
  printf( "\t-h or --help    : Display this help.\n" );
  printf( "\n" ); 
  printf( "\t-n children_nb or --nb-children children_nb : \n" );
  printf( "\t                  Use children_nb to compute Fv.\n" );
  printf( "\n" ); 
  printf( "\t-i index or --index index : \n" );
  printf( "\t                  Index of individual of whom we want information about the offsprings at each backup\n" );
  printf( "\n" ); 
  printf( "\t-r rank or --rank rank : \n" );
  printf( "\t                  Rank of individual of whom we want information about the offsprings at each backup\n" );
  printf( "\n" );
  printf( "\t-e end_gener or --end end_gener : \n" );
  printf( "\t                  Generation at which the statistics are computed\n" );
  printf( "\n" );

}

//
// Data concerning individuals wanted is printed in one or two files.
//
// The first file (output_dir/fv.out) contains information about Fv:
//
//       "r_i_fit_pr_fvexp_fvexpnob_fvth_N\n"
//
// where:
//     * _ is a blank space
//     * r:        Rank of individual in population
//     * i:        Index of individual in population
//     * fit:      FITness of the individual
//     * pr:       Probability of Reproduction
//     * fvexp:    EXPerimental estimation of Fv
//     * fvexpnob: EXPerimental estimation of Fv counting children Neutral Or Better
//     * fvth:     THeoritical estimation of Fv
//     * N:        Number of individuals in population
//      
// The second file (output_dir/replications.out), if activated (-d option),
// contains information about children obtained when experimentally estimating Fv:
//
//       "#i_gl_nfg_cl_f_afc_sdfc_psp_pbp_pnga_anga\n"
//
// where:
//     * # is # and _ is a blank space
//     * i:     Index of individual in population
//     * gl:    total Genome Length
//     * nfg:   Number of Functional Genes
//     * cl:    total Coding Length
//     * f:     Fitness of individual
//     * afc:   Average Fitness of the Children
//     * sdfc:  Standard Deviation of Fitness values of the Children
//     * psp:   Proportion of children with Same fitness as Parent
//     * pbp:   Proportion of children Better than Parent
//     * pnga:  Proportion of children with No Gene Affected
//     * anga:  Average Number of Gene Affected by replication
//
// One additional line per child is printed with its fitness and the number of affected 
// genes ("fc_nga\n").
//
//
// Examples :
//
// For generation 20000, compute statistics for all the
// individuals and print them in directory out_020000 :
//
//    robustness -f backup/gen_020000.ae -o out_020000
//
// For generation 20000, write the best individual's statistics in
// out_020000_best with details about replications :
//
//    robustness -d -r 1 -f backup/gen_020000.ae -o out_020000_best
//

