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


//#define FV_FILE "fv.out"
//#define REP_FILE "replications.out"

int main( int argc, char* argv[] ) 
{
  // ----------------------------------------
  //     command-line option parsing
  // ----------------------------------------
  char* env_file_name = NULL;
  char* pop_file_name = NULL;
  char* output_dir    = NULL;
  int nb_children     = 1000;
  int wanted_rank     = -1;
  int wanted_index    = -1;
  bool details        = false;
  int32_t num_gener   = 100;  
  
  char* exp_setup_file_name = new char[63];
  char* out_prof_file_name  = new char[63];
  strcpy( exp_setup_file_name,  "exp_setup.ae" );
  strcpy( out_prof_file_name,   "output_profile.ae" );
  char* sp_struct_file_name = NULL;

  const char * options_list = "h:e:o:n:r:i:d"; 
  static struct option long_options_list[] = {
    {"help",          no_argument,        NULL, 'h'},
    {"end",           required_argument,  NULL, 'e' }, 
    {"output",        1,                  NULL, 'o'},
    {"nb-children",   1,                  NULL, 'n'},
    {"rank",          1,                  NULL, 'r'},
    {"index",         1,                  NULL, 'i'},
    {"with-details",  0,                  NULL, 'd'},
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
        
        env_file_name       = new char[255];
        pop_file_name       = new char[255];
        sp_struct_file_name = new char[255];
        
        sprintf( env_file_name,       ENV_FNAME_FORMAT,       num_gener );
        sprintf( pop_file_name,       POP_FNAME_FORMAT,       num_gener );
        sprintf( sp_struct_file_name, SP_STRUCT_FNAME_FORMAT, num_gener );
		  
        // Check existence of optional files in file system.
        // Missing files will cause the corresponding file_name variable to be nullified
        struct stat stat_buf;
        if ( stat( sp_struct_file_name, &stat_buf ) == -1 )
        {
          if ( errno == ENOENT )
          {
            delete [] sp_struct_file_name;
            sp_struct_file_name = NULL;
          }
          else
          {
            printf( "%s:%d: error: unknown error.\n", __FILE__, __LINE__ );
            exit( EXIT_FAILURE );
          }
        }
        
        break;
      }
    case 'o' : 
      output_dir = new char[strlen(optarg) + 1];
      sprintf( output_dir, "%s", optarg );
      break;
    case 'n' :
      nb_children = atoi(optarg);
      break;  
    case 'r' :
      wanted_rank = atoi(optarg);
      wanted_index = -1;
      break;  
    case 'i' :
      wanted_index = atoi(optarg);
      wanted_rank = -1;
      break;  
    case 'd' :
      details=true;
      break;  
    }
  }
  
  if ( env_file_name == NULL || pop_file_name == NULL )
  {
    printf( "%s: error: You must provide a generation number.\n", argv[0] );
    exit( EXIT_FAILURE );
  }
  
  analysis_type type = ROBUSTNESS;

  population_statistics* population_statistics_compute = new population_statistics(type, nb_children, output_dir, wanted_rank, wanted_index, details);
  
  // Load simulation  
  #ifndef __NO_X
    ae_exp_manager* exp_manager = new ae_exp_manager_X11();
  #else
    ae_exp_manager* exp_manager = new ae_exp_manager();
  #endif
  exp_manager->load( num_gener, exp_setup_file_name, out_prof_file_name, env_file_name, pop_file_name, sp_struct_file_name, true );
  
  population_statistics_compute->compute_population_f_nu(exp_manager);

  delete exp_manager;
  delete population_statistics_compute;
  delete [] output_dir;
  delete [] pop_file_name;
  delete [] env_file_name;
  delete [] exp_setup_file_name;
  delete [] out_prof_file_name;
  
  // TODO: update
  // rename dstory.bak.gz
  //rename("kept_dstory.bak.gz","dstory.bak.gz");
  //if(generation%param::back_step == 0) delete po

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
#ifdef __REGUL
  printf( "Usage : rrobustness -h\n");
  printf( "or :    rrobustness -e end_gener [-o output_dir] [-n children_nb] [-r rank | -i index] [-d boolean]\n");
#else
  printf( "Usage : robustness -h\n");
  printf( "or :    robustness -e end_gener [-o output_dir] [-n children_nb] [-r rank | -i index] [-d boolean]\n" );
#endif
  printf( "\n" ); 
  printf( "This program takes individuals and computes f_nu\n" );
  printf( "\n" ); 
  printf( "WARNING: This program should not be used for simulations run with lateral\n" ); 
  printf( "transfer. When an individual has more than one parent, the notion of lineage\n" ); 
  printf( "used here is not relevant.\n" );
  printf( "\n" );  
  printf( "\t-h or --help    : Display this help.\n" );
  printf( "\n" ); 
  printf( "\t-o output_dir or --output output_dir : \n" );
  printf( "\t                  Save the statistics inside output_dir.\n" );
  printf( "\n" ); 
  printf( "\t-n children_nb or --nb-children children_nb : \n" );
  printf( "\t                  Use children_nb to compute Fv.\n" );
  printf( "\n" ); 
  printf( "\t-i index or --index index : \n" );
  printf( "\t                  Get individual with index index (default -1 -> all individuals).\n" );
  printf( "\n" ); 
  printf( "\t-r rank or --rank rank : \n" );
  printf( "\t                  Get individual with rank rank (default -1 -> all individuals).\n" );
  printf( "\n" );
  printf( "\t-e end_gener or --end end_gener : \n" );
  printf( "\t                  Compute f_nu at end_gener\n" );
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

