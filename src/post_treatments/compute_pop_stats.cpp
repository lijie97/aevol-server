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

#include <ae_macros.h>
#include <ae_utils.h>
#include <population_statistics.h>
//#include <ae_common.h>
#ifndef __NO_X
  #include <ae_exp_manager_X11.h>
#else
  #include <ae_exp_manager.h>
#endif
//#include <ae_param_loader.h>


// =======================================================================
//                       Secondary Functions
// =======================================================================



// TODO: update this function...
// reconstruct final individual from backup and lineage
// ae_individual * get_final_individual_using_dstory();     

void print_help( void );

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
  char*   output_dir        = NULL;
  int     nb_children       = 1000;
  int     backup_step       = 0;
  int     generation_number = -1;

  const char * options_list = "he:o:n:"; 
  static struct option long_options_list[] = {
  	{"help",      no_argument,        NULL, 'h'},
    {"end",       required_argument,  NULL, 'e' },
    {"output", 1, NULL, 'o'},
    {"nb-children", 1, NULL, 'n'},
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
        
        generation_number = atol( optarg );
        
        break;
      }
    case 'o' : 
      output_dir = new char[strlen(optarg) + 1];
      sprintf( output_dir, "%s", optarg );
      break;
    case 'n' :
      nb_children = atoi(optarg);
      break;  
    }
  }
  
  if ( generation_number == -1 )
  {
    printf( "%s: error: You must provide a generation number.\n", argv[0] );
    exit( EXIT_FAILURE );
  }
  
  analysis_type type = EVOLVABILITY;

  population_statistics* population_statistics_compute = new population_statistics(type, nb_children, output_dir);
  
  
  // Load simulation  
  #ifndef __NO_X
    ae_exp_manager* exp_manager = new ae_exp_manager_X11();
  #else
    ae_exp_manager* exp_manager = new ae_exp_manager();
  #endif
  exp_manager->load( generation_number, false, true );
  
  backup_step = exp_manager->get_backup_step();
  
  
  
  for ( int32_t i = 0 ; i <= generation_number ; i += backup_step )
  {
  	printf("\n\n Generation : %d\n\n", i);
    
    delete exp_manager;
    #ifndef __NO_X
    	exp_manager = new ae_exp_manager_X11();
  	#else
    	exp_manager = new ae_exp_manager();
  	#endif
    exp_manager->load( generation_number, false, true );
    
    population_statistics_compute->compute_population_f_nu(exp_manager);
    population_statistics_compute->compute_evolvability_stats(i);
  }

  delete exp_manager;
  //delete [] pop_file_name;
  //delete [] env_file_name;
  delete population_statistics_compute;
  delete [] output_dir;
  
  
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
  printf( "*    Population statistics computation post-treatment program          * \n" );
  printf( "*                                                                      * \n" );
  printf( "************************************************************************ \n" );
  printf( "\n\n" ); 
  printf( "This program is Free Software. No Warranty.\n" );
  printf( "Copyright (C) 2009  LIRIS.\n" );
  printf( "\n" ); 
#ifdef __REGUL
  printf( "Usage : rcomputate_pop_stats -h\n");
  printf( "or :    rcomputate_pop_stats -e end_gener [-o output_dir] [-n children_nb]\n" );
#else
  printf( "Usage : computate_pop_stats -h\n");
  printf( "or :    computate_pop_stats -e end_gener [-o output_dir] [-n children_nb]\n" );
#endif
  printf( "\n" ); 
  printf( "This program computes some population statistics at each available backup until end_gener\n" );
  printf( "and save this statistics in files inside output_dir. The children_nb is used to compute Fv.\n");
  printf( "The population backup files and environment backup file are required.\n" );
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
  printf( "\t-e end_gener or --end end_gener : \n" );
  printf( "\t                  Retrieve the lineage of the individual of end_gener \n" );
  printf( "\n" );

}

