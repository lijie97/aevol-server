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
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <sys/stat.h>  // for the permission symbols used with mkdir



// =================================================================
//                            Project Files
// =================================================================
#include <ae_common.h>
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_environment.h>
#include <ae_list.h>
#include <ae_simulation.h>
#include <ae_param_overloader.h>


void print_help( void );





int main( int argc, char* argv[] )
{
  // =================================================================
  //                      Get command-line options
  // =================================================================
  //
  // 1) Initialize command-line option variables with default values
  char* backup_file_name = NULL;
  long int ndiv = 100000;
  ae_param_overloader* param_overloader = new ae_param_overloader();
  char* output_file_name = NULL;
  
  // 2) Define allowed options
  const char * options_list = "hp:n:f:o:";
  static struct option long_options_list[] = {
    { "number", 1, NULL, 'n' },
    { "file", 1, NULL, 'f' },
    { "output", 1, NULL, 'o' },
    { "param",    required_argument, NULL, 'p' },
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
        print_help();
        exit( EXIT_SUCCESS );
        break;
      }
      case 'n' :
      {
        ndiv = strtol( optarg, NULL, 10 );
        break;
      }
      case 'p' :
      {
        if( strcmp( optarg, "") == 0)
        {
          printf( "ERROR : Option -p or --param : missing agument.\n" );
          exit( EXIT_FAILURE );
        }
        
        param_overloader->store_overload( optarg );
        
        break;
      }  
      case 'f' :
      {
        backup_file_name = new char[strlen(optarg) + 1];
        sprintf( backup_file_name, "%s", optarg );
        break;
      }
      case 'o' :
      {
        output_file_name = new char[strlen(optarg) + 1];
        sprintf( output_file_name, "%s", optarg );
        break;
      }
    }
  }
  
  
  // =================================================================
  //                       Read the backup file
  // =================================================================
  
  ae_individual*  best_indiv;
  ae_environment* env;
  int32_t         num_gener;
  
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
    ae_common::sim  = new ae_simulation ( backup_file_name, false, param_overloader );
    best_indiv      = ae_common::sim->get_pop()->get_best();
    env             = ae_common::sim->get_env();
    num_gener       = ae_common::sim->get_num_gener();
    printf("done\n");

  }
  
  delete [] backup_file_name;
  
  
  
  // The constructor of the ae_simulation has read the genomes of the individuals
  // and located their promoters, but has not performed the translation nor the
  // phenotype computation. We must do it now.
  // However, as the individuals in the backups are sorted, we don't need to evaluate
  // all the individuals, only those we are interested in (here only the best one)
      
  best_indiv->evaluate( env );
  
  
  
  // =================================================================
  //                     Open output file
  // =================================================================

  FILE * output = fopen(output_file_name, "w");
  if ( output == NULL )
  {
    fprintf( stderr, "ERROR : Could not create the output file %s\n", output_file_name );
    exit( EXIT_FAILURE );
  }

  // Write the header
  fprintf( output, "# #################################################################\n" );
  fprintf( output, "#              Mutations produced by mutagenesissecretion\n" );
  fprintf( output, "# #################################################################\n" );
  fprintf( output, "#\n" );
  fprintf( output, "impact on metabolism ; impact on secretion\n" );
  fprintf( output, "#\n" );
  


  
  // =================================================================
  //                     Proceed to mutagenesis
  // =================================================================
  // Set some data
  ae_individual*  initial_indiv = best_indiv;

  double          initial_metabolic_error = initial_indiv->get_dist_to_target_by_feature( METABOLISM );
  double          initial_secretion_error = initial_indiv->get_dist_to_target_by_feature( SECRETION );
  
  double final_metabolic_error      = 0.0;
  double impact_on_metabolic_error  = 0.0;
  
  double final_secretion_error      = 0.0;
  double impact_on_secretion_error  = 0.0;
  
  ae_individual*  indiv       = NULL;
  ae_vis_a_vis*   alignment_1 = NULL;
  ae_vis_a_vis*   alignment_2 = NULL;
  int32_t         nb_pairs;
  
  int16_t i;
  bool    rear_done;

  // Perform ndiv reproductions with mutations
  
  for ( i = 0; i < ndiv; i++ )
  {
    if ( i % 100 == 0 )
    {
      printf( "*" );
      fflush(stdout);
    }
    
    indiv = ae_common::sim->get_pop()->do_replication( best_indiv, -1 );
    
    indiv->reevaluate(env);
    final_metabolic_error     = indiv->get_dist_to_target_by_feature( METABOLISM );
    impact_on_metabolic_error = final_metabolic_error - initial_metabolic_error;
    
    final_secretion_error     = indiv->get_dist_to_target_by_feature( SECRETION );
    impact_on_secretion_error = final_secretion_error - initial_secretion_error;

    fprintf( output, "%+.15f %+.15f \n",impact_on_metabolic_error, impact_on_secretion_error );
    
    delete indiv;

  }
  

  delete ae_common::sim;

  fclose(output);
  delete [] output_file_name;
  
  return EXIT_SUCCESS;
}




void print_help( void ) 
{
}




