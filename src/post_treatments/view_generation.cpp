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
#include <ae_experiment_X11.h>





void print_help( char* prog_name );
void open_protein_list_file( void );
void write_protein_list_file( ae_population * pop );
FILE* protein_list_file     = NULL;
int32_t end_gener;

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
      ae_common::sim = ae_common::sim_display = new ae_experiment_X11();
      ae_common::sim->load_backup( backup_file_name, false, NULL );
      printf("done\n");
    }



  // =================================================================
  //                      Evaluate the individuals
  // =================================================================


  // The constructor of the ae_experiment has read the genomes of the individuals
  // and located their promoters, but has not performed the translation nor the
  // phenotype computation. We must do it now.
  
  printf( "Computing phenotypes... " );
  fflush( stdout );

  ae_list_node  * indiv_node = ae_common::pop->get_indivs()->get_first();
  ae_individual * indiv      = NULL;
  ae_population * pop = ae_common::pop;
  
  end_gener = 100;
  int32_t size = 100;
  char size_file_name[60];
  FILE* size_file     = NULL;
  char codant_file_name[60];
  FILE* codant_file     = NULL;
  
  end_gener = atoi((const char*)&backup_file_name[11]);
  open_protein_list_file();
  snprintf( size_file_name, 60, "genome-size-e%06"PRId32, end_gener );
  size_file = fopen( size_file_name, "w" );
  
  snprintf( codant_file_name, 60, "codant-size-e%06"PRId32, end_gener );
  codant_file = fopen( codant_file_name, "w" );


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
    
    size = indiv->get_total_genome_size();
    fprintf( size_file, "%"PRId32"\n", size );
    
    size = indiv->get_overall_size_functional_genes();
    fprintf( codant_file, "%"PRId32"\n", size );

    indiv_node = indiv_node->get_next();
  }
  write_protein_list_file(pop);
  fclose(size_file);
  fclose(codant_file);
  fclose(protein_list_file);
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

void open_protein_list_file( void )
{
  // Open file
  char prot_list_file_name[60];

  snprintf( prot_list_file_name, 60, "view_generation_best-e%06"PRId32"_protein_list.out", end_gener );
  protein_list_file = fopen( prot_list_file_name, "w" );
  
  // Write headers
  fprintf( protein_list_file, "# Proteins of the best individual of generation %"PRId32"\n", end_gener );
  fprintf( protein_list_file, "# 1:  Strand on which is the protein (0 = chromosome) \n" );
  fprintf( protein_list_file, "# 2:  Shine-Dalgarno position \n" );
  fprintf( protein_list_file, "# 3:  Mean \n" );
  fprintf( protein_list_file, "# 4:  Half-width \n" );
  fprintf( protein_list_file, "# 5:  Height \n" );
  fprintf( protein_list_file, "# 6:  Length (number of codons) \n" );
  fprintf( protein_list_file, "# 7:  Concentration \n" );
  fprintf( protein_list_file, "# 8:  Whether the protein is functional (0 = False, 1 = True)\n" );
  fprintf( protein_list_file, "# \n" );
}

void write_protein_list_file( ae_population * pop )
{
  ae_individual* indiv1     = pop->get_best();
  // For each protein : Strand, Shine-Dal pos, mean, width, height, length, concentration, is-functionnal
  ae_list_node * protnode  = indiv1->get_protein_list()->get_first();
  ae_protein *  prot       = NULL;
  while ( protnode != NULL )
  {
    prot = (ae_protein *) protnode->get_obj();
    fprintf( protein_list_file, "%"PRId32" %d %f %f %f %d %f %d \n", prot->get_strand(), 
									  prot->get_shine_dal_pos(), 
									  prot->get_mean() ,
									  prot->get_width(),
									  prot->get_height(),
									  prot->get_length(),
									  prot->get_concentration(),
									  prot->get_is_functional() );
    protnode = protnode->get_next();
  }
  
  fprintf( protein_list_file, "\n" );
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
