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
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>




// =================================================================
//                            Project Files
// =================================================================
#include <ae_utils.h>
#include <ae_simulation.h>
#include <ae_individual.h>
#include <ae_genetic_unit.h>
#include <ae_list.h>
#include <ae_tree.h>
#include <ae_replication_report.h>
#include <ae_dna_replic_report.h>
#include <ae_mutation.h>
#include <ae_param_loader.h>




void print_help( void );


int main(int argc, char** argv) 
{
  // The output file (lineage.ae or lineage.rae) contains the following informations:
  //
  // - common data                                                (ae_common::write_to_backup) 
  // - begin gener                                                (int32_t)
  // - end gener                                                  (int32_t)
  // - final individual index                                     (int32_t)
  // - initial genome size                                        (int32_t)
  // - initial ancestor (nb genetic units + sequences)            (ae_individual::write_to_backup)
  // - replication report of ancestor at generation begin_gener+1 (ae_replic_report::write_to_backup)    
  // - replication report of ancestor at generation begin_gener+2 (ae_replic_report::write_to_backup)    
  // - replication report of ancestor at generation begin_gener+3 (ae_replic_report::write_to_backup)    
  // - ...
  // - replication report of ancestor at generation end_gener     (ae_replic_report::write_to_backup)    


  //  TODO : trouver la valeur du parametre ae_common::tree_mode
  // OK si NORMAL, sortir si LIGHT


  printf("\n  WARNING : Parameters' change in the middle of a generation is not managed.n");
 

  // =====================
  //  Parse command line
  // =====================

  // Default values
  bool    check_genome    = true;
  bool    verbose         = false;
  int32_t begin_gener = 0; 
  int32_t end_gener   = 100;
  int32_t final_indiv = -1; 

  char* backup_file_name = NULL;
  char tree_file_name[50];

  const char * short_options = "hvnb:i:f:"; 
  static struct option long_options[] = {
    {"help",    no_argument,       NULL, 'h'},
    {"verbose", no_argument,       NULL, 'v'},
    {"nocheck", no_argument,       NULL, 'n'},
    {"begin",   required_argument, NULL, 'b'},
    {"indiv",   required_argument, NULL, 'i'},
    {"file",    required_argument, NULL, 'f'},
    {0, 0, 0, 0}
  };

  int option;
  while( (option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1 ) 
  {
    switch( option ) 
    {
      case 'h' : print_help(); exit(EXIT_SUCCESS); break;
      case 'v' : verbose = true;                   break;
      case 'n' : check_genome = false;             break;
      case 'b' : begin_gener  = atol(optarg);      break;
      case 'i' : final_indiv  = atol(optarg);      break;
      case 'f' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "ERROR : Option -f or --file : missing argument.\n" );
          exit( EXIT_FAILURE );
        }
       
        backup_file_name = new char[strlen(optarg) + 1];
        sprintf( backup_file_name, "%s", optarg );

        // Updating the value of end_gener with the common options just charged
        end_gener = atoi((const char*)&backup_file_name[11]);
       
        break;      
      }
    }
  }

  // Check mandatory arguments
  if ( backup_file_name == NULL )
  {
    print_help();
    exit( EXIT_FAILURE );
  }


  // =========================
  //  Load the last tree file
  // =========================

  if ( verbose )
  {
    printf( "\n\n" );
    printf( "====================================\n" );
    printf( " Loading the last tree file ... " );
    fflush( stdout );
  }


  // Example for ae_common::tree_step == 100 :
  //
  // tree_000100.ae ==>  generations   1 to 100.
  // tree_000200.ae ==>  generations 101 to 200.
  // tree_000300.ae ==>  generations 201 to 300.
  // etc.
  // 
  // Thus, the informations for generation end_gener are located 
  // in the file called (end_gener/ae_common::tree_step + 1) * ae_common::tree_step,
  // except if end_gener%ae_common::tree_step==0.

  #ifdef __REGUL
    sprintf( tree_file_name,"tree/tree_%06"PRId32".rae", end_gener ); 
  #else
    sprintf( tree_file_name,"tree/tree_%06"PRId32".ae", end_gener ); 
  #endif

  ae_tree * tree = new ae_tree( backup_file_name, tree_file_name );

  if ( verbose )
  {
    printf("OK\n");
    printf("====================================\n");
  }


  // ========================================
  //  Find the index of the final individual
  // ========================================


  if ( final_indiv == -1 )
  {
    // No index was given in the command line.
    // By default, we construct the lineage of the best individual.
    // Since individuals are indexed according to their ranks in the
    // population, we know that the best individual has the index
    // N-1, if N is the size of the population.
    final_indiv = tree->get_nb_indivs( end_gener ) - 1;
    if ( verbose ) printf( "The final individual has the index %"PRId32"\n", final_indiv );
  }


  // =======================
  //  Open the output file
  // =======================

  char output_file_name[101];
#ifdef __REGUL
  snprintf(output_file_name, 100, "lineage-b%06"PRId32"-e%06"PRId32"-i%"PRId32".rae", begin_gener, end_gener, final_indiv);
#else 
  snprintf(output_file_name, 100, "lineage-b%06"PRId32"-e%06"PRId32"-i%"PRId32".ae", begin_gener, end_gener, final_indiv);
#endif

  gzFile * lineage_file = (gzFile*) gzopen(output_file_name, "w");
  if(lineage_file == NULL)
  {
    fprintf(stderr, "File %s could not be created, exiting.\n", output_file_name);
    fprintf(stderr, "Please check your permissions in this directory.\n");
    exit(EXIT_FAILURE);
  }




  // ===================================================
  //  Retrieve the replication reports of the ancestors
  // ===================================================

  if ( verbose )
  {
    printf("\n\n\n");
    printf("================================================================================\n");
    printf(" Parsing tree files to retrieve the ancestors and their replication reports... \n");
    printf("================================================================================\n");
  }

  // Indices and replication reports of the individuals in the lineage
  int32_t *                 indices = new int32_t[end_gener - begin_gener + 1];
  ae_replication_report **  reports = new ae_replication_report*[end_gener - begin_gener];
  // NB: we do not need the report of the ancestor at generation begin_gener
  // (it might be the generation 0, for which we have no reports)
  // reports[0] = how ancestor at generation begin_gener + 1 was created 
  // reports[i] = how ancestor at generation begin_gener + i + 1 was created
  // reports[end_gener - begin_gener - 1] = how the final individual was created
  
  // Copy the replication report of the individual whose lineage we are computing
  indices[end_gener - begin_gener] = final_indiv;

  
  int32_t i, num_gener;
 
  for ( i = end_gener - begin_gener - 1 ; i >= 0 ; i-- )
  {
    // We want to fill reports[i], that is to say, how the ancestor
    // at generation begin_gener + i + 1  was created

    num_gener = begin_gener + i + 1;
    if ( verbose ) printf( "Getting the replication report for the ancestor at generation %"PRId32"\n", num_gener );

    
    if ( utils::mod( num_gener, ae_common::tree_step ) == 0 ) 
    {
      // Change the tree file 
      delete tree;

      #ifdef __REGUL
        sprintf( tree_file_name,"tree/tree_%06"PRId32".rae", num_gener ); 
      #else
        sprintf( tree_file_name,"tree/tree_%06"PRId32".ae",  num_gener ); 
      #endif
      
      tree = new ae_tree( backup_file_name, tree_file_name );
    }
    

    // Copy the replication report of the ancestor
    reports[i] = new ae_replication_report( *(tree->get_report(num_gener, indices[i + 1])) );

    // Retreive the index of the ancestor from the report
    indices[i] = reports[i]->get_parent_index();
  }
  
  for ( i = 0 ; i < end_gener - begin_gener ; i++ )
  {
    printf( "  indices[%"PRId32"] : %"PRId32"\n", i, indices[i] );
  }
  

  if ( verbose )  printf("OK\n");


  // =============================================================================
  //  Get the intial genome from the backup file, and write it in the output file
  // =============================================================================

  if ( verbose )
  {
    printf("\n\n\n");
    printf("=============================================== \n");
    printf(" Getting the initial genome sequence... ");
    fflush(NULL);
  }

  char genomes_file_name[50];
#ifdef __REGUL
  sprintf( genomes_file_name,"backup/gen_%06"PRId32".rae", begin_gener );
#else
  sprintf( genomes_file_name,"backup/gen_%06"PRId32".ae",  begin_gener );
#endif


  ae_simulation * sim = new ae_simulation ( genomes_file_name, false );

  // Copy the initial ancestor
  // NB : The list of individuals is sorted according to the index
  ae_individual * initial_ancestor_tmp  = (ae_individual *)sim->get_pop()->get_indivs()->get_object( indices[0] );
  ae_individual * initial_ancestor      = new ae_individual( *initial_ancestor_tmp );
  
  ae_common::write_to_backup( lineage_file );
  
  delete sim; // we don't need the other individuals
    
    
    
  gzwrite( lineage_file, &begin_gener,  sizeof(begin_gener) );
  gzwrite( lineage_file, &end_gener,    sizeof(end_gener)   );
  gzwrite( lineage_file, &final_indiv,  sizeof(final_indiv) );
  
  initial_ancestor->write_to_backup( lineage_file );


  if ( verbose )
  {
    printf("OK\n");      
    printf("=============================================== \n");              
  }


  // ===============================================================================
  //  Write the replication reports of the successive ancestors in the output file
  //  (and, optionally, check that the rebuilt genome is correct each time a backup
  //  is available)
  // ===============================================================================

  if ( verbose )
  {
    printf("\n\n\n");
    printf("============================================================ \n");              
    printf(" Write the replication reports in the output file... \n");
    printf("============================================================ \n"); 
  }


  ae_list_node * report_node = NULL;
  ae_dna_replic_report * rep = NULL;
  ae_list_node * mut_node = NULL;
  ae_mutation *  mut = NULL;
  ae_list_node*   gen_unit_node = NULL;
  ae_genetic_unit *  unit = NULL;

  ae_individual *   stored_indiv          = NULL;
  ae_list_node *    stored_gen_unit_node  = NULL;
  ae_genetic_unit * stored_gen_unit       = NULL;

  

  // NB: I must keep the genome encapsulated inside an ae_individual, because
  // replaying the mutations has side effects on the list of promoters,
  // which is stored in the individual

  for ( i = 0 ; i < end_gener - begin_gener ; i++ )
  {
    num_gener = begin_gener + i + 1; // where are we in time...
    
    if ( verbose )
    {
      printf( "Writing the replication report for generation %"PRId32" (built from indiv %"PRId32" at generation %"PRId32")\n",
              num_gener, indices[i], num_gener-1 );
    }
    reports[i]->write_to_backup( lineage_file );
    if ( verbose ) printf( " OK\n" );
    

    if( check_genome )
    {
      if ( utils::mod( num_gener, ae_common::backup_step ) == 0 )
      {
#ifdef __REGUL
        sprintf( genomes_file_name,"backup/gen_%06"PRId32".rae", num_gener );
#else
        sprintf( genomes_file_name,"backup/gen_%06"PRId32".ae", num_gener );
#endif
        if ( verbose )
        {
          printf( "Loading the data stored in backup file %s\n",genomes_file_name );
        }
        
        // Copy the ancestor from the backup
        sim = new ae_simulation ( genomes_file_name, false );
        ae_individual * stored_indiv_tmp = (ae_individual *) sim->get_pop()->get_indivs()->get_object( indices[i+1] );
        stored_indiv = new ae_individual( *stored_indiv_tmp );
        
        delete sim;
      }
          
  
      // Warning: this portion of code won't work if the number of units changes
      // during the evolution
      
      report_node   = reports[i]->get_dna_replic_reports()->get_first();
      gen_unit_node = initial_ancestor->get_genetic_unit_list()->get_first();
      if ( utils::mod( num_gener, ae_common::backup_step ) == 0 )
      {
        stored_gen_unit_node = stored_indiv->get_genetic_unit_list()->get_first();
      }
      
      while ( report_node != NULL )
      {
        assert( gen_unit_node != NULL );
        if ( utils::mod( num_gener, ae_common::backup_step ) == 0 )
        {
          assert( stored_gen_unit_node != NULL );
        }
        
        rep = (ae_dna_replic_report *) report_node->get_obj();
        unit = (ae_genetic_unit *) gen_unit_node->get_obj();
        if ( utils::mod( num_gener, ae_common::backup_step ) == 0 )
        {
          stored_gen_unit = (ae_genetic_unit *) stored_gen_unit_node->get_obj();
        }
        
        mut_node = (rep->get_rearrangements())->get_first();              
        while ( mut_node != NULL )
        {
          mut = (ae_mutation *) mut_node->get_obj();
          (unit->get_dna())->undergo_this_mutation(mut);
          mut_node = mut_node->get_next();
        }
        
        mut_node = (rep->get_mutations())->get_first();              
        while (mut_node != NULL)
        {
          mut = (ae_mutation *) mut_node->get_obj();
          (unit->get_dna())->undergo_this_mutation(mut); 
          mut_node = mut_node->get_next();
        }

        if ( utils::mod( num_gener, ae_common::backup_step ) == 0 )
        {
          if ( verbose )
          {
            printf( "Checking the sequence of the unit..." );
            fflush( stdout );
          }

          char * str1 = new char[(unit->get_dna())->get_length() + 1];
          memcpy(str1, (unit->get_dna())->get_data(), \
                 (unit->get_dna())->get_length()*sizeof(char));
          str1[(unit->get_dna())->get_length()] = '\0';
          
          char * str2 = new char[(stored_gen_unit->get_dna())->get_length() + 1];
          memcpy(str2, (stored_gen_unit->get_dna())->get_data(), (stored_gen_unit->get_dna())->get_length()*sizeof(char));
          str2[(stored_gen_unit->get_dna())->get_length()] = '\0';
          
          if ( strncmp( str1, str2, (stored_gen_unit->get_dna())->get_length() ) == 0 )
          {
            if ( verbose ) printf( " OK\n" );
          }
          else
          {
            if ( verbose ) printf( " ERROR !\n" );
            fprintf( stderr, "Error: the rebuilt unit is not the same as \n");
            fprintf( stderr, "the one stored in backup file %s\n", genomes_file_name);
            fprintf( stderr, "Rebuilt unit : %"PRId32" bp\n %s\n", (int32_t)strlen(str1), str1 );
            fprintf( stderr, "Stored unit  : %"PRId32" bp\n %s\n", (int32_t)strlen(str2), str2 );
            delete [] str1;
            delete [] str2;
            gzclose(lineage_file);
            delete initial_ancestor;
            delete stored_indiv;
            delete [] reports;
            exit(EXIT_FAILURE);
          }
          
          delete [] str1;
          delete [] str2;
        }

    
        report_node = report_node->get_next();
        gen_unit_node = gen_unit_node->get_next();
        if ( utils::mod(num_gener, ae_common::backup_step) == 0 )
        {
          stored_gen_unit_node = stored_gen_unit_node->get_next();
        }
      }
      
      assert(gen_unit_node == NULL);
      if ( utils::mod( num_gener, ae_common::backup_step ) == 0 )
      {
        assert(stored_gen_unit_node == NULL);
        delete stored_indiv;
      }

    }
    //      initial_ancestor->evaluate();

  }

  //  ae_common::clean();

  gzclose(lineage_file);
  delete initial_ancestor;
  delete [] reports;


  exit(EXIT_SUCCESS);
  
}



void print_help( void )
{
  // default values :
  // begin_gener = 0 
  // indiv  = best individual at generation end_gener
  
  // there must be a genome backup file for begin_gener 

  // not relevant if crossover

  printf( "\n" ); 
  printf( "*********************** aevol - Artificial Evolution ******************* \n" );
  printf( "*                                                                      * \n" );
  printf( "*                      Lineage post-treatment program                  * \n" );
  printf( "*                                                                      * \n" );
  printf( "************************************************************************ \n" );
  printf( "\n\n" ); 
  printf( "This program is Free Software. No Warranty.\n" );
  printf( "Copyright (C) 2009  LIRIS.\n" );
  printf( "\n" ); 
#ifdef __REGUL
  printf( "Usage : rlineage -h\n");
  printf( "or :    rlineage [-vn] [-i index] [-b gener1] -f backup_file \n" );
#else
  printf( "Usage : lineage -h\n");
  printf( "or :    lineage [-vn] [-i index] [-b gener1] -f backup_file \n" );
#endif
  printf( "\n" ); 
#ifdef __REGUL
  printf( "This program retrieves the ancestral lineage of an individual and writes \n" );
  printf( "it in an output file called lineage.rae. Specifically, it retrieves the \n");
  printf( "lineage of the individual of backup_file whose index is index, going \n" );
  printf( "back in time up to gener1. This program requires at least one genome backup\n" );
  printf( "file (for the generation gener1) and all tree files for generations gener1\n" );
  printf( "to that corresponding to backup_file.\n" );
#else
  printf( "This program retrieves the ancestral lineage of an individual and writes \n" );
  printf( "it in an output file called lineage.rae. Specifically, it retrieves the \n");
  printf( "lineage of the individual of backup_file whose index is index, going \n" );
  printf( "back in time up to gener1. This program requires at least one genome backup\n" );
  printf( "file (for the generation gener1) and all tree files for generations gener1\n" );
  printf( "to that corresponding to backup_file.\n" );
#endif
  printf( "\n" ); 
  printf( "WARNING: This program should not be used for simulations run with lateral\n" ); 
  printf( "transfer. When an individual has more than one parent, the notion of lineage\n" ); 
  printf( "used here is not relevant.\n" );
  printf( "\n" );  
  printf( "\t-h or --help    : Display this screen.\n" );
  printf( "\n" ); 
  printf( "\t-v or --verbose : Be verbose, listing generations as they are \n" );
  printf( "\t                  treated.\n" );
  printf( "\n" ); 
  printf( "\t-n or --nocheck : Disable genome sequence checking. Makes the \n"); 
  printf( "\t                  program faster, but it is not recommended. \n");
  printf( "\t                  It is better to let the program check that \n");
  printf( "\t                  when we rebuild the genomes of the ancestors\n");
  printf( "\t                  from the tree files, we get the same sequences\n");
  printf( "\t                  as those stored in the backup files.\n" );
  printf( "\n" ); 
  printf( "\t-i index or --indiv index : \n" );
  printf( "\t                  Retrieve the lineage of the individual whose\n" );
  printf( "\t                  index is index. The index must be comprised \n" );
  printf( "\t                  between 0 and N-1, with N the size of the \n" );
  printf( "\t                  population at generation gener2. If not \n" );
  printf( "\t                  specified, the program retrieves the lineage \n" );
  printf( "\t                  of the best individual of generation gener2.\n");
  printf( "\n" ); 
  printf( "\t-b gener1 or --begin gener1 : \n" );
  printf( "\t                  Retrieve the lineage up to generation gener1.\n" );
  printf( "\t                  There must be a genome backup file for this\n" );
  printf( "\t                  generation. If not specified, the program \n" );
  printf( "\t                  retrieves the lineage up to generation 0.\n"); 
  printf( "\n" ); 
  printf( "\t-f backup_file or --file backup_file : \n" );
  printf( "\t                  Retrieve the lineage of an individual that \n" );
  printf( "\t                  belongs to the generation of the backup_file.\n" );
  printf( "\n" ); 

}
