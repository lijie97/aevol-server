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

//debug
//#include <ae_gaussian.h>


void print_help( void );


int main(int argc, char** argv) 
{
  // The input file (lineage.ae or lineage.rae) must contain the following informations:
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



 
  // =====================
  //  Parse command line
  // =====================

  // Default values
  char*   lineage_file_name   = NULL;
  bool    verbose             = false;
  bool    check               = true;

  const char * short_options = "hvnf:"; 
  static struct option long_options[] =
  {
    {"help",    no_argument,       NULL, 'h'},
    {"verbose", no_argument,       NULL, 'v'},
    {"nocheck", no_argument,       NULL, 'n'},
    {"file",    required_argument, NULL, 'f'},
    {0, 0, 0, 0}
  };

  int option;
  while( (option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1 ) 
  {
    switch( option ) 
    {
      case 'h' : print_help(); exit(EXIT_SUCCESS);  break;
      case 'v' : verbose = true;                    break;
      case 'n' : check = false;                     break;
      case 'f' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          fprintf( stderr, "ERROR : Option -f or --file : missing argument.\n" );
          exit( EXIT_FAILURE );
        }
        
        lineage_file_name = new char[strlen(optarg) + 1];
        sprintf( lineage_file_name, "%s", optarg );
        break;      
      }
      default :
      {
        fprintf( stderr, "ERROR : Unknown option, check your syntax.\n" );
        print_help();
        exit( EXIT_FAILURE );
      }
    }
  }
  
  if (lineage_file_name == NULL )
  {
    fprintf( stderr, "ERROR : Option -f or --file missing. \n" );
    exit( EXIT_FAILURE );
  }


  // =======================
  //  Open the lineage file
  // =======================

  gzFile * lineage_file = (gzFile*) gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL)
  {
    fprintf( stderr, "ERROR : Could not read the lineage file %s\n", lineage_file_name );
    exit( EXIT_FAILURE );
  }

  ae_common::read_from_backup( lineage_file );

  int32_t begin_gener, end_gener, final_index;
  gzread(lineage_file, &begin_gener, sizeof(begin_gener));
  gzread(lineage_file, &end_gener, sizeof(end_gener));
  gzread(lineage_file, &final_index, sizeof(final_index) );

  if ( verbose )
  {
    printf("\n\n");
    printf("================================================================================\n");
    printf(" Statistics of the ancestors of indiv. #%"PRId32" (t=%"PRId32" to %"PRId32")\n", final_index, begin_gener, end_gener);
    printf("================================================================================\n");
  }


  // =======================
  //  Open the output file
  // =======================

  char output_file_name[60];
  snprintf(output_file_name, 60, "fixedmut-b%06"PRId32"-e%06"PRId32"-i%"PRId32".out", begin_gener, end_gener, final_index);

  FILE * output = fopen(output_file_name, "w");
  if (output == NULL)
  {
    fprintf( stderr, "ERROR : Could not create the output file %s\n", output_file_name );
    exit( EXIT_FAILURE );
  }


  // Write the header
  fprintf( output, "# Mutations in the lineage of indiv #%"PRId32" bw t=%"PRId32" and t=%"PRId32"\n", final_index, begin_gener, end_gener );
  fprintf( output, "# 1:  Generation (mut. occurred when producing the indiv. of this generation) \n" );
  fprintf( output, "# 2:  Genetic unit which underwent the mutation (0 = chromosome) \n" );
  fprintf( output, "# 3:  Mutation type (0: switch, 1: smallins, 2: smalldel, 3:duplication, 4: deletion, 5:transloc, 6:inversion) \n" );
  fprintf( output, "# 4:  pos[0] (position for the small events, begin_segment for the rearrangements) \n" );
  fprintf( output, "# 5:  pos[1] (-1 for the small events, end_segment for the rearrangements) \n" );
  fprintf( output, "# 6:  pos[2] (reinsertion point for duplic., cutting point in segment for transloc., -1 for other events)\n" );
  fprintf( output, "# 7:  pos[3] (reinsertion point for transloc., -1 for other events)\n" );
  fprintf( output, "# 8:  For transloc., was the segment inverted (0/1)? (-1 for other events)\n" );
  fprintf( output, "# 9:  Length of the segment\n" );
  fprintf( output, "# 10: Length of the genetic unit before the event\n" );
  fprintf( output, "# 11: Impact of the mutation on the metabolic error (negative value = smaller gap after = beneficial mutation) \n" );

  // =========================
  //  Prepare the environment
  // =========================

  if ( verbose )
  {
    printf("Preparing the environment... ");  
    fflush(NULL);
  }


  // Initialize the environment according the ae_common::... values
  // This instruction also creates the random generator and sets its
  // seed, according to the ae_common::env_seed value.

  ae_environment * env = new ae_environment();

  int32_t t = 0;
  for (t = 0; t < begin_gener; t++)
  {
    env->apply_variation();
  }


  if ( verbose ) printf("OK\n");

  char backup_file_name[50];
  if ( check )
  {
    // check that the environment is now identical to the one stored
    // in the backup file of generation begin_gener
    
#ifdef __REGUL
    sprintf(backup_file_name,"backup/gen_%06"PRId32".rae", begin_gener);
#else
    sprintf(backup_file_name,"backup/gen_%06"PRId32".ae",  begin_gener);
#endif
    if ( verbose )
    {
      printf("Comparing the environment with the one in %s... ", backup_file_name);  
      fflush(NULL);
    }
    
    ae_simulation * sim = new ae_simulation ( backup_file_name, false );
    
    ae_environment * stored_env = sim->get_env();
    
    // debug
    //       ae_list_node * gaussnode = env->get_gaussians()->get_first();
    //       ae_list_node * storedgaussnode = stored_env->get_gaussians()->get_first();
    //       ae_gaussian * gauss = NULL;
    //       ae_gaussian * storedgauss = NULL;
    //       while(gaussnode != NULL)
    //       {
    //         assert(storedgaussnode != NULL);
    //         gauss = (ae_gaussian *) gaussnode->get_obj();
    //         storedgauss = (ae_gaussian *) storedgaussnode->get_obj();
            
    //         printf("gaussianA : m %f, w %f, h %f\n", gauss->get_mean(), gauss->get_width(), gauss->get_height());
    //         printf("gaussianB : m %f, w %f, h %f\n", storedgauss->get_mean(), storedgauss->get_width(), storedgauss->get_height());
            
    //         gaussnode = gaussnode->get_next();
    //         storedgaussnode = storedgaussnode->get_next();
    //       }
    //       assert(storedgaussnode == NULL);
    // end debug

    if ( ! (env->is_identical_to(stored_env)) )
    {
      fprintf(stderr, "ERROR: The replayed environment is not the same\n");
      fprintf(stderr, "       as the one in %s\n", backup_file_name);
      exit(EXIT_FAILURE);
    }
    
    if ( verbose ) printf("OK\n");
    delete sim; 
  }

 


  // ==============================
  //  Prepare the initial ancestor
  // ==============================

  ae_individual * indiv = new ae_individual(lineage_file);
  indiv->evaluate( env );
  indiv->compute_statistical_data();
  if ( verbose )
    {
      printf("Initial fitness     = %f\n", indiv->get_fitness());
      printf("Initial genome size = %"PRId32"\n", indiv->get_total_genome_size());
    }



  // ===============================================================================
  //  Replay the mutation to get the successive ancestors and analyze them
  //  (and, optionally, check that the rebuilt envir and genome are correct each 
  //   time a backup is available)
  // ===============================================================================


  ae_replication_report * rep = NULL;
  ae_list_node * dnarepnode = NULL;
  ae_dna_replic_report * dnarep = NULL;

  ae_list_node * mnode = NULL;
  ae_mutation *  mut = NULL;

  ae_list_node*   unitnode = NULL;
  ae_genetic_unit *  unit = NULL;

  ae_individual * stored_indiv = NULL;
  ae_list_node*   storedunitnode = NULL;
  ae_genetic_unit *  storedunit = NULL;

  int32_t i, index, genetic_unit_number, unitlen_before, seglen;
  double metabolic_error_before, metabolic_error_after, impact_on_metabolic_error;
  char mut_descr_string[80];

  for ( i = 0; i < end_gener - begin_gener; i++ )
  {
    t = begin_gener + i + 1;  // where are we in time...

    rep = new ae_replication_report(lineage_file);
    index = rep->get_index(); // who are we building...
    indiv->set_replication_report(rep);


    if ( verbose ) printf("Rebuilding ancestor at generation %"PRId32" (index %"PRId32")...", t, index); 

    env->apply_variation();
    
    if ( check && utils::mod(t, ae_common::backup_step) == 0 )
    {
      // check that the environment is now identical to the one stored
      // in the backup file of generation begin_gener
      
      #ifdef __REGUL
        sprintf(backup_file_name,"backup/gen_%06"PRId32".rae", t);
      #else
        sprintf(backup_file_name,"backup/gen_%06"PRId32".ae",  t);
      #endif
      
      if ( verbose )
      {
        printf("Comparing the environment with the one in %s... ", backup_file_name);  
        fflush(NULL);
      }
      
      ae_simulation * sim = new ae_simulation ( backup_file_name, false );
      stored_indiv = new ae_individual( (ae_individual *) \
                     ((sim->get_pop())->get_indivs())->get_object(index) ); // copy
      

      ae_environment * stored_env = sim->get_env();
      
      // debug
//               ae_list_node * gaussnode = env->get_gaussians()->get_first();
//               ae_list_node * storedgaussnode = stored_env->get_gaussians()->get_first();
//               ae_gaussian * gauss = NULL;
//               ae_gaussian * storedgauss = NULL;
//               while(gaussnode != NULL)
//                 {
//                   assert(storedgaussnode != NULL);
//                   gauss = (ae_gaussian *) gaussnode->get_obj();
//                   storedgauss = (ae_gaussian *) storedgaussnode->get_obj();
          
//                   printf("gaussianA : m %f, w %f, h %f\n", gauss->get_mean(), gauss->get_width(), gauss->get_height());
//                   printf("gaussianB : m %f, w %f, h %f\n", storedgauss->get_mean(), storedgauss->get_width(), storedgauss->get_height());
          
//                   gaussnode = gaussnode->get_next();
//                   storedgaussnode = storedgaussnode->get_next();
//                 }
//               assert(storedgaussnode == NULL);
      // end debug


      if ( ! (env->is_identical_to(stored_env)) )
        {
          fprintf(stderr, "ERROR: The replayed environment is not the same\n");
          fprintf(stderr, "       as the one in %s\n", backup_file_name);
          exit(EXIT_FAILURE);
        }
      
      if ( verbose ) printf("OK\n");
      delete sim;
    }

    // Warning: this portion of code won't work if the number of units changes
    // during the evolution
    
    genetic_unit_number = 0;
    dnarepnode = (rep->get_dna_replic_reports())->get_first();
    unitnode = (indiv->get_genetic_unit_list())->get_first();
    
    if ( check && utils::mod(t, ae_common::backup_step) == 0 )
    {
      storedunitnode = stored_indiv->get_genetic_unit_list()->get_first();
    }
    
    while ( dnarepnode != NULL )
    {
      assert( unitnode != NULL );

      dnarep = (ae_dna_replic_report *) dnarepnode->get_obj();
      unit = (ae_genetic_unit *) unitnode->get_obj();
      unit->get_dna()->set_replic_report( dnarep );

      mnode = dnarep->get_rearrangements()->get_first();              
      while ( mnode != NULL )
      {
        mut = (ae_mutation *) mnode->get_obj();
        
        metabolic_error_before = indiv->get_dist_to_target_by_feature( METABOLISM );
        unitlen_before = unit->get_dna()->get_length();
        seglen = mut->segment_length( unitlen_before );
        // TODO : number of affected genes

        unit->get_dna()->undergo_this_mutation( mut ); 

        indiv->reevaluate(env);
        metabolic_error_after = indiv->get_dist_to_target_by_feature( METABOLISM );
        impact_on_metabolic_error = metabolic_error_after - metabolic_error_before;

        mut->get_generic_description_string( mut_descr_string );
        fprintf( output, "%"PRId32" %"PRId32" %s %"PRId32" %"PRId32" %.15f \n",\
                 t, genetic_unit_number, \
                 mut_descr_string, seglen, \
                 unitlen_before, \
                 impact_on_metabolic_error );


        mnode = mnode->get_next();
      }
      
      mnode = dnarep->get_mutations()->get_first();              
      while ( mnode != NULL )
      {
        mut = (ae_mutation *) mnode->get_obj();

        metabolic_error_before = indiv->get_dist_to_target_by_feature( METABOLISM );
        unitlen_before = unit->get_dna()->get_length();
        seglen = mut->segment_length( unitlen_before );
        // TODO : number of affected genes

        unit->get_dna()->undergo_this_mutation( mut ); 

        indiv->reevaluate(env);
        metabolic_error_after = indiv->get_dist_to_target_by_feature( METABOLISM );
        impact_on_metabolic_error = metabolic_error_after - metabolic_error_before;

        mut->get_generic_description_string( mut_descr_string );
        fprintf( output, "%"PRId32" %"PRId32" %s %"PRId32" %"PRId32" %.15f \n",\
                 t, genetic_unit_number, \
                 mut_descr_string, seglen, \
                 unitlen_before, \
                 impact_on_metabolic_error );

        mnode = mnode->get_next();
      }

      if ( check && utils::mod(t, ae_common::backup_step) == 0 )
      {
        if ( verbose )
        {
          printf("Checking the sequence of the unit...");
          fflush(NULL);
        }
        
        assert( storedunitnode != NULL );
        storedunit = (ae_genetic_unit *) storedunitnode->get_obj();
        
        char * str1 = new char[unit->get_dna()->get_length() + 1];
        memcpy(str1, unit->get_dna()->get_data(), \
               unit->get_dna()->get_length()*sizeof(char));
        str1[unit->get_dna()->get_length()] = '\0';
        
        char * str2 = new char[(storedunit->get_dna())->get_length() + 1];
        memcpy(str2, (storedunit->get_dna())->get_data(), (storedunit->get_dna())->get_length()*sizeof(char));
        str2[(storedunit->get_dna())->get_length()] = '\0';
        
        if(strncmp(str1,str2, (storedunit->get_dna())->get_length())==0)
        {
          if ( verbose ) printf(" OK\n");
        }
        else
        {
          if ( verbose ) printf( " ERROR !\n" );
          fprintf( stderr, "Error: the rebuilt unit is not the same as \n");
          fprintf( stderr, "the one stored in backup file %s\n", backup_file_name);
          fprintf( stderr, "Rebuilt unit : %ld bp\n %s\n", strlen(str1), str1 );
          fprintf( stderr, "Stored unit  : %ld bp\n %s\n", strlen(str2), str2 );
          delete [] str1;
          delete [] str2;
          gzclose(lineage_file);
          delete indiv;
          delete stored_indiv;
          delete env;
          ae_common::clean();
          exit(EXIT_FAILURE);
        }
        
        delete [] str1;
        delete [] str2;
        
        storedunitnode = storedunitnode->get_next();
      }
      
      
      dnarepnode = dnarepnode->get_next();
      unitnode = unitnode->get_next();
      genetic_unit_number ++;
    }
        
    assert( unitnode == NULL );


    if ( verbose ) printf(" OK\n");

    delete rep;
    if ( utils::mod(t, ae_common::backup_step) == 0 )
    {
      assert(storedunitnode == NULL);
      delete stored_indiv;
    }
  }

  
  ae_common::clean();

  gzclose(lineage_file);
  fclose(output);
  delete indiv;
  delete env;


  exit(EXIT_SUCCESS);
  
}



void print_help( void )
{
  // TODO
}
