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
#include <ae_exp_manager.h>
#include <ae_individual.h>
#include <ae_genetic_unit.h>
#include <ae_list.h>
#include <ae_tree.h>
#include <ae_replication_report.h>
#include <ae_dna_replic_report.h>
#include <ae_mutation.h>
//#include <ae_param_loader.h>




enum check_type
{
  FULL_CHECK  = 0,
  LIGHT_CHECK = 1,
  NO_CHECK    = 2
};




void print_help( void );

double* dist_to_target_segment;

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
  check_type  check               = LIGHT_CHECK;
  char*       lineage_file_name   = NULL;
  bool        verbose             = false;
  bool        log                 = false;
  
  char* exp_setup_file_name = new char[63];
  char* out_prof_file_name  = new char[63];
  strcpy( exp_setup_file_name,  "exp_setup.ae" );
  strcpy( out_prof_file_name,   "output_profile.ae" );
  char* sp_struct_file_name = NULL;

  const char * short_options = "hvncf:l"; 
  static struct option long_options[] =
  {
    {"help",      no_argument,       NULL, 'h'},
    {"verbose",   no_argument,       NULL, 'v'},
    {"nocheck",   no_argument,       NULL, 'n'},
    {"fullcheck", no_argument,       NULL, 'c'},
    {"file",      required_argument, NULL, 'f'},
    {"log",         no_argument,       NULL, 'l'},
    {0, 0, 0, 0}
  };

  int option;
  while( (option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1 ) 
  {
    switch( option ) 
    {
      case 'h' : print_help(); exit(EXIT_SUCCESS);  break;
      case 'v' : verbose = true;                    break;
      case 'n' : check = NO_CHECK;                  break;
      case 'c' : check = FULL_CHECK;                break;
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
      case 'l' : log = true;                        break;
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

  // ================================
  //  Open the log file for overload
  // ================================
  
  /*ae_param_loader* log_overload = NULL;
  int32_t num_generation_overload = -10;
  
  fflush(stdout);
  
  if (log == true) 
  {
    if ( verbose)
    {
      printf( "Loading the log file from backup \n" );
    }
    
    log_overload = new ae_param_loader("log_load_from_backup.out");
    
    f_line* line;

    line = log_overload->get_line();
  
    while ( (line != NULL) && (strcmp( line->words[0], "GENERATION_OVERLOAD") != 0) ) 
    {
      line = log_overload->get_line();
    }
      
    if( line != NULL ) num_generation_overload =  atol(line->words[1]) ;
    
    delete line;
  }
  else
  {
    printf("\n");
    printf( "WARNING : Parameter change during simulation is not managed (consider -l option)\n" );
    printf("\n");
  }*/
  
  printf("\n");
  printf( "WARNING : Parameter change during simulation is not managed (consider -l option)\n" );
  printf("\n");
  
  // =======================
  //  Open the lineage file
  // =======================

  gzFile * lineage_file = (gzFile*) gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL)
  {
    fprintf( stderr, "ERROR : Could not read the lineage file %s\n", lineage_file_name );
    exit( EXIT_FAILURE );
  }

  //ae_common::read_from_backup( lineage_file );

  int32_t begin_gener, end_gener, final_index, num_gener, final_indiv_rank;
  gzread(lineage_file, &begin_gener, sizeof(begin_gener));
  gzread(lineage_file, &end_gener, sizeof(end_gener));
  gzread(lineage_file, &final_index, sizeof(final_index) );
  gzread( lineage_file, &final_indiv_rank,   sizeof(final_indiv_rank) );

  if ( verbose )
  {
    printf("\n\n");
    printf("================================================================================\n");
    printf(" Statistics of the ancestors of indiv. #%"PRId32" (t=%"PRId32" to %"PRId32")\n", final_index, begin_gener, end_gener);
    printf("================================================================================\n");
  }


  // =========================
  //  Open the output file(s)
  // =========================

  char output_file_name[60];
  snprintf(output_file_name, 60, "stats/fixedmut-b%06"PRId32"-e%06"PRId32"-i%"PRId32".out", begin_gener, end_gener, final_index);

  FILE * output = fopen(output_file_name, "w");
  if ( output == NULL )
  {
    fprintf( stderr, "ERROR : Could not create the output file %s\n", output_file_name );
    exit( EXIT_FAILURE );
  }


  // Write the header
  fprintf( output, "# #################################################################\n" );
  fprintf( output, "#  Mutations in the lineage of the best indiv at generation %"PRId32"\n", end_gener );
  fprintf( output, "# #################################################################\n" );
  fprintf( output, "#  1.  Generation       (mut. occurred when producing the indiv. of this generation)\n" );
  fprintf( output, "#  2.  Genetic unit     (which underwent the mutation, 0 = chromosome) \n" );
  fprintf( output, "#  3.  Mutation type    (0: switch, 1: smallins, 2: smalldel, 3:dupl, 4: del, 5:trans, 6:inv) \n" );
  fprintf( output, "#  4.  pos_0            (position for the small events, begin_segment for the rearrangements) \n" );
  fprintf( output, "#  5.  pos_1            (-1 for the small events, end_segment for the rearrangements) \n" );
  fprintf( output, "#  6.  pos_2            (reinsertion point for duplic., cutting point in segment for transloc., -1 for other events)\n" );
  fprintf( output, "#  7.  pos_3            (reinsertion point for transloc., -1 for other events)\n" );
  fprintf( output, "#  8.  invert           (transloc only, was the segment inverted (0/1)? (-1 for other events))\n" );
  fprintf( output, "#  9.  align_score      (score that was needed for this rearrangement to occur)\n" );
  fprintf( output, "#  10. align_score2     (transloc only, score for the reinsertion)\n" );
  fprintf( output, "#  11. segment_length   \n" );
  fprintf( output, "#  12. GU_length        (before the event)\n" );
  fprintf( output, "#  13. Impact of the mutation on the metabolic error (negative value = smaller gap after = beneficial mutation) \n" );
  fprintf( output, "####################################################################################################################\n" );
  fprintf( output, "#\n" );
  fprintf( output, "# Header for R\n" );
  fprintf( output, "gener gen_unit mut_type pos_0 pos_1 pos_2 pos_3 invert align_score align_score_2 seg_len GU_len impact\n" );
  fprintf( output, "#\n" );

  // =========================
  //  Prepare the environment
  // =========================

  if ( verbose )
  {
    printf("Preparing the environment... ");  
    fflush(NULL);
  }

  //  Open the experience manager  
  char environment_file_name[50];
  #ifdef __REGUL
  	sprintf( environment_file_name,"environment/env_%06"PRId32".rae", begin_gener );
  #else
  	sprintf( environment_file_name,"environment/env_%06"PRId32".ae",  begin_gener );
  #endif
  
  char genomes_file_name[50];
  #ifdef __REGUL
  	sprintf( genomes_file_name,"populations/pop_%06"PRId32".rae", begin_gener );
  #else
  	sprintf( genomes_file_name,"populations/pop_%06"PRId32".ae",  begin_gener );
  #endif
  
  // Load the simulation
  #ifndef __NO_X
    ae_exp_manager* exp_manager = new ae_exp_manager_X11();
  #else
    ae_exp_manager* exp_manager = new ae_exp_manager();
  #endif
  exp_manager->load_experiment( exp_setup_file_name, out_prof_file_name, environment_file_name, genomes_file_name, sp_struct_file_name, true );
  ae_environment* env = exp_manager->get_env();
  
  int32_t backup_step = exp_manager->get_backup_step();

  /*// Initialize the environment according the ae_common::... values
  // This instruction also creates the random generator and sets its
  // seed, according to the ae_common::env_seed value.

  ae_environment * env = new ae_environment();
  
  for ( num_gener = 0 ; num_gener < begin_gener ; num_gener++ )
  {
    env->apply_variation();
  }


  if ( verbose ) printf("OK\n");

  char backup_file_name[50];
  if ( check != NO_CHECK )
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
      printf( "Comparing the environment with the one in %s... ", backup_file_name );
      fflush( stdout );
    }
  
    ae_experiment * sim_backup = new ae_experiment();
    sim_backup->load_backup( backup_file_name, false, NULL );
    
    ae_environment* env_backup = sim_backup->get_env();

    if ( ! env->is_identical_to(env_backup) )
    {
      fprintf(stderr, "ERROR: The replayed environment is not the same\n");
      fprintf(stderr, "       as the one in %s\n", backup_file_name);
      exit(EXIT_FAILURE);
    }
    
    if ( verbose ) printf("OK\n");
    delete sim_backup; 
  }*/

 


  // ==============================
  //  Prepare the initial ancestor
  // ==============================

  ae_individual * indiv = new ae_individual(exp_manager, lineage_file );
  indiv->evaluate( env );
  indiv->compute_statistical_data();
  
  if ( verbose )
  {
    printf("Initial fitness     = %f\n", indiv->get_fitness());
    printf("Initial genome size = %"PRId32"\n", indiv->get_total_genome_size());
  }

  delete exp_manager;


  // ===============================================================================
  //  Replay the mutation to get the successive ancestors and analyze them
  //  (and, optionally, check that the rebuilt envir and genome are correct each 
  //   time a backup is available)
  // ===============================================================================
  ae_replication_report * rep   = NULL;
  ae_list_node * dnarepnode     = NULL;
  ae_dna_replic_report * dnarep = NULL;

  ae_list_node * mnode  = NULL;
  ae_mutation *  mut    = NULL;

  ae_list_node*   unitnode  = NULL;
  ae_genetic_unit *  unit   = NULL;

  ae_individual *     stored_indiv    = NULL;
  ae_list_node*       stored_unit_node  = NULL;
  ae_genetic_unit *   stored_unit      = NULL;

  int32_t i, index, genetic_unit_number, unitlen_before, seglen;
  double metabolic_error_before, metabolic_error_after, impact_on_metabolic_error;
  char mut_descr_string[80];
  
  bool check_now = false;

  for ( i = 0; i < end_gener - begin_gener; i++ )
  {
    num_gener = begin_gener + i + 1;  // where are we in time...
       
    /*if (log == true)
    {
      // overload of environment variations
      while (num_gener == num_generation_overload)
      {
        f_line* line;
        int16_t nb_param_overloaded;
        
        line = log_overload->get_line();
        if (strcmp( line->words[0], "NB_PARAM_OVERLOADED") == 0)
        {
          nb_param_overloaded = atol( line->words[1] );
        }
        else
        {
          printf( "ERROR in log file for overload : unknown number of parameters overloaded\n");
          exit( EXIT_FAILURE );
        }
	
	
        for ( int16_t i = 0 ; i < nb_param_overloaded ; i++ )
        {
          line = log_overload->get_line();
          if (strcmp( line->words[0], "ENV_VARIATION") == 0)
          {
            if ( strcmp( line->words[1], "none" ) == 0 )
            {
              env->set_variation_method( NONE );
            }
            else if ( strcmp( line->words[1], "autoregressive_mean_variation" ) == 0 )
            {
              env->set_variation_method( AUTOREGRESSIVE_MEAN_VAR );
              env->set_var_sigma_tau( atof( line->words[2] ), atol( line->words[3] ) );
            }
            else if ( strcmp( line->words[1], "add_local_gaussians" ) == 0 )
            {
              env->set_variation_method( LOCAL_GAUSSIANS_VAR );
            }
            else
            {
              printf( "ERROR in log file for overload : unknown environment variation method\n" );
              exit( EXIT_FAILURE );
            }
          }
          else if ( strcmp( line->words[0], "ENV_AXIS_SEGMENTS") == 0 || strcmp( line ->words[0], "ENV_AXIS_FEATURES") == 0 )
          {
            // Static variables for env axis segmentation
            static int16_t              env_axis_nb_segments        = -1;
            static double*              env_axis_segment_boundaries = NULL;
            static ae_env_axis_feature* env_axis_features           = NULL;
            static bool                 env_axis_separate_segments  = false;
            
            // Initialize static variables for env axis segmentation
            if ( env_axis_segment_boundaries == NULL && strcmp( line->words[0], "ENV_AXIS_SEGMENTS") == 0 )
            {
              if ( env_axis_features != NULL )
              {
                if ( line->nb_words != env_axis_nb_segments + 1 )
                {
                  printf( "ERROR : Number of segments defined by ENV_AXIS_SEGMENTS and ENV_AXIS_FEATURES do not match. \"%s\".\n", line->words[i] );
                }
              }
              else
              {
                env_axis_nb_segments = line->nb_words - 1;
              }
              
              env_axis_segment_boundaries = new double[env_axis_nb_segments + 1];
            }
            else if ( env_axis_features == NULL && strcmp( line ->words[0], "ENV_AXIS_FEATURES") == 0 )
            {
              if ( env_axis_segment_boundaries != NULL )
              {
                if ( line->nb_words != env_axis_nb_segments )
                {
                  printf( "ERROR : Number of segments defined by ENV_AXIS_SEGMENTS and ENV_AXIS_FEATURES do not match. \"%s\".\n", line->words[i] );
                }
              }
              else
              {
                env_axis_nb_segments = line->nb_words;
              }
              
              env_axis_features = new ae_env_axis_feature[env_axis_nb_segments];
            }
            
            // Set temp data (static variables) according to values in param file
            if ( strcmp( line->words[0], "ENV_AXIS_SEGMENTS") == 0 )
            {
              env_axis_segment_boundaries[0] = MIN_X;
              for ( int16_t i = 1 ; i < env_axis_nb_segments ; i++ )
              {
                env_axis_segment_boundaries[i] = atof( line->words[i] );
              }
              env_axis_segment_boundaries[env_axis_nb_segments] = MAX_X;
            }
            else // ( strcmp( line ->words[0], "ENV_AXIS_FEATURES") == 0 )
            {
              for ( int16_t i = 0 ; i < env_axis_nb_segments ; i++ )
              {
                if ( strcmp( line->words[(i+1)], "NEUTRAL" ) == 0 )
                {
                  env_axis_features[i] = NEUTRAL;
                }
                else if ( strcmp( line->words[(i+1)], "METABOLISM" ) == 0 )
                {
                  env_axis_features[i] = METABOLISM;
                }
                else if ( strcmp( line->words[(i+1)], "SECRETION" ) == 0 )
                {
                  ae_common::params->set_use_secretion( true );
                  env_axis_features[i] = SECRETION;
                }
                else if ( strcmp( line->words[(i+1)], "TRANSFER" ) == 0 )
                {
                  env_axis_features[i] = TRANSFER;
                }
                else
                {
                  printf( "ERROR : unknown axis feature \"%s\".\n",line->words[i] );
                  exit( EXIT_FAILURE );
                }
              }
            }
            
            // If both boundaries and features have been defined, flush into environment
            if ( env_axis_segment_boundaries != NULL && env_axis_features != NULL )
            {
              env->set_segmentation( env_axis_nb_segments, env_axis_segment_boundaries, env_axis_features, env_axis_separate_segments );
              // do not delete env_axis_segment_boundaries and env_axis_features, no in-depth copy !
              
              env_axis_segment_boundaries = NULL;
              env_axis_features = NULL;
            }
          }
          else if (strcmp(line ->words[0], "ENV_ADD_GAUSSIAN") == 0)
          {
             env->add_gaussian( atof( line->words[1] ),
                                atof( line->words[2] ),
                                atof( line->words[3] ) );
          }
          
        }
	
        while ( (line != NULL) && (strcmp( line->words[0], "GENERATION_OVERLOAD") != 0) ) 
        {
          line = log_overload->get_line();
        }
          
         
        if( line != NULL )
        {
          num_generation_overload =  atol(line->words[1]) ;
          if( num_generation_overload < num_gener)
          {
              printf( "ERROR in log file for overload : overload of an anterior generation\n" );
              exit( EXIT_FAILURE );
          }
        }
        else
        {
          num_generation_overload =  -10 ;
        }
	
        // If there was a modification in the segmentation of the environment
        // TODO : Check that (it wasn't correctly tested and might not be functional
        //~ if ( ae_common::env_axis_is_segmented )
        //~ {
          //~ dist_to_target_segment = new double [ae_common::init_params->get_env_axis_nb_segments()];
          
          //~ for ( int16_t i = 0 ; i < ae_common::init_params->get_env_axis_nb_segments() ; i++ )
          //~ {
            //~ dist_to_target_segment[i] = 0;
          //~ }
          //~ indiv->reset_dist_to_target_segment( dist_to_target_segment );
          
          //~ indiv->renew_dist_to_target_by_feature();
          //~ indiv->renew_fitness_by_feature();
        //~ }
        
        delete line; 
      }
    }*/
    
    env->build();

    rep = new ae_replication_report( lineage_file, indiv );
    index = rep->get_id(); // who are we building...
    indiv->set_replication_report(rep);
    
    // Check now?
    check_now = ( ( check == FULL_CHECK && ae_utils::mod( num_gener, backup_step ) == 0 ) ||
                  ( check == LIGHT_CHECK && num_gener == end_gener ) );


    if ( verbose ) printf("Rebuilding ancestor at generation %"PRId32" (index %"PRId32")...", num_gener, index); 

    env->apply_variation();
    
    if ( check_now && ae_utils::mod(num_gener, backup_step) == 0 )
    {
      // check that the environment is now identical to the one stored
      // in the backup file of generation begin_gener
      
      char environment_file_name[50];
      #ifdef __REGUL
  		sprintf( environment_file_name,"environment/env_%06"PRId32".rae", begin_gener );
      #else
  	  	sprintf( environment_file_name,"environment/env_%06"PRId32".ae",  begin_gener );
      #endif
      
      char genomes_file_name[50];
      #ifdef __REGUL
  		sprintf( genomes_file_name,"populations/pop_%06"PRId32".rae", begin_gener );
  	  #else
  		sprintf( genomes_file_name,"populations/pop_%06"PRId32".ae",  begin_gener );
      #endif
      
      // Load the simulation
      #ifndef __NO_X
      	exp_manager = new ae_exp_manager_X11();
      #else
      	exp_manager = new ae_exp_manager();
      #endif
      exp_manager->load_experiment( exp_setup_file_name, out_prof_file_name, environment_file_name, genomes_file_name, sp_struct_file_name, true );
      ae_environment* backup_env = exp_manager->get_env();
      stored_indiv = new ae_individual( * (ae_individual *)exp_manager->get_indiv_by_id( index ) );
      delete exp_manager;
      
      if ( verbose )
      {
        printf("Comparing the environment with the one in %s... \n", environment_file_name);  
        fflush(NULL);
      }
      
      if ( ! env->is_identical_to(backup_env) )
      {
        fprintf(stderr, "ERROR: The replayed environment is not the same\n");
        fprintf(stderr, "       as the one in %s\n", environment_file_name);
        exit(EXIT_FAILURE);
      }
      
      if ( verbose ) printf("OK\n");
    }

    // Warning: this portion of code won't work if the number of units changes
    // during the evolution
    
    genetic_unit_number = 0;
    dnarepnode = (rep->get_dna_replic_reports())->get_first();
    unitnode   = (indiv->get_genetic_unit_list())->get_first();
    
    if ( check_now && ae_utils::mod(num_gener, backup_step) == 0)
    {
      stored_unit_node = stored_indiv->get_genetic_unit_list()->get_first();
    }
    
    while ( dnarepnode != NULL )
    {
      assert( unitnode != NULL );

      dnarep = (ae_dna_replic_report *) dnarepnode->get_obj();
      unit   = (ae_genetic_unit *) unitnode->get_obj();
      
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
                 num_gener, genetic_unit_number, \
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
                 num_gener, genetic_unit_number, \
                 mut_descr_string, seglen, \
                 unitlen_before, \
                 impact_on_metabolic_error );

        mnode = mnode->get_next();
      }

      if ( check_now && ae_utils::mod(num_gener, backup_step) == 0)
      {
        if ( verbose )
        {
          printf("Checking the sequence of the unit...");
          fflush(NULL);
        }
        
        assert( stored_unit_node != NULL );
        stored_unit = (ae_genetic_unit *) stored_unit_node->get_obj();
        
        char * str1 = new char[unit->get_dna()->get_length() + 1];
        memcpy(str1, unit->get_dna()->get_data(), \
               unit->get_dna()->get_length()*sizeof(char));
        str1[unit->get_dna()->get_length()] = '\0';
        
        char * str2 = new char[(stored_unit->get_dna())->get_length() + 1];
        memcpy(str2, (stored_unit->get_dna())->get_data(), (stored_unit->get_dna())->get_length()*sizeof(char));
        str2[(stored_unit->get_dna())->get_length()] = '\0';
        
        if(strncmp(str1,str2, (stored_unit->get_dna())->get_length())==0)
        {
          if ( verbose ) printf(" OK\n");
        }
        else
        {
          if ( verbose ) printf( " ERROR !\n" );
          fprintf( stderr, "Error: the rebuilt unit is not the same as \n");
          fprintf( stderr, "the one stored in backup file %s\n", genomes_file_name);
          fprintf( stderr, "Rebuilt unit : %zu bp\n %s\n", strlen(str1), str1 );
          fprintf( stderr, "Stored unit  : %zu bp\n %s\n", strlen(str2), str2 );
          delete [] str1;
          delete [] str2;
          gzclose(lineage_file);
          delete indiv;
          delete stored_indiv;
          delete env;
          //ae_common::clean();
          exit(EXIT_FAILURE);
        }
        
        delete [] str1;
        delete [] str2;
        
        stored_unit_node = stored_unit_node->get_next();
      }
      
      
      dnarepnode = dnarepnode->get_next();
      unitnode = unitnode->get_next();
      genetic_unit_number ++;
    }
        
    assert( unitnode == NULL );


    if ( verbose ) printf(" OK\n");

    delete rep;
    
    if ( ae_utils::mod(num_gener, backup_step) == 0 )
    {
      assert(stored_unit_node == NULL);
      delete stored_indiv;
    }
  }

  
  //ae_common::clean();

  gzclose(lineage_file);
  fclose(output);
  delete indiv;
  delete env;

  delete [] exp_setup_file_name;
  delete [] out_prof_file_name;

  exit(EXIT_SUCCESS);
  
}



void print_help( void )
{
  printf( "\n" ); 
  printf( "*********************** aevol - Artificial Evolution ******************* \n" );
  printf( "*                                                                      * \n" );
  printf( "*               Fixed mutations post-treatment program                 * \n" );
  printf( "*                                                                      * \n" );
  printf( "************************************************************************ \n" );
  printf( "\n\n" ); 
  printf( "This program is Free Software. No Warranty.\n" );
  printf( "Copyright (C) 2009  LIRIS.\n" );
  printf( "\n" ); 
#ifdef __REGUL
  printf( "Usage : rfixed_mutations -h\n");
  printf( "or :    rfixed_mutations [-vn] -f lineage_file \n" );
#else
  printf( "Usage : fixed_mutations -h\n");
  printf( "or :    fixed_mutations [-vn] -f lineage_file \n" );
#endif
  printf( "\n" ); 
  printf( "This program computes the fixed mutations of the individuals within the lineage\n" );
  printf( "of lineage_file\n" );
  printf( "\n" ); 
  printf( "WARNING: This program should not be used for simulations run with lateral\n" ); 
  printf( "transfer. When an individual has more than one parent, the notion of lineage\n" ); 
  printf( "used here is not relevant.\n" );
  printf( "\n" );  
  printf( "\t-h or --help       : Display this help.\n" );
  printf( "\n" ); 
  printf( "\t-v or --verbose    : Be verbose, listing generations as they are \n" );
  printf( "\t                       treated.\n" );
  printf( "\n" );
  printf( "\t-n or --nocheck    : Disable genome sequence checking. Makes the \n"); 
  printf( "\t                       program faster, but it is not recommended. \n");
  printf( "\t                       It is better to let the program check that \n");
  printf( "\t                       when we rebuild the genomes of the ancestors\n");
  printf( "\t                       from the lineage file, we get the same sequences\n");
  printf( "\t                       as those stored in the backup files.\n" );
  printf( "\n" );
  printf( "\t-c or --fullcheck  : Will perform the genome and environment checks every\n" );
  printf( "\t                       <BACKUP_STEP> generations. Default behaviour is\n" );
  printf( "\t                       lighter as it only perform sthese checks at the\n" );
  printf( "\t                       ending generation.\n" );
  printf( "\n" ); 
  printf( "\t-f lineage_file or --file lineage_file : \n" );
  printf( "\t                       Compute the fixed mutations of the individuals within lineage_file.\n" );
  printf( "\n" );
  printf( "\t-l or --log        : Will take on account the parameter change during\n");
  printf( "\t                       the simulation (rerun from backup) by loading the \n" );
  printf( "\t                       file log_load_from_backup.out, generated with the option \n" );
  printf( "\t                       log = load in param.in\n" );
  printf( "\n" );
}
