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
#include <ae_protein.h>
#include <ae_rna.h>
#include <ae_list.h>
#include <ae_experiment.h>
#include <ae_align.h>



// ========
//  TO DO 
// ========
// 
//  * option --color ?
//  * Raevol-specific output (EPS file with the network) ?





void print_help( void );





int main( int argc, char* argv[] )
{
  // =================================================================
  //                      Get command-line options
  // =================================================================
  //
  // 1) Initialize command-line option variables with default values
  bool  verbose               = false;
  bool  use_single_indiv_file = false;
  char* backup_file_name      = NULL;
  
  // 2) Define allowed options
  const char * options_list = "hvf:";
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
        print_help();
        exit( EXIT_SUCCESS );
        break;
      case 'v' :
        verbose = true;
        break;
      case 'f' :
        if ( strstr( optarg, "best" ) != NULL )
        {
          use_single_indiv_file = true;
        }
        else
        {
          printf( "%s|\n", optarg );
        }
        backup_file_name = new char[strlen(optarg) + 1];
        sprintf( backup_file_name, "%s", optarg );
        break;
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
    if ( use_single_indiv_file )
    {
      // TODO : best* backups don't look right...
      printf( "Reading single individual backup file <%s>... ", backup_file_name );
      gzFile* backup_file = (gzFile*) gzopen( backup_file_name, "r" );
      ae_common::read_from_backup( backup_file, verbose );
      env = new ae_environment(); // Uses the ae_common data
      best_indiv = new ae_individual( backup_file );
      
      num_gener = -1; // TODO!!!
      printf("done\n");
    }
    else
    {
      printf( "Reading backup file <%s>... ", backup_file_name );
      fflush(stdout);

      // Load simulation from backup
      ae_common::sim = new ae_experiment();
      ae_common::sim->load_backup( backup_file_name, false, NULL );
      
      best_indiv      = ae_common::pop->get_best();
      env             = ae_common::sim->get_env();
      num_gener       = ae_common::sim->get_num_gener();
      printf("done\n");
    }
  }
  
  delete [] backup_file_name;
  
  
  
  // The constructor of the ae_experiment has read the genomes of the individuals
  // and located their promoters, but has not performed the translation nor the
  // phenotype computation. We must do it now.
  // However, as the individuals in the backups are sorted, we don't need to evaluate
  // all the individuals, only those we are interested in (here only the best one)
      
  best_indiv->evaluate( env );
  
  
  
  
  // =================================================================
  //                     Open output file
  // =================================================================

  char output_file_name[60];
  snprintf( output_file_name, 60, "mutagenesis-g%06"PRId32".out", num_gener );

  FILE * output = fopen(output_file_name, "w");
  if ( output == NULL )
  {
    fprintf( stderr, "ERROR : Could not create the output file %s\n", output_file_name );
    exit( EXIT_FAILURE );
  }


  // Write the header
  fprintf( output, "# #################################################################\n" );
  fprintf( output, "#              Mutations produced by mutagenesis\n" );
  fprintf( output, "# #################################################################\n" );
  fprintf( output, "#  1.  ---              Padded line (generation in similar files)\n" );
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
  


  
  // =================================================================
  //                     Proceed to mutagenesis
  // =================================================================
  // Set some data
  ae_individual*  initial_indiv = best_indiv;
  ae_dna*         initial_dna   = initial_indiv->get_genetic_unit(0)->get_dna();
  int32_t         initial_len   = initial_dna->get_length();
  double          initial_metabolic_error = initial_indiv->get_dist_to_target_by_feature( METABOLISM );
  
  ae_mutation* mut = NULL;
  char mut_descr_string[80];
  
  double final_metabolic_error      = 0.0;
  double impact_on_metabolic_error  = 0.0;
  
  ae_individual*  indiv       = NULL;
  ae_vis_a_vis*   alignment_1 = NULL;
  ae_vis_a_vis*   alignment_2 = NULL;
  int32_t         nb_pairs;
  
  int16_t i;
  bool    rear_done;
  
  
  // Perform 100 duplications (reinitializing the genome after each of them)
  for ( i = 0 ; i < 1000 ; i++ )
  {
    if ( i % 100 == 0 )
    {
      printf( "*" );
      fflush(stdout);
    }
    rear_done   = false;
    
    indiv = new ae_individual( *best_indiv );
    assert( indiv->get_dist_to_target_by_feature( METABOLISM ) == best_indiv->get_dist_to_target_by_feature( METABOLISM ) );
    
    do
    {
      nb_pairs = initial_len;
      alignment_1 = ae_dna::search_alignment( initial_dna, initial_dna, nb_pairs, DIRECT );
    }
    while ( alignment_1 == NULL );
    
    
    // Remember the length of the segment to be duplicated and of the former genome
    int32_t segment_length = ae_utils::mod( alignment_1->get_i_2() - alignment_1->get_i_1(), initial_len );
    
    // Perform in situ (tandem) DUPLICATION
    rear_done = indiv->get_genetic_unit(0)->get_dna()->do_duplication( alignment_1->get_i_1(), alignment_1->get_i_2(), alignment_1->get_i_2() );
    
    if ( rear_done )
    {
      // Create a temporary report for the duplication
      mut = new ae_mutation();
      mut->report_duplication( alignment_1->get_i_1(), alignment_1->get_i_2(), alignment_1->get_i_2(), segment_length, alignment_1->get_score() );
      
      // Evaluate the mutated individual and write a line in the output file
      indiv->reevaluate(env);
      final_metabolic_error     = indiv->get_dist_to_target_by_feature( METABOLISM );
      impact_on_metabolic_error = final_metabolic_error - initial_metabolic_error;

      mut->get_generic_description_string( mut_descr_string );
      fprintf( output, "%"PRId32" %"PRId32" %s %"PRId32" %"PRId32" %.15f \n",\
               -1, 0, mut_descr_string, segment_length, \
               initial_len, impact_on_metabolic_error );
    }
    
    delete indiv;
  }
  printf( "\n" );
  
  // Perform 100 deletions (reinitializing the genome after each of them)
  for ( i = 0 ; i < 1000 ; i++ )
  {
    if ( i % 100 == 0 )
    {
      printf( "*" );
      fflush(stdout);
    }
    rear_done = false;
    
    indiv = new ae_individual( *best_indiv );
    assert( indiv->get_dist_to_target_by_feature( METABOLISM ) == best_indiv->get_dist_to_target_by_feature( METABOLISM ) );
    
    do
    {
      nb_pairs = initial_len;
      alignment_1 = ae_dna::search_alignment( initial_dna, initial_dna, nb_pairs, DIRECT );
    }
    while ( alignment_1 == NULL );
    
    
    // Remember the length of the segment to be duplicated and of the former genome
    int32_t segment_length = ae_utils::mod( alignment_1->get_i_2() - alignment_1->get_i_1(), initial_len );
    
    // Perform DELETION
    rear_done = indiv->get_genetic_unit(0)->get_dna()->do_deletion( alignment_1->get_i_1(), alignment_1->get_i_2() );
    
    if ( rear_done )
    {
      // Create a temporary report for the duplication
      mut = new ae_mutation();
      mut->report_deletion( alignment_1->get_i_1(), alignment_1->get_i_2(), segment_length, alignment_1->get_score() );
      
      // Evaluate the mutated individual and write a line in the output file
      indiv->reevaluate(env);
      final_metabolic_error     = indiv->get_dist_to_target_by_feature( METABOLISM );
      impact_on_metabolic_error = final_metabolic_error - initial_metabolic_error;

      mut->get_generic_description_string( mut_descr_string );
      fprintf( output, "%"PRId32" %"PRId32" %s %"PRId32" %"PRId32" %.15f \n",\
               -1, 0, mut_descr_string, segment_length, \
               initial_len, impact_on_metabolic_error );
      
      delete mut;
    }
    
    delete indiv;
  }
  printf( "\n" );
  
  
  // Perform 100 inversions (reinitializing the genome after each of them)
  for ( i = 0 ; i < 1000 ; i++ )
  {
    if ( i % 100 == 0 )
    {
      printf( "*" );
      fflush(stdout);
    }
    rear_done = false;
    
    indiv = new ae_individual( *best_indiv );
    assert( indiv->get_dist_to_target_by_feature( METABOLISM ) == best_indiv->get_dist_to_target_by_feature( METABOLISM ) );
    
    do
    {
      nb_pairs = initial_len;
      alignment_1 = ae_dna::search_alignment( initial_dna, initial_dna, nb_pairs, INDIRECT );
    } while ( alignment_1 == NULL );
    
    // Make sure the segment to be inverted doesn't contain OriC
    if ( alignment_1->get_i_1() > alignment_1->get_i_2() )
    {
      alignment_1->swap();
    }
    
    // Remember the length of the segment to be duplicated and of the former genome
    int32_t segment_length = ae_utils::mod( alignment_1->get_i_2() - alignment_1->get_i_1(), initial_len );
    
    // Perform INVERSION
    rear_done = indiv->get_genetic_unit(0)->get_dna()->do_inversion( alignment_1->get_i_1(), alignment_1->get_i_2() );
    
    if ( rear_done )
    {
      // Create a temporary report for the duplication
      mut = new ae_mutation();
      mut->report_inversion( alignment_1->get_i_1(), alignment_1->get_i_2(), segment_length, alignment_1->get_score() );
      
      // Evaluate the mutated individual and write a line in the output file
      indiv->reevaluate(env);
      final_metabolic_error     = indiv->get_dist_to_target_by_feature( METABOLISM );
      impact_on_metabolic_error = final_metabolic_error - initial_metabolic_error;

      mut->get_generic_description_string( mut_descr_string );
      fprintf( output, "%"PRId32" %"PRId32" %s %"PRId32" %"PRId32" %.15f \n",\
               -1, 0, mut_descr_string, segment_length, \
               initial_len, impact_on_metabolic_error );
      
      delete mut;
    }
    
    delete indiv;
  }
  printf( "\n" );
  
  
  // Perform 100 translocations (reinitializing the genome after each of them)
  for ( i = 0 ; i < 1000 ; i++ )
  {
    if ( i % 100 == 0 )
    {
      printf( "*" );
      fflush(stdout);
    }
    rear_done   = false;
    
    // Make a working copy of the initial individual
    indiv = new ae_individual( *best_indiv );
    assert( indiv->get_dist_to_target_by_feature( METABOLISM ) == best_indiv->get_dist_to_target_by_feature( METABOLISM ) );
    
    // Look for an alignment
    do
    {
      nb_pairs = initial_len;
      alignment_1 = ae_dna::search_alignment( initial_dna, initial_dna, nb_pairs, DIRECT );
    }
    while ( alignment_1 == NULL );
    
    // Remember the length of the segment to be duplicated and of the former genome
    int32_t segment_length = ae_utils::mod( alignment_1->get_i_2() - alignment_1->get_i_1(), initial_len );
    
    // Extract the segment to be translocated
    ae_genetic_unit* tmp_segment = indiv->get_genetic_unit(0)->get_dna()->extract_into_new_GU( alignment_1->get_i_1(), alignment_1->get_i_2() );
    
    // Look for a "new" alignments between this segment and the remaining of the chromosome
    do
    {
      nb_pairs = initial_len;
      alignment_2 = ae_dna::search_alignment( tmp_segment->get_dna(), indiv->get_genetic_unit(0)->get_dna(), nb_pairs, BOTH_SENSES );
    }
    while ( alignment_2 == NULL );
    
    // Reinsert the segment into the chromosome
    indiv->get_genetic_unit(0)->get_dna()->insert_GU( tmp_segment, alignment_2->get_i_2(), alignment_2->get_i_1(), (alignment_2->get_sense() == INDIRECT) );
    rear_done = true;
    
    if ( rear_done )
    {
      // Create a temporary report for the duplication
      mut = new ae_mutation();
      mut->report_translocation( alignment_1->get_i_1(), alignment_1->get_i_2(), alignment_2->get_i_1(), alignment_2->get_i_2(),
                                 segment_length, alignment_1->get_score(), alignment_2->get_score() );
      
      // Evaluate the mutated individual and write a line in the output file
      indiv->reevaluate(env);
      final_metabolic_error     = indiv->get_dist_to_target_by_feature( METABOLISM );
      impact_on_metabolic_error = final_metabolic_error - initial_metabolic_error;

      mut->get_generic_description_string( mut_descr_string );
      fprintf( output, "%"PRId32" %"PRId32" %s %"PRId32" %"PRId32" %.15f \n",\
               -1, 0, mut_descr_string, segment_length, \
               initial_len, impact_on_metabolic_error );
      
      delete mut;
    }
    
    delete indiv;
  }
  printf( "\n" );


  
  if ( use_single_indiv_file )
  {
    delete best_indiv;
    delete env;
  }
  else
  {
    delete ae_common::sim;
  }

  return EXIT_SUCCESS;
}




void print_help( void ) 
{
}




